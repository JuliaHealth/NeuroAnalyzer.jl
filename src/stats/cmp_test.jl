export cmp_test

"""
    cmp_test(seg1, seg2; <keyword arguments>)

Compare two vectors using an automatically selected or explicitly requested statistical test.

Test selection (when `type = :auto`):

1. Jarque–Bera test is applied to the pooled data to assess normality.
2. If normal: variances are compared with an F-test; an equal- or unequal-variance t-test (or a paired one-sample t-test) is used.
3. If non-normal: Wilcoxon signed-rank test (paired) or Mann–Whitney U test (unpaired) is used.

# Arguments

- `s1::AbstractVector`: first signal vector
- `s2::AbstractVector`: second signal vector; must have the same length as `s1` when `paired = true`
- `paired::Bool`" if `true`, treat the vectors as paired observations"
- `alpha::Float64=0.05`: significance level; must be in `(0, 1)`
- `type::Symbol=:auto`: test selection strategy:
    - `:auto`: choose automatically based on normality and variance tests
    - `:perm`: permutation-based test
    - `:p`: force parametric test
    - `:np`: force non-parametric test
- `exact::Bool=false`: if `true`, use the exact Wilcoxon signed-rank test (only applies when `type ∈ {:auto, :np}` and `paired = true`)
- `nperm::Int64=1000`: number of permutations for `type = :perm`; must be ≥ 1
- `verbose::Bool=true`: if `true`, print test selection and result details

# Returns

Permutation test (`type = :perm`), named tuple:

- `t::NamedTuple{perm_diff, obs_diff}`: permutation null distribution and observed difference
- `p1::Float64`: one-tailed p-value
- `p2::Float64`: two-tailed p-value (`2 × p1`, clamped to 1.0)

All other tests, named tuple:

- `t`: the HypothesisTests.jl test object
- `ts::Tuple{Float64, String}`: test statistic value and its name
- `tc`: confidence interval of the test statistic (tuple or scalar NaN for non-parametric)
- `df::Float64`: degrees of freedom
- `p::Float64`: two-tailed p-value (clamped to `eps()` if below machine epsilon)

# Throws

- `ArgumentError`: if `type` is invalid, `alpha ∉ (0, 1)`, `nperm < 1`, or `paired = true` with unequal-length vectors
"""
function cmp_test(
    s1::AbstractVector,
    s2::AbstractVector;
    paired::Bool,
    alpha::Float64 = 0.05,
    type::Symbol = :auto,
    exact::Bool = false,
    nperm::Int64 = 1000,
    verbose::Bool = true,
)::Union{
    @NamedTuple{t::OneSampleTTest, ts::Tuple{Float64, String}, tc::Tuple{Float64, Float64}, df::Float64, p::Float64},
    @NamedTuple{t::EqualVarianceTTest, ts::Tuple{Float64, String}, tc::Tuple{Float64, Float64}, df::Float64, p::Float64},
    @NamedTuple{t::UnequalVarianceTTest, ts::Tuple{Float64, String}, tc::Tuple{Float64, Float64}, df::Float64, p::Float64},
    @NamedTuple{t::ExactSignedRankTest{Float64}, ts::Tuple{Float64, String}, tc::Float64, df::Float64, p::Float64},
    @NamedTuple{t::ApproximateSignedRankTest{Float64}, ts::Tuple{Float64, String}, tc::Float64, df::Float64, p::Float64},
    @NamedTuple{t::ExactMannWhitneyUTest{Float64}, ts::Tuple{Float64, String}, tc::Float64, df::Float64, p::Float64},
    @NamedTuple{t::ApproximateMannWhitneyUTest{Float64}, ts::Tuple{Float64, String}, tc::Float64, df::Float64, p::Float64},
    @NamedTuple{t::@NamedTuple{perm_diff::Vector{Float64}, obs_diff::Float64}, p1::Float64, p2::Float64},
}

    _check_var(type, [:auto, :perm, :p, :np], "type")
    @assert alpha > 0.0  "alpha must be > 0."
    @assert alpha < 1.0  "alpha must be < 1."
    @assert nperm >= 1   "nperm must be ≥ 1."
    paired && @assert length(s1) == length(s2) "Paired test requires equal-length vectors."

    # --- normality test (Jarque–Bera on pooled data) ---
    jb = JarqueBeraTest([s1; s2])
    pjb = round(pvalue(jb), digits=3)
    jb_stat = round(jb.JB, digits=3)
    if verbose
        dist_label = pjb < alpha ? "non-normal" : "normal"
        println("Distribution: $dist_label (Jarque-Bera test: JB=$jb_stat, p=$pjb)")
    end

    if type !== :perm

        # --- parametric or non-parametric branch ---
        use_parametric = (pjb >= alpha && type === :auto) || type === :p

        if use_parametric

            # parametric branch
            if paired

                verbose && println("Using one-sample T-test (paired data)")
                t  = OneSampleTTest(s1, s2)

            else

                pf = pvalue(VarianceFTest(s1, s2))
                if pf >= alpha
                    verbose && println("Using equal-variance two-sample T-test")
                    t = EqualVarianceTTest(s1, s2)
                else
                    verbose && println("Using unequal-variance two-sample T-test (Welch)")
                    t = UnequalVarianceTTest(s1, s2)
                end

            end
            df = t.df
            ts = t.t
            tc = confint(t; level=(1 - alpha))
            tn = "t"

        else

            # non-parametric branch
            if paired

                if exact
                    verbose && println("Using exact Wilcoxon signed-rank test (paired data)")
                    t = ExactSignedRankTest(s1, s2)
                else
                    verbose && println("Using approximate Wilcoxon signed-rank test (paired data)")
                    t = SignedRankTest(s1, s2)
                end
                ts = t.W
                df = t.n - 1
                tn = "W"

            else

                verbose && println("Using Mann-Whitney U test")
                t  = MannWhitneyUTest(s1, s2)
                ts = t.U
                df = length(s1) + length(s2) - 2
                tn = "U"

            end
            tc = NaN

        end

        p = pvalue(t)
        p < eps() && (p = eps())
        return (
            t = t,
            ts = (round(ts, digits=3), tn),
            tc = round.(tc, digits=3),
            df = round(df, digits=3),
            p = p
        )

    else

        # --- permutation test branch ---
        n1 = length(s1)
        n2 = length(s2)
        g  = [s1; s2]

        verbose && println("Group 1 mean: $(round(mean(s1), digits=3))")
        verbose && println("Group 2 mean: $(round(mean(s2), digits=3))")

        # observed difference (group 1 − group 2)
        observed_diff = round(mean(s1) - mean(s2), digits=3)

        # build null distribution by repeatedly shuffling group labels
        perm_diff = zeros(nperm)
        progbar   = Progress(nperm, dt=1, barlen=20, color=:white, enabled=progress_bar)

        @inbounds for idx in 1:nperm
            # random label assignment: first n1 positions → group 0, rest → group 1
            shuffle_idx = randperm(n1 + n2)
            perm_diff[idx] = mean(g[shuffle_idx[1:n1]]) - mean(g[shuffle_idx[(n1 + 1):end]])
            progress_bar && next!(progbar)
        end

        # one-tailed p-value: proportion of permuted differences ≥ observed
        if pjb < alpha
            # non-normal: use empirical proportion directly
            verbose && println("H0 distribution is non-normal (Jarque-Bera p=$pjb); using empirical p-value.")
            p1 = sum(perm_diff .>= observed_diff) / nperm
        else
            # normal null: standardise and use Z-score
            z  = round((observed_diff - mean(perm_diff)) / std(perm_diff), digits=3)
            verbose && println("Permutation Z-score: $z")
            p1 = 1 - cdf(Distributions.Normal(), abs(z))
        end

        # two-tailed p-value = 2 × one-tailed, clamped to 1.0
        p2 = min(1.0, 2 * p1)

        return (
            t  = (perm_diff=perm_diff, obs_diff=observed_diff),
            p1 = p1,
            p2 = p2,
        )

    end

end
