export cmp_test

"""
    cmp_test(seg1, seg2, paired, alpha, type, exact)

Compare two vectors; Kruskall-Wallis test is used first, next t-test (paired on non-paired) or non-parametric test (paired: Wilcoxon signed rank, non-paired: Mann-Whitney U test) is applied.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`
- `paired::Bool`
- `alpha::Float64=0.05`: confidence level
- `type::Symbol=:auto`: choose test automatically (`:auto`), permutation-based (`:perm`), parametric (`:p`) or non-parametric (`:np`)
- `exact::Bool=false`: if true, use exact Wilcoxon test
- `nperm::Int64=1000`: number of permutation for `:perm` method

# Returns

Named tuple containing for type !== `:perm`:
- `t`: test results
- `ts::Tuple{Float64, String}`: test statistics
- `tc::Tuple{Float64, Float64}`: test statistics confidence interval
- `df::Float64`: degrees of freedom
- `p::Float64`: p-value

Named tuple containing for type === `:perm`:
- `t::Tuple{perm_diff::Vector{Float64}, obs_diff::Float64}`: test results: (permutation difference, observed difference)
- `p1::Float64`: one-sided p-value
- `p2::Float64`: two-sided p-value
"""
function cmp_test(s1::AbstractVector, s2::AbstractVector; paired::Bool, alpha::Float64=0.05, type::Symbol=:auto, exact::Bool=false, nperm::Int64=1000)::Union{@NamedTuple{t::OneSampleTTest, ts::Tuple{Float64, String}, tc::Tuple{Float64, Float64}, df::Float64, p::Float64}, @NamedTuple{t::EqualVarianceTTest, ts::Tuple{Float64, String}, tc::Tuple{Float64, Float64}, df::Float64, p::Float64}, @NamedTuple{t::ApproximateMannWhitneyUTest{Float64}, ts::Tuple{Float64, String}, tc::Float64, df::Float64, p::Float64}, @NamedTuple{t::@NamedTuple{perm_diff::Vector{Float64}, obs_diff::Float64}, p1::Float64, p2::Float64}}

    _check_var(type, [:auto, :perm, :p, :np], "type")
    paired && length(s1) != length(s2) && @error "For paired test both vectors must have the same size."

    ks = ApproximateTwoSampleKSTest(s1, s2)
    pks = pvalue(ks)
    if type !== :perm
        if (pks > alpha && type === :auto) || type === :p
            if paired
                _info("Using one sample T-test")
                t = OneSampleTTest(s1, s2)
            else
                pf = pvalue(VarianceFTest(s1, s2))
                if pf > alpha
                    _info("Using equal variance two samples T-test")
                    t = EqualVarianceTTest(s1, s2)
                else
                    _info("Using unequal variance two samples T-test")
                    t = UnequalVarianceTTest(s1, s2)
                end
            end
            df = t.df
            ts = t.t
            tc = confint(t, level=(1 - alpha))
            tn = "t"
        elseif (pks <= alpha && type === :auto) || type === :np
            if paired
                if exact
                    _info("Using exact signed rank (Wilcoxon) test")
                    t = ExactSignedRankTest(s1, s2)
                else
                    _info("Using signed rank (Wilcoxon) test")
                    t = SignedRankTest(s1, s2)
                end
                ts = t.W
                df = t.n - 1
                tn = "W"
            else
                _info("Using Mann-Whitney U test")
                t = MannWhitneyUTest(s1, s2)
                ts = t.U
                df = length(s1) + length(s2) - 2
                tn = "U"
            end
            tc = NaN
        end
    else
        # combine groups
        g = [s1; s2]
        n1 = length(s1)
        n2 = length(s2)
        g_idx = [repeat([0], n1); repeat([1], n2)]
        # permute signals
        p = randperm(n1 + n2)
        g = g[p]
        g_idx = g_idx[p]
        _info("Group 1 mean: $(round(mean(s1), digits=3))")
        _info("Group 2 mean: $(round(mean(s2), digits=3))")
        perm_diff = zeros(nperm)

        # initialize progress bar
        progress_bar && (progbar = Progress(nperm, dt=1, barlen=20, color=:white))

        @inbounds for idx in 1:nperm
            f_idx = randperm(n1 + n2)
            f_idx[f_idx .< n1 + 1] .= 0
            f_idx[f_idx .> 0 ] .= 1
            perm_diff[idx] = mean(g[f_idx .== 0]) - mean(g[f_idx .== 1])

            # update progress bar
            progress_bar && next!(progbar)
        end
        observed_difference = mean(g[g_idx .== 0]) - mean(g[g_idx .== 1])
        observed_difference = round(observed_difference, digits=3)
    end

    if type !== :perm
        p = pvalue(t)
        p < eps() && (p = eps())
        return (t=t, ts=(round(ts, digits=4), tn), tc=round.(tc, digits=4), df=round(df, digits=4), p=p)
    else
        if pks < alpha
            _info("H0 has non-normal distribution; p-value for KS-test: $pks")
            p_one_tailed = sum(perm_diff .> observed_difference) / nperm
        else
            z = (observed_difference - mean(perm_diff)) / std(perm_diff)
            z = round(z, digits=3)
            p_one_tailed = 1 - cdf(Distributions.Normal(), abs(z))
        end
        return (t=(perm_diff=perm_diff, obs_diff=observed_difference), p1=p_one_tailed, p2=(p_one_tailed / 2))
    end

end
