export cmp_test

"""
    cmp_test(seg1, seg2; <keyword arguments>)

Compare two vectors. Jarque–Bera is used first to test normality of distribution, next variances are tested using F-test. Finally, paired or non-paired t-test or non-parametric test (Wilcoxon signed rank for paired data; Mann-Whitney U test for non-paired data) is applied. Permutation-based testing is also 

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`
- `paired::Bool`
- `alpha::Float64=0.05`: confidence level
- `type::Symbol=:auto`
    - `:auto`: choose test automatically
    - `:perm`: permutation-based
    - `:p`: parametric
    - `:np`: non-parametric
- `exact::Bool=false`: if true, use exact Wilcoxon test
- `nperm::Int64=1000`: number of permutation for `:perm` method
- `verbose::Bool=true`: print detailed output

# Returns

For permutation-based test, named tuple containing:
- `t::Tuple{perm_diff::Vector{Float64}, obs_diff::Float64}`: test result (permutation difference, observed difference)
- `p1::Float64`: one-tiled p value
- `p2::Float64`: two-tiled p value

Otherwise, named tuple containing:
- `t::Union{OneSampleTTest, EqualVarianceTTest, UnequalVarianceTTest, ExactSignedRankTest, ApproximateSignedRankTest, ExactMannWhitneyUTest, ApproximateMannWhitneyUTest}`: statistical test
- `ts::Tuple{Float64, String}`: test statistic value and name
- `tc::Tuple{Float64, Float64}`: test statistic confidence interval
- `df::Float64`: degrees of freedom
- `p::Float64`: p value
"""
function cmp_test(s1::AbstractVector, s2::AbstractVector; paired::Bool, alpha::Float64=0.05, type::Symbol=:auto, exact::Bool=false, nperm::Int64=1000, verbose::Bool=true)::Union{@NamedTuple{t::OneSampleTTest, ts::Tuple{Float64, String}, tc::Tuple{Float64, Float64}, df::Float64, p::Float64}, @NamedTuple{t::EqualVarianceTTest, ts::Tuple{Float64, String}, tc::Tuple{Float64, Float64}, df::Float64, p::Float64}, @NamedTuple{t::UnequalVarianceTTest, ts::Tuple{Float64, String}, tc::Tuple{Float64, Float64}, df::Float64, p::Float64}, @NamedTuple{t::ExactSignedRankTest{Float64}, ts::Tuple{Float64, String}, tc::Float64, df::Float64, p::Float64}, @NamedTuple{t::ApproximateSignedRankTest{Float64}, ts::Tuple{Float64, String}, tc::Float64, df::Float64, p::Float64}, @NamedTuple{t::ExactMannWhitneyUTest{Float64}, ts::Tuple{Float64, String}, tc::Float64, df::Float64, p::Float64}, @NamedTuple{t::ApproximateMannWhitneyUTest{Float64}, ts::Tuple{Float64, String}, tc::Float64, df::Float64, p::Float64}, @NamedTuple{t::@NamedTuple{perm_diff::Vector{Float64}, obs_diff::Float64}, p1::Float64, p2::Float64}}

    _check_var(type, [:auto, :perm, :p, :np], "type")
    paired && @assert length(s1) == length(s2) "For paired test both vectors must have the same size."
    @assert alpha < 1.0 "alpha must be < 1.0."
    @assert alpha > 0.0 "alpha must be > 0.0."
    @assert nperm > 0 "nperm must be > 0."

    # ks = ApproximateTwoSampleKSTest(s1, s2)
    jb = JarqueBeraTest([s1; s2])
    pjb = round(pvalue(jb), digits=3)
    jb = round(jb.JB, digits=3)
    pjb < 0.05 ? (verbose && println("Distribution: non-normal (Jarque-Bera test: JB=$jb, p=$pjb)")) : (verbose && println("Distribution: normal (Jarque-Bera test: JB=$jb, p=$pjb)"))
    if type !== :perm
        if (pjb >= alpha && type === :auto) || type === :p
            if paired
                verbose && println("Using one sample T-test (paired data)")
                t = OneSampleTTest(s1, s2)
            else
                pf = pvalue(VarianceFTest(s1, s2))
                if pf >= alpha
                    verbose && println("Using equal variance two samples T-test")
                    t = EqualVarianceTTest(s1, s2)
                else
                    verbose && println("Using unequal variance two samples T-test")
                    t = UnequalVarianceTTest(s1, s2)
                end
            end
            df = t.df
            ts = t.t
            tc = confint(t, level=(1 - alpha))
            tn = "t"
        elseif (pjb < alpha && type === :auto) || type === :np
            if paired
                if exact
                    verbose && println("Using exact signed rank (Wilcoxon) test (paired data)")
                    t = ExactSignedRankTest(s1, s2)
                else
                    verbose && println("Using signed rank (Wilcoxon) test (paired data)")
                    t = SignedRankTest(s1, s2)
                end
                ts = t.W
                df = t.n - 1
                tn = "W"
            else
                verbose && println("Using Mann-Whitney U test")
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
        verbose && println("Group 1 mean: $(round(mean(s1), digits=3))")
        verbose && println("Group 2 mean: $(round(mean(s2), digits=3))")
        perm_diff = zeros(nperm)

        # initialize progress bar
        progress_bar && (progbar = Progress(nperm, dt=1, barlen=20, color=:white))

        @inbounds for idx in 1:nperm
            f_idx = randperm(n1 + n2)
            f_idx[f_idx .< n1 + 1] .= 0
            f_idx[f_idx .> 0 ] .= 1
            perm_diff[idx] = mean(g[f_idx .== 0]) - mean(g[f_idx .== 1])

            # update progress bar
            next!(progbar)
        end
        observed_difference = mean(g[g_idx .== 0]) - mean(g[g_idx .== 1])
        observed_difference = round(observed_difference, digits=3)
    end

    if type !== :perm
        p = pvalue(t)
        p < eps() && (p = eps())
        return (t=t, ts=(round(ts, digits=3), tn), tc=round.(tc, digits=3), df=round(df, digits=3), p=p)
    else
        if pjb < alpha
            verbose && println("H0 has non-normal distribution; p value for JB-test: $pjb")
            p_one_tailed = sum(perm_diff .> observed_difference) / nperm
        else
            z = (observed_difference - mean(perm_diff)) / std(perm_diff)
            z = round(z, digits=3)
            p_one_tailed = 1 - cdf(Distributions.Normal(), abs(z))
        end
        return (t=(perm_diff=perm_diff, obs_diff=observed_difference), p1=p_one_tailed, p2=(p_one_tailed / 2))
    end

end
