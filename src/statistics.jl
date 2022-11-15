"""
    hildebrand_rule(x)

Calculate Hildebrand rule for vector `x`.
If H < 0.2 then the vector `x` is symmetrical.

# Arguments

- `x::AbstractVector`

# Returns

- `h::Float64`
"""
function hildebrand_rule(x::AbstractVector)
    return (mean(x) - median(x)) ./ std(x)
end

"""
    jaccard_similarity(x, y)

Calculate Jaccard similarity between two vectors `x` and `y`.

# Arguments

- `x::AbstractVector`
- `y::AbstractVector`

# Returns

- `j::Float64`
"""
function jaccard_similarity(x::AbstractVector, y::AbstractVector)

    i = length(intersect(x, y))
    u = length(x) + length(y) - i
    j = i / u

    return j
end

"""
    z_score(x)

Calculate Z-scores for each value of the vector `x`.

# Arguments

- `x::AbstractVector`

# Returns

- `z_score::Vector{Float64}`
"""
function z_score(x::AbstractVector)

    return (x .- mean(x)) ./ std(x)
end

"""
    k_categories(n)

Calculate number of categories for a given sample size `n`.

# Arguments

- `n::Int64`

# Returns

Named tuple containing:
- `k1::Float64`: sqrt(n)
- `k2::Float64`: 1 + 3.222 * log10(n)
"""
function k_categories(n::Int64)
    return (k1=sqrt(n), k2=(1 + 3.222 * log10(n)))
end

"""
    effsize(x1, x2)

Calculate Cohen's d and Hedges g effect sizes.

# Arguments

- `x1::AbstractVector`
- `x2::AbstractVector`

# Returns

Named tuple containing:
- `d::Float64`: Cohen's d
- `g::Float64`: Hedges g
"""
function effsize(x1::AbstractVector, x2::AbstractVector)
    d = (mean(x2) - mean(x1)) / sqrt((std(x1)^2 + std(x2)^2) / 2)
    g = (mean(x2) - mean(x1)) / sqrt((((length(x1) - 1) * (std(x1)^2)) + ((length(x2) - 1) * (std(x2)^2))) / (length(x1) + length(x2) - 2))
    return (cohen=d, hedges=g)
end

"""
    infcrit(m)

Calculate Akaike’s Information Criterion (AIC) and Bayesian Information Criterion (BIC) for a linear regression model `m`.

# Arguments

- `m::StatsModels.TableRegressionModel`

# Returns

Named tuple containing:
- `aic::Float64`
- `bic::Float64`
"""
function infcrit(m::T) where {T<:StatsModels.TableRegressionModel}

    k = length(coef(m)) - 1
    n = length(MultivariateStats.predict(m))
    aic = 2 * k - 2 * log(r2(m))
    bic = k * log(n) - 2 * log(r2(m))
    
    return (aic=aic, bic=bic)
end

"""
    grubbs(x; alpha, t)

Perform Grubbs test for outlier in vector `x`.

# Arguments

- `x::AbstractVector`
- `alpha::Float64=0.95`
- `t::Int64=0`: test type: -1 test whether the minimum value is an outlier; 0 two-sided test; 1 test whether the maximum value is an outlier

# Returns

Named tuple containing:
- `g::Bool`: true: outlier exists, false: there is no outlier
"""
function grubbs(x::AbstractVector; alpha::Float64=0.95, t::Int64=0)

    std(x) == 0 && throw(ArgumentError("Standard deviation of the input vector must not be 0."))

    n = length(x)
    df = n - 2

    if t == 0
        two_sided = true
        g = maximum(abs.(x .- mean(x))) / std(x)
    elseif t == -1
        two_sided = false
        g = (mean(x) - minimum(x)) / std(x)
    elseif t == 1
        two_sided = false
        g = (maximum(x) - mean(x)) / std(x)
    else
        throw(ArgumentError("type must be -1, 0 or 1."))
    end

    p = two_sided == true ? (1 - alpha) / (2 * n) : (1 - alpha) / n
    t_critical = quantile(TDist(df), 1 - p)
    h = (n - 1) * t_critical / sqrt(n * (df + t_critical^2))

    return g < h ? false : true
end

"""
    outlier_detect(x; method)

Detect outliers in `x`.

# Arguments

- `x::AbstractVector`
- `method::Symbol=iqr`: methods: `:iqr` (interquartile range), `:z` (z-score) or `:g` (Grubbs test)

# Returns

- `o::Vector{Bool}`: index of outliers
"""
function outlier_detect(x::AbstractVector; method::Symbol=:iqr)
    _check_var(method, [:iqr, :z, :g], "method")
    o = zeros(Bool, length(x))
    
    if method === :iqr
        m1 = quantile(x, 0.25) - 1.5 * iqr(x)
        m2 = quantile(x, 0.75) + 1.5 * iqr(x)
        o[x .< m1] .= true
        o[x .> m2] .= true
    elseif method === :z
        z = z_score(x)
        o[z .< -3] .= true
        o[z .> 3] .= true
    else
        length(x) > 6 || throw(ArgumentError("For :g method length(x) must be > 6."))
        x_tmp = deepcopy(x)
        for idx in length(x_tmp):-1:6
            _, m_idx = findmax(x_tmp)
            if grubbs(x_tmp, t=1) == true
                o[m_idx] = true
                deleteat!(x_tmp, m_idx)
            end
        end
        x_tmp = deepcopy(x)
        for idx in length(x_tmp):-1:6
            _, m_idx = findmin(x_tmp)
            if grubbs(x_tmp, t=-1) == true
                o[m_idx] = true
                deleteat!(x_tmp, m_idx)
            end
        end
    end
    
    return o
end

"""
    seg_mean(seg)

Calculate mean of a segment (e.g. spectrogram).

# Arguments

- `seg::AbstractArray`

# Returns

- `seg::Vector{Float64}`: averaged segment
"""
function seg_mean(seg::AbstractArray)

    return reshape(mean(mean(seg, dims=1), dims=2), size(seg, 3))
end

"""
    seg2_mean(seg1, seg2)

Calculate mean of two segments (e.g. spectrograms).

# Arguments

- `seg1::AbstractArray`
- `seg2::AbstractArray`

# Returns

Named tuple containing:
- `seg1::Vector{Float64}`: averaged segment 1
- `seg2::Vector{Float64}`: averaged segment 2
"""
function seg2_mean(seg1::AbstractArray, seg2::AbstractArray)

    seg1_avg = reshape(mean(mean(seg1, dims=1), dims=2), size(seg1, 3))
    seg2_avg = reshape(mean(mean(seg2, dims=1), dims=2), size(seg2, 3))

    return (seg1=seg1_avg, seg2=seg2_avg)
end

"""
    binom_prob(p, r, n)

Calculate probability of exactly `r` successes in `n` trials.

# Arguments

- `p::Float64`: proportion of successes
- `r::Int64`: number of successes
- `n::Int64`: number of trials

# Returns

- `binomp::Float64`: probability
"""
function binom_prob(p::Float64, r::Int64, n::Int64)
    return binomial(n, r) * (p^r) * (1 - p)^(n - r)
end

"""
    binom_stat(p, n)

Calculate mean and standard deviation for probability `p`.

# Arguments

- `p::Float64`: proportion of successes
- `n::Int64`: number of trials

# Returns

- `mean::Float64`
- `std::Float64`
"""
function binom_stat(p::Float64, n::Int64)
    return n * p, sqrt(n * p * (1 - p))
end

"""
    cvar_mean(x)

Calculate coefficient of variation for a mean.

# Arguments

- `x::AbstractVector`

# Returns

- `cvar::Float64`
"""
function cvar_mean(x::AbstractVector)
    return std(x) / mean(x)
end

"""
    cvar_median(x)

Calculate coefficient of variation for a median.

# Arguments

- `x::AbstractVector`

# Returns

- `cvar::Float64`
"""
function cvar_median(x::AbstractVector)
    return ((quantile(x, 0.75) - quantile(x, 0.25)) / 2) / median(x)
end

"""
    cvar(se, s)

Calculate coefficient of variation for statistic `s`.

# Arguments

- `se::Real`: standard error
- `s::Real`: statistics, e.g. mean value

# Returns

- `cvar::Float64`
"""
function cvar(se::Real, s::Real)
    return 100 * (se / s)
end

"""
    effsize(p1, p2)

Calculate effect size for two proportions `p1` and `p2`.

# Arguments

- `p1::Float64`: 1st proportion, e.g. 0.7
- `p2::Float64`: 2nd proportion, e.g. 0.3

# Returns

- `e::Float64`
"""
function effsize(p1::Float64, p2::Float64)
    p1 + p2 == 1.0 || throw(ArgumentError("Proportions must add to 1.0."))
    return 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))
end

"""
    meang(x)

Calculate geometric mean.

# Arguments

- `x::AbstractVector`

# Returns

- `m::Float64`
"""
function meang(x::AbstractVector)
    return exp(mean(log.(x[x .> 0])))
end

"""
    meanh(x)

Calculate harmonic mean.

# Arguments

- `x::AbstractVector`

# Returns

- `m::Float64`
"""
function meanh(x::AbstractVector)
    return length(x) / sum(1 ./ x)
end

"""
    meanw(x, w)

Calculate weighted mean.

# Arguments

- `x::AbstractVector`
- `w::AbstractVector`: weights

# Returns

- `m::Float64`
"""
function meanw(x::AbstractVector, w::AbstractVector)
    length(x) == length(w) || throw(ArgumentError("Weights and values vectors must have the same length."))
    return length(x) / sum(1 ./ x)
end

"""
    moe(n)

Calculate margin of error for given sample size `n`.

# Arguments

- `n::Int64`

# Returns

- `moe::Float64`
"""
function moe(n::Int64)
    return 1 / sqrt(n)
end

"""
    rng(x)

Calculate range.

# Arguments

- `x::AbstractVector`

# Returns

- `r::Float64`
"""
function rng(x::AbstractVector)
    return maximum(x) - minimum(x)
end

"""
    se(x)

Calculate standard error.

# Arguments

- `x::AbstractVector`

# Returns

- `se::Float64`
"""
function se(x::AbstractVector)
    return std(x) / sqrt(length(x))
end

"""
    pred_int(n)

Calculates the prediction interval (95% CI adjusted for sample size)

# Arguments

- `n::Int64`: sample size

# Returns

- `pred_int::Float64`
"""
function pred_int(n::Int64)
    n < 1 && throw(ArgumentError("n must be ≥ 1."))
    (n > 0 && n < 21) && return [NaN, 15.56, 4.97, 3.56, 3.04, 2.78, 2.62, 2.51, 2.43, 2.37, 2.33, 2.29, 2.26, 2.24, 2.22, 2.18, 2.17, 2.16, 2.10][n]
    @warn "Result may not be accurate."
    n > 20 && n <= 25 && return 2.10
    n > 25 && n <= 30 && return 2.08
    n > 31 && n <= 35 && return 2.06
    n > 35 && n <= 40 && return 2.05
    n > 41 && n <= 50 && return 2.03
    n > 51 && n <= 60 && return 2.02
    n > 61 && n <= 70 && return 2.01
    n > 71 && n <= 80 && return 2.00
    n > 81 && n <= 90 && return 2.00
    n > 91 && n <= 100 && return 1.99
    n > 100 && return 1.98
end

"""
    sem_diff(x::AbstractVector, y::AbstractVector)

Calculate SEM (standard error of the mean) for the difference of two means.

# Arguments

- `x::AbstractVector`
- `y::AbstractVector`

# Returns

- `sem::Float64`
"""
function sem_diff(x::AbstractVector, y::AbstractVector)
    length(x) == length(y) || throw(ArgumentError("Both vectors must have the same length."))
    return sqrt((std(x)^2 / sqrt(length(x))) + (std(y)^2 / sqrt(length(y))))
end

"""
    prank(x)

Calculate percentile rank.

# Arguments

- `x::AbstractVector`: the vector to analyze

# Returns

- `prnk::Vector{Float64}`
"""
function prank(x::AbstractVector)
    xorder = sortperm(x)
    x = sort(x)
    prnk = zeros(length(x))
    for idx in eachindex(x)
        percentile = length(x[x .< x[idx]]) / length(x) * 100
        prnk[idx] = percentile / (100 * (length(x) + 1))
    end
    return prnk[xorder]
end

"""
    linreg(x, y)

Linear regression between `x` and `y`.

# Arguments

- `x::AbstractVector`
- `y::AbstractVector`

# Notes

To predict, use: `new_x = DataFrame(x = [3.5, 7]); predict(lr, new_x)

# Returns

Named tuple containing:
- `lr::StatsModels.TableRegressionModel`: model
- `radj::Flpoat64`: R^2
- `c::Vector{Float64}`: coefficients
- `se::Vector{Float64}`: standard error for coefficients
- `aic::Float64`:: Akaike’s Information Criterion (AIC)
- `bic::Float64`:: Bayesian Information Criterion (BIC)
- `lf::Vector{Float64}`: linear fit (plot(x, lf))
"""
function linreg(x::AbstractVector, y::AbstractVector)

    df = DataFrame(:x => x, :y => y)
    lr = lm(@formula(y ~ x), df)
    radj = r2(lr)
    c = coef(lr)
    se = stderror(lr)
    aic, bic = infcrit(lr)
    lf = MultivariateStats.predict(lr)

    return (lr=lr, radj=radj, c=c, se=se, aic=aic, bic=bic, lf=lf)
end

"""
    s2_cmp(seg1, seg2, paired, alpha, type, exact)

Compare two vectors; Kruskall-Wallis test is used first, next t-test (paired on non-paired) or non-parametric test (paired: Wilcoxon signed rank, non-paired: Mann-Whitney U test) is applied.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`
- `paired::Bool`
- `alpha::Float64=0.05`: confidence level
- `type::Symbol=:auto`: choose test automatically (:auto), parametric (:p) or non-parametric (:np)
- `exact::Bool=false`: if true, use exact Wilcoxon test

# Returns

Named tuple containing:
- `t`: test results
- `ts::Tuple{Float64, String}`: test statistics
- `tc::Tuple{Float64, Float64}`: test statistics confidence interval
- `df::Int64`: degrees of freedom
- `p::Float64`: p-value
"""
function s2_cmp(s1::AbstractVector, s2::AbstractVector; paired::Bool, alpha::Float64=0.05, type::Symbol=:auto, exact::Bool=false)

    _check_var(type, [:auto, :p, :np], "type")
    paired == true && size(s1) != size(s2) && throw(ArgumentError("For paired test both segments must have the same size."))

    ks = ApproximateTwoSampleKSTest(s1, s2)
    pks = pvalue(ks)
    if (pks < alpha && type === :auto) || type === :p
        if paired == true
            verbose == true && @info "Using one sample T-test."
            t = OneSampleTTest(s1, s2)
        else
            pf = pvalue(VarianceFTest(s1, s2))
            if pf < alpha
                verbose == true && @info "Using equal variance two samples T-test."
                t = EqualVarianceTTest(s1, s2)
            else
                verbose == true && @info "Using unequal variance two samples T-test."
                t = UnequalVarianceTTest(s1, s2)
            end
        end
        df = t.df
        ts = t.t
        tc = confint(t, level=(1 - alpha))
        tn = "t"
    elseif (pks >= alpha && type === :auto) || type === :np
        if paired == true
            if exact == false
                verbose == true && @info "Using signed rank (Wilcoxon) test."
                t = SignedRankTest(s1, s2)
            else
                verbose == true && @info "Using exact signed rank (Wilcoxon) test."
                t = ExactSignedRankTest(s1, s2)
            end
            ts = t.W
            df = t.n - 1
            tn = "W"
        else
            verbose == true && @info "Using Mann-Whitney U test."
            t = MannWhitneyUTest(s1, s2)
            ts = t.U
            df = length(s1) + length(s2) - 2
            tn = "U"
        end
        tc = NaN
    end

    p = pvalue(t)
    p < eps() && (p = eps())

    return (t=t, ts=(ts, tn), tc=tc, df=df, p=p)
end

"""
    s2_cor(seg1, seg2)

Calculate correlation between two vectors.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`

# Returns

Named tuple containing:
- `t::CorrelationTest{Float64}`
- `r::Float64`: correlation coefficient
- `rc::Tuple{Float64, Float64}`: correlation coefficient confidence interval
- `tt::Tuple{Float64, String}`: t-statistics
- `df::Int64`: degrees of freedom
- `p::Float64`: p-value
"""
function s2_cor(s1::AbstractVector, s2::AbstractVector)

    length(s1) == length(s2) || throw(ArgumentError("Both vectors must have the same length."))
    t = CorrelationTest(s1, s2)
    p = pvalue(t)
    p < eps() && (p = eps())
    df = length(s1) + length(s2) - 2

    return (t=t, r=t.r, rc=confint(t), ts=(t.t, "t"), df=df, p=p)
end
