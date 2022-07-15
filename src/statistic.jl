"""
    hildebrand_rule(x)

Calculate Hildebrand rule for vector `x`.
If H < 0.2 then the vector `x` is symmetrical.

# Arguments

- `x::Vector{<:Real}`

# Returns

- `h::Float64`
"""
function hildebrand_rule(x::Vector{<:Real})

    return (mean(x) - median(x)) ./ std(x)
end

"""
    jaccard_similarity(x, y)

Calculate Jaccard similarity between two vectors `x` and `y`.

# Arguments

- `n::Int64`

# Returns

- `j::Float64`
"""
function jaccard_similarity(x::Vector{<:Real}, y::Vector{<:Real})

    i = length(intersect(x, y))
    u = length(x) + length(y) - i
    j = i / u

    return j
end

"""
    z_score(x)

Calculate Z-scores for each value of the vector `x`.

# Arguments

- `x::Vector{<:Real}`

# Returns

- `z_score::Vector{Float64}`
"""
function z_score(x::Vector{<:Real})

    return (x .- mean(x)) ./ std(x)
end

"""
    k_categories(n)

Calculate number of categories for a given sample size `n`.

# Arguments

- `n::Int64`

# Returns

- `k::Float64`
"""
function k_categories(n::Int64)

    return (sqrt(n), (1 + 3.222 * log10(n)))
end

"""
    effsize(x1, x2)

Calculate Cohen's d and Hedges g effect sizes.

# Arguments

- `x1::Vector{Float64}`
- `x2::Vector{Float64}`

# Returns

Named tuple containing:
- `d::Float64`: Cohen's d
- `g::Float64`: Hedges g
"""
function effsize(x1::Vector{<:Real}, x2::Vector{<:Real})
    d = (mean(x2) - mean(x1)) / sqrt((std(x1)^2 + std(x2)^2) / 2)
    g = (mean(x2) - mean(x1)) / sqrt((((length(x1) - 1) * (std(x1)^2)) + ((length(x2) - 1) * (std(x2)^2))) / (length(x1) + length(x2) - 2))
    return (cohen=d, hedges=g)
end

"""
    infcrit(m)

Calculate Akaike’s Information Criterion (AIC) and Bayesian Information Criterion (BIC) for a linear regression `model`.

# Arguments

- `m::StatsModels.TableRegressionModel`

# Returns

Named tuple containing:
- `aic::Float64`
- `bic::Float64`
"""
function infcrit(m)

    typeof(m) <: StatsModels.TableRegressionModel || throw(ArgumentError("Argument must be a regression model."))

    k = length(coef(m)) - 1
    n = length(MultivariateStats.predict(m))
    aic = 2 * k - 2 * log(r2(m))
    bic = k * log(n) - 2 * log(r2(m))
    
    return (aic=aic, bic=bic)
end

"""
    grubbsx; alpha, t)

Perform Grubbs test for outlier in vector `x`.

# Arguments

- `x::Vector{<:Real}`
- `alpha::Float64=0.95`
- `t::Int64=0`: test type: -1 test whether the minimum value is an outlier; 0 two-sided test; 1 test whether the maximum value is an outlier

# Returns

Named tuple containing:
- `g::Bool`: true: outlier exists, false: there is no outlier
"""
function grubbs(x::Vector{<:Real}; alpha::Float64=0.95, t::Int64=0)

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

- `x::Vector{<:Real}`
- `method::Symbol=iqr`: methods: `:iqr` (interquartile range), `:z` (z-score) or `:g` (Grubbs test)

# Returns

- `o::Vector{Bool}`: index of outliers
"""
function outlier_detect(x::Vector{<:Real}; method::Symbol=:iqr)
    method in [:iqr, :z, :g] || throw(ArgumentError("method must be :iqr, :z or :g."))

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
    seg_tcmp(seg1, seg2, paired)

Compare two segments; Kruskall-Wallis test is used first, next t-test (paired on non-paired) or non-parametric test (paired: Wilcoxon signed rank, non-paired: Mann-Whitney U test) is applied.

# Arguments

- `seg1::Array{Float64, 3}`
- `seg2::Array{Float64, 3}`
- `paired::Bool`
- `alpha::Float64=0.05`: confidence level
- `type::Symbol=:auto`: choose test automatically (:auto, :p for parametric and :np for non-parametric)

# Returns

Named tuple containing:
- `tt`: test results
- `t::Tuple{Float64, String}`: test value and name
- `c::Tuple{Float64, Float64}`: test value confidence interval
- `df::Int64`: degrees of freedom
- `p::Float64`: p-value
- `seg1::Vector{Float64}`: averaged segment 1
- `seg2::Vector{Float64}`: averaged segment 2
"""
function seg_cmp(seg1::Array{Float64, 3}, seg2::Array{Float64, 3}; paired::Bool, alpha::Float64=0.05, type::Symbol=:auto)

    type in [:auto, :p, :np] || throw(ArgumentError("type must be :auto, :p or :np."))
    paired == true && size(seg1) != size(seg2) && throw(ArgumentError("For paired test both segments must have the same size."))

    seg1_avg = reshape(mean(mean(seg1, dims=1), dims=2), size(seg1, 3))
    seg2_avg = reshape(mean(mean(seg2, dims=1), dims=2), size(seg2, 3))

    ks = ApproximateTwoSampleKSTest(seg1_avg, seg2_avg)
    pks = pvalue(ks)
    if (pks < alpha && type === :auto) || type === :p
        if paired == true
            tt = OneSampleTTest(seg1_avg, seg2_avg)
        else
            pf = pvalue(VarianceFTest(seg1_avg, seg2_avg))
            if pf < alpha
                tt = EqualVarianceTTest(seg1_avg, seg2_avg)
            else
                tt = UnequalVarianceTTest(seg1_avg, seg2_avg)
            end
        end
        df = tt.df
        t = round(tt.t, digits=2)
        c = round.(confint(tt, level=(1 - alpha)), digits=2)
        tn = "t"
    elseif (pks >= alpha && type === :auto) || type === :np
        if paired == true
            tt = SignedRankTest(seg1_avg, seg2_avg)
            t = round(tt.W, digits=2)
            df = tt.n - 1
            tn = "W"
        else
            tt = MannWhitneyUTest(seg1_avg, seg2_avg)
            t = round(tt.U, digits=2)
            df = 2 * size(seg1, 3) - 2
            tn = "U"
        end
        c = NaN
    end

    p = pvalue(tt)
    p < eps() && (p = 0.0001)
    p = round(p, digits=4)

    return (tt=tt, t=(t, tn), c=c, df=df, p=p, seg1=seg1_avg, seg2=seg2_avg)
end