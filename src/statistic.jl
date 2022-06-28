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
    grubbs(x; alpha, t)

Perform Grubbs test for outlier in vector `x`.

# Arguments

- `x::Vector{<:Real}`
- `alpha::Float64=0.95`
- `t::Int64=0`: test type: -1 test whether the minimum value is an outlier; 0 two-sided test; 1 test whether the maximum value is an outlier

# Returns

- `g::Bool`: true: outlier exists, false: there is no outlier
"""
function grubbs(x::Vector{<:Real}; alpha::Float64=0.95, t::Int64=0)

    n = length(x)
    d = n - 2

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

    if two_sided == true
        h = (n - 1) / sqrt(n) * sqrt(quantile(TDist(d), 1 - (alpha / (2 * n)))^2 / (n - 2 + quantile(TDist(d), 1 - (alpha / (2 * n)))^2))
    else
        h = (n - 1) / sqrt(n) * sqrt(quantile(TDist(d), 1 - (alpha / n))^2 / (n - 2 + quantile(TDist(d), 1 - (alpha / n))^2))
    end

    g < h && return false
    g > h && return true

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

