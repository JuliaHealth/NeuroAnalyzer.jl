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
    grubbs(signal)

Perform Grubbs test for outlier in `signal`.

# Arguments

- `signal::Vector{<:Real}`

# Returns

- `g::Bool`: true: outlier exists, false: there is no outlier
"""
function grubbs(signal::Vector{<:Real})
    n = length(signal)
    d = n - 2
    alpha = 0.95
    g = maximum(abs.(signal .- mean(signal))) / std(signal)
    h = (n - 1) / sqrt(n) * sqrt(quantile(TDist(d), 1 - (alpha  / (2 * n)))^2 / (n - 2 + quantile(TDist(d), 1 - (alpha / (2 * n)))^2))
    if g > h
        g = true
    else
        g = false
    end

    return g
end