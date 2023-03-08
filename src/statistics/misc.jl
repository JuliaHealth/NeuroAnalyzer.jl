export z_score
export k_categories
export se
export rng
export moe
export cvar
export effsize
export binom_prob
export binom_stat
export cvar_mean
export cvar_median

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
