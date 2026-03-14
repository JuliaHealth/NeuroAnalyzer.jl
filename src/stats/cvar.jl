export cvm
export cvmd
export fano

"""
    cvm(x)

Calculate the coefficient of variation (CV) for the mean.

Computed as `σ / μ`, where `σ = std(x)` and `μ = mean(x)`. Expresses the standard deviation as a fraction of the mean; meaningful only when `mean(x) ≠ 0`.

# Arguments

- `x::AbstractVector`: input vector; must contain at least 2 elements and have a non-zero mean

# Returns

- `Float64`: coefficient of variation (dimensionless ratio)

# Throws

- `ArgumentError`: if `length(x) < 2` or `mean(x) == 0`

# See also

[`cvmd`](@ref), [`fano`](@ref)
"""
function cvm(x::AbstractVector)::Float64

    !(length(x) >= 2) && throw(ArgumentError("x must contain at least 2 elements."))
    m = mean(x)
    !(m != 0) && throw(ArgumentError("mean(x) must not be zero (division by zero)."))

    return std(x) / m

end

"""
    cvmd(x)

Calculate the coefficient of variation (CV) for the median.

Uses the robust formula `(Q3 − Q1) / 2 / median(x)`, where `Q1` and `Q3` are the 25th and 75th percentiles (the semi-interquartile range normalised by the median). Meaningful only when `median(x) ≠ 0`.

# Arguments

- `x::AbstractVector`: input vector; must contain at least 2 elements and have a non-zero median

# Returns

- `Float64`: robust coefficient of variation (dimensionless ratio)

# Throws

- `ArgumentError`: if `length(x) < 2` or `median(x) == 0`

# See also

[`cvm`](@ref), [`fano`](@ref)
"""
function cvmd(x::AbstractVector)::Float64

    !(length(x) >= 2) && throw(ArgumentError("x must contain at least 2 elements."))
    md = median(x)
    !(md != 0) && throw(ArgumentError("median(x) must not be zero (division by zero)."))

    return (quantile(x, 0.75) - quantile(x, 0.25)) / 2 / md

end

"""
    fano(x)

Calculate the Fano factor.

Computed as `σ² / μ`, where `σ² = var(x)` and `μ = mean(x)`. The Fano factor is a measure of the dispersion of a probability distribution relative to its mean; a value of 1 corresponds to a Poisson process.

# Arguments

- `x::AbstractVector`: input vector; must contain at least 2 elements and have a non-zero mean

# Returns

- `Float64`: Fano factor (units of `x`)

# Throws

- `ArgumentError`: if `length(x) < 2` or `mean(x) == 0`

# See also

[`cvm`](@ref), [`cvmd`](@ref)
"""
function fano(x::AbstractVector)::Float64

    !(length(x) >= 2) && throw(ArgumentError("x must contain at least 2 elements."))
    m = mean(x)
    !(m != 0) && throw(ArgumentError("mean(x) must not be zero (division by zero)."))

    return var(x) / m

end
