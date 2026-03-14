export meanp
export meanc
export meang
export meanh
export meanw
export meancirc
export meant

"""
    meanp(p, n)

Calculate the expected count (mean) of a proportion.

Computed as `n × p`.

# Arguments

- `p::Float64`: proportion; must be in `[0, 1]`
- `n::Int64`: number of observations; must be ≥ 1

# Returns

- `Float64`: expected count `n × p`

# Throws

- `ArgumentError`: if `p ∉ [0, 1]` or `n < 1`

# See also

[`meanc`](@ref)
"""
function meanp(p::Float64, n::Int64)::Float64

    _in(p, (0.0, 1.0), "p")
    !(n >= 1) && throw(ArgumentError("n must be ≥ 1."))

    return n * p

end

"""
    meanc(g, x)

Calculate the weighted mean of categorical data.

Computed as `Σ(g × x) / Σx`.

# Arguments

- `g::Vector{Int64}`: group labels (e.g. `[0, 1, 2, 3]`)
- `x::Vector{Int64}`: subject counts per group (e.g. `[2, 8, 27, 45]`); must be the same length as `g`, non-empty, and sum to ≥ 1

# Returns

- `Float64`: weighted categorical mean

# Throws

- `ArgumentError`: if `g` is empty, lengths differ, or `sum(x) == 0`

# See also

[`meanp`](@ref)
"""
function meanc(g::Vector{Int64}, x::Vector{Int64})::Float64

    !(length(g) > 0) && throw(ArgumentError("g must not be empty."))
    !(length(g) == length(x)) && throw(ArgumentError("g and x must have the same length."))
    !(sum(x) > 0) && throw(ArgumentError("sum(x) must be > 0 (division by zero)."))

    return sum(g .* x) / sum(x)

end

"""
    meang(x)

Calculate the geometric mean.

Computed as `exp(mean(log.(x)))`, which is numerically stable for large vectors (avoids overflow from `prod(x)`). All elements must be strictly positive.

# Arguments

- `x::AbstractVector`: input vector; must be non-empty and contain only positive values

# Returns

- `Float64`: geometric mean

# Throws

- `ArgumentError`: if `x` is empty or contains any non-positive value

# See also

[`meanh`](@ref), [`meanw`](@ref)
"""
function meang(x::AbstractVector)::Float64

    !(length(x) > 0) && throw(ArgumentError("x must not be empty."))
    !(all(>(0), x)) && throw(ArgumentError("All elements of x must be > 0."))

    # use log-sum-exp form for numerical stability (avoids overflow from prod)
    return exp(mean(log.(x)))

end

"""
    meanh(x)

Calculate the harmonic mean.

Computed as `n / Σ(1/xᵢ)`. All elements must be non-zero.

# Arguments

- `x::AbstractVector`: input vector; must be non-empty and contain no zero values

# Returns

- `Float64`: harmonic mean

# Throws

- `ArgumentError`: if `x` is empty or contains any zero

# See also

[`meang`](@ref), [`meanw`](@ref)
"""
function meanh(x::AbstractVector)::Float64

    !(length(x) > 0) && throw(ArgumentError("x must not be empty."))
    !(!any(iszero, x)) && throw(ArgumentError("x must not contain zeros."))

    return length(x) / sum(1 ./ x)

end

"""
    meanw(x, w)

Calculate the weighted mean.

Computed as `Σ(xᵢ × wᵢ) / Σwᵢ`.

# Arguments

- `x::AbstractVector`: values vector; must be non-empty
- `w::AbstractVector`: weights vector; must have the same length as `x` and sum to a non-zero value

# Returns

- `Float64`: weighted mean

# Throws

- `ArgumentError`: if `x` is empty, lengths differ, or `sum(w) == 0`

# See also

[`meang`](@ref), [`meanh`](@ref)
"""
function meanw(x::AbstractVector, w::AbstractVector)::Float64

    !(length(x) > 0) && throw(ArgumentError("x must not be empty."))
    !(length(x) == length(w)) && throw(ArgumentError("x and w must have the same length."))
    !(sum(w) != 0) && throw(ArgumentError("sum(w) must not be zero (division by zero)."))

    return sum(x .* w) / sum(w)

end

"""
    meancirc(x; <keyword arguments>)

Calculate the circular (angular) mean.

Uses the two-argument arctangent of the mean sine and cosine components.

# Arguments

- `x::AbstractVector`: angles; must be non-empty
- `rad::Bool=false`: if `true`, `x` is in radians; if `false`, `x` is in degrees

# Returns

- `Float64`: circular mean in the same unit as the input (radians or degrees)

# Throws

- `ArgumentError`: if `x` is empty

# See also

[`meant`](@ref)
"""
function meancirc(x::AbstractVector; rad::Bool = false)::Float64

    !(length(x) > 0) && throw(ArgumentError("x must not be empty."))

    if rad
        return atan(sum(sin.(x)), sum(cos.(x)))
    else
        return rad2deg(atan(sum(sind.(x)), sum(cosd.(x))))
    end

end

"""
    meant(x; <keyword arguments>)

Calculate the trimmed mean.

Sorts `x` and removes the bottom and top `n × 100 %` of values before computing the mean.

# Arguments

- `x::AbstractVector`: input vector; must contain enough elements so that trimming leaves at least one value
- `n::Float64=0.1`: fraction of values to trim from each tail; must be in `(0, 0.5)` to ensure a non-empty trimmed set

# Returns

- `Float64`: trimmed mean

# Throws

- `ArgumentError`: if `n ∉ (0, 0.5)` or trimming leaves no observations

# See also

[`meancirc`](@ref)
"""
function meant(x::AbstractVector; n::Float64 = 0.1)::Float64

    !(n > 0.0 && n < 0.5) && throw(ArgumentError("n must be in (0, 0.5)."))
    xs  = sort(x)
    xn  = round(Int64, length(xs) * n)
    !(xn + 1 <= length(xs) - xn) && throw(ArgumentError("n is too large: no observations remain after trimming."))

    return mean(xs[(xn + 1):(end - xn)])

end
