export zscore

"""
    zscore(x)

Calculate Z-scores for each element of a vector.

Computed as `(xᵢ − mean(x)) / std(x)`.

# Arguments

- `x::AbstractVector`: input vector; must contain at least 2 elements and have non-zero standard deviation

# Returns

- `Vector{Float64}`: Z-scores with mean ≈ 0 and SD ≈ 1

# Throws

- `ArgumentError`: if `length(x) < 2` or `std(x) == 0`

# See also

[`zscore(::Real, ::Real, ::Real)`](@ref)
"""
function zscore(x::AbstractVector)::Vector{Float64}

    @assert length(x) >= 2 "x must contain at least 2 elements."
    m = mean(x)
    s = std(x)
    @assert s != 0 "std(x) must not be zero."

    return (x .- m) ./ s

end

"""
    zscore(x, m, sd)

Calculate the Z-score of a single value given a known mean and standard deviation.

Computed as `(x − m) / sd`.

# Arguments

- `x::Real`: observed value
- `m::Real`: population or sample mean
- `sd::Real`: population or sample standard deviation; must not be zero

# Returns

- `Float64`: Z-score

# Throws

- `ArgumentError`: if `sd == 0`

# See also

[`zscore(::AbstractVector)`](@ref)
"""
function zscore(x::Real, m::Real, sd::Real)::Float64

    @assert sd != 0 "sd must not be zero."

    return (x - m) / sd

end
