export vsearch
export vsplit

"""
    vsearch(y, x; return_distance)

Return the positions of the `y` value in the vector `x`.

# Arguments

- `y::Real`: value of interest
- `x::AbstractVector`: vector to search within
- `acc::Bool=false`: if true, return the difference between `y` and `x[y_idx]`

# Returns

- `y_idx::Int64`
- `y_dist::Real`: the difference between `y` and `x[y_idx]`
"""
function vsearch(y::Real, x::AbstractVector; acc::Bool=false)
    y_dist, y_idx = findmin(abs.(x .- y))
    return acc == true ? (y_idx, y_dist) : y_idx
end

"""
    vsearch(y, x; acc)

Return the positions of the `y` vector in the vector `x`.

# Arguments

- `y::AbstractVector`: vector of interest
- `x::AbstractVector`: vector to search within
- `acc::Bool=false`: if true, return the difference between `y` and `x[y_idx:y_idx + length(y)]`

# Returns

- `y_idx::Int64`
- `y_dist::Real`: the difference between `y` and `x[y_idx:y_idx + length(y)]`
"""
function vsearch(y::AbstractVector, x::AbstractVector; acc::Bool=false)

    length(y) > length(x) && throw(ArgumentError("Length of 'y' cannot be larger than length 'x'"))

    y_idx = zeros(length(y))
    y_dist = zeros(length(y))

    @inbounds @simd for idx in eachindex(y)
        y_dist[idx], y_idx[idx] = findmin(abs.(x .- y[idx]))
    end

    return acc == true ? (convert.(Int64, y_idx), y_dist) : y_idx
end

"""
    vsplit(x, n)

Splits the vector `x` into `n`-long pieces.

# Argument

- `x::AbstractVector`
- `n::Int64`

# Returns

- `x::Vector{AbstractVector}`
"""
function vsplit(x::AbstractVector, n::Int64=1)

    n < 1 && throw(ArgumentError("n must be ≥ 1."))
    length(x) % n == 0 || throw(ArgumentError("Length of x must be a multiple of n."))

    x_m = reshape(x, length(x) ÷ n, n)
    result = [x_m[1, :]]
    @inbounds @simd for idx in 2:size(x_m, 1)
        result = vcat(result, [x_m[idx, :]])
    end

    return result
end
