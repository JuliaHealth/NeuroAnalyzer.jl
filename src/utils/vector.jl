export vsearch
export vsplit

"""
    vsearch(y, x; acc)

Return the positions of the `y` value in the vector `x`.

# Arguments

- `y::T`: value of interest
- `x::AbstractVector`: vector to search within
- `acc::Bool=false`: if true, return the difference between `y` and `x[idx]`

# Returns

- `idx::Int64`
- `d::Real`: the difference between `y` and `x[idx]`
"""
function vsearch(y::T, x::AbstractVector; acc::Bool=false) where {T<:Real}

    d, idx = findmin(abs.(x .- y))

    return acc == true ? (idx, d) : idx

end

"""
    vsearch(y, x; acc)

Return the positions of the `y` vector in the vector `x`.

# Arguments

- `y::AbstractVector`: vector of interest
- `x::AbstractVector`: vector to search within
- `acc::Bool=false`: if true, return the difference between `y` and `x[idx:idx + length(y)]`

# Returns

- `idx::Int64`
- `d::Real`: the difference between `y` and `x[idx:idx + length(y)]`
"""
function vsearch(y::AbstractVector, x::AbstractVector; acc::Bool=false)

    length(y) > length(x) && throw(ArgumentError("Length of 'y' cannot be larger than length 'x'"))

    idx = zeros(length(y))
    d = zeros(length(y))

    @inbounds @simd for y_idx in 1:length(y)
        d[y_idx], idx[y_idx] = findmin(abs.(x .- y[y_idx]))
    end

    return acc == true ? (convert.(Int64, idx), d) : idx
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