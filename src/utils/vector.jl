export vsearch
export vsplit
export minat
export maxat
export vreduce

"""
    vsearch(y, x)

Return the index of the first occurrence of string `y` in vector `x`, or `nothing` if `y` is not present.

# Arguments

- `y::String`: value to search for
- `x::Vector{String}`: vector to search within

# Returns

- `Union{Int64, Nothing}`: index of the first match, or `nothing` if not found

# See also

[`vsearch(::Real, ::AbstractVector)`](@ref)
"""
function vsearch(
    y::String,
    x::Vector{String}
)::Union{
    Int64,
    Nothing
}

    # findfirst already returns nothing on no match
    return findfirst(isequal(y), x)

end

"""
    vsearch(y, x; <keyword arguments>)

Return the index of the element in `x` nearest to scalar `y`.

# Arguments

- `y::Real`: value of interest
- `x::AbstractVector`: vector to search within; must not be empty
- `acc::Bool=false`: if `true`, also return the absolute difference between `y` and `x[idx]`

# Returns

- `Int64`: index of the nearest element
- `Tuple{Int64, Real}`: `(index, |y − x[index]|)` when `acc=true`

# Throws

- `ArgumentError`: if `x` is empty

# See also

[`vsearch(::AbstractVector, ::AbstractVector)`](@ref)
"""
function vsearch(
    y::Real,
    x::AbstractVector;
    acc::Bool = false
)::Union{
    Int64,
    Tuple{Int64, Real}
}

    !(length(x) > 0) && throw(ArgumentError("x must not be empty."))
    d, idx = findmin(abs.(x .- y))

    return acc ? (idx, d) : idx

end

"""
    vsearch(y, x; <keyword arguments>)

Return the indices of the elements in `x` nearest to each element of `y`.

# Arguments

- `y::AbstractVector`: vector of interest
- `x::AbstractVector`: vector to search within; must not be empty and must be at least as long as `y`.
- `acc::Bool=false`: if `true`, also return the absolute differences between each `y[i]` and its nearest match in `x`

# Returns

- `Vector{Int64}`: indices of nearest elements (one per entry of `y`)
- `Tuple{Vector{Int64}, Vector{Real}}`: `(indices, differences)` when `acc=true`

# Throws

- `ArgumentError`: if `x` is empty or `length(y) > length(x)`

# See also

[`vsearch(::Real, ::AbstractVector)`](@ref)
"""
function vsearch(
    y::AbstractVector,
    x::AbstractVector;
    acc::Bool = false
)::Union{
    AbstractVector,
    Tuple{AbstractVector, AbstractVector}
}

    !(length(x) > 0) && throw(ArgumentError("x must not be empty."))
    !(length(y) <= length(x)) && throw(ArgumentError("length(y) must be ≤ length(x)."))

    idx = Vector{Int64}(undef, length(y))
    d = Vector{Float64}(undef, length(y))
    @inbounds for i in eachindex(y)
        d[i], idx[i] = findmin(abs.(x .- y[i]))
    end

    return acc ? (idx, d) : idx

end

"""
    vsplit(x, n)

Split a vector into contiguous pieces of equal length `n`.

# Argument

- `x::AbstractVector`: input vector; must not be empty and its length must be divisible by `n`
- `n::Int64`: length of each piece; must be ≥ 1

# Returns

- `Vector{AbstractVector}`: vector of `length(x) ÷ n` sub-vectors, each of length `n`

# Throws

- `ArgumentError`: if `x` is empty, `n < 1`, or `length(x)` is not a multiple of `n`
"""
function vsplit(x::AbstractVector, n::Int64 = 1)::Vector{AbstractVector}

    !(length(x) > 0) && throw(ArgumentError("x must not be empty."))
    !(n >= 1) && throw(ArgumentError("n must be ≥ 1."))
    !(length(x) % n == 0) && throw(ArgumentError("length(x) must be a multiple of n."))

    n_pieces = length(x) ÷ n

    # pre-allocate and fill
    return [x[(i - 1) * n + 1 : i * n] for i in 1:n_pieces]

end

"""
    minat(x, y)

Find the minimum value of `x` and return the corresponding value from `y` at that index.

# Argument

- `x::AbstractVector`: vector to find the minimum of; must not be empty
- `y::AbstractVector`: vector to read the value from; must be the same length as `x`

# Returns

- `Real`: `y[idx]` where `idx = argmin(x)`
- `Int64`: index of the minimum value in `x`

# Throws

- `ArgumentError`: if `x` or `y` is empty, or their lengths differ

# See also

[`maxat`](@ref)
"""
function minat(x::AbstractVector, y::AbstractVector)::Tuple{Real, Int64}

    !(length(x) > 0) && throw(ArgumentError("x must not be empty."))
    !(length(y) > 0) && throw(ArgumentError("y must not be empty."))
    !(length(x) == length(y)) && throw(ArgumentError("x and y must have the same length."))

    idx = vsearch(minimum(x), x)

    return y[idx], idx

end

"""
    maxat(x, y)

Find the maximum value of `x` and return the corresponding value from `y` at that index.

# Argument

- `x::AbstractVector`: vector to find the maximum of; must not be empty
- `y::AbstractVector`: vector to read the value from; must be the same length as `x`

# Returns

- `Real`: `y[idx]` where `idx = argmin(x)`
- `Int64`: index of the maximum value in `x`

# Throws

- `ArgumentError`: if `x` or `y` is empty, or their lengths differ

# See also

[`minat`](@ref)
"""
function maxat(x::AbstractVector, y::AbstractVector)::Tuple{Real, Int64}

    !(length(x) > 0) && throw(ArgumentError("x must not be empty."))
    !(length(y) > 0) && throw(ArgumentError("y must not be empty."))
    !(length(x) == length(y)) && throw(ArgumentError("x and y must have the same length."))

    idx = vsearch(maximum(x), x)

    return y[idx], idx

end

"""
    vreduce(x, f; <keyword arguments>)

Reduce two same-length vectors by resampling `f` at regular intervals of `n` and returning the corresponding values from `x`.

Useful for downsampling a frequency axis (and its associated data) when the number of frequency bins is large. The nearest existing entry in `f` is used for each target frequency (via `vsearch`); no interpolation is performed.

# Arguments

- `x::AbstractVector`: data vector (e.g. spectral power); must not be empty
- `f::AbstractVector`: index vector (e.g. frequencies); must have the same length as `x`
- `n::Float64=0.5`: step size for the reduced grid (in the same units as `f`)

# Returns

- `AbstractVector`: reduced data values
- `AbstractVector`: reduced frequency grid

# Throws

- `ArgumentError`: if `x` or `f` is empty, or their lengths differ

# See also

[`vsearch`](@ref)
"""
function vreduce(
    x::AbstractVector,
    f::AbstractVector;
    n::Float64 = 0.5
)::Tuple{AbstractVector, AbstractVector}

    !(length(x) > 0) && throw(ArgumentError("x must not be empty."))
    !(length(f) > 0) && throw(ArgumentError("f must not be empty."))
    !(length(x) == length(f)) && throw(ArgumentError("x and f must have the same length."))

    # build the reduced frequency grid from rounded min/max frequencies
    f1 = round(f[vsearch(round(f[1]), f)])
    f2 = round(f[vsearch(round(f[end]), f)])
    f_new = collect(f1:n:f2)

    x_new = [x[vsearch(freq, f)] for freq in f_new]

    return x_new, f_new

end
