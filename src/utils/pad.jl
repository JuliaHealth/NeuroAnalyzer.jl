export pad0
export pad2

"""
    pad0(x, n)

Pad vector / rows of matrix / array with zeros. Works with 1-, 2- and 3-dimensional arrays.

# Arguments

- `x::Union{AbstractVector, AbstractArray`
- `n::Int64`: number of zeros to add

# Returns

- `pad0::Union{AbstractVector, AbstractArray`
"""
function pad0(x::Union{AbstractVector, AbstractArray}, n::Int64)

    n < 0 && throw(ArgumentError("n must be ≥ 0."))
    
    ndims(x) == 1 && return vcat(x, zeros(eltype(x), n))
    ndims(x) == 2 && return hcat(x, zeros(eltype(x), size(x, 1), n))
    ndims(x) == 3 && return hcat(x, zeros(eltype(x), size(x, 1), n, size(x, 3)))
    ndims(x) > 3 && throw(ArgumentError("pad0() works only for 1-, 2- or 3-dimension array."))

end

"""
    pad2(x)

Pad vector / rows of matrix / array with zeros to the nearest power of 2 length.

# Arguments

- `x::Union{AbstractVector, AbstractArray`

# Returns

- `pad2::Union{AbstractVector, AbstractArray`
"""
function pad2(x::Union{AbstractVector, AbstractArray})

    ndims(x) == 1 && return pad0(x, nextpow2(length(x)) - length(x))
    ndims(x) == 2 && return hcat(x, zeros(eltype(x), size(x, 1), nextpow2(size(x, 2)) - size(x, 2)))
    ndims(x) == 3 && return hcat(x, zeros(eltype(x), size(x, 1), nextpow2(size(x, 2)) - size(x, 2), size(x, 3)))
    ndims(x) > 3 && throw(ArgumentError("pad2() works only for 1-, 2- or 3-dimension array."))
    
end
