export pad0
export pad2
export padm

"""
    pad0(x, n)

Append `n` zeros along the second axis (or at the end for 1-D input). Works with 1-, 2-, and 3-dimensional arrays. The zero values match the element type of `x`.

# Arguments

- `x::Union{AbstractVector, AbstractArray}`: input array; must be 1-, 2-, or 3-dimensional
- `n::Int64`: number of zeros to append; must be ≥ 0

# Returns

- `Union{AbstractVector, AbstractArray}`: padded array with the same number of dimensions as `x`

# Throws

- `ArgumentError`: if `n < 0` or `ndims(x) > 3`

# See also

[`pad2`](@ref), [`padm`](@ref)
"""
function pad0(
    x::Union{AbstractVector, AbstractArray},
    n::Int64
)::Union{AbstractVector, AbstractArray}

    @assert n >= 0 "n must be ≥ 0."
    @assert ndims(x) <= 3 "pad0() supports only 1-, 2-, or 3-dimensional arrays."

    # fast path: nothing to append
    n == 0 && return x
    if ndims(x) == 1
        return vcat(x, zeros(eltype(x), n))
    elseif ndims(x) == 2
        return hcat(x, zeros(eltype(x), size(x, 1), n))
    else
        # cat along dim 2 is safer than hcat for 3-D arrays
        return cat(x, zeros(eltype(x), size(x, 1), n, size(x, 3)); dims=2)
    end

end

"""
    pad2(x)

Pad an array with zeros along its second axis to the next power-of-2 length. Works with 1-, 2-, and 3-dimensional arrays. If the relevant dimension is already a power of 2, the array is returned unchanged.

# Arguments

- `x::Union{AbstractVector, AbstractArray}`: input array; must be 1-, 2-, or 3-dimensional

# Returns

- `Union{AbstractVector, AbstractArray}`: padded array whose second dimension (or length, for 1-D) is a power of 2.

# Throws

- `ArgumentError`: if `ndims(x) > 3`

# See also

[`pad0`](@ref), [`padm`](@ref)
"""
function pad2(
    x::Union{AbstractVector, AbstractArray}
)::Union{AbstractVector, AbstractArray}

    @assert ndims(x) <= 3 "pad2() supports only 1-, 2-, or 3-dimensional arrays."

    if ndims(x) == 1
        n = nextpow2(length(x)) - length(x)
        return pad0(x, n)
    elseif ndims(x) == 2
        n = nextpow2(size(x, 2)) - size(x, 2)
        return n == 0 ? x : hcat(x, zeros(eltype(x), size(x, 1), n))
    else
        n = nextpow2(size(x, 2)) - size(x, 2)
        # cat along dim 2 is safer than hcat for 3-D arrays
        return n == 0 ? x : cat(x, zeros(eltype(x), size(x, 1), n, size(x, 3)); dims=2)
    end

end

"""
    padm(x, n)

Pad an array with mean values along its second axis (or at the end for 1-D input). Works with 1-, 2-, and 3-dimensional arrays.

# Arguments

- `x::Union{AbstractVector, AbstractArray}`: input array; must be 1-, 2-, or 3-dimensional
- `n::Int64`: number of mean-value samples to append; must be ≥ 0
- `mode::Symbol=:row`: how the mean is computed for 2-D and 3-D arrays:
    - `:all`: single mean of all elements
    - `:row`: separate mean per row (broadcasts across the padded columns)

# Returns

- `Union{AbstractVector, AbstractArray}`: padded array with the same number of dimensions as `x`

# Throws

- `ArgumentError`: if `n < 0`, `ndims(x) > 3`, or `mode` is not `:all`/`:row`

# See also

[`pad0`](@ref), [`pad2`](@ref)
"""
function padm(
    x::Union{AbstractVector, AbstractArray},
    n::Int64;
    mode::Symbol = :all
)::Union{AbstractVector, AbstractArray}

    _check_var(mode, [:all, :row], "mode")
    @assert n >= 0 "n must be ≥ 0."
    @assert ndims(x) <= 3 "padm() supports only 1-, 2-, or 3-dimensional arrays."

    # fast path: nothing to append
    n == 0 && return x
    if ndims(x) == 1
        # mode has no effect for 1-D: there is only one "row"
        return vcat(x, fill(mean(x), n))
    elseif ndims(x) == 2
        m = mode === :all ? mean(x) : mean(x; dims=2)
        return hcat(x, m .* ones(eltype(x), size(x, 1), n))
    else
        m = mode === :all ? mean(x) : mean(x; dims=2)
        # cat along dim 2 is safer than hcat for 3-D arrays
        return cat(x, m .* ones(eltype(x), size(x, 1), n, size(x, 3)); dims=2)
    end

end
