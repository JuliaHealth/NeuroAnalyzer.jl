export pad0
export pad2
export padm

"""
    pad0(x, n)

Pad vector / rows of matrix / array with zeros. Works with 1-, 2- and 3-dimensional arrays.

# Arguments

- `x::Union{AbstractVector, AbstractArray}`
- `n::Int64`: number of zeros to add

# Returns

- `pad0::Union{AbstractVector, AbstractArray`
"""
function pad0(x::Union{AbstractVector, AbstractArray}, n::Int64)

    @assert n >= 0 "n must be ≥ 0."
    
    ndims(x) == 1 && return vcat(x, zeros(eltype(x), n))
    ndims(x) == 2 && return hcat(x, zeros(eltype(x), size(x, 1), n))
    ndims(x) == 3 && return hcat(x, zeros(eltype(x), size(x, 1), n, size(x, 3)))
    @assert ndims(x) <= 3 "pad0() works only for 1-, 2- or 3-dimension array."

end

"""
    pad2(x)

Pad vector / rows of matrix / array with zeros to the nearest power of 2 length.

# Arguments

- `x::Union{AbstractVector, AbstractArray}`

# Returns

- `pad2::Union{AbstractVector, AbstractArray`
"""
function pad2(x::Union{AbstractVector, AbstractArray})

    ndims(x) == 1 && return pad0(x, nextpow2(length(x)) - length(x))
    ndims(x) == 2 && return hcat(x, zeros(eltype(x), size(x, 1), nextpow2(size(x, 2)) - size(x, 2)))
    ndims(x) == 3 && return hcat(x, zeros(eltype(x), size(x, 1), nextpow2(size(x, 2)) - size(x, 2), size(x, 3)))
    @assert ndims(x) <= 3 "pad2() works only for 1-, 2- or 3-dimension array."
    
end

"""
    padm(x, n)

Pad vector / rows of matrix / array with mean value. Works with 1-, 2- and 3-dimensional arrays.

# Arguments

- `x::Union{AbstractVector, AbstractArray}`
- `n::Int64`: number of values to add
- `mode::Symbol=:row`: how the mean is calculated:
    - `:all`: mean of all rows
    - `:row`: separate mean per each row

# Returns

- `padm::Union{AbstractVector, AbstractArray`
"""
function padm(x::Union{AbstractVector, AbstractArray}, n::Int64; mode::Symbol=:all)

    _check_var(mode, [:all, :row], "mode")
    @assert n >= 0 "n must be ≥ 0."
    @assert ndims(x) <= 3 "pad0() works only for 1-, 2- or 3-dimension array."
    
    if ndims(x) == 1
        m = mean(x)
        return vcat(x, m .* ones(eltype(x), n))
    elseif ndims(x) == 2
        mode === :all && (m = mean(x))
        mode === :row && (m = mean(x, dims=2))
        return hcat(x, m .* ones(eltype(x), size(x, 1), n))
    elseif ndims(x) == 3
        mode === :all && (m = mean(x))
        mode === :row && (m = mean(x, dims=2))
        return hcat(x, m .* ones(eltype(x), size(x, 1), n, size(x, 3)))
    end

end