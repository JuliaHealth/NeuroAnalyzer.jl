export gradient

"""
    gradient(x; <keyword arguments>)

Calculate the gradient of a 1-dimensional scalar field.

Each element of the returned vector field `g` is a 1-element `Vector{Float64}` giving the gradient direction at that position. `g_len` contains the corresponding gradient magnitudes.

# Arguments

- `x::AbstractVector`: 1-D scalar field
- `rev::Bool=false`: if `false` (default), the gradient direction points toward the maximum value; if `true`, it points toward the minimum value

# Returns

Named tuple:

- `g::Vector{Vector{Float64}}`: vector field of gradients (one gradient vector per element)
- `g_len::Vector{Float64}`: scalar field of gradient magnitudes

# See also

[`gradient(::AbstractMatrix)`](@ref), [`gradient(::AbstractArray)`](@ref)
"""
function gradient(
    x::AbstractVector;
    rev::Bool = false
)::@NamedTuple{
    g::Vector{Vector{Float64}},
    g_len::Vector{Float64}
}

    g_tmp, g_len = _gradient(x, rev = rev)

    g = Vector{Vector{Float64}}(undef, length(g_tmp))
    # copy each gradient vector from the internal representation into g
    copyto!(g, g_tmp)

    return (g=g, g_len=g_len)

end

"""
    gradient(x; <keyword arguments>)

Calculate the gradient of a 2-dimensional scalar field.

Each element of the returned matrix field `g` is a 2-element `Vector{Float64}` giving the gradient direction (row, column) at that position. `g_len` contains the corresponding gradient magnitudes.

# Arguments

- `x::AbstractMatrix`: 2-D scalar field.
- `rev::Bool=false`: if `false` (default), the gradient direction points toward the maximum value; if `true`, it points toward the minimum value

# Returns

Named tuple:

- `g::Matrix{Vector{Float64}}`: vector field of gradients (one gradient vector per element)
- `g_len::Matrix{Float64}`: scalar field of gradient magnitudes

# See also

[`gradient(::AbstractVector)`](@ref), [`gradient(::AbstractArray)`](@ref)
"""
function gradient(
    x::AbstractMatrix;
    rev::Bool = false
)::@NamedTuple{
    g::Matrix{Vector{Float64}},
    g_len::Matrix{Float64}
}

    g_tmp, g_len = _gradient(x, rev = rev)

    g = Matrix{Vector{Float64}}(undef, size(g_tmp))
    # copy each gradient vector from the internal representation into g
    copyto!(g, g_tmp)
    return (g=g, g_len=g_len)

end

"""
    gradient(x; <keyword arguments>)

Calculate the gradient of a 3-or-higher-dimensional scalar field.

Dispatches when `x` is neither a vector nor a matrix (i.e. `ndims(x) ≥ 3`). Each element of the returned array field `g` is an `ndims(x)`-element `Vector{Float64}` giving the gradient direction at that position. `g_len` contains the corresponding gradient magnitudes.

# Arguments

- `x::AbstractArray`: scalar field with `ndims(x) ≥ 3`
- `rev::Bool=false`: if `false` (default), the gradient direction points toward the maximum value; if `true`, it points toward the minimum value


# Returns

Named tuple:

- `g::Array{Vector{Float64}, 3}`: vector field of gradients
- `g_len::Array{Float64, 3}`: scalar field of gradient magnitudes

# See also

[`gradient(::AbstractVector)`](@ref), [`gradient(::AbstractMatrix)`](@ref)
"""
function gradient(
    x::AbstractArray;
    rev::Bool = false
)::@NamedTuple{
    g::Array{Vector{Float64}, 3},
    g_len::Array{Float64, 3}
}

    g_tmp, g_len = _gradient(x; rev = rev)

    g = Array{Vector{Float64}}(undef, size(g_tmp))
    # copy each gradient vector from the internal representation into g
    copyto!(g, g_tmp)

    return (g=g, g_len=g_len)

end
