export gradient

"""
    gradient(x; <keyword arguments>)

Calculate the gradient of a 1-dimensional scalar field.

Each element of the returned vector field `grad_vf` is a 1-element `Vector{Float64}` giving the gradient direction at that position. `grad_mag` contains the corresponding gradient magnitudes.

# Arguments

- `x::AbstractVector`: 1-D scalar field
- `rev::Bool=false`: if `false` (default), the gradient direction points toward the maximum value; if `true`, it points toward the minimum value

# Returns

Named tuple:

- `grad_vf::Vector{Vector{Float64}}`: vector field of gradients (one gradient vector per element)
- `grad_mag::Vector{Float64}`: scalar field of gradient magnitudes
"""
function gradient(
        x::AbstractVector;
        rev::Bool = false,
    )::@NamedTuple{
        grad_vf::Vector{Vector{Float64}},
        grad_mag::Vector{Float64},
    }
    g_tmp, grad_mag = _gradient(x; rev = rev)

    grad_vf = Vector{Vector{Float64}}(undef, length(g_tmp))
    # copy each gradient vector from the internal representation into grad_vf
    copyto!(grad_vf, g_tmp)

    return (; grad_vf, grad_mag)
end

"""
    gradient(x; <keyword arguments>)

Calculate the gradient of a 2-dimensional scalar field.

Each element of the returned matrix field `grad_vf` is a 2-element `Vector{Float64}` giving the gradient direction (row, column) at that position. `grad_mag` contains the corresponding gradient magnitudes.

# Arguments

- `x::AbstractMatrix`: 2-D scalar field
- `rev::Bool=false`: if `false` (default), the gradient direction points toward the maximum value; if `true`, it points toward the minimum value

# Returns

Named tuple:

- `grad_vf::Matrix{Vector{Float64}}`: vector field of gradients (one gradient vector per element)
- `grad_mag::Matrix{Float64}`: scalar field of gradient magnitudes
"""
function gradient(
        x::AbstractMatrix;
        rev::Bool = false,
    )::@NamedTuple{
        grad_vf::Matrix{Vector{Float64}},
        grad_mag::Matrix{Float64},
    }
    g_tmp, grad_mag = _gradient(x; rev = rev)

    grad_vf = Matrix{Vector{Float64}}(undef, size(g_tmp))
    # copy each gradient vector from the internal representation into grad_vf
    copyto!(grad_vf, g_tmp)

    return (; grad_vf, grad_mag)
end

"""
    gradient(x; <keyword arguments>)

Calculate the gradient of a 3-or-higher-dimensional scalar field.

Dispatches when `x` is neither a vector nor a matrix (i.e. `ndims(x) ≥ 3`). Each element of the returned array field `grad_vf` is an `ndims(x)`-element `Vector{Float64}` giving the gradient direction at that position. `grad_mag` contains the corresponding gradient magnitudes.

# Arguments

- `x::AbstractArray`: scalar field with `ndims(x) ≥ 3`
- `rev::Bool=false`: if `false` (default), the gradient direction points toward the maximum value; if `true`, it points toward the minimum value


# Returns

Named tuple:

- `grad_vf::Array{Vector{Float64}, 3}`: vector field of gradients
- `grad_mag::Array{Float64, 3}`: scalar field of gradient magnitudes
"""
function gradient(
        x::AbstractArray;
        rev::Bool = false,
    )::@NamedTuple{
        grad_vf::Array{Vector{Float64}, 3},
        grad_mag::Array{Float64, 3},
    }
    g_tmp, grad_mag = _gradient(x; rev = rev)

    grad_vf = Array{Vector{Float64}}(undef, size(g_tmp))
    # copy each gradient vector from the internal representation into grad_vf
    copyto!(grad_vf, g_tmp)

    return (; grad_vf, grad_mag)
end
