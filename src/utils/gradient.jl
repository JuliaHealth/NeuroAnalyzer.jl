export gradient

"""
   gradient(x; l, demean)

Calculate gradient of a 1-dimensional scalar field.

# Arguments

- `x::AbstractVector`
- `inverse::Bool=false`: by default the direction of the gradient vector field is towards maximum value, if `inverse=true`, the direction is towards the minimum value

# Returns

Named tuple containing:
- `g::Vector{Vector{Float64}}`: vector field of gradients
- `g::Vector{Float64}`: scalar field of gradient lengths
"""
function gradient(x::AbstractVector; inverse::Bool=false)

    g_tmp, g_len = _gradient(x, inverse=inverse)

    g = Vector{Vector{Float64}}(undef, length(g_tmp))
    for idx in CartesianIndices(g_tmp)
        g[idx] = g_tmp[idx]
    end

    return (g=g, g_len=g_len)

end

"""
   gradient(x; l, demean)

Calculate gradient of a 2-dimensional scalar field.

# Arguments

- `x::AbstractMatrix`
- `inverse::Bool=false`: by default the direction of the gradient vector field is towards maximum value, if `inverse=true`, the direction is towards the minimum value

# Returns

Named tuple containing:
- `g::Matrix{Vector{Float64}}`: vector field of gradients
- `g::Matrix{Float64}`: scalar field of gradient lengths
"""
function gradient(x::AbstractMatrix; inverse::Bool=false)

    g_tmp, g_len = _gradient(x, inverse=inverse)

    g = Matrix{Vector{Float64}}(undef, size(g_tmp))
    for idx in CartesianIndices(g_tmp)
        g[idx] = g_tmp[idx]
    end

    return (g=g, g_len=g_len)

end

"""
   gradient(x; l, demean)

Calculate gradient of a ≥3-dimensional scalar field.

# Arguments

- `x::AbstractArray`
- `inverse::Bool=false`: by default the direction of the gradient vector field is towards maximum value, if `inverse=true`, the direction is towards the minimum value

# Returns

Named tuple containing:
- `g::Array{Vector{Float64}, 3}`: vector field of gradients
- `g::Array{Float64, 3}`: scalar field of gradient lengths
"""
function gradient(x::AbstractArray; inverse::Bool=false)

    g_tmp, g_len = _gradient(x, inverse=inverse)

    g = Array{Vector{Float64}}(undef, size(g_tmp))
    for idx in CartesianIndices(g_tmp)
        g[idx] = g_tmp[idx]
    end

    return (g=g, g_len=g_len)

end
