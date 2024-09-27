export gradient

"""
    gradient(x; <keyword arguments>)

Calculate gradient of a 1-dimensional scalar field.

# Arguments

- `x::AbstractVector`
- `rev::Bool=false`: by default the direction of the gradient vector field is towards maximum value, if `rev=true`, the direction is towards the minimum value

# Returns

Named tuple containing:
- `g::Vector{Vector{Float64}}`: vector field of gradients
- `g_len::Vector{Float64}`: scalar field of gradient lengths
"""
function gradient(x::AbstractVector; rev::Bool=false)::NamedTuple{(:g, :g_len), Tuple{Vector{Vector{Float64}}, Vector{Float64}}}

    g_tmp, g_len = _gradient(x, rev=rev)

    g = Vector{Vector{Float64}}(undef, length(g_tmp))
    for idx in CartesianIndices(g_tmp)
        g[idx] = g_tmp[idx]
    end

    return (g=g, g_len=g_len)

end

"""
    gradient(x; <keyword arguments>)

Calculate gradient of a 2-dimensional scalar field.

# Arguments

- `x::AbstractMatrix`
- `rev::Bool=false`: by default the direction of the gradient vector field is towards maximum value, if `rev=true`, the direction is towards the minimum value

# Returns

Named tuple containing:
- `g::Matrix{Vector{Float64}}`: vector field of gradients
- `g_len::Matrix{Float64}`: scalar field of gradient lengths
"""
function gradient(x::AbstractMatrix; rev::Bool=false)::NamedTuple{(:g, :g_len), Tuple{Matrix{Vector{Float64}}, Matrix{Float64}}}

    g_tmp, g_len = _gradient(x, rev=rev)

    g = Matrix{Vector{Float64}}(undef, size(g_tmp))
    for idx in CartesianIndices(g_tmp)
        g[idx] = g_tmp[idx]
    end

    return (g=g, g_len=g_len)

end

"""
    gradient(x; <keyword arguments>)

Calculate gradient of a ≥3-dimensional scalar field.

# Arguments

- `x::AbstractArray`
- `rev::Bool=false`: by default the direction of the gradient vector field is towards maximum value, if `rev=true`, the direction is towards the minimum value

# Returns

Named tuple containing:
- `g::Array{Vector{Float64}, 3}`: vector field of gradients
- `g_len::Array{Float64, 3}`: scalar field of gradient lengths
"""
function gradient(x::AbstractArray; rev::Bool=false)::NamedTuple{(:g, :g_len), Tuple{Array{Vector{Float64}, 3}, Array{Float64}}}


    g_tmp, g_len = _gradient(x, rev=rev)

    g = Array{Vector{Float64}}(undef, size(g_tmp))
    for idx in CartesianIndices(g_tmp)
        g[idx] = g_tmp[idx]
    end

    return (g=g, g_len=g_len)

end
