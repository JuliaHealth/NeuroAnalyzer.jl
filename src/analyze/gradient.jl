export gradient

"""
   gradient(z; l, demean)

Calculate auto-correlation

# Arguments

- `z::AbstractMatrix`
- `inverse::Bool=false`: by default the direction of the gradient vector field is towards maximum value, if `inverse=true`, the direction is towards the minimum value

# Returns

Named tuple containing:
- `g::Matrix{Vector{Float64}}`: vector field of gradients
- `g::Matrix{Float64}`: scalar field of gradient lengths
"""
function gradient(z::AbstractMatrix; inverse::Bool=false)

    itp = inverse == true ? linear_interpolation(axes(z), -z) : linear_interpolation(axes(z), -z)
    g_tmp = [Interpolations.gradient(itp, idx...) for idx in knots(itp)]
    g_len = norm.(g_tmp)

    g = Matrix{Vector{Float64}}(undef, size(g_tmp))
    for idx in CartesianIndices(g_tmp)
        g[idx] = g_tmp[idx]
    end

    return (g=g, g_len=g_len)

end
