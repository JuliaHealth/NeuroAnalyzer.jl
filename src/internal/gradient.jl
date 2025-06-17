function _gradient(x::Union{AbstractVector, AbstractMatrix, AbstractArray}; rev::Bool=false)::Tuple{Union{Vector{Interpolations.SVector{1, Float64}}, Matrix{Interpolations.SVector{2, Float64}}, Array{Interpolations.SVector{3, Float64}}}, Union{Vector{Float64}, Matrix{Float64}, Array{Float64, 3}}}
    itp = rev ? linear_interpolation(axes(x), -x) : linear_interpolation(axes(x), x)
    g_tmp = [Interpolations.gradient(itp, idx...) for idx in knots(itp)]
    g_len = norm.(g_tmp)
    return g_tmp, g_len
end
