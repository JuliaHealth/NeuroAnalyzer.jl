function _gradient(x::Union{AbstractVector, AbstractMatrix, AbstractArray}; rev::Bool=false)::Tuple{Union{Vector{StaticArraysCore.SVector{1, Float64}}, Matrix{StaticArraysCore.SVector{2, Float64}}, Array{StaticArraysCore.SVector{3, Float64}, 3}}, Union{Vector{Float64}, Matrix{Float64}, Array{Float64, 3}}}
    itp = rev ? LinearInterpolation(axes(x), -x) : LinearInterpolation(axes(x), x)
    g_tmp = [Interpolations.gradient(itp, idx...) for idx in knots(itp)]
    g_len = norm.(g_tmp)
    return g_tmp, g_len
end
