function _gradient(x::Union{AbstractVector, AbstractMatrix, AbstractArray}; inverse::Bool=false)
    itp = inverse == true ? linear_interpolation(axes(x), -x) : linear_interpolation(axes(x), x)
    g_tmp = [Interpolations.gradient(itp, idx...) for idx in knots(itp)]
    g_len = norm.(g_tmp)
    return g_tmp, g_len
end
