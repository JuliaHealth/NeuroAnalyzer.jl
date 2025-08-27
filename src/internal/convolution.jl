function _remove_kernel(s::AbstractVector, kernel::AbstractVector)::AbstractVector
    # remove in- and out- edges
    n_kernel = length(kernel)
    half_kernel = floor(Int64, n_kernel / 2)
    if mod(n_kernel, 2) == 0
        s_new = s[half_kernel:(end - half_kernel)]
    else
        s_new = s[half_kernel:(end - half_kernel - 1)]
    end
    return s_new
end

_group_delay(kernel::AbstractVector)::Int64 = round(Int64, (length(kernel) - 1) / 2)
