"""Return pinch indices in TPT signal."""
function _tpt_peaks(x::AbstractVector)::Vector{Int64}
    _wip()
    x = _detrend(x; type = :loess)
    x = _derivative(x)
    # number of windows
    # 50 is sampling rate
    # 0.5 is window length in second, 25 in samples
    wlen = 25
    # only complete windows, so the last window cannot run past the end of x
    n = length(x) ÷ wlen
    p_idx = Int64[]
    for idx = 1:n
        t1 = (idx - 1) * wlen + 1
        t2 = idx * wlen
        xw = x[t1:t2]
        p_idx_tmp, _ = findpeaks1d(xw; prominence = 4)
        # highest peak in the window, converted from a window-local to a global index
        length(p_idx_tmp) > 0 && push!(p_idx, t1 - 1 + p_idx_tmp[argmax(xw[p_idx_tmp])])
    end

    return p_idx
end
