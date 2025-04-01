function _tpt_peaks(x::AbstractVector)::Vector{Int64}
    x = detrend(x, type=:loess)
    x = derivative(x)
    # number of windows
    # 50 is sampling rate
    # 0.5 is window length in second, 25 in samples
    wlen = 25
    n = round(Int64, length(x) / (50 * 0.5))
    p_idx = Int64[]
    for idx in 1:n
        w = (idx - 1) * wlen + 1:idx  * wlen
        t1 = w[1]
        t2 = w[end]
        p_idx_tmp, _ = findpeaks1d(x[t1:t2], prominence=4)
        length(p_idx_tmp) > 0 && push!(p_idx, t1 + vsearch(maximum(x[t1:t2][p_idx_tmp]), x[t1:t2]))
    end

    return p_idx

end