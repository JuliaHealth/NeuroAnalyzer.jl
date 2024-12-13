function _tpt_peaks(x::AbstractVector, t::AbstractVector)::Vector{Int64}
    x[x .< 0] .= 0
    x = diff(x)
    t = t[1:end-1]
    x[x .< 0] .= 0
    #x = derivative(x)
    dx = env_up(x, d=4, t)
    verbose_tmp = NeuroAnalyzer.verbose
    NeuroAnalyzer.verbose = false
    f = filter_create(fprototype=:fir, ftype=:lp, cutoff=15, n=length(dx), fs=50, order=8)
    NeuroAnalyzer.verbose = verbose_tmp
    dx = filter_apply(dx, flt=f)
    p_idx, _ = findpeaks1d(dx, distance=10, height=mean(dx) + 0.25 * std(dx))
    # p_idx, _ = findpeaks1d(dx2, distance=(sr(obj) ÷ 4), height=mean(dx) + 1.5 * std(dx))
    return p_idx
end