function _tpt_peaks(x::AbstractVector, t::AbstractVector)::Vector{Int64}
    x = remove_dc(x)
    #x[x .< 0] .= 0
    #x = diff(x)
    #t = t[1:end-1]
    #x[x .< 0] .= 0
    # x = derivative(x)
    dx = env_up(x, d=4, t)
    verbose_tmp = NeuroAnalyzer.verbose
    NeuroAnalyzer.verbose = false
    f = filter_create(fprototype=:fir, ftype=:lp, cutoff=20, n=length(dx), fs=50, order=8)
    dx = filter_apply(dx, flt=f)
    NeuroAnalyzer.verbose = verbose_tmp
    # p_idx, _ = findpeaks1d(dx, distance=10, height=mean(dx) + 0.5 * std(dx))
    p_idx, _ = findpeaks1d(dx, distance=20, height=5.5)
    return p_idx
end