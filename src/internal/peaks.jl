function _tpt_peaks(x::AbstractVector, t::AbstractVector)::Vector{Int64}
    verbose_tmp = NeuroAnalyzer.verbose
    NeuroAnalyzer.verbose = false
    # f = filter_create(fprototype=:fir, ftype=:hp, cutoff=5, n=length(x), fs=50, order=8)
    # x = filter_apply(x, flt=f)
    f = filter_create(fprototype=:fir, ftype=:lp, cutoff=15, n=length(x), fs=50, order=8)
    x = filter_apply(x, flt=f)
    NeuroAnalyzer.verbose = verbose_tmp
    x[x .< 0] .= 0
    x = diff(x)
    t = t[1:end-1]
    x[x .< 0] .= 0
    dx = env_up(x, d=4, t)
    p_idx, _ = findpeaks1d(dx, distance=10, height=5)
    return p_idx
end

function _tpt_peaks(x::AbstractArray, t::AbstractVector)::Vector{Int64}
    p_idx = Vector{Vector{Int64}}()
    verbose_tmp = NeuroAnalyzer.verbose
    NeuroAnalyzer.verbose = false
    f = filter_create(fprototype=:fir, ftype=:lp, cutoff=5, n=length(x[1, :]), fs=50, order=8)
    for idx in 1:3
        tmp = diff(x[idx, :])
        tmp[tmp .< 0] .= 0
        tmp = filter_apply(tmp, flt=f)
        tmp[tmp .< 0] .= 0
        dx = env_up(tmp, d=4, t[2:end])
        idx == 1 && (h = 3)
        idx == 2 && (h = 3)
        idx == 3 && (h = 3)
        p_tmp, _ = findpeaks1d(dx, distance=50, height=h)
        push!(p_idx, p_tmp)
    end
    NeuroAnalyzer.verbose = verbose_tmp
    m = zeros(Int64, 3, length(x[1, :]))
    for idx in 1:3
        m[idx, p_idx[idx]] .= 1
    end
    p_idx = sum(m, dims=1)[:]
    p_idx, _ = findpeaks1d(p_idx, distance=10, height=1.0)
    return p_idx
end