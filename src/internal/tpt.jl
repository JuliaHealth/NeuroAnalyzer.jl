function _tpt_peaks(x::AbstractArray, t::AbstractVector)::Vector{Int64}
    x = detrend(x, type=:loess)
    x = x[:, :, 1]
    q2_x = 3.166866296141654
    q2_y = 4.0740732957574375
    q2_z = 5.614191777055647
    q2_x_acc = 1.5005706677814734
    q2_y_acc = 1.9154153635212732
    q2_z_acc = 2.620307540507893
    p_idx = Vector{Vector{Int64}}()
    verbose_tmp = NeuroAnalyzer.verbose
    NeuroAnalyzer.verbose = false
    f = filter_create(fprototype=:fir, ftype=:lp, cutoff=5, n=length(x[1, :]), fs=50, order=12)
    for idx in 1:6
        tmp = x[idx, :]
        tmp = filter_apply(tmp, flt=f)
        p, f, t = spectrogram(tmp, fs=50)
        tmp[tmp .< 0] .= 0
        idx == 1 && (h = q2_x)
        idx == 2 && (h = q2_y)
        idx == 3 && (h = q2_z)
        idx == 4 && (h = q2_x_acc)
        idx == 5 && (h = q2_y_acc)
        idx == 6 && (h = q2_z_acc)
        p_tmp, _ = findpeaks1d(tmp, distance=12, prominence=(h/2, h))
        # p_tmp, _ = findpeaks1d(tmp, distance=20)
        # p_tmp, _ = findpeaks1d(tmp, prominence=(1, h))
        push!(p_idx, p_tmp)
    end
    NeuroAnalyzer.verbose = verbose_tmp
    m = zeros(Int64, 6, length(x[1, :]))
    for idx in 1:6
        m[idx, p_idx[idx]] .= 1
    end
    p_idx = sum(m, dims=1)[:]
    p_idx, _ = findpeaks1d(p_idx, distance=25, height=0.2)
    return p_idx
end