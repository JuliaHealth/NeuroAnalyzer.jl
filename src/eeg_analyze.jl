################################
#                              #
# Low-level internal functions #
#                              #
################################

################################

"""
    eeg_total_power(eeg)

Calculate total power of the `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `stp::Matrix{Float64}`: total power for each channel per epoch
"""
function eeg_total_power(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    fs = eeg_sr(eeg)
    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)
    stp = zeros(channel_n, epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for channel in 1:channel_n
            s = @view eeg.eeg_signals[channel, :, epoch]
            stp[channel, epoch] = s_total_power(s, fs=fs)
        end
    end

    return stp
end

"""
    eeg_band_power(eeg; f)

Calculate absolute band power between frequencies `f[1]` and `f[2]` of the `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `f::Tuple(Union(Int64, Float64}, Union(Int64, Float64}}`: lower and upper frequency bounds

# Returns

- `sbp::Matrix{Float64}`: band power for each channel per epoch
"""
function eeg_band_power(eeg::NeuroJ.EEG; f::Tuple)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    fs = eeg_sr(eeg)
    length(f) != 2 && throw(ArgumentError("f must contain two frequencies."))
    f = tuple_order(f)
    f[1] <= 0 && throw(ArgumentError("Lower frequency bound must be be > 0."))
    f[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be be < $(fs / 2)."))

    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)
    sbp = zeros(channel_n, epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for channel in 1:channel_n
            s = @view eeg.eeg_signals[channel, :, epoch]
            sbp[channel, epoch] = s_band_power(s, fs=fs, f=f)
        end
    end

    return sbp
end

"""
    eeg_cov(eeg; norm=true)

Calculate covariance matrix for all channels of `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `norm::Bool`: normalize covariance

# Returns

- `cov_mat::Array{Float64, 3}`: covariance matrix for each epoch
"""
function eeg_cov(eeg::NeuroJ.EEG; norm::Bool=true)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)
    cov_mat = zeros(channel_n, channel_n, epoch_n)

    Threads.@threads for epoch in 1:epoch_n
        s = @view eeg.eeg_signals[:, :, epoch]
        cov_mat[:, :, epoch] = cov(s')
    end

    # normalize covariance matrix
    norm == true && (cov_mat = m_norm(cov_mat))

    return cov_mat
end

"""
    eeg_cor(eeg)

Calculate correlation coefficients between all channels of `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `cov_mat::Array{Float64, 3}`: correlation matrix for each epoch
"""
function eeg_cor(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)
    cor_mat = zeros(channel_n, channel_n, epoch_n)

    Threads.@threads for epoch in 1:epoch_n
        s = @view eeg.eeg_signals[:, :, epoch]
        cor_mat[:, :, epoch] = cor(s')
    end

    return cor_mat
end

"""
    eeg_crosscov(eeg; lag, demean, norm)

Calculate cross-covariance of each the `eeg` channels.

# Arguments

- `eeg::NeuroJ.EEG`
- `lag::Int64=1`: lags range is `-lag:lag`
- `demean::Bool=false`: demean signal prior to analysis
- `norm::Bool=false`: normalize cross-covariance

# Returns

Named tuple containing:
- `ccov::Matrix{Float64}`
- `lags::Vector{Float64}
"""
function eeg_crosscov(eeg::NeuroJ.EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)
    lags = collect(-lag:lag)
    ccov = zeros(channel_n^2, length(lags), epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        ccov_packed = Array{Vector{Float64}}(undef, channel_n, channel_n)
        Threads.@threads for idx1 in 1:channel_n
            for idx2 in 1:channel_n
                s1 = @view eeg.eeg_signals[idx1, :, epoch]
                s2 = @view eeg.eeg_signals[idx2, :, epoch]
                ccov_packed[idx1, idx2], lags = s_xcov(s1,
                                                       s2,
                                                       lag=lag,
                                                       demean=demean,
                                                       norm=norm)
            end
        end
        for idx in 1:channel_n^2
            ccov[idx, :, epoch] = ccov_packed[idx]
        end
    end

    lags = (eeg.eeg_time[2] - eeg.eeg_time[1]) .* collect(-lag:lag)

    return (ccov=ccov, ccov_lags=lags)
end

"""
    eeg_crosscov(eeg1, eeg2; lag, demean, norm)

Calculate cross-covariance between `eeg1` and `eeg2` channels.

# Arguments

- `eeg1::NeuroJ.EEG`
- `eeg2::NeuroJ.EEG`
- `lag::Int64=1`: lags range is `-lag:lag`
- `demean::Bool=false`: demean signal prior to analysis
- `norm::Bool=false`: normalize cross-covariance

# Returns

Named tuple containing:
- `ccov::Matrix{Float64}`
- `lags::Vector{Float64}
"""
function eeg_crosscov(eeg1::NeuroJ.EEG, eeg2::NeuroJ.EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)

    eeg_channel_n(eeg1, type=:eeg) < eeg_channel_n(eeg1, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))
    eeg_channel_n(eeg2, type=:eeg) < eeg_channel_n(eeg2, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s1 = @view eeg1.eeg_signals[idx, :, epoch]
            s2 = @view eeg2.eeg_signals[idx, :, epoch]
            ccov[idx, :, epoch], lags = s_xcov(s1,
                                               s2,
                                               lag=lag,
                                               demean=demean,
                                               norm=norm)
        end
    end

    lags = (eeg1.eeg_time[2] - eeg1.eeg_time[1]) .* collect(-lag:lag)

    return (ccov=ccov, ccov_lags=lags)
end

"""
    eeg_psd(eeg; norm)

Calculate total power for each the `eeg` channels.

# Arguments

- `eeg::NeuroJ.EEG`
- `norm::Bool=false`: normalize do dB

# Returns

Named tuple containing:
- `psd_pow::Array{Float64, 3}`:powers
- `psd_frq::Array{Float64, 3}`: frequencies
"""
function eeg_psd(eeg::NeuroJ.EEG; norm::Bool=false)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    fs = eeg_sr(eeg)
    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)
    psd_len, _ = s_psd(eeg.eeg_signals[1, :, 1], fs=fs, norm=norm)
    psd_pow = zeros(channel_n, length(psd_len), epoch_n)
    psd_frq = zeros(channel_n, length(psd_len), epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = eeg.eeg_signals[idx, :, epoch]
            psd_pow[idx, :, epoch], psd_frq[idx, :, epoch] = s_psd(s,
                                                                   fs=fs,
                                                                   norm=norm)
        end
    end

    return (psd_pow=psd_pow, psd_frq=psd_frq)
end

"""
    eeg_stationarity(eeg; window, method)

Calculate stationarity.

# Arguments

- `eeg:EEG`
- `window::Int64=10`: time window in samples
- `method::Symbol=:euclid`: stationarity method: :mean, :var, :euclid, :hilbert

# Returns

- `stationarity::Union{Matrix{Float64}, Array{Float64, 3}}`
"""
function eeg_stationarity(eeg::NeuroJ.EEG; window::Int64=10, method::Symbol=:hilbert)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))
    method in [:mean, :var, :euclid, :hilbert] || throw(ArgumentError("Method must be must be :mean, :var, :euclid or :hilbert."))
    (typeof(window) == Int64 && window < 1) && throw(ArgumentError("window must be ≥ 1."))
    (typeof(window) == Int64 && window > size(eeg.eeg_signals, 2)) && throw(ArgumentError("window must be ≤ $(size(eeg.eeg_signals, 2))."))

    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)

    if method === :mean
        s_stationarity = zeros(channel_n, window, epoch_n)
        @inbounds @simd for epoch in 1:epoch_n
            Threads.@threads for idx in 1:channel_n
                s = @view eeg.eeg_signals[idx, :, epoch]
                s_stationarity[idx, :, epoch] = s_stationarity_mean(s, window=window)
            end
        end
    end

    if method === :var
        s_stationarity = zeros(channel_n, window, epoch_n)
        @inbounds @simd for epoch in 1:epoch_n
            Threads.@threads for idx in 1:channel_n
                s = @view eeg.eeg_signals[idx, :, epoch]
                s_stationarity[idx, :, epoch] = s_stationarity_var(s, window=window)
            end
        end
    end

    if method === :hilbert
        s_stationarity = zeros(channel_n, eeg_epoch_len(eeg) - 1, epoch_n)
        @inbounds @simd for epoch in 1:epoch_n
            Threads.@threads for idx in 1:channel_n
                s = @view eeg.eeg_signals[idx, :, epoch]
                s_stationarity[idx, :, epoch] = s_stationarity_hilbert(s)
            end
        end
    end

    if method === :euclid
        # number of time windows per epoch
        window_n = eeg_epoch_len(eeg)
        cov_mat = zeros(channel_n, channel_n, window_n, epoch_n)
        s_stationarity = zeros(1 + length(2:window:window_n), epoch_n)

        @inbounds @simd for epoch in 1:epoch_n
            Threads.@threads for idx = 1:window_n
                s = @view eeg.eeg_signals[:, idx, epoch]
                cov_mat[:, :, idx, epoch] = s2_cov(s, s)
            end
        end

        @inbounds @simd for epoch in 1:epoch_n
            phase_idx = 1
            Threads.@threads for idx = 2:window:window_n
                s_stationarity[phase_idx, epoch] = euclidean(cov_mat[:, :, idx - 1, epoch],
                                                             cov_mat[:, :, idx, epoch])
                phase_idx += 1
            end
        end
    end

    return s_stationarity
end

"""
    eeg_mi(eeg)

Calculate mutual information between all channels of `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `mi::Array{Float64, 3}`
"""
function eeg_mi(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)
    mi = zeros(channel_n, channel_n, epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx1 in 1:channel_n
            for idx2 in 1:channel_n
                s = @view eeg.eeg_signals[idx2, :, epoch]
                mi[idx1, idx2, epoch] = s2_mi(s, s)
            end
        end
    end
    return mi
end

"""
    eeg_mi(eeg1, eeg2)

Calculate mutual information between all channels of `eeg1` and `eeg2`.

# Arguments

- `eeg1::NeuroJ.EEG`
- `eeg2::NeuroJ.EEG`

# Returns

- `mi::Array{Float64, 3}`
"""
function eeg_mi(eeg1::NeuroJ.EEG, eeg2::NeuroJ.EEG)

    eeg_channel_n(eeg1, type=:eeg) < eeg_channel_n(eeg1, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))
    eeg_channel_n(eeg2, type=:eeg) < eeg_channel_n(eeg2, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))
    eeg_channel_n(eeg1) == eeg_channel_n(eeg2) || throw(ArgumentError("Both signals must have the same number of channels."))
    eeg_epoch_n(eeg1) == eeg_epoch_n(eeg2) || throw(ArgumentError("Both signals must have the same number of epochs."))

    channel_n = eeg_channel_n(eeg1)
    epoch_n = eeg_epoch_n(eeg1)
    mi = zeros(channel_n, channel_n, epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx1 in 1:channel_n
            for idx2 in 1:channel_n
                s1 = eeg1.eeg_signals[idx2, :, epoch]
                s2 = eeg2.eeg_signals[idx2, :, epoch]
                mi[idx1, idx2, epoch] = s2_mi(s1, s2)
            end
        end
    end

    return mi
end

"""
    eeg_entropy(eeg)

Calculate entropy of all channels of `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `entropy::Matrix{Float64}`
"""
function eeg_entropy(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)
    s_ent = zeros(channel_n, epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view eeg.eeg_signals[idx, :, epoch]
            s_ent[idx, epoch] = s_entropy(s)
        end
    end

    return s_ent
end

"""
    eeg_band(eeg, band)

Return frequency limits for a `band` range.

# Arguments

- `eeg:EEG`
- `band::Symbol`: name of band range: :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher. If lower or upper band frequency limit exceeds Nyquist frequency of `eeg`, than bound is truncated to `eeg` range.

# Returns

- `band_frequency::Tuple{Float64, Float64}`
"""
function eeg_band(eeg; band::Symbol)

    band in [:total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher] || throw(ArgumentError("band must be: :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower or :gamma_higher."))

    band === :total && (band_frequency = (0, (eeg_sr(eeg) / 2)))
    band === :delta && (band_frequency = (0.5, 4.0))
    band === :theta && (band_frequency = (4.0, 8.0))
    band === :alpha && (band_frequency = (8.0, 13.0))
    band === :beta && (band_frequency = (14.0, 30.0))
    band === :beta_high && (band_frequency = (25.0, 30.0))
    band === :gamma && (band_frequency = (30.0, 150.0))
    band === :gamma_1 && (band_frequency = (31.0, 40.0))
    band === :gamma_2 && (band_frequency = (41.0, 50.0))
    band === :gamma_lower && (band_frequency = (30.0, 80.0))
    band === :gamma_higher && (band_frequency = (80.0, 150.0))
    
    band_frequency[1] > eeg_sr(eeg) / 2 && (band_frequency = (eeg_sr(eeg) / 2, band_frequency[2]))
    band_frequency[2] > eeg_sr(eeg) / 2 && (band_frequency = (band_frequency[1], eeg_sr(eeg) / 2))

    return band_frequency
end

"""
    eeg_coherence(eeg1, eeg2)

Calculate coherence between all channels of `eeg1` and `eeg2`.

# Arguments

- `eeg1::NeuroJ.EEG`
- `eeg2::NeuroJ.EEG`

# Returns

- `coherence::Union{Matrix{Float64}, Array{ComplexF64, 3}}`
"""
function eeg_coherence(eeg1::NeuroJ.EEG, eeg2::NeuroJ.EEG)

    eeg_channel_n(eeg1, type=:eeg) < eeg_channel_n(eeg1, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))
    eeg_channel_n(eeg2, type=:eeg) < eeg_channel_n(eeg2, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))
    size(eeg1.eeg_signals) == size(eeg2.eeg_signals) || throw(ArgumentError("Both signals must have the same size."))

    channel_n = eeg_channel_n(eeg1)
    epoch_n = eeg_epoch_n(eeg1)
    coh = zeros(ComplexF64, size(eeg1.eeg_signals))

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s1 = @view eeg1.eeg_signals[idx, :, epoch]
            s2 = @view eeg2.eeg_signals[idx, :, epoch]
            coh[idx, :, epoch] = s2_coherence(s1, s2)            
        end
    end

    return coh
end

"""
    eeg_coherence(eeg; channel1, channel2, epoch1, epoch2)

Calculate coherence between `channel1`/`epoch1` and `channel2` of `epoch2` of `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `channel1::Int64`
- `channel2::Int64`
- `epoch1::Int64`
- `epoch2::Int64`

# Returns

- `coh::Vector{ComplexF64}`
"""
function eeg_coherence(eeg::NeuroJ.EEG; channel1::Int64, channel2::Int64, epoch1::Int64, epoch2::Int64)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    (channel1 < 0 || channel2 < 0 || epoch1 < 0 || epoch2 < 0) && throw(ArgumentError("channel1/epoch1/channel2/epoch2 must be > 0."))
    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)
    (channel1 > channel_n || channel2 > channel_n) && throw(ArgumentError("channel1/channel2 must be ≤ $(channel_n)."))
    (epoch1 > epoch_n || epoch2 > epoch_n) && throw(ArgumentError("epoch1/epoch2 must be ≤ $(epoch_n)."))

    s1 = @view eeg.eeg_signals[channel1, :, epoch1]
    s2 = @view eeg.eeg_signals[channel2, :, epoch2]
    coh = s2_coherence(s1, s2)

    return coh
end

"""
    eeg_freqs(eeg)

Return vector of frequencies and Nyquist frequency for `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

Named tuple containing:
- `hz::Vector{Float64}`
- `nyquist::Float64`
"""
function eeg_freqs(eeg::NeuroJ.EEG)

    hz, nyq = s_freqs(eeg.eeg_signals[1, :, 1], eeg_sr(eeg))

    return (hz=hz, nyquist=nyq)
end

"""
    eeg_difference(eeg1, eeg2; n, method)

Calculate mean difference and its 95% CI between `eeg1` and `eeg2`.

# Arguments

- `eeg1::NeuroJ.EEG`
- `eeg2::NeuroJ.EEG`
- `n::Int64=3`: number of bootstraps
- `method::Symbol=:absdiff`
    - `:absdiff`: maximum difference
    - `:diff2int`: integrated area of the squared difference

# Returns

Named tuple containing:
- `signals_statistic::Matrix{Float64}`
- `signals_statistic_single::Vector{Float64}`
- `p::Vector{Float64}`
"""
function eeg_difference(eeg1::NeuroJ.EEG, eeg2::NeuroJ.EEG; n::Int64=3, method::Symbol=:absdiff)

    eeg_channel_n(eeg1, type=:eeg) < eeg_channel_n(eeg1, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))
    eeg_channel_n(eeg2, type=:eeg) < eeg_channel_n(eeg2, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    epoch_n = size(eeg1.eeg_signals, 3)
    signals_statistic = zeros(epoch_n, size(eeg1.eeg_signals, 1) * n)
    signals_statistic_single = zeros(epoch_n)
    p = zeros(epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        signals_statistic[epoch, :], signals_statistic_single[epoch], p[epoch] = s2_difference(eeg1.eeg_signals[:, :, epoch], eeg2.eeg_signals[:, :, epoch], n=n, method=method)
    end

    return (s_stat=signals_statistic, s_stat_single=signals_statistic_single, p=p)
end

"""
    eeg_picks(eeg; pick)

Return `pick` of electrodes for `eeg` electrodes.

# Arguments

- `pick::Vector{Symbol}`

# Returns

- `channels::Vector{Int64}`
"""
function eeg_pick(eeg::NeuroJ.EEG; pick::Union{Symbol, Vector{Symbol}})

    length(eeg_labels(eeg)) == 0 && throw(ArgumentError("EEG does not contain channel labels."))

    if typeof(pick) == Vector{Symbol}
        for idx in 1:length(pick)
            pick[idx] in [:list, :central, :c, :left, :l, :right, :r, :frontal, :f, :temporal, :t, :parietal, :p, :occipital, :o] || throw(ArgumentError("pick must be: :central, :c, :left, :l, :right, :r, :frontal, :f, :temporal, :t, :parietal, :p, :occipital, :o"))
        end

        c = Vector{Char}()
        for idx in 1:length(pick)
            (pick[idx] === :central || pick[idx] === :c) && push!(c, 'z')
            (pick[idx] === :frontal || pick[idx] === :f) && push!(c, 'F')
            (pick[idx] === :temporal || pick[idx] === :t) && push!(c, 'T')
            (pick[idx] === :parietal || pick[idx] === :p) && push!(c, 'P')
            (pick[idx] === :occipital || pick[idx] === :o) && push!(c, 'O')
        end
        
        labels = eeg_labels(eeg)
        channels = Vector{Int64}()
        for idx1 in 1:length(labels)
            for idx2 in 1:length(c)
                in(c[idx2], labels[idx1]) && push!(channels, idx1)
            end
        end

        # check for both :l and :r
        for idx1 in 1:length(pick)
            if (pick[idx1] === :left || pick[idx1] === :l)
                for idx2 in 1:length(pick)
                    if (pick[idx2] === :right || pick[idx2] === :r)
                        return channels
                    end
                end
            end
            if (pick[idx1] === :right || pick[idx1] === :r)
                for idx2 in 1:length(pick)
                    if (pick[idx2] === :left || pick[idx2] === :l)
                        return channels
                    end
                end
            end
        end

        labels = eeg_labels(eeg)
        labels = labels[channels]
        pat = nothing
        for idx in 1:length(pick)
            # for :right remove lefts
            (pick[idx] === :right || pick[idx] === :r) && (pat = r"[z13579]$")
            # for :left remove rights
            (pick[idx] === :left || pick[idx] === :l) && (pat = r"[z02468]$")
        end
        if typeof(pat) == Regex
            for idx in length(labels):-1:1
                typeof(match(pat, labels[idx])) == RegexMatch && deleteat!(channels, idx)
            end
        end

        return channels
    else
        pick in [:central, :c, :left, :l, :right, :r, :frontal, :f, :temporal, :t, :parietal, :p, :occipital, :o] || throw(ArgumentError("pick must be: :central, :c, :left, :l, :right, :r, :frontal, :f, :temporal, :t, :parietal, :p, :occipital, :o"))

        c = Vector{Char}()
        (pick === :central || pick === :c) && (c = ['z'])
        (pick === :left || pick === :l) && (c = ['1', '3', '5', '7', '9'])
        (pick === :right || pick === :r) && (c = ['2', '4', '6', '8'])
        (pick === :frontal || pick === :f) && (c = ['F'])
        (pick === :temporal || pick === :t) && (c = ['T'])
        (pick === :parietal || pick === :p) && (c = ['P'])
        (pick === :occipital || pick === :o) && (c = ['O'])

        labels = eeg_labels(eeg)
        channels = Vector{Int64}()
        for idx1 in 1:length(c)
            for idx2 in 1:length(labels)
                in(c[idx1], labels[idx2]) && push!(channels, idx2)
            end
        end

        return channels
    end
end

"""
    eeg_epochs_stats(eeg)

Calculate `eeg` epochs statistics.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

Named tuple containing:
- `e_mean::Vector(Float64)`: mean
- `e_median::Vector(Float64)`: median
- `e_std::Vector(Float64)`: standard deviation
- `e_var::Vector(Float64)`: variance
- `e_kurt::Vector(Float64)`: kurtosis
- `e_mean_diff::Vector(Float64)`: mean diff value
- `e_median_diff::Vector(Float64)`: median diff value
- `e_max_dif::Vector(Float64)`: max difference
- `e_dev_mean::Vector(Float64)`: deviation from channel mean
"""
function eeg_epochs_stats(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    epoch_n = eeg_epoch_n(eeg)
    e_mean = zeros(epoch_n)
    e_median = zeros(epoch_n)
    e_std = zeros(epoch_n)
    e_var = zeros(epoch_n)
    e_kurt = zeros(epoch_n)
    e_mean_diff = zeros(epoch_n)
    e_median_diff = zeros(epoch_n)
    e_max_dif = zeros(epoch_n)
    e_dev_mean = zeros(epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        s = @view eeg.eeg_signals[:, :, epoch]
        e_mean[epoch] = mean(s)
        e_median[epoch] = median(s)
        e_std[epoch] = std(s)
        e_var[epoch] = var(s)
        e_kurt[epoch] = kurtosis(s)
        e_mean_diff = mean(diff(s, dims=2))
        e_median_diff = median(diff(s, dims=2))
        e_max_dif = maximum(s) - minimum(s)
        e_dev_mean = abs(mean(s)) - mean(s)
    end

    return (e_mean=e_mean, e_median=e_median, e_std=e_std, e_var=e_var, e_kurt=e_kurt, e_mean_diff=e_mean_diff, e_median_diff=e_median_diff, e_max_dif=e_max_dif, e_dev_mean=e_dev_mean)
end

"""
    eeg_spectrogram(eeg; norm, demean)

Return spectrogram of `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `norm::Bool`=true: normalize powers to dB
- `demean::Bool`=true: demean signal prior to analysis

# Returns

Named tuple containing:
- `s_pow::Array{Float64, 3}`
- `s_frq::Matrix{Float64}`
- `s_t::Matrix{Float64}`
"""
function eeg_spectrogram(eeg::NeuroJ.EEG; norm::Bool=true, demean::Bool=true)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)
    fs = eeg_sr(eeg)
    p_tmp, f_tmp, t_tmp = s_spectrogram(eeg.eeg_signals[1, :, 1], fs=fs, norm=norm, demean=demean)
    s_pow = zeros(size(p_tmp, 1), size(p_tmp, 2), channel_n, epoch_n)
    s_frq = zeros(length(f_tmp), epoch_n)
    s_t = zeros(length(t_tmp), epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view eeg.eeg_signals[idx, :, epoch]
            s_pow[:, :, idx, epoch], s_frq[:, epoch], s_t[:, epoch] = s_spectrogram(s, fs=fs, norm=norm, demean=demean)
        end
    end

    return (s_pow=s_pow, s_frq=s_frq, s_t=s_t)
end

"""
    eeg_spectrum(eeg; pad)

Calculate FFT, amplitudes, powers and phases for each channel of the `eeg`. For `pad` > 0 channels are padded with 0s.

# Arguments

- `eeg::NeuroJ.EEG`
- `pad::Int64=0`: pad with `pad` zeros

# Returns

Named tuple containing:
- `fft::Array{ComplexF64, 3}`
- `amp::Array{Float64, 3}`
- `pow::Array{Float64, 3}`
- `phase::Array{Float64, 3}
"""
function eeg_spectrum(eeg::NeuroJ.EEG; pad::Int64=0)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)

    s_fft = zeros(ComplexF64, channel_n, eeg_epoch_len(eeg) + pad, epoch_n)
    s_amplitudes = zeros(channel_n, eeg_epoch_len(eeg) + pad, epoch_n)
    s_powers = zeros(channel_n, eeg_epoch_len(eeg) + pad, epoch_n)
    s_phases = zeros(channel_n, eeg_epoch_len(eeg) + pad, epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view eeg.eeg_signals[idx, :, epoch]
            s_fft[idx, :, epoch], s_amplitudes[idx, :, epoch], s_powers[idx, :, epoch], s_phases[idx, :, epoch] = s_spectrum(s, pad=pad)
        end
    end

    return (fft=s_fft, amp=s_amplitudes, pow=s_powers, phase=s_phases)
end

"""
    eeg_s2t(eeg; t)

Convert time `t` in samples to seconds using `eeg` sampling rate.

# Arguments

- `eeg::NeuroJ.EEG`
- `t::Int64`: time in samples

# Returns

- `t_s::Float64`: time in seconds
"""
function eeg_s2t(eeg::NeuroJ.EEG; t::Int64)
    t_s = round(t / eeg_sr(eeg), digits=2)
    
    return t_s
end

"""
    eeg_t2s(eeg; t)

Convert time `t` in seconds to samples using `eeg` sampling rate.

# Arguments

- `eeg::NeuroJ.EEG`
- `t::Union{Int64, Float64}`: time in seconds

# Returns

- `t_s::Float64`: time in samples
"""
function eeg_t2s(eeg::NeuroJ.EEG; t::Union{Int64, Float64})
    t_s = floor(Int64, t * eeg_sr(eeg)) + 1
    
    return t_s
end

"""
    eeg_channels_stats(eeg)

Calculate `eeg` channels statistics.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

Named tuple containing:
- `c_mean::Matrix(Float64)`: mean
- `c_median::Matrix(Float64)`: median
- `c_std::Matrix(Float64)`: standard deviation
- `c_var::Matrix(Float64)`: variance
- `c_kurt::Matrix(Float64)`: kurtosis
- `c_mean_diff::Matrix(Float64)`: mean diff value
- `c_median_diff::Matrix(Float64)`: median diff value
- `c_max_dif::Matrix(Float64)`: max difference
- `c_dev_mean::Matrix(Float64)`: deviation from channel mean
"""
function eeg_channels_stats(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)
    
    c_mean = zeros(channel_n, epoch_n)
    c_median = zeros(channel_n, epoch_n)
    c_std = zeros(channel_n, epoch_n)
    c_var = zeros(channel_n, epoch_n)
    c_kurt = zeros(channel_n, epoch_n)
    c_mean_diff = zeros(channel_n, epoch_n)
    c_median_diff = zeros(channel_n, epoch_n)
    c_max_dif = zeros(channel_n, epoch_n)
    c_dev_mean = zeros(channel_n, epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view eeg.eeg_signals[idx, :, epoch]
            c_mean[idx, epoch] = mean(s)
            c_median[idx, epoch] = median(s)
            c_std[idx, epoch] = std(s)
            c_var[idx, epoch] = var(s)
            c_kurt[idx, epoch] = kurtosis(s)
            c_mean_diff[idx, epoch] = mean(diff(s))
            c_median_diff[idx, epoch] = median(diff(s))
            c_max_dif[idx, epoch] = maximum(s) - minimum(s)
            c_dev_mean[idx, epoch] = abs(mean(s)) - mean(s)
        end
    end

    return (c_mean=c_mean, c_median=c_median, c_std=c_std, c_var=c_var, c_kurt=c_kurt, c_mean_diff=c_mean_diff, c_median_diff=c_median_diff, c_max_dif=c_max_dif, c_dev_mean=c_dev_mean)
end

"""
    eeg_snr(eeg)

Calculate SNR of `eeg` channels.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `snr::Matrix(Float64)`: SNR for each channel per epoch

# Source

D. J. Schroeder (1999). Astronomical optics (2nd ed.). Academic Press. ISBN 978-0-12-629810-9, p.278
"""
function eeg_snr(eeg::NeuroJ.EEG)

    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)
    snr = zeros(channel_n, epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view eeg.eeg_signals[idx, :, epoch]
            snr[idx, epoch] = s_snr(s)
        end
    end

    return snr
end

"""
    eeg_standardize(eeg)

Standardize `eeg` channels for ML.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `eeg_new::NeuroJ.EEG`: standardized EEG
- `scaler::Matrix{Float64}`: standardized EEG
"""
function eeg_standardize(eeg::NeuroJ.EEG)
    
    epoch_n = eeg_epoch_n(eeg)
    ss = similar(eeg.eeg_signals)
    scaler = Vector{Any}()

    @inbounds @simd for epoch in 1:epoch_n
        s = @view eeg.eeg_signals[:, :, epoch]
        push!(scaler, StatsBase.fit(ZScoreTransform, s, dims=2)) 
        ss[:,:, epoch] = StatsBase.transform(scaler[epoch], s)
    end

    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals = ss
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_standardize(EEG)")

    return eeg_new, scaler
end

"""
    eeg_standardize!(eeg)

Standardize `eeg` channels for ML.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `scaler::Matrix{Float64}`: standardized EEG
"""
function eeg_standardize!(eeg::NeuroJ.EEG)
    ss, scaler = s_standardize(eeg.eeg_signals)
    eeg.eeg_signals = ss
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_standardize!(EEG)")

    return scaler
end

"""
    eeg_fconv(eeg, kernel)

Perform convolution of all `eeg` channels in the frequency domain using `kernel`.

# Arguments

- `eeg::NeuroJ.EEG`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`: kernel for convolution

# Returns

- `s_convoluted::Union{Array{Float64, 3}, Array{ComplexF64, 3}}`: convoluted signal
"""
function eeg_fconv(eeg::NeuroJ.EEG; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)
    s_convoluted = zeros(ComplexF64, size(eeg.eeg_signals))

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view eeg.eeg_signals[idx, :, epoch]
            s_convoluted[idx, :, epoch] = s_fconv(s, kernel=kernel)
        end
    end

    return s_convoluted
end

"""
    eeg_tconv(eeg; kernel)

Perform convolution in the time domain.

# Arguments

- `eeg::NeuroJ.EEG`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`: kernel used for convolution

# Returns

- `s_convoluted::Union{Array{Float64, 3}, Array{ComplexF64, 3}}`: convoluted signal
"""
function eeg_tconv(eeg::NeuroJ.EEG; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)
    if typeof(kernel) == Vector{ComplexF64}
        s_convoluted = zeros(ComplexF64, size(eeg.eeg_signals))
    else
        s_convoluted = zeros(size(eeg.eeg_signals))
    end

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view eeg.eeg_signals[idx, :, epoch]
            s_convoluted[idx, :, epoch] = s_tconv(s, kernel=kernel)
        end
    end

    return s_convoluted
end

"""
    eeg_make_spectrum(eeg)

Returns FFT and DFT sample frequencies for a DFT for each the `eeg` channels.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

Named tuple containing:
- `fft::Array{ComplexF64, 3}`: FFT
- `sf::Array{Float64, 3}`: sample frequencies
"""
function eeg_dft(eeg::NeuroJ.EEG)

    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)
    fs = eeg_sr(eeg)
    s_fft = zeros(ComplexF64, size(eeg.eeg_signals))
    s_sf = zeros(size(eeg.eeg_signals))

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view eeg.eeg_signals[idx, :, epoch]
            s_fft[idx, :, epoch], s_sf[idx, :, epoch] = s_dft(s, fs=fs)
        end
    end

    return (fft=fft, sf=sf)
end

"""
    eeg_msci95(eeg; n::=3, method=:normal)

Calculates mean, std and 95% confidence interval for each the `eeg` channels.

# Arguments

- `eeg::NeuroJ.EEG`
- `n::Int64`: number of bootstraps
- `method::Symbol[:normal, :boot]`: use normal method or `n`-times boostrapping

# Returns

Named tuple containing:
- `s_m::Matrix{Float64}`: mean
- `s_s::Matrix{Float64}`: standard deviation
- `s_u::Matrix{Float64}`: upper 95% CI
- `s_l::Matrix{Float64}`: lower 95% CI
"""
function eeg_msci95(signal::Array{Float64, 3}; n::Int64=3, method::Symbol=:normal)

    method in [:normal, :boot] || throw(ArgumentError("method must be :normal or :boot."))
    n < 1 && throw(ArgumentError("n must be ≥ 1."))

    epoch_n = eeg_epoch_n(eeg)
    epoch_len = eeg_epoch_len(eeg)
    s_m = zeros(epoch_n, epoch_len)
    s_s = zeros(epoch_n, epoch_len)
    s_u = zeros(epoch_n, epoch_len)
    s_l = zeros(epoch_n, epoch_len)

    Threads.@threads for epoch in 1:epoch_n
        s = @view eeg.eeg_signals[:, :, epoch]
        s_m[idx, :], s_s[idx, :], s_u[idx, :], s_l[idx, :] = s_msci95(s, n=n, method=method)
    end

    return (mean=s_m, sd=s_s, upper=s_u, lower=s_l)
end

"""
    eeg_mean(eeg1, eeg2)

Calculates mean and 95% confidence interval for `eeg1` and `eeg2` channels.

# Arguments

- `eeg1::NeuroJ.EEG`
- `eeg2:NeuroJ.EEG`

# Returns

Named tuple containing:
- `s_m::Matrix{Float64}`: mean by epochs
- `s_s::Matrix{Float64}`: std by epochs
- `s_u::Matrix{Float64}`: upper 95% CI bound by epochs
- `s_l::Matrix{Float64}`: lower 95% CI bound by epochs
"""
function eeg_mean(eeg1::NeuroJ.EEG, eeg2::NeuroJ.EEG)

    size(eeg1.eeg_signals) != size(eeg2.eeg_signals) && throw(ArgumentError("Both signals must be of the same as size."))

    epoch_n = eeg_epoch_n(eeg1)
    e_len = eeg_epoch_len(eeg1)
    s_m = zeros(epoch_n, epoch_len)
    s_s = zeros(epoch_n, epoch_len)
    s_u = zeros(epoch_n, epoch_len)
    s_l = zeros(epoch_n, epoch_len)

    Threads.@threads for epoch in 1:epoch_n
        s1 = @view signal1[:, :, epoch]
        s2 = @view signal2[:, :, epoch]
        s1_mean = mean(s1, dims=1)
        s2_mean = mean(s2, dims=1)
        s_m[epoch, :] = s1_mean - s2_mean
        s1_sd = std(s1, dims=1) / sqrt(size(s1, 2))
        s2_sd = std(s2, dims=1) / sqrt(size(s2, 2))
        s_s[epoch, :] = sqrt.(s1_sd.^2 .+ s2_sd.^2)
        s_u[epoch, :] = @. s_m[epoch, :] + 1.96 * s_s[epoch, :]
        s_l[epoch, :] = @. s_m[epoch, :] - 1.96 * s_s[epoch, :]
    end

    return (mean=s_m, sd=s_s, upper=s_u, lower=s_l)
end

"""
    eeg_difference(eeg1, eeg2; n=3, method=:absdiff)

Calculates mean difference and 95% confidence interval for `eeg1` and `eeg2`.

# Arguments

- `eeg1::Array{Float64, 3}`
- `eeg2:Array{Float64, 3}`
- `n::Int64`: number of bootstraps
- `method::Symbol[:absdiff, :diff2int]`
    - `:absdiff`: maximum difference
    - `:diff2int`: integrated area of the squared difference

# Returns

Named tuple containing:
- `statistic::Matrix{Float64}`
- `statistic_single::Vector{Float64}`
- `p::Vector{Float64}`
"""
function eeg_difference(eeg1::Array{Float64, 3}, eeg2::Array{Float64, 3}; n::Int64=3, method::Symbol=:absdiff)

    size(eeg1.eeg_signals) != size(eeg2.eeg_signals) && throw(ArgumentError("Both signals must be of the same size."))
    method in [:absdiff, :diff2int] || throw(ArgumentError("method must be :absdiff or :diff2int."))

    epoch_n = eeg_epoch_n(eeg1)
    s_statistic = zeros(epoch_n, size(signal1, 1) * n)
    s_statistic_single = zeros(epoch_n)
    p = zeros(epoch_n)

    Threads.@threads for epoch in 1:epoch_n
        s1 = @view signal1[:, :, epoch]
        s2 = @view signal2[:, :, epoch]
        s_statistic[epoch, :], s_statistic_single[epoch], p[epoch] = s2_difference(s1, s2)
    end

    return (statistic=s_statistic, statsitic_single=s_statistic_single, p=p)
end

"""
   eeg_autocov(eeg; lag=1, demean=false, norm=false)

Calculate autocovariance of each the `eeg` channels.

# Arguments

- `eeg::NeuroJ.EEG`
- `lag::Int64`: lags range is `-lag:lag`
- `demean::Bool`: demean eeg prior to analysis
- `norm::Bool`: normalize autocovariance

# Returns

Named tuple containing:
- `acov::Matrix{Float64}`
- `lags::Vector{Float64}`
"""
function eeg_autocov(eeg::NeuroJ.EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)

    lag < 1 && throw(ArgumentError("lag must be ≥ 1."))

    lags = collect(-lag:lag)
    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)
    acov = zeros(channel_n, length(lags), epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view eeg.eeg_signals[idx, :, epoch]
            acov[idx, :, epoch], lags = s_acov(s,
                                               lag=lag,
                                               demean=demean,
                                               norm=norm)
        end
    end

    lags = (eeg.eeg_time[2] - eeg.eeg_time[1]) .* collect(-lag:lag)

    return (acov=acov, acov_lags=lags)
end