function _signal_total_power(signal::AbstractArray; fs::Int64)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))

    psd = welch_pgram(signal, 4*fs, fs=fs)
    # dx: frequency resolution
    dx = psd.freq[2] - psd.freq[1]
    stp = simpson(psd.power, dx=dx)

    return stp
end

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
        Threads.@threads for idx in 1:channel_n
            s = @view eeg.eeg_signals[idx, :, epoch]
            stp[idx, epoch] = _signal_total_power(s, fs=fs)
        end
    end

    return stp
end

function _signal_band_power(signal::AbstractArray; fs::Int64, f::Tuple)

    psd = welch_pgram(signal, 4*fs, fs=fs)
    psd_freq = Vector(psd.freq)
    
    f1_idx = vsearch(f[1], psd_freq)
    f2_idx = vsearch(f[2], psd_freq)
    frq_idx = [f1_idx, f2_idx]

    # dx: frequency resolution
    dx = psd_freq[2] - psd_freq[1]
    sbp = simpson(psd.power[frq_idx[1]:frq_idx[2]], psd_freq[frq_idx[1]:frq_idx[2]], dx=dx)

    return sbp
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
        Threads.@threads for idx in 1:channel_n
            s = @view eeg.eeg_signals[idx, :, epoch]
            sbp[idx, epoch] = _signal_band_power(s, fs=fs, f=f)
        end
    end

    return sbp
end

"""
    eeg_cov(eeg; norm=true)

Calculate covariance between all channels of `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `norm::Bool`: normalize covariance

# Returns

- `cov_mat::Array{Float64, 3}`
"""
function eeg_cov(eeg::NeuroJ.EEG; norm=true)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    cov_mat = signal_cov(eeg.eeg_signals, norm=norm)

    return cov_mat
end

"""
    eeg_cor(eeg)

Calculate correlation coefficients between all channels of `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `cov_mat::Array{Float64, 3}`
"""
function eeg_cor(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    cor_mat = signal_cor(eeg.eeg_signals)

    return cor_mat
end

"""
    eeg_autocov(eeg; lag, demean, norm)

Calculate autocovariance of each the `eeg` channels.

# Arguments

- `eeg::NeuroJ.EEG`
- `lag::Int64=1`: lags range is `-lag:lag`
- `demean::Bool=false`: demean signal prior to analysis
- `norm::Bool=false`: normalize autocovariance

# Returns

Named tuple containing:
- `acov::Matrix{Float64}`
- `lags::Vector{Float64}
"""
function eeg_autocov(eeg::NeuroJ.EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    acov, lags = signal_autocov(eeg.eeg_signals, lag=lag, demean=demean, norm=norm)
    lags = (eeg.eeg_time[2] - eeg.eeg_time[1]) .* collect(-lag:lag)

    return (acov=acov, acov_lags=lags)
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

    ccov, lags = signal_crosscov(eeg.eeg_signals, lag=lag, demean=demean, norm=norm)
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

    ccov, lags = signal_crosscov(eeg1.eeg_signals, eeg2.eeg_signals, lag=lag, demean=demean, norm=norm)
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
- `powers::Array{Float64, 3}`
- `frequencies::Array{Float64, 3}`
"""
function eeg_psd(eeg::NeuroJ.EEG; norm::Bool=false)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    psd_pow, psd_frq = signal_psd(eeg.eeg_signals, fs=eeg_sr(eeg), norm=norm)

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

    s_stationarity = signal_stationarity(eeg.eeg_signals, window=window, method=method)

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

    mi = signal_mi(eeg.eeg_signals)
    size(mi, 3) == 1 && (mi = reshape(mi, size(mi, 1), size(mi, 2)))

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

    mi = signal_mi(eeg1.eeg_signals, eeg2.eeg_signals)
    size(mi, 3) == 1 && (mi = reshape(mi, size(mi, 1), size(mi, 2)))

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

    ent = signal_entropy(eeg.eeg_signals)
    size(ent, 3) == 1 && (ent = reshape(ent, size(ent, 1), size(ent, 2)))

    return ent
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

    coherence = signal_coherence(eeg1.eeg_signals, eeg2.eeg_signals)
    size(coherence, 3) == 1 && (coherence = reshape(coherence, size(coherence, 1), size(coherence, 2)))

    return coherence
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

- `coherence::Vector{ComplexF64}`
"""
function eeg_coherence(eeg::NeuroJ.EEG; channel1::Int64, channel2::Int64, epoch1::Int64, epoch2::Int64)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    (channel1 < 0 || channel2 < 0 || epoch1 < 0 || epoch2 < 0) && throw(ArgumentError("channel1/epoch1/channel2/epoch2 must be > 0."))
    channel_n = eeg.eeg_header[:channel_n]
    epoch_n = eeg_epoch_n(eeg)
    (channel1 > channel_n || channel2 > channel_n) && throw(ArgumentError("channel1/channel2 must be ≤ $(channel_n)."))
    (epoch1 > epoch_n || epoch2 > epoch_n) && throw(ArgumentError("epoch1/epoch2 must be ≤ $(epoch_n)."))

    coherence = signal_coherence(eeg.eeg_signals[channel1, :, epoch1], eeg.eeg_signals[channel2, :, epoch2])

    return coherence
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

    hz, nyq = freqs(eeg.eeg_signals[1, :, 1], eeg_sr(eeg))

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
        signals_statistic[epoch, :], signals_statistic_single[epoch], p[epoch] = signal_difference(eeg1.eeg_signals[:, :, epoch], eeg2.eeg_signals[:, :, epoch], n=n, method=method)
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
- `e_mean::Matrix(Float64)`: mean
- `e_median::Matrix(Float64)`: median
- `e_std::Matrix(Float64)`: standard deviation
- `e_var::Matrix(Float64)`: variance
- `e_kurt::Matrix(Float64)`: kurtosis
- `e_mean_diff::Matrix(Float64)`: mean diff value
- `e_median_diff::Matrix(Float64)`: median diff value
- `e_max_dif::Matrix(Float64)`: max difference
- `e_dev_mean::Matrix(Float64)`: deviation from channel mean
"""
function eeg_epochs_stats(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    e_mean, e_median, e_std, e_var, e_kurt, e_mean_diff, e_median_diff, e_dev_mean, e_max_dif = signal_epochs_stats(eeg.eeg_signals)

    return (e_mean=e_mean, e_median=e_median, e_std=e_std, e_var=e_var, e_kurt=e_kurt, e_mean_diff=e_mean_diff, e_median_diff=e_median_diff, e_dev_mean=e_dev_mean, e_max_dif=e_max_dif)
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

    s_pow, s_frq, s_t = signal_spectrogram(eeg.eeg_signals, fs=eeg_sr(eeg), norm=norm, demean=demean)

    return (s_pow=s_pow, s_frq=s_frq, s_t=s_t)
end

"""
    eeg_spectrum(eeg; pad)

Calculate FFT, amplitudes, powers and phases for each channel of the `eeg`. For `pad` > 0 channels are padded with 0s.

# Arguments

- `eeg::NeuroJ.EEG`
- `pad::Int64=0`: number of 0s to pad

# Returns

- `fft::Array{ComplexF64, 3}`
- `amplitudes::Array{Float64, 3}`
- `powers::Array{Float64, 3}`
- `phases::Array{Float64, 3}
"""
function eeg_spectrum(eeg::NeuroJ.EEG; pad::Int64=0)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    s_fft, s_amp, s_pow, s_pha = signal_spectrum(eeg.eeg_signals, pad=pad)

    return s_fft, s_amp, s_pow, s_pha
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

function _signal_channels_stats(signal::AbstractArray)
    c_mean = mean(signal)
    c_median = median(signal)
    c_std = std(signal)
    c_var = var(signal)
    c_kurt = kurtosis(signal)
    c_mean_diff = mean(diff(signal))
    c_median_diff = median(diff(signal))
    c_max_dif = maximum(signal) - minimum(signal)
    c_dev_mean = abs(mean(signal)) - mean(signal)
    
    return c_mean, c_median, c_std, c_var, c_kurt, c_mean_diff, c_median_diff, c_max_dif, c_dev_mean
end

"""
    eeg_channels_stats(eeg)

Calculate `eeg` channels statistics.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

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
    c_dev_mean = zeros(channel_n, epoch_n)
    c_max_dif = zeros(channel_n, epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view eeg.eeg_signals[idx, :, epoch]
            c_mean[idx, epoch], c_median[idx, epoch], c_std[idx, epoch], c_var[idx, epoch], c_kurt[idx, epoch], c_mean_diff[idx, epoch], c_median_diff[idx, epoch], c_dev_mean[idx, epoch], c_max_dif[idx, epoch] = signal_channels_stats(s)
        end
    end

    return c_mean, c_median, c_std, c_var, c_kurt, c_mean_diff, c_median_diff, c_dev_mean, c_max_dif
end

function _signal_snr(signal::AbstractArray)

    # make signal positive
    signal .+= abs(minimum(signal))

    snr = mean(signal) / std(signal)

    return snr
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
            snr[idx, epoch] = _signal_snr(s)
        end
    end

    return snr
end

function _signal_standardize(signal::Array{Float64, 3})

    epoch_n = eeg_epoch_n(eeg)    
    ss = similar(eeg.eeg_signals)
    scaler = Vector{Any}()

    @inbounds @simd for epoch in 1:epoch_n
        s = @view signal[:, :, epoch]
        push!(scaler, StatsBase.fit(ZScoreTransform, s, dims=2)) 
        ss[:,:, epoch] = StatsBase.transform(scaler[epoch], s)
    end

    return ss, scaler
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
    ss, scaler = signal_standardize(eeg.eeg_signals)
    eeg_new = deepcopy(eeg)
    eeg.eeg_signals = ss

    push!(eeg_new.eeg_header[:history], "eeg_standardize!(EEG)")

    return eeg_new, scaler
end

"""
    eeg_fconv(eeg, kernel)

Perform convolution in the time domain.

# Arguments

- `eeg::NeuroJ.EEG`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`: kernel for convolution

# Returns

- `s_convoluted::Union{Array{Float64, 3}, Array{ComplexF64, 3}}`: convoluted signal
"""
function eeg_fconv(eeg::NeuroJ.EEG; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    s_convoluted = signal_fconv(eeg.eeg_signals, kernel=kernel)

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

    s_convoluted = signal_tconv(eeg.eeg_signals, kernel=kernel)

    return s_convoluted
end