"""
    eeg_total_power(eeg)

Calculates total power of the `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `stp::Vector{Float64}`
"""
function eeg_total_power(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    fs = eeg.eeg_header[:sampling_rate][1]
    stp = signal_total_power(eeg.eeg_signals, fs=fs)
    size(stp, 3) == 1 && (stp = reshape(stp, size(stp, 1), size(stp, 2)))

    return stp
    return
end
"""
    eeg_total_power!(eeg)

Calculates total power of the `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
"""
function eeg_total_power!(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    :total_power in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:total_power)
    fs = eeg.eeg_header[:sampling_rate][1]
    stp = signal_total_power(eeg.eeg_signals, fs=fs)
    push!(eeg.eeg_components, stp)
    push!(eeg.eeg_header[:components], :total_power)
    push!(eeg.eeg_header[:history], "eeg_total_power!(EEG)")

end
"""
    eeg_band_power(eeg; f)

Calculates absolute band power between frequencies `f[1]` and `f[2]` of the `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `f::Tuple(Union(Int64, Float64}, Union(Int64, Float64}}`: lower and upper frequency bounds

# Returns

- `sbp::Vector{Float64}`
"""
function eeg_band_power(eeg::NeuroJ.EEG; f::Tuple)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    fs = eeg.eeg_header[:sampling_rate][1]
    sbp = signal_band_power(eeg.eeg_signals, fs=fs, f=f)
    size(sbp, 3) == 1 && (sbp = reshape(sbp, size(sbp, 1), size(sbp, 2)))

    return sbp
end
"""
    eeg_cov(eeg; norm=true)

Calculates covariance between all channels of `eeg`.

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
    eeg_cov!(eeg; norm=true)

Calculates covariance between all channels of `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `norm::Bool`: normalize covariance
"""
function eeg_cov!(eeg::NeuroJ.EEG; norm=true)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    :cov_mat in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:cov_mat)
    cov_mat = eeg_cov(eeg, norm=norm)
    push!(eeg.eeg_components, cov_mat)
    push!(eeg.eeg_header[:components], :cov_mat)
    push!(eeg.eeg_header[:history], "eeg_cov!(EEG)")

    return
end
"""
    eeg_cor(eeg)

Calculates correlation coefficients between all channels of `eeg`.

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
    eeg_cor!(eeg)

Calculates correlation coefficients between all channels of `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
"""
function eeg_cor!(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    :cor_mat in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:cor_mat)
    cor_mat = eeg_cor(eeg)
    push!(eeg.eeg_components, cor_mat)
    push!(eeg.eeg_header[:components], :cor_mat)
    push!(eeg.eeg_header[:history], "eeg_cor!(EEG)")

    return
end
"""
    eeg_autocov(eeg; lag=1, demean=false, norm=false)

Calculates autocovariance of each the `eeg` channels.

# Arguments

- `eeg::NeuroJ.EEG`
- `lag::Int64`: lags range is `-lag:lag`
- `demean::Bool`: demean signal prior to analysis
- `norm::Bool`: normalize autocovariance

# Returns

- `acov::Matrix{Float64}`
- `lags::Vector{Float64}
"""
function eeg_autocov(eeg::NeuroJ.EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    acov, lags = signal_autocov(eeg.eeg_signals, lag=lag, demean=demean, norm=norm)
    size(acov, 3) == 1 && (acov = reshape(acov, size(acov, 1), size(acov, 2)))
    lags = (eeg.eeg_time[2] - eeg.eeg_time[1]) .* collect(-lag:lag)

    return acov, lags
end
"""
    eeg_autocov!(eeg; lag=1, demean=false, norm=false)

Calculates autocovariance of each the `eeg` channels.

# Arguments

- `eeg::NeuroJ.EEG`
- `lag::Int64`: lags range is `-lag:lag`
- `demean::Bool`: demean signal prior to analysis
- `norm::Bool`: normalize autocovariance
"""
function eeg_autocov!(eeg::NeuroJ.EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    :acov in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:acov)
    :acov_lags in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:acov_lags)
    acov, lags = eeg_autocov(eeg, lag=lag, demean=demean, norm=norm)
    push!(eeg.eeg_components, acov)
    push!(eeg.eeg_components, lags)
    push!(eeg.eeg_header[:components], :acov)
    push!(eeg.eeg_header[:components], :acov_lags)
    push!(eeg.eeg_header[:history], "eeg_autocov!(EEG, lag=$lag, demean=$demean, norm=$norm)")

    return
end
"""
    eeg_crosscov(eeg; lag=1, demean=false, norm=false)

Calculates cross-covariance of each the `eeg` channels.

# Arguments

- `eeg::NeuroJ.EEG`
- `lag::Int64`: lags range is `-lag:lag`
- `demean::Bool`: demean signal prior to analysis
- `norm::Bool`: normalize cross-covariance

# Returns

- `ccov::Matrix{Float64}`
- `lags::Vector{Float64}
"""
function eeg_crosscov(eeg::NeuroJ.EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    ccov, lags = signal_crosscov(eeg.eeg_signals, lag=lag, demean=demean, norm=norm)
    size(ccov, 3) == 1 && (ccov = reshape(ccov, size(ccov, 1), size(ccov, 2)))
    lags = (eeg.eeg_time[2] - eeg.eeg_time[1]) .* collect(-lag:lag)

    return ccov, lags
end
"""
    eeg_crosscov!(eeg; lag=1, demean=false, norm=false)

Calculates cross-covariance of each the `eeg` channels.

# Arguments

- `eeg::NeuroJ.EEG`
- `lag::Int64`: lags range is `-lag:lag`
- `demean::Bool`: demean signal prior to analysis
- `norm::Bool`: normalize cross-covariance
"""
function eeg_crosscov!(eeg::NeuroJ.EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    :ccov in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:ccov)
    :ccov_lags in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:ccov_lags)
    ccov, lags = eeg_crosscov(eeg, lag=lag, demean=demean, norm=norm)
    push!(eeg.eeg_components, ccov)
    push!(eeg.eeg_components, lags)
    push!(eeg.eeg_header[:components], :ccov)
    push!(eeg.eeg_header[:components], :ccov_lags)
    push!(eeg.eeg_header[:history], "eeg_crosscov!(EEG, lag=$lag, demean=$demean, norm=$norm)")

    return
end
"""
    eeg_crosscov(eeg1, eeg2; lag=1, demean=false, norm=false)

Calculates cross-covariance between same channels in `eeg1` and `eeg2`.

# Arguments

- `eeg1::NeuroJ.EEG`
- `eeg2::NeuroJ.EEG`
- `lag::Int64`: lags range is `-lag:lag`
- `demean::Bool`: demean signal prior to analysis
- `norm::Bool`: normalize cross-covariance

# Returns

- `ccov::Matrix{Float64}`
- `lags::Vector{Float64}
"""
function eeg_crosscov(eeg1::NeuroJ.EEG, eeg2::NeuroJ.EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)

    eeg_channel_n(eeg1, type=:eeg) < eeg_channel_n(eeg1, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))
    eeg_channel_n(eeg2, type=:eeg) < eeg_channel_n(eeg2, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    ccov, lags = signal_crosscov(eeg1.eeg_signals, eeg2.eeg_signals, lag=lag, demean=demean, norm=norm)
    size(ccov, 3) == 1 && (ccov = reshape(ccov, size(ccov, 1), size(ccov, 2)))
    lags = (eeg1.eeg_time[2] - eeg1.eeg_time[1]) .* collect(-lag:lag)

    return ccov, lags
end
"""
    eeg_psd(eeg; norm=false)

Calculates total power for each the `eeg` channels.

# Arguments

- `eeg::NeuroJ.EEG`
- `norm::Bool`: normalize do dB

# Returns

- `powers::Array{Float64, 3}`
- `frequencies::Array{Float64, 3}`
"""
function eeg_psd(eeg::NeuroJ.EEG; norm::Bool=false)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    s_psd_powers, s_psd_frequencies = signal_psd(eeg.eeg_signals, fs=eeg_sr(eeg), norm=norm)
    size(s_psd_powers, 3) == 1 && (s_psd_powers = reshape(s_psd_powers, size(s_psd_powers, 1), size(s_psd_powers, 2)))
    size(s_psd_frequencies, 3) == 1 && (s_psd_frequencies = reshape(s_psd_frequencies, size(s_psd_frequencies, 1), size(s_psd_frequencies, 2)))

    return s_psd_powers, s_psd_frequencies
end
"""
    eeg_psd!(eeg; norm=false)

Calculates total power for each the `eeg` channels.

# Arguments

- `eeg::NeuroJ.EEG`
- `norm::Bool`: normalize do dB
"""
function eeg_psd!(eeg::NeuroJ.EEG; norm::Bool=false)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    :psd_p in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:psd_p)
    :psd_f in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:psd_f)
    s_psd_powers, s_psd_frequencies = eeg_psd(eeg, norm=norm)
    push!(eeg.eeg_components, s_psd_powers)
    push!(eeg.eeg_components, s_psd_frequencies)
    push!(eeg.eeg_header[:components], :psd_p)
    push!(eeg.eeg_header[:components], :psd_f)
    push!(eeg.eeg_header[:history], "eeg_psd!(EEG, norm=$norm)")

    return
end
"""
    eeg_stationarity(eeg:EEG; window=10, method=:euclid)

Calculates stationarity.

# Arguments

- `eeg:EEG`
- `window::Int64`: time window in samples
- `method::Symbol[:mean, :var, :euclid, :hilbert]

# Returns

- `stationarity::Union{Matrix{Float64}, Array{Float64, 3}}`
"""
function eeg_stationarity(eeg::NeuroJ.EEG; window::Int64=10, method::Symbol=:hilbert)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    s_stationarity = signal_stationarity(eeg.eeg_signals, window=window, method=method)

    return s_stationarity
end
"""
    eeg_stationarity!(eeg:EEG; window=10, method=:euclid)

Calculates stationarity.

# Arguments

- `eeg:EEG`
- `window::Int64`: time window in samples
- `method::Symbol[:mean, :var, :euclid, :hilbert]
"""
function eeg_stationarity!(eeg::NeuroJ.EEG; window::Int64=10, method::Symbol=:hilbert)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    :stationarity in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:stationarity)
    s_stationarity = eeg_stationarity(eeg, window=window, method=method)
    push!(eeg.eeg_components, s_stationarity)
    push!(eeg.eeg_header[:components], :stationarity)
    push!(eeg.eeg_header[:history], "eeg_stationarity!(EEG, window=$window, method=$method)")

    return
end
"""
    eeg_mi(eeg)

Calculates mutual information between all channels of `eeg`.

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
    eeg_mi!(eeg)

Calculates mutual information between all channels of `eeg` and stores into :mi component.

# Arguments

- `eeg::NeuroJ.EEG`
"""
function eeg_mi!(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    :mi in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:mi)
    mi = signal_mi(eeg.eeg_signals)
    size(mi, 3) == 1 && (mi = reshape(mi, size(mi, 1), size(mi, 2)))
    push!(eeg.eeg_components, mi)
    push!(eeg.eeg_header[:components], :mi)
    push!(eeg.eeg_header[:history], "eeg_mi!(EEG)")

    return
end
"""
    eeg_mi(eeg1, eeg2)

Calculates mutual information between all channels of `eeg1` and `eeg2`.

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

Calculates entropy of all channels of `eeg`.

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
    eeg_entropy!(eeg)

Calculates entropy of all channels of `eeg` and stores into :entropy component.

# Arguments

- `eeg::NeuroJ.EEG`
"""
function eeg_entropy!(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    :entropy in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:entropy)
    ent = eeg_entropy(eeg)
    push!(eeg.eeg_components, ent)
    push!(eeg.eeg_header[:components], :entropy)
    push!(eeg.eeg_header[:history], "eeg_entropy!(EEG)")

    return
end
"""
    eeg_band(eeg, band)

Return `eeg` :band frequency limits.

# Arguments

- `eeg:EEG`
- `band::Symbol`

# Returns

- `band_frequency::Tuple{Float64, Float64}`
"""
function eeg_band(eeg; band::Symbol)

    band in [:delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher] || throw(ArgumentError("band must be: :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower or :gamma_higher."))

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

Calculates coherence between all channels of `eeg1` and `eeg2`.

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

Calculates coherence between `channel1`/`epoch1` and `channel2` of `epoch2` of `eeg`.

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

Returns vector of frequencies and Nyquist frequency for `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `hz::Vector{Float64}`
- `nyquist::Float64`
"""
function eeg_freqs(eeg::NeuroJ.EEG)

    hz, nyq = freqs(eeg.eeg_signals[1, :, 1], eeg_sr(eeg))

    return hz, nyq
end
"""
    eeg_freqs!(eeg)

Returns vector of frequencies and Nyquist frequency for `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
"""
function eeg_freqs!(eeg::NeuroJ.EEG)

    :hz in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:hz)
    :nyq in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:nyq)
    hz, nyq = eeg_freqs(eeg)
    push!(eeg.eeg_components, hz)
    push!(eeg.eeg_components, nyq)
    push!(eeg.eeg_header[:components], :hz)
    push!(eeg.eeg_header[:components], :nyq)
    push!(eeg.eeg_header[:history], "eeg_freqs!(EEG)")

    return
end
"""
    eeg_difference(eeg1, eeg2; n=3, method=:absdiff)

Calculates mean difference and 95% confidence interval for `eeg1` and `eeg2`.

# Arguments

- `eeg1::NeuroJ.EEG`
- `eeg2::NeuroJ.EEG`
- `n::Int64`: number of bootstraps
- `method::Symbol[:absdiff, :diff2int]`
    - `:absdiff`: maximum difference
    - `:diff2int`: integrated area of the squared difference

# Returns

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

    return signals_statistic, signals_statistic_single, p
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

Calculates mean, median, sd, kurtosis and variance of each `eeg` epoch.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `mean::Vector{Float64}`
- `median::Vector{Float64}`
- `sd::Vector{Float64}`
- `var::Vector{Float64}`
- `kurtosis::Vector{Float64}`
"""
function eeg_epochs_stats(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    e_mean, e_median, e_sd, e_var, e_kurt = signal_epochs_stats(eeg.eeg_signals)

    return e_mean, e_median, e_sd, e_var, e_kurt
end
"""
    eeg_epochs_stats!(eeg)

Calculates mean, sd and variance of `eeg` epochs and stores in `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
"""
function eeg_epochs_stats!(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    :epochs_mean in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:epochs_mean)
    :epochs_median in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:epochs_median)
    :epochs_sd in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:epochs_sd)
    :epochs_var in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:epochs_var)
    :epochs_kurtosis in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:epochs_kurtosis)
    e_mean, e_median, e_sd, e_var, e_kurt = signal_epochs_stats(eeg.eeg_signals)
    push!(eeg.eeg_components, e_mean)
    push!(eeg.eeg_components, e_median)
    push!(eeg.eeg_components, e_sd)
    push!(eeg.eeg_components, e_var)
    push!(eeg.eeg_components, e_kurt)
    push!(eeg.eeg_header[:components], :epochs_mean)
    push!(eeg.eeg_header[:components], :epochs_median)
    push!(eeg.eeg_header[:components], :epochs_sd)
    push!(eeg.eeg_header[:components], :epochs_var)
    push!(eeg.eeg_header[:components], :epochs_kurtosis)
    push!(eeg.eeg_header[:history], "eeg_epochs_stats!(EEG)")

    return
end
"""
    eeg_spectrogram(eeg; norm=true, demean=true)

Calculates spectrogram of `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `norm::Bool`: normalize powers to dB
- `demean::Bool`: demean signal prior to analysis

# Returns

- `spec.power::Array{Float64, 3}`
- `spec.freq::Matrix{Float64}`
- `spec.time::Matrix{Float64}`
"""
function eeg_spectrogram(eeg::NeuroJ.EEG; norm::Bool=true, demean::Bool=true)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    s_pow, s_frq, s_t = signal_spectrogram(eeg.eeg_signals, fs=eeg_sr(eeg), norm=norm, demean=demean)

    return s_pow, s_frq, s_t
end

"""
    eeg_spectrogram!(eeg, norm=true, demean=true)

Calculates spectrogram of `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `norm::Bool`: normalize powers to dB
- `demean::Bool`: demean signal prior to analysis
"""
function eeg_spectrogram!(eeg::NeuroJ.EEG; norm::Bool=true, demean::Bool=true)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    :spectrogram_pow in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:spec_pow)
    :spectrogram_frq in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:spec_frq)
    :spectrogram_t in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:spec_t)
    s_pow, s_frq, s_t = signal_spectrogram(eeg.eeg_signals, fs=eeg_sr(eeg), norm=norm, demean=demean)
    push!(eeg.eeg_components, s_pow)
    push!(eeg.eeg_components, s_frq)
    push!(eeg.eeg_components, s_t)
    push!(eeg.eeg_header[:components], :spec_pow)
    push!(eeg.eeg_header[:components], :spec_frq)
    push!(eeg.eeg_header[:components], :spec_t)
    push!(eeg.eeg_header[:history], "eeg_spectrogram!(EEG, norm=$norm, demean=$demean)")

    return
end
"""
    eeg_spectrum(eeg; pad=0)

Calculates FFT, amplitudes, powers and phases for each channel of the `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`: the signal
- `pad::Int64`: pad channels `pad` zeros

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
    eeg_spectrum!(eeg; pad=0)

Calculates FFT, amplitudes, powers and phases for each channel of the `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`: the signal
- `pad::Int64`: pad channels `pad` zeros
"""
function eeg_spectrum!(eeg::NeuroJ.EEG; pad::Int64=0)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    :spectrum_fft in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:spectrum_fft)
    :spectrum_amp in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:spectrum_amp)
    :spectrum_pow in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:spectrum_pow)
    :spectrum_phase in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:spectrum_phase)
    s_fft, s_amplitudes, s_powers, s_phases = signal_spectrum(eeg.eeg_signals)
    push!(eeg.eeg_components, s_fft)
    push!(eeg.eeg_components, s_amplitudes)
    push!(eeg.eeg_components, s_powers)
    push!(eeg.eeg_components, s_phases)
    push!(eeg.eeg_header[:components], :spectrum_fft)
    push!(eeg.eeg_header[:components], :spectrum_amp)
    push!(eeg.eeg_header[:components], :spectrum_pow)
    push!(eeg.eeg_header[:components], :spectrum_phase)
    push!(eeg.eeg_header[:history], "eeg_spectrum!(EEG, pad=$pad)")

    return
end
"""
    eeg_s2t(eeg; t)

Converts time `t` in samples to seconds using `eeg` sampling rate.

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

Converts time `t` in seconds to samples using `eeg` sampling rate.

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
