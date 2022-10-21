"""
    eeg_total_power(eeg, mt)

Calculate total power of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns
 
- `stp::Matrix{Float64}`: total power for each channel per epoch
"""
function eeg_total_power(eeg::NeuroAnalyzer.EEG, mt::Bool=false)

    fs = eeg_sr(eeg)
    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = length(channels)
    epoch_n = size(signal, 3)

    stp = zeros(channel_n, epoch_n)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            @views stp[channel_idx, epoch_idx] = s_total_power(signal[channel_idx, :, epoch_idx], fs=fs, mt=mt)
        end
    end

    return stp
end

"""
    eeg_band_power(eeg; f, mt)

Calculate absolute band power between frequencies `f[1]` and `f[2]` of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

- `sbp::Matrix{Float64}`: band power for each channel per epoch
"""
function eeg_band_power(eeg::NeuroAnalyzer.EEG; f::Tuple{Real, Real}, mt::Bool=false)

    fs = eeg_sr(eeg)
    length(f) != 2 && throw(ArgumentError("f must contain two frequencies."))
    f = tuple_order(f)
    f[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
    f[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be < $(fs / 2)."))

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = length(channels)
    epoch_n = size(signal, 3)

    sbp = zeros(channel_n, epoch_n)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            @views sbp[channel_idx, epoch_idx] = s_band_power(signal[channel_idx, :, epoch_idx], fs=fs, f=f, mt=mt)
        end
    end

    return sbp
end

"""
    eeg_cov(eeg; norm)

Calculate covariance matrix for all EEG/MEG channels of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `norm::Bool=true`: normalize matrix

# Returns

- `cov_mat::Array{Float64, 3}`: covariance matrix for each epoch
"""
function eeg_cov(eeg::NeuroAnalyzer.EEG; norm::Bool=true)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    cov_mat = zeros(channel_n, channel_n, epoch_n)
    Threads.@threads for epoch_idx in 1:epoch_n
        @views cov_mat[:, :, epoch_idx] = cov(signal[:, :, epoch_idx]')
    end

    # normalize covariance matrix
    norm == true && (cov_mat = m_norm(cov_mat))

    return cov_mat
end

"""
    eeg_cor(eeg; norm)

Calculate correlation coefficients between all channels of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `norm::Bool=true`: normalize matrix

# Returns

- `cov_mat::Array{Float64, 3}`: correlation matrix for each epoch
"""
function eeg_cor(eeg::NeuroAnalyzer.EEG; norm::Bool=true)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    cor_mat = zeros(channel_n, channel_n, epoch_n)
    Threads.@threads for epoch_idx in 1:epoch_n
        cor_mat[:, :, epoch_idx] = @views cor(signal[:, :, epoch_idx]')
    end

    # normalize covariance matrix
    norm == true && (cor_mat = m_norm(cor_mat))

    return cor_mat
end

"""
    eeg_xcov(eeg; lag, demean, norm)

Calculate cross-covariance for all `eeg` channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `lag::Int64=1`: lags range is `-lag:lag`
- `demean::Bool=false`: demean signal prior to analysis
- `norm::Bool=false`: normalize cross-covariance

# Returns

Named tuple containing:
- `xcov::Matrix{Float64}`
- `lags::Vector{Float64}`
"""
function eeg_xcov(eeg::NeuroAnalyzer.EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    lags = (eeg.eeg_time[2] - eeg.eeg_time[1]) .* collect(-lag:lag)
    xcov = zeros(channel_n^2, length(lags), epoch_n)
    @inbounds @simd for epoch_idx in 1:epoch_n
        xcov_packed = Array{Vector{Float64}}(undef, channel_n, channel_n)
        Threads.@threads for channel_idx1 in 1:channel_n
            for channel_idx2 in 1:channel_idx1
                xcov_packed[channel_idx1, channel_idx2], _ = @views s2_xcov(signal[channel_idx1, :, epoch_idx], signal[channel_idx2, :, epoch_idx], lag=lag, demean=demean, norm=norm)
            end
        end
        Threads.@threads for channel_idx1 in 1:(channel_n - 1)
            for channel_idx2 in (channel_idx1 + 1):channel_n
                xcov_packed[channel_idx1, channel_idx2] = @views xcov_packed[channel_idx2, channel_idx1]
            end
        end
        Threads.@threads for channel_idx in 1:channel_n^2
            xcov[channel_idx, :, epoch_idx] = @views xcov_packed[channel_idx]
        end
    end

    return (xcov=xcov, lags=lags)
end

"""
    eeg_xcov(eeg1, eeg2; channel1, channel2, epoch1, epoch2, lag, demean, norm)

Calculate cross-covariance between `eeg1` and `eeg2`.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2::NeuroAnalyzer.EEG`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=0`: default use all channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=0`: default use all channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=0`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=0`: default use all epochs
- `lag::Int64=1`: lags range is `-lag:lag`
- `demean::Bool=false`: demean signal prior to analysis
- `norm::Bool=false`: normalize cross-covariance

# Returns

Named tuple containing:
- `xcov::Array{Float64, 3}`
- `lags::Vector{Float64}`
"""
function eeg_xcov(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Union{Int64, Vector{Int64}, AbstractRange}=0, channel2::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch1::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch2::Union{Int64, Vector{Int64}, AbstractRange}=0, lag::Int64=1, demean::Bool=false, norm::Bool=false)

    # select channels, default is all
    channels1 = eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))
    channels2 = eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))
    channel1 == 0 && (channel1 = channels1)
    _check_channels(channels1, channel1)
    channel2 == 0 && (channel2 = channels2)
    _check_channels(channels2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    # select epochs, default is all
    epoch1 == 0 && (epoch1 = _select_epochs(eeg1, epoch1))
    _check_epochs(eeg1, epoch1)
    epoch2 == 0 && (epoch2 = _select_epochs(eeg2, epoch2))
    _check_epochs(eeg2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    eeg_epoch_len(eeg1) == eeg_epoch_len(eeg2) || throw(ArgumentError("eeg1 and eeg2 epoch lengths must be equal."))

    signal1 = @view eeg1.eeg_signals[channel1, :, epoch1]
    signal2 = @view eeg2.eeg_signals[channel2, :, epoch2]

    lags = (eeg1.eeg_time[2] - eeg1.eeg_time[1]) .* collect(-lag:lag)
    xcov = zeros(length(channel1), (2 * lag + 1), length(epoch1))
    @inbounds @simd for epoch_idx in 1:length(epoch1)
        Threads.@threads for channel_idx in 1:length(channel1)
            xcov[channel_idx, :, epoch_idx], _ = @views s2_xcov(signal1[channel_idx, :, epoch_idx],
                                                                signal2[channel_idx, :, epoch_idx],
                                                                lag=lag,
                                                                demean=demean,
                                                                norm=norm)
        end
    end

    return (xcov=xcov, lags=lags)
end

"""
    eeg_psd(eeg; norm, mt)

Calculate power spectrum density for each `eeg` channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

Named tuple containing:
- `psd_pow::Array{Float64, 3}`:powers
- `psd_frq::Vector{Float64}`: frequencies
"""
function eeg_psd(eeg::NeuroAnalyzer.EEG; norm::Bool=false, mt::Bool=false)

    fs = eeg_sr(eeg)
    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)
    
    _, psd_frq = s_psd(signal[1, :, 1], fs=fs, norm=norm, mt=mt)
    psd_pow = zeros(channel_n, length(psd_frq), epoch_n)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            psd_pow[channel_idx, :, epoch_idx], _ = s_psd(signal[channel_idx, :, epoch_idx], fs=fs, norm=norm, mt=mt)
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
- `method::Symbol=:euclid`: stationarity method: :mean, :var, :euclid, :hilbert, :adf

# Returns

- `stationarity::Union{Matrix{Float64}, Array{Float64, 3}}}`
"""
function eeg_stationarity(eeg::NeuroAnalyzer.EEG; window::Int64=10, method::Symbol=:hilbert)

    method in [:mean, :var, :euclid, :hilbert, :adf] || throw(ArgumentError("Method must be must be :mean, :var, :euclid, :hilbert or :adf."))

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    (typeof(window) == Int64 && window < 1) && throw(ArgumentError("window must be ≥ 1."))
    (typeof(window) == Int64 && window > size(signal, 2)) && throw(ArgumentError("window must be ≤ $(size(signal, 2))."))

    if method === :mean
        s_stationarity = zeros(channel_n, window, epoch_n)
        @inbounds @simd for epoch_idx in 1:epoch_n
            Threads.@threads for channel_idx in 1:channel_n
                s_stationarity[channel_idx, :, epoch_idx] = @views s_stationarity_mean(signal[channel_idx, :, epoch_idx], window=window)
            end
        end
    end

    if method === :var
        s_stationarity = zeros(channel_n, window, epoch_n)
        @inbounds @simd for epoch_idx in 1:epoch_n
            Threads.@threads for channel_idx in 1:channel_n
                s_stationarity[channel_idx, :, epoch_idx] = @views s_stationarity_var(signal[channel_idx, :, epoch_idx], window=window)
            end
        end
    end

    if method === :hilbert
        s_stationarity = zeros(channel_n, eeg_epoch_len(eeg) - 1, epoch_n)
        @inbounds @simd for epoch_idx in 1:epoch_n
            Threads.@threads for channel_idx in 1:channel_n
                s_stationarity[channel_idx, :, epoch_idx] = @views s_stationarity_hilbert(signal[channel_idx, :, epoch_idx])
            end
        end
    end

    if method === :euclid
        # number of time windows per epoch
        window_n = eeg_epoch_len(eeg)
        cov_mat = zeros(channel_n, channel_n, window_n, epoch_n)
        s_stationarity = zeros(1 + length(2:window:window_n), epoch_n)

        @inbounds @simd for epoch_idx in 1:epoch_n
            Threads.@threads for window_idx = 1:window_n
                cov_mat[:, :, window_idx, epoch_idx] = @views s2_cov(signal[:, window_idx, epoch_idx], signal[:, window_idx, epoch_idx])
            end
        end

        @inbounds @simd for epoch_idx in 1:epoch_n
            phase_idx = 1
            Threads.@threads for window_idx = 2:window:window_n
                s_stationarity[phase_idx, epoch_idx] = @views euclidean(cov_mat[:, :, window_idx - 1, epoch_idx], cov_mat[:, :, window_idx, epoch_idx])
                phase_idx += 1
            end
        end
    end

    if method === :adf
        s_stationarity = zeros(channel_n, 2, epoch_n)
        @inbounds @simd for epoch_idx in 1:epoch_n
            Threads.@threads for channel_idx = 1:channel_n
                adf = @views ADFTest(signal[channel_idx, :, epoch_idx], :constant, window)
                a = adf.stat
                p = pvalue(adf)
                p < eps() && (p = 0.0001)
                a = round(a, digits=2)
                p = round(p, digits=4)
                p == 0.0 && (p = 0.0001)
                s_stationarity[channel_idx, :, epoch_idx] = [a, p]
            end
        end
    end

    return s_stationarity
end

"""
    eeg_mi(eeg)

Calculate mutual information between all channels of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `mi::Array{Float64, 3}`
"""
function eeg_mi(eeg::NeuroAnalyzer.EEG)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    mi = zeros(channel_n, channel_n, epoch_n)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx1 in 1:channel_n
            for channel_idx2 in 1:channel_idx1
                mi[channel_idx1, channel_idx2, epoch_idx] = @views s2_mi(signal[channel_idx1, :, epoch_idx], signal[channel_idx2, :, epoch_idx])
            end
        end
        Threads.@threads for channel_idx1 in 1:(channel_n - 1)
            for channel_idx2 in (channel_idx1 + 1):channel_n
                mi[channel_idx1, channel_idx2, epoch_idx] = @views mi[channel_idx2, channel_idx1, epoch_idx]
            end
        end
    end

    return mi
end

"""
    eeg_mi(eeg1, eeg2)

Calculate mutual information between all channels of `eeg1` and `eeg2`.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2::NeuroAnalyzer.EEG`

# Returns

- `mi::Array{Float64, 3}`
"""
function eeg_mi(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG)

    eeg_channel_n(eeg1) == eeg_channel_n(eeg2) || throw(ArgumentError("Both EEG objects must have the same number of channels."))
    eeg_epoch_n(eeg1) == eeg_epoch_n(eeg2) || throw(ArgumentError("Both EEG objects must have the same number of epochs."))

    channels1 = eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))
    channels2 = eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))
    signal1 = @view eeg1.eeg_signals[channels1, :, :]
    signal2 = @view eeg2.eeg_signals[channels2, :, :]
    channel_n = size(signal1, 1)
    epoch_n = size(signal1, 3)

    mi = zeros(channel_n, channel_n, epoch_n)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx1 in 1:channel_n
            for channel_idx2 in 1:channel_n
                mi[channel_idx1, channel_idx2, epoch_idx] = @views s2_mi(signal1[channel_idx1, :, epoch_idx], signal2[channel_idx2, :, epoch_idx])
            end
        end
    end

    return mi
end

"""
    eeg_entropy(eeg)

Calculate entropy of all channels of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

Named tuple containing:
- `ent::Array{Float64, 2}`
- `sent::Array{Float64, 2}`: Shanon entropy
- `leent::Array{Float64, 2}`: log energy entropy
"""
function eeg_entropy(eeg::NeuroAnalyzer.EEG)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    ent = zeros(channel_n, epoch_n)
    sent = zeros(channel_n, epoch_n)
    leent = zeros(channel_n, epoch_n)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            ent[channel_idx, epoch_idx], sent[channel_idx, epoch_idx], leent[channel_idx, epoch_idx] = @views s_entropy(signal[channel_idx, :, epoch_idx])
        end
    end

    return (ent=ent, sent=sent, leent=leent)
end

"""
    eeg_negentropy(eeg)

Calculate negentropy of all channels of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `ne::Matrix{Float64}`
"""
function eeg_negentropy(eeg::NeuroAnalyzer.EEG)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    ne = zeros(channel_n, epoch_n)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            ne[channel_idx, epoch_idx] = @views s_negentropy(signal[channel_idx, :, epoch_idx])
        end
    end

    return ne
end

"""
    eeg_band(eeg, band)

Return frequency limits for a `band` range.

# Arguments

- `eeg:EEG`
- `band::Symbol`: name of band range: :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher. If lower or upper band frequency limit exceeds Nyquist frequency of `eeg`, than bound is truncated to `eeg` range.

# Returns

- `band_frequency::Tuple{Real, Real}`
"""
function eeg_band(eeg::NeuroAnalyzer.EEG; band::Symbol)

    band in [:total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher] || throw(ArgumentError("band must be: :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower or :gamma_higher."))

    band === :total && (band_frequency = (0.1, round(eeg_sr(eeg) / 2, digits=1)))
    band === :delta && (band_frequency = (0.1, 4.0))
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
    band_frequency[2] > eeg_sr(eeg) / 2 && (band_frequency = (band_frequency[1], (eeg_sr(eeg) / 2) - 0.1))

    return band_frequency
end

"""
    eeg_tcoherence(eeg1, eeg2; channel1, channel2, epoch1, epoch2)

Calculate coherence (mean over time) and MSC (magnitude-squared coherence) between `eeg1` and `eeg2`.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2::NeuroAnalyzer.EEG`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=0`: default use all channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=0`: default use all channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=0`: default use all channels
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=0`: default use all channels

# Returns

Named tuple containing:
- `c::Array{Float64, 3}`: coherence
- `msc::Array{Float64, 3}`: MSC
- `ic::Array{Float64, 3}`: imaginary part of coherence
"""
function eeg_tcoherence(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Union{Int64, Vector{Int64}, AbstractRange}=0, channel2::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch1::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch2::Union{Int64, Vector{Int64}, AbstractRange}=0)

    # select channels, default is all
    channels1 = eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))
    channels2 = eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))
    channel1 == 0 && (channel1 = channels1)
    _check_channels(channels1, channel1)
    channel2 == 0 && (channel2 = channels2)
    _check_channels(channels2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))

    # select epochs, default is all
    epoch1 == 0 && (epoch1 = _select_epochs(eeg1, epoch1))
    _check_epochs(eeg1, epoch1)
    epoch2 == 0 && (epoch2 = _select_epochs(eeg2, epoch2))
    _check_epochs(eeg2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    eeg_epoch_len(eeg1) == eeg_epoch_len(eeg2) || throw(ArgumentError("eeg1 and eeg2 epoch lengths must be equal."))

    c = zeros(length(channel1), eeg_epoch_len(eeg1), length(epoch1))
    msc = zeros(length(channel1), eeg_epoch_len(eeg1), length(epoch1))
    ic = zeros(length(channel1), eeg_epoch_len(eeg1), length(epoch1))

    @inbounds @simd for epoch_idx in 1:length(epoch1)
        Threads.@threads for channel_idx in 1:length(channel1)
            c[channel_idx, :, epoch_idx], msc[channel_idx, :, epoch_idx], ic[channel_idx, :, epoch_idx] = @views s2_tcoherence(eeg1.eeg_signals[channel1[channel_idx], :, epoch1[epoch_idx]], eeg2.eeg_signals[channel2[channel_idx], :, epoch2[epoch_idx]])
        end
    end

    return (c=c, msc=msc, ic=ic)
end

"""
    eeg_freqs(eeg)

Return vector of frequencies and Nyquist frequency for `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

Named tuple containing:
- `hz::Vector{Float64}`
- `nyquist::Float64`
"""
function eeg_freqs(eeg::NeuroAnalyzer.EEG)

    hz, nyq = s_freqs(eeg.eeg_signals[1, :, 1], eeg_sr(eeg))

    return (hz=hz, nyquist=nyq)
end

"""
    eeg_difference(eeg1, eeg2; n, method)

Calculate mean difference and its 95% CI between `eeg1` and `eeg2`.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2::NeuroAnalyzer.EEG`
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
function eeg_difference(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; n::Int64=3, method::Symbol=:absdiff)

    channels1 = eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))
    channels2 = eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))
    signal1 = @view eeg1.eeg_signals[channels1, :, :]
    signal2 = @view eeg2.eeg_signals[channels2, :, :]
    epoch_n = size(signal1, 3)

    s_stat = zeros(epoch_n, size(signal1, 1) * n)
    s_stat_single = zeros(epoch_n)
    p = zeros(epoch_n)

    @inbounds @simd for epoch_idx in 1:epoch_n
        s_stat[epoch_idx, :], s_stat_single[epoch_idx], p[epoch_idx] = s2_difference(signal1[:, :, epoch_idx], signal2[:, :, epoch_idx], n=n, method=method)
    end

    return (s_stat=s_stat, s_stat_single=s_stat_single, p=p)
end

"""
    eeg_picks(eeg; pick)

Return `pick` of electrodes for `eeg` electrodes.

# Arguments

- `pick::Vector{Symbol}`

# Returns

- `channels::Vector{Int64}`
"""
function eeg_pick(eeg::NeuroAnalyzer.EEG; pick::Union{Symbol, Vector{Symbol}})

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

- `eeg::NeuroAnalyzer.EEG`

# Returns

Named tuple containing:
- `e_mean::Vector(Float64)`: mean
- `e_median::Vector(Float64)`: median
- `e_std::Vector(Float64)`: standard deviation
- `e_var::Vector(Float64)`: variance
- `e_kurt::Vector(Float64)`: kurtosis
- `e_skew::Vector(Float64)`: skewness
- `e_mean_diff::Vector(Float64)`: mean diff value
- `e_median_diff::Vector(Float64)`: median diff value
- `e_max_dif::Vector(Float64)`: max difference
- `e_dev_mean::Vector(Float64)`: deviation from channel mean
"""
function eeg_epochs_stats(eeg::NeuroAnalyzer.EEG)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    epoch_n = size(signal, 3)

    e_mean = zeros(epoch_n)
    e_median = zeros(epoch_n)
    e_std = zeros(epoch_n)
    e_var = zeros(epoch_n)
    e_kurt = zeros(epoch_n)
    e_skew = zeros(epoch_n)
    e_mean_diff = zeros(epoch_n)
    e_median_diff = zeros(epoch_n)
    e_max_dif = zeros(epoch_n)
    e_dev_mean = zeros(epoch_n)

    @inbounds @simd for epoch_idx in 1:epoch_n
        e_mean[epoch_idx] = @views mean(signal[:, :, epoch_idx])
        e_median[epoch_idx] = @views median(signal[:, :, epoch_idx])
        e_std[epoch_idx] = @views std(signal[:, :, epoch_idx])
        e_var[epoch_idx] = @views var(signal[:, :, epoch_idx])
        e_kurt[epoch_idx] = @views kurtosis(signal[:, :, epoch_idx])
        e_skew[epoch_idx] = @views skewness(signal[:, :, epoch_idx])
        e_mean_diff = @views mean(diff(signal[:, :, epoch_idx], dims=2))
        e_median_diff = @views median(diff(signal[:, :, epoch_idx], dims=2))
        e_max_dif = @views maximum(signal[:, :, epoch_idx]) - minimum(signal[:, :, epoch_idx])
        e_dev_mean = @views abs(mean(signal[:, :, epoch_idx])) - mean(signal[:, :, epoch_idx])
    end

    return (e_mean=e_mean, e_median=e_median, e_std=e_std, e_var=e_var, e_kurt=e_kurt, e_skew=e_skew, e_mean_diff=e_mean_diff, e_median_diff=e_median_diff, e_max_dif=e_max_dif, e_dev_mean=e_dev_mean)
end

"""
    eeg_spectrogram(eeg; norm, mt, st, demean)

Return spectrogram of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `norm::Bool=true`: normalize powers to dB
- `mt::Bool=false`: if true use multi-tapered spectrogram
- `st::Bool=false`: if true use short time Fourier transform
- `demean::Bool=true`: demean signal prior to analysis

# Returns

Named tuple containing:
- `s_pow::Array{Float64, 3}`
- `s_frq::Vector{Float64}`
- `s_t::Vector{Float64}`
"""
function eeg_spectrogram(eeg::NeuroAnalyzer.EEG; norm::Bool=true, mt::Bool=false, st::Bool=false, demean::Bool=true)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)
    fs = eeg_sr(eeg)
    p_tmp, s_frq, s_t = @views s_spectrogram(signal[1, :, 1], fs=fs, norm=norm, mt=mt, st=st, demean=demean)
    s_pow = zeros(size(p_tmp, 1), size(p_tmp, 2), channel_n, epoch_n)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            s_pow[:, :, channel_idx, epoch_idx], _, _ = @views s_spectrogram(signal[channel_idx, :, epoch_idx], fs=fs, norm=norm, mt=mt, st=st, demean=demean)
        end
    end

    s_frq = round.(s_frq, digits=2)
    s_t = round.(s_t, digits=2)
    s_t .+= eeg.eeg_epochs_time[1]

    return (s_pow=s_pow, s_frq=s_frq, s_t=s_t)
end

"""
    eeg_spectrum(eeg; pad, h)

Calculate FFT, amplitudes, powers and phases for each channel of `eeg`. For `pad` > 0 channels are padded with 0s.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `pad::Int64=0`: pad with `pad` zeros
- `h::Bool=false`: use Hilbert transform for calculations instead of FFT

# Returns

Named tuple containing:
- `fft::Array{ComplexF64, 3}`: Fourier or Hilbert components
- `amp::Array{Float64, 3}`: amplitudes
- `pow::Array{Float64, 3}`: powers
- `phase::Array{Float64, 3}: phase angles
"""
function eeg_spectrum(eeg::NeuroAnalyzer.EEG; pad::Int64=0, h::Bool=false)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    s_fft = zeros(ComplexF64, channel_n, eeg_epoch_len(eeg) + pad, epoch_n)
    s_amplitudes = zeros(channel_n, eeg_epoch_len(eeg) + pad, epoch_n)
    s_powers = zeros(channel_n, eeg_epoch_len(eeg) + pad, epoch_n)
    s_phases = zeros(channel_n, eeg_epoch_len(eeg) + pad, epoch_n)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            if h == false
                s_fft[channel_idx, :, epoch_idx], s_amplitudes[channel_idx, :, epoch_idx], s_powers[channel_idx, :, epoch_idx], s_phases[channel_idx, :, epoch_idx] = @views s_spectrum(signal[channel_idx, :, epoch_idx], pad=pad)
            else
                s_fft[channel_idx, :, epoch_idx], s_amplitudes[channel_idx, :, epoch_idx], s_powers[channel_idx, :, epoch_idx], s_phases[channel_idx, :, epoch_idx] = @views s_hspectrum(signal[channel_idx, :, epoch_idx], pad=pad)
            end
        end
    end

    return (fft=s_fft, amp=s_amplitudes, pow=s_powers, phase=s_phases)
end

"""
    eeg_s2t(eeg; t)

Convert time `t` in samples to seconds using `eeg` sampling rate.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `t::Int64`: time in samples

# Returns

- `t_s::Float64`: time in seconds
"""
function eeg_s2t(eeg::NeuroAnalyzer.EEG; t::Int64)
    return round(t / eeg_sr(eeg), digits=2)
end

"""
    eeg_t2s(eeg; t)

Convert time `t` in seconds to samples using `eeg` sampling rate.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `t::Real`: time in seconds

# Returns

- `t_s::Int64`: time in samples
"""
function eeg_t2s(eeg::NeuroAnalyzer.EEG; t::Real)
    return floor(Int64, t * eeg_sr(eeg)) + 1
end

"""
    eeg_channels_stats(eeg)

Calculate `eeg` channels statistics.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

Named tuple containing:
- `c_mean::Matrix(Float64)`: mean
- `c_median::Matrix(Float64)`: median
- `c_std::Matrix(Float64)`: standard deviation
- `c_var::Matrix(Float64)`: variance
- `c_kurt::Matrix(Float64)`: kurtosis
- `c_skew::Matrix(Float64)`: skewness
- `c_mean_diff::Matrix(Float64)`: mean diff value
- `c_median_diff::Matrix(Float64)`: median diff value
- `c_max_dif::Matrix(Float64)`: max difference
- `c_dev_mean::Matrix(Float64)`: deviation from channel mean
"""
function eeg_channels_stats(eeg::NeuroAnalyzer.EEG)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)
    
    c_mean = zeros(channel_n, epoch_n)
    c_median = zeros(channel_n, epoch_n)
    c_std = zeros(channel_n, epoch_n)
    c_var = zeros(channel_n, epoch_n)
    c_kurt = zeros(channel_n, epoch_n)
    c_skew = zeros(channel_n, epoch_n)
    c_mean_diff = zeros(channel_n, epoch_n)
    c_median_diff = zeros(channel_n, epoch_n)
    c_max_dif = zeros(channel_n, epoch_n)
    c_dev_mean = zeros(channel_n, epoch_n)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            c_mean[channel_idx, epoch_idx] = @views mean(signal[channel_idx, :, epoch_idx])
            c_median[channel_idx, epoch_idx] = @views median(signal[channel_idx, :, epoch_idx])
            c_std[channel_idx, epoch_idx] = @views std(signal[channel_idx, :, epoch_idx])
            c_var[channel_idx, epoch_idx] = @views var(signal[channel_idx, :, epoch_idx])
            c_kurt[channel_idx, epoch_idx] = @views kurtosis(signal[channel_idx, :, epoch_idx])
            c_skew[channel_idx, epoch_idx] = @views skewness(signal[channel_idx, :, epoch_idx])
            c_mean_diff[channel_idx, epoch_idx] = @views mean(diff(signal[channel_idx, :, epoch_idx]))
            c_median_diff[channel_idx, epoch_idx] = @views median(diff(signal[channel_idx, :, epoch_idx]))
            c_max_dif[channel_idx, epoch_idx] = @views maximum(signal[channel_idx, :, epoch_idx]) - minimum(signal[channel_idx, :, epoch_idx])
            c_dev_mean[channel_idx, epoch_idx] = @views abs(mean(signal[channel_idx, :, epoch_idx])) - mean(signal[channel_idx, :, epoch_idx])
        end
    end

    return (c_mean=c_mean, c_median=c_median, c_std=c_std, c_var=c_var, c_kurt=c_kurt, c_skew=c_skew, c_mean_diff=c_mean_diff, c_median_diff=c_median_diff, c_max_dif=c_max_dif, c_dev_mean=c_dev_mean)
end

"""
    eeg_snr(eeg)

Calculate SNR of `eeg` channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `snr::Matrix(Float64)`: SNR for each channel per epoch
"""
function eeg_snr(eeg::NeuroAnalyzer.EEG)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    snr = zeros(channel_n, epoch_n)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            snr[channel_idx, epoch_idx] = @views s_snr(signal[channel_idx, :, epoch_idx])
        end
    end

    return snr
end

"""
    eeg_standardize(eeg)

Standardize `eeg` channels for ML.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `eeg_new::NeuroAnalyzer.EEG`: standardized EEG
- `scaler::Matrix{Float64}`: standardizing matrix
"""
function eeg_standardize(eeg::NeuroAnalyzer.EEG)
    
    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    epoch_n = size(signal, 3)

    ss = similar(signal)
    scaler = Vector{Any}()

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:epoch_n
        @views push!(scaler, StatsBase.fit(ZScoreTransform, signal[:, :, epoch_idx], dims=2)) 
        @views eeg_new.eeg_signals[channels,:, epoch_idx] = StatsBase.transform(scaler[epoch_idx], signal[:, :, epoch_idx])
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_standardize(EEG)")

    return eeg_new, scaler
end

"""
    eeg_standardize!(eeg)

Standardize `eeg` channels for ML.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `scaler::Matrix{Float64}`: standardizing matrix
"""
function eeg_standardize!(eeg::NeuroAnalyzer.EEG)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    epoch_n = size(signal, 3)

    ss = similar(signal)
    scaler = Vector{Any}()

    @inbounds @simd for epoch_idx in 1:epoch_n
        @views push!(scaler, StatsBase.fit(ZScoreTransform, signal[:, :, epoch_idx], dims=2)) 
        @views eeg.eeg_signals[channels,:, epoch_idx] = StatsBase.transform(scaler[epoch_idx], signal[:, :, epoch_idx])
    end

    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_standardize!(EEG)")

    return scaler
end

"""
    eeg_fconv(eeg, kernel, norm)

Perform convolution of all `eeg` channels in the frequency domain using `kernel`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`: kernel for convolution
- `norm::Bool=false`: normalize kernel

# Returns

- `s_convoluted::Union{Array{Float64, 3}, Array{ComplexF64, 3}}`: convoluted signal
"""
function eeg_fconv(eeg::NeuroAnalyzer.EEG; kernel::Union{Vector{<:Real}, Vector{ComplexF64}}, norm::Bool=false)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    verbose == true && @info "This will take a while.."
    s_convoluted = zeros(ComplexF64, size(signal))
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            s_convoluted[channel_idx, :, epoch_idx] = @views s_fconv(signal[channel_idx, :, epoch_idx], kernel=kernel, norm=norm)
        end
    end

    return s_convoluted
end

"""
    eeg_tconv(eeg; kernel)

Perform convolution in the time domain.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`: kernel used for convolution

# Returns

- `s_convoluted::Union{Array{Float64, 3}, Array{ComplexF64, 3}}`: convoluted signal
"""
function eeg_tconv(eeg::NeuroAnalyzer.EEG; kernel::Union{Vector{<:Real}, Vector{ComplexF64}})

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    if typeof(kernel) == Vector{ComplexF64}
        s_convoluted = zeros(ComplexF64, size(signal))
    else
        s_convoluted = zeros(size(signal))
    end

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            s_convoluted[channel_idx, :, epoch_idx] = @views s_tconv(signal[channel_idx, :, epoch_idx], kernel=kernel)
        end
    end

    return s_convoluted
end

"""
    eeg_dft(eeg)

Returns FFT and DFT sample frequencies for a DFT for each `eeg` channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

Named tuple containing:
- `sfft::Array{ComplexF64, 3}`: FFT
- `sf::Vector{Float64}`: sample frequencies
"""
function eeg_dft(eeg::NeuroAnalyzer.EEG)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    fs = eeg_sr(eeg)
    sfft = zeros(ComplexF64, size(signal))
    sf = nothing

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            sfft[channel_idx, :, epoch_idx], sf = @views s_dft(signal[channel_idx, :, epoch_idx], fs=fs)
        end
    end

    return (sfft=sfft, sf=sf)
end

"""
    eeg_msci95(eeg; n::=3, method=:normal)

Calculates mean, std and 95% confidence interval for each `eeg` channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `n::Int64`: number of bootstraps
- `method::Symbol[:normal, :boot]`: use normal method or `n`-times boostrapping

# Returns

Named tuple containing:
- `s_m::Matrix{Float64}`: mean
- `s_s::Matrix{Float64}`: standard deviation
- `s_u::Matrix{Float64}`: upper 95% CI
- `s_l::Matrix{Float64}`: lower 95% CI
"""
function eeg_msci95(eeg::NeuroAnalyzer.EEG; n::Int64=3, method::Symbol=:normal)

    method in [:normal, :boot] || throw(ArgumentError("method must be :normal or :boot."))
    n < 1 && throw(ArgumentError("n must be ≥ 1."))

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_len = size(signal, 2)
    epoch_n = size(signal, 3)

    s_m = zeros(epoch_n, epoch_len)
    s_s = zeros(epoch_n, epoch_len)
    s_u = zeros(epoch_n, epoch_len)
    s_l = zeros(epoch_n, epoch_len)

    Threads.@threads for epoch_idx in 1:epoch_n
        s_m[epoch_idx, :], s_s[epoch_idx, :], s_u[epoch_idx, :], s_l[epoch_idx, :] = @views s_msci95(signal[:, :, epoch_idx], n=n, method=method)
    end

    return (mean=s_m, sd=s_s, upper=s_u, lower=s_l)
end

"""
    eeg_mean(eeg1, eeg2)

Calculates mean and 95% confidence interval for `eeg1` and `eeg2` channels.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2:NeuroAnalyzer.EEG`

# Returns

Named tuple containing:
- `s_m::Matrix{Float64}`: mean by epochs
- `s_s::Matrix{Float64}`: std by epochs
- `s_u::Matrix{Float64}`: upper 95% CI bound by epochs
- `s_l::Matrix{Float64}`: lower 95% CI bound by epochs
"""
function eeg_mean(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG)

    channels1 = eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))
    channels2 = eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))
    signal1 = @view eeg1.eeg_signals[channels1, :, :]
    signal2 = @view eeg2.eeg_signals[channels2, :, :]
    size(signal1) == size(signal1) || throw(ArgumentError("Both EEG signals must have the same size."))
    epoch_n = size(signal1, 3)
    epoch_len = size(signal1, 2)

    s_m = zeros(epoch_n, epoch_len)
    s_s = zeros(epoch_n, epoch_len)
    s_u = zeros(epoch_n, epoch_len)
    s_l = zeros(epoch_n, epoch_len)

    Threads.@threads for epoch_idx in 1:epoch_n
        s1_mean = @views mean(signal1[:, :, epoch_idx], dims=1)
        s2_mean = @views mean(signal2[:, :, epoch_idx], dims=1)
        s_m[epoch_idx, :] = s1_mean - s2_mean
        s1_sd = @views std(signal1[:, :, epoch_idx], dims=1) / sqrt(size(signal1[:, :, epoch_idx], 2))
        s2_sd = @views std(signal2[:, :, epoch_idx], dims=1) / sqrt(size(signal2[:, :, epoch_idx], 2))
        s_s[epoch_idx, :] = sqrt.(s1_sd.^2 .+ s2_sd.^2)
        s_u[epoch_idx, :] = @. s_m[epoch_idx, :] + 1.96 * s_s[epoch_idx, :]
        s_l[epoch_idx, :] = @. s_m[epoch_idx, :] - 1.96 * s_s[epoch_idx, :]
    end

    return (s_m=s_m, s_s=s_s, s_u=s_u, s_l=s_l)
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
- `s_stat::Matrix{Float64}`
- `s_stat_single::Vector{Float64}`
- `p::Vector{Float64}`
"""
function eeg_difference(eeg1::Array{Float64, 3}, eeg2::Array{Float64, 3}; n::Int64=3, method::Symbol=:absdiff)

    method in [:absdiff, :diff2int] || throw(ArgumentError("method must be :absdiff or :diff2int."))
    channels1 = eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))
    channels2 = eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))
    signal1 = @view eeg1.eeg_signals[channels1, :, :]
    signal2 = @view eeg2.eeg_signals[channels2, :, :]
    size(signal1) == size(signal1) || throw(ArgumentError("Both EEG signals must have the same size."))
    epoch_n = size(signal1, 3)

    s_stat = zeros(epoch_n, size(signal1, 1) * n)
    s_stat_single = zeros(epoch_n)
    p = zeros(epoch_n)

    Threads.@threads for epoch_idx in 1:epoch_n
        s_stat[epoch_idx, :], s_stat_single[epoch_idx], p[epoch_idx] = @views s2_difference(signal1[:, :, epoch_idx], signal2[:, :, epoch_idx])
    end

    return (s_stat=s_stat, statsitic_single=s_stat_single, p=p)
end

"""
   eeg_acov(eeg; lag=1, demean=false, norm=false)

Calculate autocovariance of each `eeg` channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `lag::Int64`: lags range is `-lag:lag`
- `demean::Bool`: demean eeg prior to analysis
- `norm::Bool`: normalize autocovariance

# Returns

Named tuple containing:
- `acov::Matrix{Float64}`
- `lags::Vector{Float64}`
"""
function eeg_acov(eeg::NeuroAnalyzer.EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)

    lag < 1 && throw(ArgumentError("lag must be ≥ 1."))

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    acov = zeros(channel_n, length(-lag:lag), epoch_n)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            acov[channel_idx, :, epoch_idx], _ = @views s_acov(signal[channel_idx, :, epoch_idx], lag=lag, demean=demean, norm=norm)
        end
    end

    lags = (eeg.eeg_time[2] - eeg.eeg_time[1]) .* collect(-lag:lag)

    return (acov=acov, acov_lags=lags)
end

"""
    eeg_tenv(eeg; d)

Calculate temporal envelope of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `t_env::Array{Float64, 3}`: temporal envelope
- `s_t::Vector{Float64}`: signal time
"""
function eeg_tenv(eeg::NeuroAnalyzer.EEG; d::Int64=32)
    
    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)
    t_env = similar(signal)
    s_t = eeg.eeg_epochs_time

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            s = @view signal[channel_idx, :, epoch_idx]
            p_idx = s_findpeaks(s, d=d)
            pushfirst!(p_idx, 1)
            push!(p_idx, length(s))
            if length(p_idx) > 4
                model = CubicSpline(s_t[p_idx], s[p_idx])
                try
                    t_env[channel_idx, :, epoch_idx] = model(s_t)
                catch
                    @error "CubicSpline error."
                end
            else
                verbose == true && @info "Less than 5 peaks detected, using Loess."
                model = loess(s_t[p_idx], s[p_idx], span=0.5)
                t_env[channel_idx, :, epoch_idx] = Loess.predict(model, s_t)
            end
            t_env[channel_idx, 1, epoch_idx] = t_env[channel_idx, 2, epoch_idx]
        end
    end
    
    return (t_env=t_env, s_t=s_t)
end

"""
    eeg_tenv_mean(eeg; dims, d)

Calculate temporal envelope of `eeg`: mean and 95% CI.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `t_env_m::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: mean
- `t_env_u::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI upper bound
- `t_env_l::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function eeg_tenv_mean(eeg::NeuroAnalyzer.EEG; dims::Int64, d::Int64=32)
    
    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)

    (channel_n == 1 || epoch_n == 1) && throw(ArgumentError("Number of channels and/or number of epochs must be ≥ 2."))

    s_a, s_t = eeg_tenv(eeg, d=d)

    if dims == 1
        t_env_m = zeros(length(s_t), epoch_n)
        t_env_u = zeros(length(s_t), epoch_n)
        t_env_l = zeros(length(s_t), epoch_n)

        @inbounds @simd for epoch_idx in 1:epoch_n
            t_env_m[:, epoch_idx] = mean(s_a[:, :, epoch_idx], dims=1)
            s = std(t_env_m[:, epoch_idx]) / sqrt(length(t_env_m[:, epoch_idx]))
            t_env_u[:, epoch_idx] = @. t_env_m[:, epoch_idx] + 1.96 * s
            t_env_l[:, epoch_idx] = @. t_env_m[:, epoch_idx] - 1.96 * s
        end
    elseif dims == 2
        t_env_m = zeros(length(s_t), channel_n)
        t_env_u = zeros(length(s_t), channel_n)
        t_env_l = zeros(length(s_t), channel_n)

        @inbounds @simd for channel_idx in 1:channel_n
            t_env_m[:, channel_idx] = mean(s_a[channel_idx, :, :], dims=2)
            s = std(t_env_m[:, channel_idx]) / sqrt(length(t_env_m[:, channel_idx]))
            t_env_u[:, channel_idx] = @. t_env_m[:, channel_idx] + 1.96 * s
            t_env_l[:, channel_idx] = @. t_env_m[:, channel_idx] - 1.96 * s
        end
    else
        t_env_m, t_env_u, t_env_l, _ = eeg_tenv_mean(eeg, dims=1, d=d)
        t_env_m = mean(t_env_m, dims=2)
        t_env_u = mean(t_env_u, dims=2)
        t_env_l = mean(t_env_l, dims=2)
        t_env_m = reshape(t_env_m, size(t_env_m, 1))
        t_env_u = reshape(t_env_u, size(t_env_u, 1))
        t_env_l = reshape(t_env_l, size(t_env_l, 1))
    end

    return (t_env_m=t_env_m, t_env_u=t_env_u, t_env_l=t_env_l, s_t=s_t)
end

"""
    eeg_tenv_median(eeg; dims, d)

Calculate temporal envelope of `eeg`: median and 95% CI.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `t_env_m::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: median
- `t_env_u::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI upper bound
- `t_env_l::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function eeg_tenv_median(eeg::NeuroAnalyzer.EEG; dims::Int64, d::Int64=32)
    
    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)

    (channel_n == 1 || epoch_n == 1) && throw(ArgumentError("Number of channels and/or number of epochs must be ≥ 2."))

    s_a, s_t = eeg_tenv(eeg, d=d)

    if dims == 1
        t_env_m = zeros(length(s_t), epoch_n)
        t_env_u = zeros(length(s_t), epoch_n)
        t_env_l = zeros(length(s_t), epoch_n)

        @inbounds @simd for epoch_idx in 1:epoch_n
            t_env_m[:, epoch_idx] = median(s_a[:, :, epoch_idx], dims=1)
            t_idx = s_findpeaks(t_env_m[:, epoch_idx], d=d)
            pushfirst!(t_idx, 1)
            push!(t_idx, length(t_env_m[:, epoch_idx]))
            if length(t_idx) > 4
                model = CubicSpline(s_t[t_idx], t_env_m[t_idx])
                try
                    t_env_m[:, epoch_idx] = model(s_t)
                catch
                    verbose == true && @info "CubicSpline could not be calculated, using non-smoothed variant instead."
                end
            end
            s = iqr(t_env_m[:, epoch_idx]) / sqrt(length(t_env_m[:, epoch_idx]))
            t_env_u[:, epoch_idx] = @. t_env_m[:, epoch_idx] + 1.96 * s
            t_env_l[:, epoch_idx] = @. t_env_m[:, epoch_idx] - 1.96 * s
        end
    elseif dims == 2
        t_env_m = zeros(length(s_t), channel_n)
        t_env_u = zeros(length(s_t), channel_n)
        t_env_l = zeros(length(s_t), channel_n)

        @inbounds @simd for channel_idx in 1:channel_n
            t_env_m[:, idx] = median(s_a[channel_idx, :, :], dims=2)
            t_idx = s_findpeaks(t_env_m[:, channel_idx], d=d)
            pushfirst!(t_idx, 1)
            push!(t_idx, length(t_env_m[:, channel_idx]))
            if length(t_idx) > 4
                model = CubicSpline(s_t[t_idx], t_env_m[t_idx])
                try
                    t_env_m[:, channel_idx] = model(s_t)
                catch
                    verbose == true && @info "CubicSpline could not be calculated, using non-smoothed variant instead."
                end
            end
            s = iqr(t_env_m[:, channel_idx]) / sqrt(length(t_env_m[:, channel_idx]))
            t_env_u[:, channel_idx] = @. t_env_m[:, channel_idx] + 1.96 * s
            t_env_l[:, channel_idx] = @. t_env_m[:, channel_idx] - 1.96 * s
        end
    else
        t_env_m, t_env_u, t_env_l, _ = eeg_tenv_median(eeg, dims=1, d=d)
        t_env_m = median(t_env_m, dims=2)
        t_env_u = median(t_env_u, dims=2)
        t_env_l = median(t_env_l, dims=2)
        t_env_m = reshape(t_env_m, size(t_env_m, 1))
        t_env_u = reshape(t_env_u, size(t_env_u, 1))
        t_env_l = reshape(t_env_l, size(t_env_l, 1))
    end

    return (t_env_m=t_env_m, t_env_u=t_env_u, t_env_l=t_env_l, s_t=s_t)
end

"""
    eeg_penv(eeg; d)

Calculate power (in dB) envelope of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

Named tuple containing:
- `p_env::Array{Float64, 3}`: power spectrum envelope
- `p_env_frq::Vector{Float64}`: frequencies for each envelope
"""
function eeg_penv(eeg::NeuroAnalyzer.EEG; d::Int64=8, mt::Bool=false)
    
    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)
    p_env = similar(signal)
    fs = eeg_sr(eeg)

    psd_tmp, frq = s_psd(signal[1, :, 1], fs=fs, mt=mt)
    p_env = zeros(channel_n, length(psd_tmp), epoch_n)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            psd_pow, _ = s_psd(signal[channel_idx, :, epoch_idx], fs=fs, mt=mt, norm=true)
            p_idx = s_findpeaks(psd_pow, d=d)
            pushfirst!(p_idx, 1)
            push!(p_idx, length(psd_pow))
            if length(p_idx) > 4
                model = CubicSpline(frq[p_idx], psd_pow[p_idx])
                try
                    p_env[channel_idx, :, epoch_idx] = model(frq)
                catch
                    verbose == true && @info "CubicSpline could not be calculated, using non-smoothed variant instead."
                end
            else
                p_env[channel_idx, :, epoch_idx] = psd_pow
            end
            p_env[channel_idx, 1, epoch_idx] = p_env[channel_idx, 2, epoch_idx]
        end
    end
    
    return (p_env=p_env, p_env_frq=frq)
end

"""
    eeg_penv_mean(eeg; dims, d)

Calculate power (in dB) envelope of `eeg`: mean and 95% CI.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

Named tuple containing:
- `p_env_m::Array{Float64, 3}`: power spectrum envelope: mean
- `p_env_u::Array{Float64, 3}`: power spectrum envelope: 95% CI upper bound
- `p_env_l::Array{Float64, 3}`: power spectrum envelope: 95% CI lower bound
- `p_env_frq::Vector{Float64}`: power spectrum envelope (useful for plotting over PSD)
"""
function eeg_penv_mean(eeg::NeuroAnalyzer.EEG; dims::Int64, d::Int64=8, mt::Bool=false)
    
    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)

    (channel_n == 1 || epoch_n == 1) && throw(ArgumentError("Number of channels and/or number of epochs must be ≥ 2."))

    s_p, s_f = eeg_psd(eeg, norm=true, mt=mt)

    if dims == 1
        p_env_m = zeros(length(s_f), epoch_n)
        p_env_u = zeros(length(s_f), epoch_n)
        p_env_l = zeros(length(s_f), epoch_n)

        @inbounds @simd for epoch_idx in 1:epoch_n
            p_env_m[:, epoch_idx] = mean(s_p[:, :, epoch_idx], dims=1)

            p_idx = s_findpeaks(p_env_m[:, epoch_idx], d=d)
            pushfirst!(p_idx, 1)
            push!(p_idx, length(p_env_m[:, epoch_idx]))
            if length(p_idx) > 4
                model = CubicSpline(s_f[p_idx], p_env_m[p_idx])
                try
                    p_env_m[:, epoch_idx] = model(s_f)
                catch
                    verbose == true && @info "CubicSpline could not be calculated, using non-smoothed variant instead."
                end
            end
            s = std(p_env_m[:, epoch_idx]) / sqrt(length(p_env_m[:, epoch_idx]))
            p_env_u[:, epoch_idx] = @. p_env_m[:, epoch_idx] + 1.96 * s
            p_env_l[:, epoch_idx] = @. p_env_m[:, epoch_idx] - 1.96 * s
        end
    elseif dims == 2
        p_env_m = zeros(length(s_f), channel_n)
        p_env_u = zeros(length(s_f), channel_n)
        p_env_l = zeros(length(s_f), channel_n)

        @inbounds @simd for channel_idx in 1:channel_n
            p_env_m[:, channel_idx] = mean(s_p[channel_idx, :, :], dims=2)

            p_idx = s_findpeaks(p_env_m[:, channel_idx], d=d)
            pushfirst!(p_idx, 1)
            push!(p_idx, length(p_env_m[:, channel_idx]))
            if length(p_idx) > 4
                model = CubicSpline(s_f[p_idx], p_env_m[p_idx])
                try
                    p_env_m[:, channel_idx] = model(s_f)
                catch
                    verbose == true && @info "CubicSpline could not be calculated, using non-smoothed variant instead."
                end
            end
            s = std(p_env_m[:, channel_idx]) / sqrt(length(p_env_m[:, channel_idx]))
            p_env_u[:, channel_idx] = @. p_env_m[:, channel_idx] + 1.96 * s
            p_env_l[:, channel_idx] = @. p_env_m[:, channel_idx] - 1.96 * s
        end
    else
        p_env_m, p_env_u, p_env_l, _ = eeg_penv_mean(eeg, dims=1, d=d)
        p_env_m = mean(p_env_m, dims=2)
        p_env_u = mean(p_env_u, dims=2)
        p_env_l = mean(p_env_l, dims=2)
        p_env_m = reshape(p_env_m, size(p_env_m, 1))
        p_env_u = reshape(p_env_u, size(p_env_u, 1))
        p_env_l = reshape(p_env_l, size(p_env_l, 1))
    end
    
    return (p_env_m=p_env_m, p_env_u=p_env_u, p_env_l=p_env_l, p_env_frq=s_f)
end

"""
    eeg_penv_median(eeg; dims, d)

Calculate power (in dB) envelope of `eeg`: median and 95% CI.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `dims::Int64`: median over channels (dims = 1) or epochs (dims = 2)
- `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

Named tuple containing:
- `p_env_m::Array{Float64, 3}`: power spectrum envelope: median
- `p_env_u::Array{Float64, 3}`: power spectrum envelope: 95% CI upper bound
- `p_env_l::Array{Float64, 3}`: power spectrum envelope: 95% CI lower bound
- `p_env_frq::Vector{Float64}`: power spectrum envelope (useful for plotting over PSD)
"""
function eeg_penv_median(eeg::NeuroAnalyzer.EEG; dims::Int64, d::Int64=8, mt::Bool=false)
    
    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)

    (channel_n == 1 || epoch_n == 1) && throw(ArgumentError("Number of channels and/or number of epochs must be ≥ 2."))

    s_p, s_f = eeg_psd(eeg, norm=true, mt=mt)

    if dims == 1
        p_env_m = zeros(length(s_f), epoch_n)
        p_env_u = zeros(length(s_f), epoch_n)
        p_env_l = zeros(length(s_f), epoch_n)

        @inbounds @simd for epoch_idx in 1:epoch_n
            p_env_m[:, epoch_idx] = median(s_p[:, :, epoch_idx], dims=1)

            p_idx = s_findpeaks(p_env_m[:, epoch_idx], d=d)
            pushfirst!(p_idx, 1)
            push!(p_idx, length(p_env_m[:, epoch_idx]))
            if length(p_idx) > 4
                model = CubicSpline(s_f[p_idx], p_env_m[p_idx])
                try
                    p_env_m[:, epoch_idx] = model(s_f)
                catch
                    verbose == true && @info "CubicSpline could not be calculated, using non-smoothed variant instead."
                end
            end
            s = iqr(p_env_m[:, epoch_idx]) / sqrt(length(p_env_m[:, epoch_idx]))
            p_env_u[:, epoch_idx] = @. p_env_m[:, epoch_idx] + 1.96 * s
            p_env_l[:, epoch_idx] = @. p_env_m[:, epoch_idx] - 1.96 * s
        end
    elseif dims == 2
        p_env_m = zeros(length(s_f), channel_n)
        p_env_u = zeros(length(s_f), channel_n)
        p_env_l = zeros(length(s_f), channel_n)

        @inbounds @simd for channel_idx in 1:channel_n
            p_env_m[:, channel_idx] = median(s_p[channel_idx, :, :], dims=2)

            p_idx = s_findpeaks(p_env_m[:, channel_idx], d=d)
            pushfirst!(p_idx, 1)
            push!(p_idx, length(p_env_m[:, channel_idx]))
            if length(p_idx) > 4
                model = CubicSpline(s_f[p_idx], p_env_m[p_idx])
                try
                    p_env_m[:, channel_idx] = model(s_f)
                catch
                    verbose == true && @info "CubicSpline could not be calculated, using non-smoothed variant instead."
                end
            end
            s = iqr(p_env_m[:, channel_idx]) / sqrt(length(p_env_m[:, channel_idx]))
            p_env_u[:, channel_idx] = @. p_env_m[:, channel_idx] + 1.96 * s
            p_env_l[:, channel_idx] = @. p_env_m[:, channel_idx] - 1.96 * s
        end
    else
        p_env_m, p_env_u, p_env_l, _ = eeg_penv_median(eeg, dims=1, d=d)
        p_env_m = median(p_env_m, dims=2)
        p_env_u = median(p_env_u, dims=2)
        p_env_l = median(p_env_l, dims=2)
        p_env_m = reshape(p_env_m, size(p_env_m, 1))
        p_env_u = reshape(p_env_u, size(p_env_u, 1))
        p_env_l = reshape(p_env_l, size(p_env_l, 1))
    end
    
    return (p_env_m=p_env_m, p_env_u=p_env_u, p_env_l=p_env_l, p_env_frq=s_f)
end

"""
    eeg_senv(eeg; d, mt, t)

Calculate spectral envelope of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
- `mt::Bool=false`: if true use multi-tapered spectrogram
- `t::Union{Real, Nothing}=nothing`: spectrogram threshold (maximize all powers > t)

# Returns

Named tuple containing:
- `s_env::Array{Float64, 3}`: spectral envelope
- `s_env_t::Vector{Float64}`: spectrogram time
"""
function eeg_senv(eeg::NeuroAnalyzer.EEG; d::Int64=2, mt::Bool=false, t::Union{Real, Nothing}=nothing)
    
    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)
    fs = eeg_sr(eeg)
    s_env_s = similar(signal)
    
    s_tmp = @view signal[1, :, 1]
    interval = fs
    overlap = round(Int64, fs * 0.85)
    length(s_tmp) < 4 * fs && (mt = true)
    if mt == false
        spec_tmp = spectrogram(s_tmp, interval, overlap, nfft=length(s_tmp), fs=fs, window=hanning)
    else
        spec_tmp = mt_spectrogram(s_tmp, fs=fs)
    end
    sp_t = collect(spec_tmp.time)
    sp_t .+= eeg.eeg_epochs_time[1]

    s_env = zeros(channel_n, length(sp_t), epoch_n)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            if mt == false
                spec = @views spectrogram(signal[channel_idx, :, epoch_idx], interval, overlap, nfft=length(s_tmp), fs=fs, window=hanning)
            else
                spec = @views mt_spectrogram(signal[channel_idx, :, epoch_idx], fs=fs)
            end

            s_frq = Vector(spec.freq)
            s_p = pow2db.(spec.power)
            if t !== nothing
                s_p[s_p .> t] .= 0
                reverse!(s_p)
                reverse!(s_frq)
            end
            
            f_idx = zeros(length(spec.time))
            m = maximum(s_p, dims=1)
            for idx2 in 1:length(m)
                f_idx[idx2] = s_frq[vsearch(m[idx2], s_p[:, idx2])]
            end
            p_idx = s_findpeaks(f_idx, d=d)
            pushfirst!(p_idx, 1)
            push!(p_idx, length(spec.time))
            if length(p_idx) > 4
                model = CubicSpline(sp_t[p_idx], f_idx[p_idx])
                try
                    s_env[channel_idx, :, epoch_idx] = model(sp_t)
                catch
                    verbose == true && @info "CubicSpline could not be calculated, using non-smoothed variant instead."
                end
            else
                s_env[channel_idx, :, epoch_idx] = f_idx
            end
            s_env[channel_idx, 1, epoch_idx] = s_env[channel_idx, 2, epoch_idx]
        end
    end
    
    return (s_env=s_env, senv_t=sp_t)
end

"""
    eeg_senv_mean(eeg; dims, d, mt, t)

Calculate spectral envelope of `eeg`: mean and 95% CI.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
- `mt::Bool=false`: if true use multi-tapered spectrogram
- `t::Union{Real, Nothing}=nothing`: spectrogram threshold (maximize all powers > t)

# Returns

Named tuple containing:
- `s_env_m::Array{Float64, 3}`: spectral envelope: mean
- `s_env_u::Array{Float64, 3}`: spectral envelope: 95% CI upper bound
- `s_env_l::Array{Float64, 3}`: spectral envelope: 95% CI lower bound
- `s_env_t::Vector{Float64}`: spectral envelope (useful for plotting over spectrogram)
"""
function eeg_senv_mean(eeg::NeuroAnalyzer.EEG; dims::Int64, d::Int64=2, mt::Bool=false, t::Union{Real, Nothing}=nothing)

    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)

    (channel_n == 1 || epoch_n == 1) && throw(ArgumentError("Number of channels and/or number of epochs must be ≥ 2."))

    s_p, s_t = eeg_senv(eeg, d=d, mt=mt, t=t)

    if dims == 1
        s_env_m = zeros(length(s_t), epoch_n)
        s_env_u = zeros(length(s_t), epoch_n)
        s_env_l = zeros(length(s_t), epoch_n)

        @inbounds @simd for epoch_idx in 1:epoch_n
            s_env_m[:, epoch_idx] = mean(s_p[:, :, epoch_idx], dims=1)

            s_idx = s_findpeaks(s_env_m[:, epoch_idx], d=d)
            pushfirst!(s_idx, 1)
            push!(s_idx, length(s_env_m[:, epoch_idx]))
            if length(s_idx) > 4
                model = CubicSpline(s_t[s_idx], s_env_m[s_idx])
                try
                    s_env_m[:, epoch_idx] = model(s_t)
                catch
                    verbose == true && @info "CubicSpline could not be calculated, using non-smoothed variant instead."
                end
            end
            s = std(s_env_m[:, epoch_idx]) / sqrt(length(s_env_m[:, epoch_idx]))
            s_env_u[:, epoch_idx] = @. s_env_m[:, epoch_idx] + 1.96 * s
            s_env_l[:, epoch_idx] = @. s_env_m[:, epoch_idx] - 1.96 * s
        end
    elseif dims == 2
        s_env_m = zeros(length(s_t), channel_n)
        s_env_u = zeros(length(s_t), channel_n)
        s_env_l = zeros(length(s_t), channel_n)

        @inbounds @simd for channel_idx in 1:channel_n
            s_env_m[:, channel_idx] = mean(s_p[channel_idx, :, :], dims=2)

            s_idx = s_findpeaks(s_env_m[:, channel_idx], d=d)
            pushfirst!(s_idx, 1)
            push!(s_idx, length(s_env_m[:, channel_idx]))
            if length(s_idx) > 4
                model = CubicSpline(s_t[s_idx], s_env_m[s_idx])
                try
                    s_env_m[:, channel_idx] = model(s_t)
                catch
                    verbose == true && @info "CubicSpline could not be calculated, using non-smoothed variant instead."
                end
            end
            s = std(s_env_m[:, channel_idx]) / sqrt(length(s_env_m[:, channel_idx]))
            s_env_u[:, channel_idx] = @. s_env_m[:, channel_idx] + 1.96 * s
            s_env_l[:, channel_idx] = @. s_env_m[:, channel_idx] - 1.96 * s
        end
    else
        s_env_m, s_env_u, s_env_l, _ = eeg_senv_mean(eeg, dims=1, d=d, mt=mt)
        s_env_m = mean(s_env_m, dims=2)
        s_env_u = mean(s_env_u, dims=2)
        s_env_l = mean(s_env_l, dims=2)
        s_env_m = reshape(s_env_m, size(s_env_m, 1))
        s_env_u = reshape(s_env_u, size(s_env_u, 1))
        s_env_l = reshape(s_env_l, size(s_env_l, 1))
    end
    
    return (s_env_m=s_env_m, s_env_u=s_env_u, s_env_l=s_env_l, s_env_t=s_t)
end

"""
    eeg_senv_median(eeg; dims, d, mt)

Calculate spectral envelope of `eeg`: median and 95% CI.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `dims::Int64`: mean over chan (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
- `mt::Bool=false`: if true use multi-tapered spectrogram
- `t::Union{Real, Nothing}=nothing`: spectrogram threshold (maximize all powers > t)

# Returns

Named tuple containing:
- `s_env_m::Array{Float64, 3}`: spectral envelope: median
- `s_env_u::Array{Float64, 3}`: spectral envelope: 95% CI upper bound
- `s_env_l::Array{Float64, 3}`: spectral envelope: 95% CI lower bound
- `s_env_t::Vector{Float64}`: spectral envelope (useful for plotting over spectrogram)
"""
function eeg_senv_median(eeg::NeuroAnalyzer.EEG; dims::Int64, d::Int64=2, mt::Bool=false, t::Union{Real, Nothing}=nothing)
    
    channel_n = eeg_channel_n(eeg)
    epoch_n = eeg_epoch_n(eeg)

    (channel_n == 1 || epoch_n == 1) && throw(ArgumentError("Number of channels and/or number of epochs must be ≥ 2."))

    s_p, s_t = eeg_senv(eeg, d=d, mt=mt, t=t)

    if dims == 1
        s_env_m = zeros(length(s_t), epoch_n)
        s_env_u = zeros(length(s_t), epoch_n)
        s_env_l = zeros(length(s_t), epoch_n)

        @inbounds @simd for epoch_idx in 1:epoch_n
            s_env_m[:, epoch_idx] = median(s_p[:, :, epoch_idx], dims=1)

            s_idx = s_findpeaks(s_env_m[:, epoch_idx], d=d)
            pushfirst!(s_idx, 1)
            push!(s_idx, length(s_env_m[:, epoch_idx]))
            if length(s_idx) > 4
                model = CubicSpline(s_t[s_idx], s_env_m[s_idx])
                try
                    s_env_m[:, epoch_idx] = model(s_t)
                catch
                    verbose == true && @info "CubicSpline could not be calculated, using non-smoothed variant instead."
                end
            end
            s = iqr(s_env_m[:, epoch_idx]) / sqrt(length(s_env_m[:, epoch_idx]))
            s_env_u[:, epoch_idx] = @. s_env_m[:, epoch_idx] + 1.96 * s
            s_env_l[:, epoch_idx] = @. s_env_m[:, epoch_idx] - 1.96 * s
        end
    elseif dims == 2
        s_env_m = zeros(length(s_t), channel_n)
        s_env_u = zeros(length(s_t), channel_n)
        s_env_l = zeros(length(s_t), channel_n)

        @inbounds @simd for channel_idx in 1:channel_n
            s_env_m[:, channel_idx] = median(s_p[channel_idx, :, :], dims=2)

            s_idx = s_findpeaks(s_env_m[:, channel_idx], d=d)
            pushfirst!(s_idx, 1)
            push!(s_idx, length(s_env_m[:, channel_idx]))
            if length(s_idx) > 4
                model = CubicSpline(s_t[s_idx], s_env_m[s_idx])
                try
                    s_env_m[:, channel_idx] = model(s_t)
                catch
                    verbose == true && @info "CubicSpline could not be calculated, using non-smoothed variant instead."
                end
            end
            s = iqr(s_env_m[:, channel_idx]) / sqrt(length(s_env_m[:, channel_idx]))
            s_env_u[:, channel_idx] = @. s_env_m[:, channel_idx] + 1.96 * s
            s_env_l[:, channel_idx] = @. s_env_m[:, channel_idx] - 1.96 * s
        end
    else
        s_env_m, s_env_u, s_env_l, _ = eeg_senv_median(eeg, dims=1, d=d, mt=mt)
        s_env_m = median(s_env_m, dims=2)
        s_env_u = median(s_env_u, dims=2)
        s_env_l = median(s_env_l, dims=2)
        s_env_m = reshape(s_env_m, size(s_env_m, 1))
        s_env_u = reshape(s_env_u, size(s_env_u, 1))
        s_env_l = reshape(s_env_l, size(s_env_l, 1))
    end
    
    return (s_env_m=s_env_m, s_env_u=s_env_u, s_env_l=s_env_l, s_env_t=s_t)
end

"""
    eeg_ispc(eeg1, eeg2; channel1, channel2, epoch1, epoch2)

Calculate ISPC (Inter-Site-Phase Clustering) between `channel1`/`epoch1` and `channel2` of `epoch2` of `eeg`.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2::NeuroAnalyzer.EEG`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=0`: default use all channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=0`: default use all channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=0`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=0`: default use all epochs

# Returns

Named tuple containing:
- `ispc::Array{Float64, 2}`: ISPC value
- `ispc_angle::Array{Float64, 2}`: ISPC angle
- `signal_diff::Array{Float64, 3}`: signal difference (signal2 - signal1)
- `phase_diff::Array{Float64, 3}`: phase difference (signal2 - signal1)
- `s1_phase::Array{Float64, 3}`: signal 1 phase
- `s2_phase::Array{Float64, 3}`: signal 2 phase
"""
function eeg_ispc(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Union{Int64, Vector{Int64}, AbstractRange}=0, channel2::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch1::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch2::Union{Int64, Vector{Int64}, AbstractRange}=0)

    # select channels, default is all
    channels1 = eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))
    channels2 = eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))
    channel1 == 0 && (channel1 = channels1)
    _check_channels(channels1, channel1)
    channel2 == 0 && (channel2 = channels2)
    _check_channels(channels2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    # select epochs, default is all
    epoch1 == 0 && (epoch1 = _select_epochs(eeg1, epoch1))
    _check_epochs(eeg1, epoch1)
    epoch2 == 0 && (epoch2 = _select_epochs(eeg2, epoch2))
    _check_epochs(eeg2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    eeg_epoch_len(eeg1) == eeg_epoch_len(eeg2) || throw(ArgumentError("eeg1 and eeg2 epoch lengths must be equal."))

    ispc = zeros(length(channel1), length(epoch1))
    ispc_angle = zeros(length(channel1), length(epoch1))
    signal_diff = zeros(length(channel1), eeg_epoch_len(eeg1), length(epoch1))
    phase_diff = zeros(length(channel1), eeg_epoch_len(eeg1), length(epoch1))
    s1_phase = zeros(length(channel1), eeg_epoch_len(eeg1), length(epoch1))
    s2_phase = zeros(length(channel1), eeg_epoch_len(eeg1), length(epoch1))

    @inbounds @simd for epoch_idx in 1:length(epoch1)
        Threads.@threads for channel_idx in 1:length(channel1)
            ispc[channel_idx, epoch_idx], ispc_angle[channel_idx, epoch_idx], signal_diff[channel_idx, :, epoch_idx], phase_diff[channel_idx, :, epoch_idx], s1_phase[channel_idx, :, epoch_idx], s2_phase[channel_idx, :, epoch_idx] = @views s2_ispc(eeg1.eeg_signals[channel1[channel_idx], :, epoch1[epoch_idx]], eeg2.eeg_signals[channel2[channel_idx], :, epoch2[epoch_idx]])
        end
    end

    return (ispc=ispc, ispc_angle=ispc_angle, signal_diff=signal_diff, phase_diff=phase_diff, s1_phase=s1_phase, s2_phase=s2_phase)
end

"""
    eeg_itpc(eeg; channel, t, w)

Calculate ITPC (Inter-Trial-Phase Clustering) at time `t` over epochs/trials of `channel` of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64`
- `t::Int64`: time point (sample number) at which ITPC is calculated
- `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation

# Returns

Named tuple containing:
- `itpc::Float64`: ITPC or wITPC value
- `itpcz::Float64`: Rayleigh's ITPC Z value
- `itpc_angle::Float64`: ITPC angle
- `phase_diff::Array{Float64, 3}`: phase difference (channel2 - channel1)
"""
function eeg_itpc(eeg::NeuroAnalyzer.EEG; channel::Int64, t::Int64, w::Union{Vector{<:Real}, Nothing}=nothing)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]

    channel in channels || throw(ArgumentError("channel does not mach signal."))
    
    itpc, itpcz, itpc_angle, itpc_phases = @views s_itpc(reshape(signal[channel, :, :], 1, size(signal[channel, :, :], 1), size(signal[channel, :, :], 2)), t=t, w=w)

    return (itpc=itpc, itpcz=itpcz, itpc_angle=itpc_angle, itpc_phases=itpc_phases)
end

"""
    eeg_pli(eeg1, eeg2; channel1, channel2, epoch1, epoch2)

Calculate PLI (Phase Lag Index) between `eeg1` and `eeg2`.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2::NeuroAnalyzer.EEG`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=0`: default use all channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=0`: default use all channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=0`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=0`: default use all epochs

# Returns

Named tuple containing:
- `pli::Array{Float64, 2}`: PLI value
- `signal_diff::Array{Float64, 3}`: signal difference (signal2 - signal1)
- `phase_diff::Array{Float64, 3}`: phase difference (signal2 - signal1)
- `s1_phase::Array{Float64, 3}`: signal 1 phase
- `s2_phase::Array{Float64, 3}`: signal 2 phase
"""
function eeg_pli(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Union{Int64, Vector{Int64}, AbstractRange}=0, channel2::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch1::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch2::Union{Int64, Vector{Int64}, AbstractRange}=0)

    # select channels, default is all
    channels1 = eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))
    channels2 = eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))
    channel1 == 0 && (channel1 = channels1)
    _check_channels(channels1, channel1)
    channel2 == 0 && (channel2 = channels2)
    _check_channels(channels2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    # select epochs, default is all
    epoch1 == 0 && (epoch1 = _select_epochs(eeg1, epoch1))
    _check_epochs(eeg1, epoch1)
    epoch2 == 0 && (epoch2 = _select_epochs(eeg2, epoch2))
    _check_epochs(eeg2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    eeg_epoch_len(eeg1) == eeg_epoch_len(eeg2) || throw(ArgumentError("eeg1 and eeg2 epoch lengths must be equal."))

    pli = zeros(length(channel1), length(epoch1))
    signal_diff = zeros(length(channel1), eeg_epoch_len(eeg1), length(epoch1))
    phase_diff = zeros(length(channel1), eeg_epoch_len(eeg1), length(epoch1))
    s1_phase = zeros(length(channel1), eeg_epoch_len(eeg1), length(epoch1))
    s2_phase = zeros(length(channel1), eeg_epoch_len(eeg1), length(epoch1))

    @inbounds @simd for epoch_idx in 1:length(epoch1)
        Threads.@threads for channel_idx in 1:length(channel1)
            pli[channel_idx, epoch_idx], signal_diff[channel_idx, :, epoch_idx], phase_diff[channel_idx, :, epoch_idx], s1_phase[channel_idx, :, epoch_idx], s2_phase[channel_idx, :, epoch_idx] = @views s2_pli(eeg1.eeg_signals[channel1[channel_idx], :, epoch1[epoch_idx]], eeg2.eeg_signals[channel2[channel_idx], :, epoch2[epoch_idx]])
        end
    end

    return (pli=pli, signal_diff=signal_diff, phase_dif=phase_diff, s1_phase=s1_phase, s2_phase=s2_phase)
end

"""
    eeg_pli(eeg)

Calculate PLIs (Phase Lag Index) between all channels of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `pli_m::Array{Float64, 3}`: PLI value matrices over epochs
"""
function eeg_pli(eeg::NeuroAnalyzer.EEG)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    pli_m = zeros(channel_n, channel_n, epoch_n)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx1 in 1:channel_n
            for channel_idx2 in 1:channel_idx1
                pli, _, _, _, _ = @views s2_pli(signal[channel_idx1, :, epoch_idx], signal[channel_idx2, :, epoch_idx])
                pli_m[channel_idx1, channel_idx2, epoch_idx] = round(pli, digits=4)
            end
        end
        Threads.@threads for channel_idx1 in 1:(channel_n - 1)
            for channel_idx2 in (channel_idx1 + 1):channel_n
                pli_m[channel_idx1, channel_idx2, epoch_idx] = @views pli_m[channel_idx2, channel_idx1, epoch_idx]
            end
        end
    end

    return pli_m
end

"""
    eeg_ispc(eeg)

Calculate ISPCs (Inter-Site-Phase Clustering) between all channels of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `ispc_m::Array{Float64, 3}`: ISPC value matrices over epochs
"""
function eeg_ispc(eeg::NeuroAnalyzer.EEG)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    ispc_m = zeros(channel_n, channel_n, epoch_n)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx1 in 1:channel_n
            for channel_idx2 in 1:channel_idx1
                ispc, _, _, _, _, _ = @views s2_ispc(signal[channel_idx1, :, epoch_idx], signal[channel_idx2, :, epoch_idx])
                # idx1 == idx2 && (ispc = 0)
                ispc_m[channel_idx1, channel_idx2, epoch_idx] = round(ispc, digits=4)
            end
        end
        Threads.@threads for channel_idx1 in 1:(channel_n - 1)
            for channel_idx2 in (channel_idx1 + 1):channel_n
                ispc_m[channel_idx1, channel_idx2, epoch_idx] = @views ispc_m[channel_idx2, channel_idx1, epoch_idx]
            end
        end
    end

    return ispc_m
end

"""
    eeg_aec(eeg1, eeg2; channel1, channel2, epoch1, epoch2)

Calculate amplitude envelope correlation between `channel1`/`epoch1` and `channel2` of `epoch2` of `eeg`.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2::NeuroAnalyzer.EEG`
- `channel1::Int64`
- `channel2::Int64`
- `epoch1::Int64`
- `epoch2::Int64`

# Returns

Named tuple containing:
- `aec::Float64`: power correlation value
- `aec_p::Float64`: power correlation p-value
"""
function eeg_aec(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Int64, channel2::Int64, epoch1::Int64, epoch2::Int64)

    eeg_epoch_len(eeg1) == eeg_epoch_len(eeg2) || throw(ArgumentError("eeg1 and eeg2 must have the same epoch length."))
    (channel1 < 0 || channel2 < 0 || epoch1 < 0 || epoch2 < 0) && throw(ArgumentError("channel1/epoch1/channel2/epoch2 must be > 0."))

    channels1 = eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))
    signal1 = @view eeg1.eeg_signals[channels1, :, :]
    channel_n1 = size(signal1, 1)
    epoch_n1 = size(signal1, 3)

    channels2 = eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))
    signal2 = @view eeg2.eeg_signals[channels2, :, :]
    channel_n2 = size(signal2, 1)
    epoch_n2 = size(signal2, 3)

    (channel1 > channel_n1) && throw(ArgumentError("channel1 must be ≤ $(channel_n1)."))
    (epoch1 > epoch_n1) && throw(ArgumentError("epoch1 must be ≤ $(epoch_n1)."))
    (channel2 > channel_n2) && throw(ArgumentError("channel2 must be ≤ $(channel_n2)."))
    (epoch2 > epoch_n2) && throw(ArgumentError("epoch2 must be ≤ $(epoch_n2)."))

    if epoch_n1 > 1
        e1 = eeg_keep_epoch(eeg1, epoch=epoch1)
    else
        e1 = deepcopy(eeg1)
    end
    if epoch_n2 > 1
        e2 = eeg_keep_epoch(eeg2, epoch=epoch1)
    else
        e2 = deepcopy(eeg2)
    end

    eeg_keep_channel!(e1, channel=channel1)
    eeg_keep_channel!(e2, channel=channel2)
    s1, _ = eeg_tenv(e1)
    s2, _ = eeg_tenv(e2)
    aec = CorrelationTest(vec(s1), vec(s2))
    aec_r = aec.r
    aec_p = pvalue(aec)

    return (aec=aec_r, aec_p=aec_p)
end

"""
    eeg_ged(eeg1, eeg2)

Perform generalized eigendecomposition between `eeg1` and `eeg2`.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`: signal data to be analyzed
- `eeg2::NeuroAnalyzer.EEG`: original signal data

# Returns

- `sged::Array{Float64, 3}`
- `ress::Matrix{Float64}`
- `ress_normalized::Matrix{Float64}`
"""
function eeg_ged(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG)

    size(eeg1.eeg_signals) == size(eeg2.eeg_signals) || throw(ArgumentError("eeg1 and eeg2 signal data must have the same size."))

    channels1 = eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))
    signal1 = @view eeg1.eeg_signals[channels1, :, :]
    channel_n = size(signal1, 1)
    epoch_n = size(signal1, 3)
    channels2 = eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))
    signal2 = @view eeg2.eeg_signals[channels2, :, :]

    sged = similar(eeg1.eeg_signals)
    ress = zeros(channel_n, epoch_n)
    ress_normalized = zeros(channel_n, epoch_n)

    Threads.@threads for epoch_idx in 1:epoch_n
        sged[:, :, epoch_idx], ress[:, epoch_idx], ress_normalized[:, epoch_idx] = @views s2_ged(signal1[:, :, epoch_idx], signal2[:, :, epoch_idx])
    end

    return (sged=sged, ress=ress, ress_normalized=ress_normalized)
end

"""
    eeg_frqinst(eeg)

Calculate instantaneous frequency of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `frqinst::Array{Float64, 3}`
"""
function eeg_frqinst(eeg::NeuroAnalyzer.EEG)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    frqinst = similar(signal)
    fs = eeg_sr(eeg)

    verbose == true && @info "eeg_frqinst() uses Hilbert transform, the signal should be narrowband for best results."

    Threads.@threads for epoch_idx in 1:epoch_n
       @inbounds @simd  for channel_idx in 1:channel_n
            frqinst[channel_idx, :, epoch_idx] = @views s_frqinst(signal[channel_idx, :, epoch_idx], fs=fs)
        end
    end

    return frqinst
end

"""
    eeg_itpc_s(eeg; <keyword arguments>)

Calculate spectrogram of ITPC (Inter-Trial-Phase Clustering) for `channel` of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64`
- `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
- `frq_n::Int64`: number of frequencies
- `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
- `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation

# Returns

Named tuple containing:
- `itpc_s::Array{Float64, 3}`: spectrogram of ITPC values
- `itpc_z_s::Array{Float64, 3}`: spectrogram ITPCz values
- `itpc_frq::Vector{Float64}`: frequencies list
"""
function eeg_itpc_s(eeg::NeuroAnalyzer.EEG; channel::Int64, frq_lim::Tuple{Real, Real}, frq_n::Int64, frq::Symbol=:log, w::Union{Vector{<:Real}, Nothing}=nothing)

    frq in [:log, :lin] || throw(ArgumentError("frq must be :log or :lin."))
    frq_lim = tuple_order(frq_lim)
    frq_lim[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
    frq_lim[2] > eeg_sr(eeg) ÷ 2 && throw(ArgumentError("Upper frequency bound must be ≤ $(eeg_sr(eeg) ÷ 2)."))
    frq_n < 2 && throw(ArgumentError("frq_n frequency bound must be ≥ 2."))
    if frq === :log
        frq_lim[1] == 0 && (frq_lim = (0.01, frq_lim[2]))
        frq_lim = (frq_lim[1], frq_lim[2])
        frq_list = logspace(log10(frq_lim[1]), log10(frq_lim[2]), frq_n)
    else
        frq_list = linspace(frq_lim[1], frq_lim[2], frq_n)
    end

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    _check_channels(channels, channel)
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)
    epoch_n < 2 && throw(ArgumentError("eeg must contain ≥ 2 epochs."))
    epoch_len = size(signal, 2)

    itpc_s = zeros(frq_n, epoch_len)
    itpc_z_s = zeros(frq_n, epoch_len)

    # initialize progress bar
    progress_bar == true && (p = Progress(frq_n, 1))

    Threads.@threads for frq_idx in 1:frq_n
        kernel = generate_morlet(eeg_sr(eeg), frq_list[frq_idx], 1, ncyc=10)
        half_kernel = floor(Int64, length(kernel) / 2) + 1
        s_conv = zeros(Float32, 1, epoch_len, epoch_n)
        @inbounds @simd for epoch_idx in 1:epoch_n
            s_conv[1, :, epoch_idx] = @views conv(signal[channel, :, epoch_idx], kernel)[(half_kernel - 1):(end - half_kernel)]
        end

        @inbounds @simd for t_idx in 1:epoch_len
            itpc, itpc_z, _, _ = s_itpc(s_conv, t=t_idx, w=w)
            itpc_s[frq_idx, t_idx] = itpc
            itpc_z_s[frq_idx, t_idx] = itpc_z
        end

        # update progress bar
        progress_bar == true && next!(p)
    end

    return (itpc_s=itpc_s, itpc_z_s=itpc_z_s, itpc_f=frq_list)
end

"""
    eeg_wspectrogram(eeg; pad, norm, frq_lim, frq_n, frq, ncyc, demean)

Return spectrogram of `eeg` using Morlet wavelet convolution.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `pad::Int64`: pad the `signal` with `pad` zeros
- `norm::Bool`=true: normalize powers to dB
- `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
- `frq_n::Int64`: number of frequencies
- `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
- `fs::Int64`: sampling rate
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n) for frq === :log or ncyc = linspace(ncyc[1], ncyc[2], frq_n) for frq === :lin
- `demean::Bool`=true: demean signal prior to analysis

# Returns

Named tuple containing:
- `w_pow::Array{Float64, 4}`
- `w_frq::Matrix{Float64}`
- `w_t::Matrix{Float64}`
"""
function eeg_wspectrogram(eeg::NeuroAnalyzer.EEG; pad::Int64=0, norm::Bool=true, frq_lim::Tuple{Real, Real}, frq_n::Int64, frq::Symbol=:lin, ncyc::Union{Int64, Tuple{Int64, Int64}}=6, demean::Bool=true)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    fs = eeg_sr(eeg)
    _, p_tmp, _, w_frq = s_wspectrogram(signal[1, :, 1], pad=pad, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc, demean=demean)
    w_pow = zeros(size(p_tmp, 1), size(p_tmp, 2), channel_n, epoch_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(epoch_n, 1))

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            _, w_pow[:, :, channel_idx, epoch_idx], _, _ = @views s_wspectrogram(signal[channel_idx, :, epoch_idx], pad=pad, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc, demean=demean)
        end

        # update progress bar
        progress_bar == true && next!(p)
    end

    return (w_pow=w_pow, w_frq=round.(w_frq, digits=2), w_t=eeg.eeg_epochs_time)
end

"""
    eeg_tkeo(eeg)

Calculate Teager-Kaiser energy-tracking operator: y(t) = x(t)^2 - x(t-1) × x(t+1)

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `tkeo::Array{Float64, 3}`
"""
function eeg_tkeo(eeg::NeuroAnalyzer.EEG)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    tkeo = similar(signal)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            tkeo[channel_idx, :, epoch_idx] = @views s_tkeo(signal[channel_idx, :, epoch_idx])
        end
    end

    return tkeo
end

"""
    eeg_wspectrum(eeg; pad, norm, frq_lim, frq_n, frq, ncyc)

Return power spectrogrum of `eeg` using Morlet wavelet convolution.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `pad::Int64`: pad the `signal` with `pad` zeros
- `norm::Bool`=true: normalize powers to dB
- `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
- `frq_n::Int64`: number of frequencies
- `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
- `fs::Int64`: sampling rate
- `ncyc::Int64=6`: number of cycles for Morlet wavelet

# Returns

Named tuple containing:
- `w_pow::Array{Float64, 4}`
- `w_frq::Matrix{Float64}`
"""
function eeg_wspectrum(eeg::NeuroAnalyzer.EEG; pad::Int64=0, norm::Bool=true, frq_lim::Tuple{Real, Real}, frq_n::Int64, frq::Symbol=:lin, ncyc::Int64=6)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    fs = eeg_sr(eeg)
    p_tmp, f_tmp = @views s_wspectrum(signal[1, :, 1], fs=fs, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)
    w_pow = zeros(length(p_tmp), channel_n, epoch_n)
    w_frq = zeros(length(f_tmp), epoch_n)
    verbose == true && @info "This will take a while.."
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            w_pow[:, channel_idx, epoch_idx], w_frq[:, epoch_idx] = @views s_wspectrum(signal[channel_idx, :, epoch_idx], pad=pad, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)
        end
    end

    return (w_pow=w_pow, w_frq=w_frq)
end

"""
    eeg_fcoherence(eeg1, eeg2; channel1, channel2, epoch1, epoch2, frq_lim)

Calculate coherence (mean over frequencies) and MSC (magnitude-squared coherence) between `eeg1` and `eeg2`.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2::NeuroAnalyzer.EEG`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=0`: default use all channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=0`: default use all channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=0`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=0`: default use all epochs
- `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: return coherence only for the given frequency range

# Returns

Named tuple containing:
- `c::Array{Float64, 3}`: coherence
- `msc::Array{Float64, 3}`: MSC
- `f::Vector{Float64}`: frequencies
"""
function eeg_fcoherence(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Union{Int64, Vector{Int64}, AbstractRange}=0, channel2::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch1::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch2::Union{Int64, Vector{Int64}, AbstractRange}=0, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)

    # select channels, default is all
    channels1 = eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))
    channels2 = eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))
    channel1 == 0 && (channel1 = channels1)
    _check_channels(channels1, channel1)
    channel2 == 0 && (channel2 = channels2)
    _check_channels(channels2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    # select epochs, default is all
    epoch1 == 0 && (epoch1 = _select_epochs(eeg1, epoch1))
    _check_epochs(eeg1, epoch1)
    epoch2 == 0 && (epoch2 = _select_epochs(eeg2, epoch2))
    _check_epochs(eeg2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    eeg_epoch_len(eeg1) == eeg_epoch_len(eeg2) || throw(ArgumentError("eeg1 and eeg2 epoch lengths must be equal."))

    eeg_sr(eeg1) == eeg_sr(eeg2) || throw(ArgumentError("EEG1 and EEG2 must have the same sampling rate."))

    c_tmp, _, f = @views s2_fcoherence(eeg1.eeg_signals[1, :, 1], eeg1.eeg_signals[1, :, 1], fs=eeg_sr(eeg1), frq_lim=frq_lim)
    c = zeros(length(channel1), length(c_tmp), length(epoch1))
    msc = zeros(length(channel1), length(c_tmp), length(epoch1))
    f = zeros(length(channel1), length(c_tmp), length(epoch1))
    @inbounds @simd for epoch_idx in 1:length(epoch1)
        Threads.@threads for channel_idx in 1:length(channel1)
            c[channel_idx, :, epoch_idx], msc[channel_idx, :, epoch_idx], _ = @views s2_fcoherence(eeg1.eeg_signals[channel1[channel_idx], :, epoch1[epoch_idx]], eeg2.eeg_signals[channel2[channel_idx], :, epoch2[epoch_idx]], fs=eeg_sr(eeg1), frq_lim=frq_lim)
        end
    end

    return (c=c, msc=msc, f=f)
end

"""
    eeg_vartest(eeg)

Calculate variance F-test for all channels of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

Named tuple containing:
- `f::Array{Float64, 3}`
- `p::Array{Float64, 3}`
"""
function eeg_vartest(eeg::NeuroAnalyzer.EEG)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    f = zeros(channel_n, channel_n, epoch_n)
    p = zeros(channel_n, channel_n, epoch_n)
    @inbounds @simd for epoch_idx in 1:epoch_n
       Threads.@threads for channel_idx1 in 1:channel_n
           for channel_idx2 in 1:channel_idx1
                ftest = @views VarianceFTest(signal[channel_idx1, :, epoch_idx], signal[channel_idx2, :, epoch_idx])
                f[channel_idx1, channel_idx2, epoch_idx] = ftest.F
                p[channel_idx1, channel_idx2, epoch_idx] = pvalue(ftest)
            end
        end
        Threads.@threads for channel_idx1 in 1:(channel_n - 1)
            for channel_idx2 in (channel_idx1 + 1):channel_n
                f[channel_idx1, channel_idx2, epoch_idx] = @views f[channel_idx2, channel_idx1, epoch_idx]
                p[channel_idx1, channel_idx2, epoch_idx] = @views p[channel_idx2, channel_idx1, epoch_idx]
            end
        end
    end

    return (f=f, p=p)
end

"""
    eeg_vartest(eeg1, eeg2)

Calculate variance F-test for all channels of `eeg1` and `eeg2`.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2::NeuroAnalyzer.EEG`

# Returns

Named tuple containing:
- `f::Array{Float64, 3}`
- `p::Array{Float64, 3}`
"""
function eeg_vartest(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG)

    channels1 = eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))
    signal1 = @view eeg1.eeg_signals[channels1, :, :]
    channel_n1 = size(signal1, 1)
    epoch_n1 = size(signal1, 3)

    channels2 = eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))
    signal2 = @view eeg2.eeg_signals[channels2, :, :]
    channel_n2 = size(signal2, 1)
    epoch_n2 = size(signal2, 3)

    channel_n1 == channel_n2 || throw(ArgumentError("Both EEG objects must have the same number of channels."))
    epoch_n1 == epoch_n2 || throw(ArgumentError("Both EEG objects must have the same number of epochs."))

    f = zeros(channel_n1, channel_n1, epoch_n1)
    p = zeros(channel_n1, channel_n1, epoch_n1)

    Threads.@threads for epoch_idx in 1:epoch_n1
       @inbounds @simd for channel_idx1 in 1:channel_n1
           for channel_idx2 in 1:channel_n1
                ftest = @views VarianceFTest(signal1[channel_idx1, :, epoch_idx], signal2[channel_idx2, :, epoch_idx])
                f[channel_idx1, channel_idx2, epoch_idx] = ftest.F
                p[channel_idx1, channel_idx2, epoch_idx] = pvalue(ftest)
            end
        end
    end

    return (f=f, p=p)
end

"""
    eeg_band_mpower(eeg; f, mt)

Calculate mean and maximum band power and its frequency.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

Named tuple containing:
- `mbp::Matrix{Float64}`: mean band power [μV^2/Hz] per channel per epoch
- `maxfrq::Matrix{Float64}`: frequency of maximum band power [Hz] per channel per epoch
- `maxbp::Matrix{Float64}`: power at maximum band frequency [μV^2/Hz] per channel per epoch
"""
function eeg_band_mpower(eeg::NeuroAnalyzer.EEG; f::Tuple{Real, Real}, mt::Bool=false)

    fs = eeg_sr(eeg)
    length(f) != 2 && throw(ArgumentError("f must contain two frequencies."))
    f = tuple_order(f)
    f[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
    f[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be < $(fs / 2)."))

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    mbp = zeros(channel_n, epoch_n)
    maxfrq = zeros(channel_n, epoch_n)
    maxbp = zeros(channel_n, epoch_n)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            mbp[channel_idx, epoch_idx], maxfrq[channel_idx, epoch_idx], maxbp[channel_idx, epoch_idx] = @views s_band_mpower(signal[channel_idx, :, epoch_idx], fs=fs, f=f, mt=mt)
        end
    end

    return (mbp=mbp, maxfrq=maxfrq, maxbp=maxbp)
end

"""
    eeg_rel_psd(eeg; norm, mt, f)

Calculate relative power spectrum density for each `eeg` channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `f::Union{Tuple{Real, Real}, Nothing}=nothing`: calculate power relative to frequency range or total power

# Returns

Named tuple containing:
- `psd_pow::Array{Float64, 3}`:powers
- `psd_frq::Array{Float64, 3}`: frequencies
"""
function eeg_rel_psd(eeg::NeuroAnalyzer.EEG; norm::Bool=false, mt::Bool=false, f::Union{Tuple{Real, Real}, Nothing}=nothing)

    fs = eeg_sr(eeg)
    if f !== nothing
        length(f) != 2 && throw(ArgumentError("f must contain two frequencies."))
        f = tuple_order(f)
        f[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
        f[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be < $(fs / 2)."))
    end

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    psd_tmp, psd_frq = @views s_rel_psd(signal[1, :, 1], fs=fs, norm=norm, mt=mt, f=f)
    psd_pow = zeros(channel_n, length(psd_tmp), epoch_n)
    @inbounds @simd for epoch_idx in 1:epoch_n
        for channel_idx in 1:channel_n
            psd_pow[channel_idx, :, epoch_idx], _ = @views s_rel_psd(signal[channel_idx, :, epoch_idx], fs=fs, norm=norm, mt=mt, f=f)
        end
    end

    return (psd_pow=psd_pow, psd_frq=psd_frq)
end

"""
    eeg_fbsplit(eeg; order)

Split EEG signal into frequency bands.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `order::Int64=8`: bandpass filter order

# Returns

Named tuple containing:
- `band_names::Vector{Symbol}`
- `band_frq::Vector{Tuple{Real, Real}}`
- `signal_split::Array{Float64, 4}`
"""
function eeg_fbsplit(eeg::NeuroAnalyzer.EEG; order::Int64=8)
    
    band = [:delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher]

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    fs = eeg_sr(eeg)
    signal_split = zeros(length(band), channel_n, size(signal, 2), epoch_n)
    band_frq = Vector{Tuple{Real, Real}}()
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for band_idx in 1:length(band)
            band_f = eeg_band(eeg, band=band[band_idx])
            push!(band_frq, band_f)
            for channel_idx in 1:channel_n
                signal_split[band_idx, channel_idx, :, epoch_idx] = @views s_filter(signal[channel_idx, :, epoch_idx], fs=fs, fprototype=:butterworth, ftype=:bp, cutoff=band_f, order=8)
            end
        end
    end

    return (band_names=band, band_frq=band_frq, signal_split=signal_split)
end

"""
    eeg_chdiff(eeg1, eeg2; channel1, channel2)

Calculate difference between `channel1` of `eeg1` and `channel2` of `eeg2`.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2::NeuroAnalyzer.EEG`
- `channel1::Int64`
- `channel2::Int64`

# Returns

- `ch_diff::Matrix{Float64}`
"""
function eeg_chdiff(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Int64, channel2::Int64)

    channel1 < 0 || channel2 < 0 && throw(ArgumentError("channel1 and channel2 must be > 0."))
    channel1 > eeg_channel_n(eeg1) && throw(ArgumentError("channel1 must be ≤ $(eeg_channel_n(eeg1))."))
    channel2 > eeg_channel_n(eeg2) && throw(ArgumentError("channel2 must be ≤ $(eeg_channel_n(eeg2))."))
    size(eeg1.eeg_signals[channel1, :, :]) == size(eeg2.eeg_signals[channel2, :, :]) || throw(ArgumentError("Both EEG channels must have the same epoch length."))

    ch_diff = zeros(size(eeg1.eeg_signals[channel1, :, :]))
    @inbounds @simd for epoch_idx in 1:eeg_epoch_n(eeg1)
        ch_diff[:, epoch_idx] = @views eeg1.eeg_signals[channel1, :, epoch_idx] .- eeg2.eeg_signals[channel2, :, epoch_idx]
    end

    return ch_diff
end

"""
    eeg_cps(eeg; norm)

Calculate cross power spectrum between all channels of `eeg`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `norm::Bool=true`: normalize do dB

# Returns

Named tuple containing:
- `cps_pw::Array{Float64, 4}`: cross power spectrum power
- `cps_ph::Array{Float64, 4}`: cross power spectrum phase (in radians)
- `cps_fq::Vector{Float64}`: cross power spectrum frequencies
"""
function eeg_cps(eeg::NeuroAnalyzer.EEG; norm::Bool=true)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    fs = eeg_sr(eeg)
    
    cps_pw_tmp, cps_ph_tmp, cps_fq = @views s2_cps(signal[1, :, 1], signal[1, :, 1], fs=fs)
    cps_pw = zeros(channel_n, channel_n, length(cps_pw_tmp), epoch_n)
    cps_ph = zeros(channel_n, channel_n, length(cps_ph_tmp), epoch_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(epoch_n, 1))

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx1 in 1:channel_n
           for channel_idx2 in 1:channel_idx1
                cps_pw[channel_idx1, channel_idx2, :, epoch_idx], cps_ph[channel_idx1, channel_idx2, :, epoch_idx], _ = @views s2_cps(signal[channel_idx1, :, epoch_idx], signal[channel_idx2, :, epoch_idx], fs=fs, norm=norm)
            end
        end
        # update progress bar
        progress_bar == true && next!(p)
    end
    @inbounds @simd for time_idx in 1:size(cps_pw, 3)
        Threads.@threads for epoch_idx in 1:epoch_n
            for channel_idx1 in 1:(channel_n - 1)
                for channel_idx2 in (channel_idx1 + 1):channel_n
                    cps_pw[channel_idx1, channel_idx2, time_idx, epoch_idx] = @views cps_pw[channel_idx2, channel_idx1, time_idx, epoch_idx]
                    cps_ph[channel_idx1, channel_idx2, time_idx, epoch_idx] = @views cps_ph[channel_idx2, channel_idx1, time_idx, epoch_idx]
                end
            end
        end
    end

    return (cps_pw=cps_pw, cps_ph=cps_ph, cps_fq=cps_fq)
end

"""
    eeg_cps(eeg1, eeg2; channel1, channel2, epoch1, epoch2, norm)

Calculate cross power spectrum between `eeg` channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel1::Int64`
- `channel2::Int64`
- `epoch1::Int64`
- `epoch2::Int64`
- `norm::Bool=true`: normalize do dB

# Returns

Named tuple containing:
- `cps_pw::Vector{Float64}`: cross power spectrum power
- `cps_ph::Vector{Float64}`: cross power spectrum phase (in radians)
- `cps_fq::Vector{Float64}`: cross power spectrum frequencies
"""
function eeg_cps(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Int64, channel2::Int64, epoch1::Int64, epoch2::Int64, norm::Bool=true)

    eeg_sr(eeg1) == eeg_sr(eeg2) || throw(ArgumentError("EEG1 and EEG2 must have the same sampling rate."))
    (channel1 < 0 || channel2 < 0 || epoch1 < 0 || epoch2 < 0) && throw(ArgumentError("channel1/epoch1/channel2/epoch2 must be > 0."))

    channels1 = eeg_channel_idx(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))
    signal1 = @view eeg1.eeg_signals[channels1, :, :]
    channel_n1 = size(signal1, 1)
    epoch_n1 = size(signal1, 3)

    channels2 = eeg_channel_idx(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))
    signal2 = @view eeg2.eeg_signals[channels2, :, :]
    channel_n2 = size(signal2, 1)
    epoch_n2 = size(signal2, 3)

    (channel1 > channel_n1) && throw(ArgumentError("channel1 must be ≤ $(channel_n1)."))
    (epoch1 > epoch_n1) && throw(ArgumentError("epoch1 must be ≤ $(epoch_n1)."))
    (channel2 > channel_n2) && throw(ArgumentError("channel2 must be ≤ $(channel_n2)."))
    (epoch2 > epoch_n2) && throw(ArgumentError("epoch2 must be ≤ $(epoch_n2)."))
    fs = eeg_sr(eeg1)
    cps_pw, cps_ph, cps_fq = @views s2_cps(signal1[channel1, :, epoch1], signal2[channel2, :, epoch2], fs=fs, norm=norm)

    return (cps_pw=cps_pw, cps_ph=cps_ph, cps_fq=cps_fq)
end

"""
    eeg_phdiff(eeg; channel, pad, h)

Calculate phase difference between each `eeg` channel and mean phase of `channel`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: reference channels, default is all channels except the analyzed one
- `avg::Symbol=:phase`: method of averaging: `:phase` or `:signal`; for :signal `channel` signals are averaged prior to phase calculation; for :phase phase is calculated for each reference channel separately and then averaged
- `pad::Int64=0`: pad signals with 0s
- `h::Bool=false`: use FFT or Hilbert transformation (if h=true)

# Returns
 
- `ph_diff::Array{Float64, 3}`
"""
function eeg_phdiff(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=0, avg::Symbol=:phase, pad::Int64=0, h::Bool=false)

    avg in [:phase, :signal] || throw(ArgumentError("avg must be :phase or :signal."))

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    # select channels, default is all
    channel == 0 && (channel = channels)
    _check_channels(channels, channel)

    ph_diff = similar(signal)
    if avg === :phase
        @inbounds @simd for epoch_idx in 1:epoch_n
            Threads.@threads for channel_idx in 1:channel_n
                ref_channels = setdiff(channel, channel_idx)
                    ph_ref = zeros(length(ref_channels), size(signal, 2))
                    for ref_idx in 1:length(ref_channels)
                        if h
                            _, _, _, ph = @views s_hspectrum(signal[ref_channels[ref_idx], :, epoch_idx], pad=pad)
                        else
                            _, _, _, ph = @views s_spectrum(signal[ref_channels[ref_idx], :, epoch_idx], pad=pad)
                        end
                        ph_ref[ref_idx, :] = ph
                    end
                    ph_ref = vec(mean(ph_ref, dims=1))
                    if h
                        _, _, _, ph = @views s_hspectrum(signal[channel_idx, :, epoch_idx], pad=pad)
                    else
                        _, _, _, ph = @views s_spectrum(signal[channel_idx, :, epoch_idx], pad=pad)
                    end
                    ph_diff[channel_idx, :, epoch_idx] = ph - ph_ref
            end
        end
    else
        @inbounds @simd for epoch_idx in 1:epoch_n
            Threads.@threads for channel_idx in 1:channel_n
                ref_channels = setdiff(channel, channel_idx)
                signal_m = @views vec(mean(signal[ref_channels, :, epoch_idx], dims=1))
                ph_diff[channel_idx, :, epoch_idx] = @views s_phdiff(signal[channel_idx, :, epoch_idx], signal_m)
            end
        end
    end

    return ph_diff
end

"""
    eeg_ampdiff(eeg; channel)

Calculate amplitude difference between each `eeg` channel and mean amplitude of `channel`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=0`: reference channels, default is all channels except the analyzed one

# Returns
 
- `amp_diff::Array{Float64, 3}`
"""
function eeg_ampdiff(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=0)

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    # select channels, default is all
    channel == 0 && (channel = channels)
    _check_channels(channels, channel)

    amp_diff = similar(signal)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            ref_channels = setdiff(channel, channel_idx)
            amp_ref = @views vec(mean(signal[ref_channels, :, epoch_idx], dims=1))
            amp_diff[channel_idx, :, epoch_idx] = @views signal[channel_idx, :, epoch_idx] - amp_ref
        end
    end

    return amp_diff
end

"""
    eeg_dwt(eeg; wt, type, l)

Perform discrete wavelet transformation (DWT) of each `eeg` channel.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: transformation type: Stationary Wavelet Transforms (:sdwt) or Autocorrelation Wavelet Transforms (:acdwt)
- `l::Int64=0`: number of levels, default maximum number of levels available or total transformation

# Returns
 
- `dwt_c::Array{Float64, 4}`: DWT coefficients cAl, cD1, ..., cDl (by rows)
"""
function eeg_dwt(eeg::NeuroAnalyzer.EEG; wt::T, type::Symbol, l::Int64=0) where {T <: DiscreteWavelet}

    if l == 0
        l = maxtransformlevels(eeg.eeg_signals[1, :, 1])
        verbose == true && @info "Calculating DWT using maximum level: $l."
    end

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)
    dwt_c = zeros(channel_n, (l + 1), size(signal, 2), epoch_n)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            dwt_c[channel_idx, :, :, epoch_idx] = @views s_dwt(signal[channel_idx, :, epoch_idx], wt=wt, type=type, l=l)
        end
    end

    return dwt_c
end

"""
    eeg_cwt(eeg; wt)

Perform continuous wavelet transformation (CWT) of each `eeg` channel.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns
 
- `cwt_c::Array{Float64, 4}`: CWT coefficients (by rows)
"""
function eeg_cwt(eeg::NeuroAnalyzer.EEG; wt::T) where {T <: CWT}

    l = size(ContinuousWavelets.cwt(eeg.eeg_signals[1, :, 1], wt), 2)
    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)
    cwt_c = zeros(channel_n, l, size(signal, 2), epoch_n)
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            cwt_c[channel_idx, :, :, epoch_idx] = @views s_cwt(signal[channel_idx, :, epoch_idx], wt=wt)
        end
    end

    return cwt_c
end

"""
    eeg_psdslope(eeg; f, norm, mt)

Calculate PSD linear fit and slope.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `f::Union{Real, Real}=(0, 0)`: calculate slope of the total power (default) or frequency range f[1] to f[2]
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

Named tuple containing:
- `lf::Array{Float64, 3}`: linear fit for each channel and epoch
- `psd_slope::Array{Float64, 2}`: slopes of each linear fit
- `frq::Vector{Float64}`: range of frequencies for the linear fits
"""
function eeg_psdslope(eeg::NeuroAnalyzer.EEG; f::Tuple{Real, Real}=(0, 0), norm::Bool=false, mt::Bool=false)

    fs = eeg_sr(eeg)
    length(f) != 2 && throw(ArgumentError("f must contain two frequencies."))
    f == (0, 0) && (f = (0, fs/2))
    f = tuple_order(f)
    f[1] < 0 && throw(ArgumentError("Lower frequency bound must be be ≥ 0."))
    f[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be be < $(fs / 2)."))

    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))
    signal = @view eeg.eeg_signals[channels, :, :]
    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)

    _, frq = s_psd(signal[1, :, 1], fs=fs, norm=norm, mt=mt)
    f1_idx = vsearch(f[1], frq)
    f2_idx = vsearch(f[2], frq)
    lf = zeros(channel_n, length(frq[f1_idx:f2_idx]), epoch_n)
    psd_slope = zeros(channel_n, epoch_n)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            pow, _ = s_psd(signal[channel_idx, :, epoch_idx], fs=fs, norm=norm, mt=mt)
            _, _, _, _, _, _, lf[channel_idx, :, epoch_idx] = @views linreg(frq[f1_idx:f2_idx], pow[f1_idx:f2_idx])
            psd_slope[channel_idx, epoch_idx] = lf[channel_idx, 2, epoch_idx] - lf[channel_idx, 1, epoch_idx]
        end
    end

    return (lf=lf, psd_slope=psd_slope, frq=frq[f1_idx:f2_idx])
end