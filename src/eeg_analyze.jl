"""
    eeg_xcov(eeg; channel, lag, demean, norm)

Calculate cross-covariance.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `lag::Int64=1`: lags range is `-lag:lag`
- `demean::Bool=false`: demean signal prior to analysis
- `norm::Bool=false`: normalize cross-covariance

# Returns

Named tuple containing:
- `xcov::Matrix{Float64}`: ch1-ch1, ch1-ch2, ch1-ch3, etc.
- `lags::Vector{Float64}`: lags in ms
"""
function eeg_xcov(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), lag::Int64=1, demean::Bool=false, norm::Bool=false)

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    # create vector of lags
    lags = 1/eeg_sr(eeg) .* collect(-lag:lag) .* 1000

    xcov = zeros(ch_n^2, length(lags), ep_n)
    @inbounds @simd for epoch_idx in 1:ep_n
        
        # create half of the covariance matrix
        xcov_packed = Array{Vector{Float64}}(undef, ch_n, ch_n)
        Threads.@threads for channel_idx1 in 1:ch_n
            for channel_idx2 in 1:channel_idx1
                xcov_packed[channel_idx1, channel_idx2], _ = @views s2_xcov(eeg.eeg_signals[channel[channel_idx1], :, epoch_idx], eeg.eeg_signals[channel[channel_idx2], :, epoch_idx], lag=lag, demean=demean, norm=norm)
            end
        end
        
        # copy to the other half
        Threads.@threads for channel_idx1 in 1:(ch_n - 1)
            for channel_idx2 in (channel_idx1 + 1):ch_n
                xcov_packed[channel_idx1, channel_idx2] = @views xcov_packed[channel_idx2, channel_idx1]
            end
        end

        # unpack by channels
        Threads.@threads for channel_idx in 1:ch_n^2
            xcov[channel_idx, :, epoch_idx] = @views xcov_packed[channel_idx]
        end
    end

    return (xcov=xcov, lags=lags)
end

"""
    eeg_xcov(eeg1, eeg2; channel1, channel2, epoch1, epoch2, lag, demean, norm)

Calculate cross-covariance between two EEG signals.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2::NeuroAnalyzer.EEG`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs
- `lag::Int64=1`: lags range is `-lag:lag`
- `demean::Bool=false`: demean signal prior to analysis
- `norm::Bool=false`: normalize cross-covariance

# Returns

Named tuple containing:
- `xcov::Array{Float64, 3}`
- `lags::Vector{Float64}`: lags in ms
"""
function eeg_xcov(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2)), lag::Int64=1, demean::Bool=false, norm::Bool=false)

    # check channels
    _check_channels(eeg1, channel1)
    _check_channels(eeg2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    # check epochs
    _check_epochs(eeg1, epoch1)
    _check_epochs(eeg2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    eeg_epoch_len(eeg1) == eeg_epoch_len(eeg2) || throw(ArgumentError("eeg1 and eeg2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)

    # create vector of lags
    lags = 1/eeg_sr(eeg1) .* collect(-lag:lag) .* 1000

    xcov = zeros(length(channel1), (2 * lag + 1), length(epoch1))
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            xcov[channel_idx, :, epoch_idx], _ = @views s2_xcov(eeg1.eeg_signals[channel1[channel_idx], :, epoch1[epoch_idx]], eeg2.eeg_signals[channel2[channel_idx], :, epoch2[epoch_idx]], lag=lag, demean=demean, norm=norm)
        end
    end

    return (xcov=xcov, lags=lags)
end

"""
    eeg_psd(eeg; channel, norm, mt)

Calculate power spectrum density.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns

Named tuple containing:
- `psd_pow::Array{Float64, 3}`:powers
- `psd_frq::Vector{Float64}`: frequencies
"""
function eeg_psd(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), norm::Bool=false, mt::Bool=false, nt::Int64=8)

    fs = eeg_sr(eeg)

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)
    
    _, psd_frq = s_psd(eeg.eeg_signals[1, :, 1], fs=fs, norm=norm, mt=mt, nt=nt)
    psd_pow = zeros(length(channel), length(psd_frq), ep_n)
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            psd_pow[channel_idx, :, epoch_idx], _ = s_psd(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], fs=fs, norm=norm, mt=mt, nt=nt)
        end
    end

    return (psd_pow=psd_pow, psd_frq=psd_frq)
end

"""
    eeg_stationarity(eeg; channel, window, method)

Calculate stationarity.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `window::Int64=10`: time window in samples
- `method::Symbol=:euclid`: stationarity method:
    - `:mean`: mean across `window`-long windows
    - `:var`: variance across `window`-long windows
    - `:cov`: covariance stationarity based on Euclidean distance between covariance matrix of adjacent time windows
    - `:hilbert`: phase stationarity using Hilbert transformation
    - `:adf`: Augmented Dickey–Fuller test; returns ADF-test value and p-value (H0: signal is non-stationary; p-value < alpha means that signal is stationary)

# Returns

- `stationarity::Array{Float64, 3}`
"""
function eeg_stationarity(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), window::Int64=10, method::Symbol=:hilbert)

    _check_var(method, [:mean, :var, :cov, :hilbert, :adf], "method")
    window < 1 && throw(ArgumentError("window must be ≥ 1."))
    window > eeg_epoch_len(eeg) && throw(ArgumentError("window must be ≤ $(eeg_epoch_len(eeg))."))

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    if method === :mean
        s_stationarity = zeros(ch_n, window, ep_n)
        @inbounds @simd for epoch_idx in 1:ep_n
            Threads.@threads for channel_idx in 1:ch_n
                s_stationarity[channel_idx, :, epoch_idx] = @views s_stationarity_mean(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], window=window)
            end
        end
    end

    if method === :var
        s_stationarity = zeros(ch_n, window, ep_n)
        @inbounds @simd for epoch_idx in 1:ep_n
            Threads.@threads for channel_idx in 1:ch_n
                s_stationarity[channel_idx, :, epoch_idx] = @views s_stationarity_var(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], window=window)
            end
        end
    end

    if method === :hilbert
        s_stationarity = zeros(ch_n, eeg_epoch_len(eeg) - 1, ep_n)
        @inbounds @simd for epoch_idx in 1:ep_n
            Threads.@threads for channel_idx in 1:ch_n
                s_stationarity[channel_idx, :, epoch_idx] = @views s_stationarity_hilbert(eeg.eeg_signals[channel[channel_idx], :, epoch_idx])
            end
        end
    end

    if method === :cov
        # number of time windows per epoch
        window_n = eeg_epoch_len(eeg)
        cov_mat = zeros(ch_n, ch_n, window_n, ep_n)
        s_stationarity = zeros(1 + length(2:window:window_n), ep_n)
        ch_n == 1 && throw(ArgumentError("For :cov method, number of channels must be ≥ 2."))

        # create covariance matrices per each window
        @inbounds @simd for epoch_idx in 1:ep_n
            Threads.@threads for window_idx = 1:window_n
                cov_mat[:, :, window_idx, epoch_idx] = @views s2_cov(eeg.eeg_signals[channel, window_idx, epoch_idx], eeg.eeg_signals[channel, window_idx, epoch_idx])
            end
        end

        # calculate Euclidean distance between adjacent matrices
        @inbounds @simd for epoch_idx in 1:ep_n
            w_idx = 1
            Threads.@threads for window_idx = 2:window:window_n
                s_stationarity[w_idx, epoch_idx] = @views euclidean(cov_mat[:, :, window_idx - 1, epoch_idx], cov_mat[:, :, window_idx, epoch_idx])
                w_idx += 1
            end
        end
    end

    if method === :adf
        s_stationarity = zeros(ch_n, 2, ep_n)

        # initialize progress bar
        progress_bar == true && (pb = Progress(ep_n * ch_n, 1))

        # perform Augmented Dickey–Fuller test
        @inbounds @simd for epoch_idx in 1:ep_n
            Threads.@threads for channel_idx = 1:ch_n
                adf = @views ADFTest(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], :constant, window)
                a = adf.stat
                p = pvalue(adf)
                p < eps() && (p = 0.0001)
                a = round(a, digits=2)
                p = round(p, digits=4)
                p == 0.0 && (p = 0.0001)
                s_stationarity[channel_idx, :, epoch_idx] = [a, p]

                # update progress bar
                progress_bar == true && next!(pb)
            end
        end
    end

    return s_stationarity
end

"""
    eeg_mi(eeg; channel)

Calculate mutual information between EEG channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels

# Returns

- `mi::Array{Float64, 3}`
"""
function eeg_mi(eeg::NeuroAnalyzer.EEG; channel::Union{Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg))

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    mi = zeros(ch_n, ch_n, ep_n)
    @inbounds @simd for epoch_idx in 1:ep_n
        
        # create half of the matrix
        Threads.@threads for channel_idx1 in 1:ch_n
            for channel_idx2 in 1:channel_idx1
                mi[channel_idx1, channel_idx2, epoch_idx] = @views s2_mi(eeg.eeg_signals[channel[channel_idx1], :, epoch_idx], eeg.eeg_signals[channel[channel_idx2], :, epoch_idx])
            end
        end

        # copy to the other half
        Threads.@threads for channel_idx1 in 1:(ch_n - 1)
            for channel_idx2 in (channel_idx1 + 1):ch_n
                mi[channel_idx1, channel_idx2, epoch_idx] = @views mi[channel_idx2, channel_idx1, epoch_idx]
            end
        end
    end

    return mi
end

"""
    eeg_mi(eeg1, eeg2; channel1, channel2, epoch1, epoch2)

Calculate mutual information between two EEG channels.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2::NeuroAnalyzer.EEG`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs

# Returns

- `mi::Array{Float64, 3}`
"""
function eeg_mi(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2)))

    # check channels
    _check_channels(eeg1, channel1)
    _check_channels(eeg2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    # check epochs
    _check_epochs(eeg1, epoch1)
    _check_epochs(eeg2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    eeg_epoch_len(eeg1) == eeg_epoch_len(eeg2) || throw(ArgumentError("eeg1 and eeg2 epoch lengths must be equal."))

    ch_n = length(channel1)
    ep_n = length(epoch1)

    mi = zeros(ch_n, ch_n, ep_n)
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx1 in 1:ch_n
            for channel_idx2 in 1:ch_n
                mi[channel_idx1, channel_idx2, epoch_idx] = @views s2_mi(eeg1.eeg_signals[channel1[channel_idx1], :, epoch1[epoch_idx]], eeg2.eeg_signals[channel2[channel_idx2], :, epoch2[epoch_idx]])
            end
        end
    end

    return mi
end

"""
    eeg_entropy(eeg; channel)

Calculate entropy.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels

# Returns

Named tuple containing:
- `ent::Array{Float64, 2}`
- `s_ent::Array{Float64, 2}`: Shanon entropy
- `le_ent::Array{Float64, 2}`: log energy entropy
"""
function eeg_entropy(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg))

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    ent = zeros(ch_n, ep_n)
    sent = zeros(ch_n, ep_n)
    leent = zeros(ch_n, ep_n)
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            ent[channel_idx, epoch_idx], sent[channel_idx, epoch_idx], leent[channel_idx, epoch_idx] = @views s_entropy(eeg.eeg_signals[channel[channel_idx], :, epoch_idx])
        end
    end

    return (ent=ent, s_ent=sent, le_ent=leent)
end

"""
    eeg_negentropy(eeg; channel)

Calculate negentropy.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels

# Returns

- `ne::Matrix{Float64}`
"""
function eeg_negentropy(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg))

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    ne = zeros(ch_n, ep_n)
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            ne[channel_idx, epoch_idx] = @views s_negentropy(eeg.eeg_signals[channel[channel_idx], :, epoch_idx])
        end
    end

    return ne
end

"""
    eeg_band(eeg, band)

Return frequency limits for a `band`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `band::Symbol`: band range name:
    - `:list`
    - `:total`
    - `:delta`
    - `:theta`
    - `:alpha`
    - `:beta`
    - `:beta_high`
    - `:gamma`
    - `:gamma_1`
    - `:gamma_2`
    - `:gamma_lower`
    - `:gamma_higher`.
# Returns

- `band_frequency::Tuple{Real, Real}`
"""
function eeg_band(eeg::NeuroAnalyzer.EEG; band::Symbol)

    bands = [:list, :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher]
    _check_var(band, [:list, :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher], "band")
    if band === :list
        print("Available band names: ")
        [print(":$(bands[x]), ") for x in 2:(length(bands) - 1)]
        println(":$(bands[end])")
        return nothing
    end

    band === :total && (band_frq = (0.1, round(eeg_sr(eeg) / 2, digits=1)))
    band === :delta && (band_frq = (0.1, 4.0))
    band === :theta && (band_frq = (4.0, 8.0))
    band === :alpha && (band_frq = (8.0, 13.0))
    band === :beta && (band_frq = (14.0, 30.0))
    band === :beta_high && (band_frq = (25.0, 30.0))
    band === :gamma && (band_frq = (30.0, 150.0))
    band === :gamma_1 && (band_frq = (31.0, 40.0))
    band === :gamma_2 && (band_frq = (41.0, 50.0))
    band === :gamma_lower && (band_frq = (30.0, 80.0))
    band === :gamma_higher && (band_frq = (80.0, 150.0))
    
    if band_frq[1] > eeg_sr(eeg) / 2
        _info("Nyquist frequency based on EEG sampling rate ($(eeg_sr(eeg) / 2)) is lower than $band range: $band_frq, band frequency truncated to: ($(eeg_sr(eeg) / 2 - 0.2), $(eeg_sr(eeg) / 2 - 0.1))")
        band_frq = (eeg_sr(eeg) / 2 - 0.2, eeg_sr(eeg) / 2 - 0.1)
    end
    if band_frq[2] > eeg_sr(eeg) / 2
        _info("Nyquist frequency based on EEG sampling rate ($(eeg_sr(eeg) / 2)) is lower than $band range: $band_frq, band frequency truncated to: ($(band_frq[1]), $(eeg_sr(eeg) / 2 - 0.1))")
        band_frq = (band_frq[1], eeg_sr(eeg) / 2 - 0.1)
    end

    return band_frq
end

"""
    eeg_band(fs, band)

Return frequency limits of a `band`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `band::Symbol`: band range name:
    - `:list`
    - `:total`
    - `:delta`
    - `:theta`
    - `:alpha`
    - `:beta`
    - `:beta_high`
    - `:gamma`
    - `:gamma_1`
    - `:gamma_2`
    - `:gamma_lower`
    - `:gamma_higher`.

# Returns

- `band_frequency::Tuple{Real, Real}`
"""
function eeg_band(fs::Int64; band::Symbol)

    bands = [:list, :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher]
    _check_var(band, [:list, :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher], "band")
    if band === :list
        print("Available band names: ")
        [print(":$(bands[x]), ") for x in 2:(length(bands) - 1)]
        println(":$(bands[end])")
        return nothing
    end

    band === :total && (band_frq = (0.1, round(fs / 2, digits=1)))
    band === :delta && (band_frq = (0.1, 4.0))
    band === :theta && (band_frq = (4.0, 8.0))
    band === :alpha && (band_frq = (8.0, 13.0))
    band === :beta && (band_frq = (14.0, 30.0))
    band === :beta_high && (band_frq = (25.0, 30.0))
    band === :gamma && (band_frq = (30.0, 150.0))
    band === :gamma_1 && (band_frq = (31.0, 40.0))
    band === :gamma_2 && (band_frq = (41.0, 50.0))
    band === :gamma_lower && (band_frq = (30.0, 80.0))
    band === :gamma_higher && (band_frq = (80.0, 150.0))
    
    if band_frq[1] > fs / 2
        _info("Nyquist frequency based on EEG sampling rate ($(fs / 2)) is lower than $band range: $band_frq, band frequency truncated to: ($(fs / 2 - 0.2), $(fs / 2 - 0.1))")
        band_frq = (fs / 2 - 0.2, fs / 2 - 0.1)
    end
    if band_frq[2] > fs / 2
        _info("Nyquist frequency based on EEG sampling rate ($(fs / 2)) is lower than $band range: $band_frq, band frequency truncated to: ($(band_frq[1]), $(fs / 2 - 0.1))")
        band_frq = (band_frq[1], fs / 2 - 0.1)
    end

    return band_frq
end

"""
    eeg_tcoherence(eeg1, eeg2; channel1, channel2, epoch1, epoch2)

Calculate coherence (mean over time) and MSC (magnitude-squared coherence).

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2::NeuroAnalyzer.EEG`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs
- `pad::Int64=0`: number of zeros to add signal for FFT

# Returns

Named tuple containing:
- `c::Array{Float64, 3}`: coherence
- `msc::Array{Float64, 3}`: MSC
- `ic::Array{Float64, 3}`: imaginary part of coherence
"""
function eeg_tcoherence(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2)), pad::Int64=0)

    _check_channels(eeg1, channel1)
    _check_channels(eeg2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))

    _check_epochs(eeg1, epoch1)
    _check_epochs(eeg2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    eeg_epoch_len(eeg1) == eeg_epoch_len(eeg2) || throw(ArgumentError("eeg1 and eeg2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)

    c = zeros(length(channel1), eeg_epoch_len(eeg1), length(epoch1))
    msc = zeros(length(channel1), eeg_epoch_len(eeg1), length(epoch1))
    ic = zeros(length(channel1), eeg_epoch_len(eeg1), length(epoch1))

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            c[channel_idx, :, epoch_idx], msc[channel_idx, :, epoch_idx], ic[channel_idx, :, epoch_idx] = @views s2_tcoherence(eeg1.eeg_signals[channel1[channel_idx], :, epoch1[epoch_idx]], eeg2.eeg_signals[channel2[channel_idx], :, epoch2[epoch_idx]], pad=pad)
        end
    end

    return (c=c, msc=msc, ic=ic)
end

"""
    eeg_freqs(eeg)

Return vector of frequencies and Nyquist frequency.

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
    eeg_difference(eeg; channel, n, method)

Calculate mean difference and its 95% CI between EEG channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
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
function eeg_difference(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), n::Int64=3, method::Symbol=:absdiff)

    ep_n = eeg_epoch_n(eeg)
    _check_channels(eeg, channel)

    s_stat = zeros(ep_n, length(channel) * n)
    s_stat_single = zeros(ep_n)
    p = zeros(ep_n)

    @inbounds @simd for epoch_idx in 1:ep_n
        s_stat[epoch_idx, :], s_stat_single[epoch_idx], p[epoch_idx] = s2_difference(eeg.eeg_signals[channel, :, epoch_idx], eeg.eeg_signals[channel, :, epoch_idx], n=n, method=method)
    end

    return (s_stat=s_stat, s_stat_single=s_stat_single, p=p)
end

"""
    eeg_channel_pick(eeg; pick)

Return set of channel indices corresponding with `pick` of electrodes

# Arguments

- `pick::Vector{Symbol}`: pick of electrodes; picks may be combined, e.g. `[:left, :frontal]`
    - `:list`
    - `:central` (or `:c`)
    - `:left` (or `:l`)
    - `:right` (or `:r`)
    - `:frontal` (or `:f`)
    - `:temporal` (or `:t`)
    - `:parietal` (or `:p`)
    - `:occipital` (or `:o`)

# Returns

- `channels::Vector{Int64}`: channel numbers
"""
function eeg_channel_pick(eeg::NeuroAnalyzer.EEG; pick::Union{Symbol, Vector{Symbol}})

    length(eeg_labels(eeg)) == 0 && throw(ArgumentError("EEG does not contain channel labels."))

    if typeof(pick) == Vector{Symbol}
        for idx in pick
            _check_var(idx, [:list, :central, :c, :left, :l, :right, :r, :frontal, :f, :temporal, :t, :parietal, :p, :occipital, :o], "pick")
        end

        # convert picks to channel labels
        c = Vector{Char}()
        for idx in pick
            (idx === :central || idx === :c) && push!(c, 'z')
            (idx === :frontal || idx === :f) && push!(c, 'F')
            (idx === :temporal || idx === :t) && push!(c, 'T')
            (idx === :parietal || idx === :p) && push!(c, 'P')
            (idx === :occipital || idx === :o) && push!(c, 'O')
        end
        
        # check which channels are in the picks list
        labels = eeg_labels(eeg)
        channels = Vector{Int64}()
        for idx1 in eachindex(labels)
            for idx2 in eachindex(c)
                in(c[idx2], labels[idx1]) && push!(channels, idx1)
            end
        end

        # check for both :l and :r
        for idx1 in eachindex(pick)
            if (pick[idx1] === :left || pick[idx1] === :l)
                for idx2 in eachindex(pick)
                    if (pick[idx2] === :right || pick[idx2] === :r)
                        return channels
                    end
                end
            end
            if (pick[idx1] === :right || pick[idx1] === :r)
                for idx2 in eachindex(pick)
                    if (pick[idx2] === :left || pick[idx2] === :l)
                        return channels
                    end
                end
            end
        end

        labels = eeg_labels(eeg)
        labels = labels[channels]
        pat = nothing
        for idx in pick
            # for :right remove lefts
            (idx === :right || idx === :r) && (pat = r"[z13579]$")
            # for :left remove rights
            (idx === :left || idx === :l) && (pat = r"[z02468]$")
        end
        if typeof(pat) == Regex
            for idx in length(labels):-1:1
                typeof(match(pat, labels[idx])) == RegexMatch && deleteat!(channels, idx)
            end
        end

        return channels
    else
        _check_var(pick, [:central, :c, :left, :l, :right, :r, :frontal, :f, :temporal, :t, :parietal, :p, :occipital, :o], "pick")

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
        for idx1 in eachindex(c)
            for idx2 in eachindex(labels)
                in(c[idx1], labels[idx2]) && push!(channels, idx2)
            end
        end

        return channels
    end
end

"""
    eeg_epoch_stats(eeg)

Calculate epochs statistics.

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
function eeg_epoch_stats(eeg::NeuroAnalyzer.EEG)

    ep_n = eeg_epoch_n(eeg)

    e_mean = zeros(ep_n)
    e_median = zeros(ep_n)
    e_std = zeros(ep_n)
    e_var = zeros(ep_n)
    e_kurt = zeros(ep_n)
    e_skew = zeros(ep_n)
    e_mean_diff = zeros(ep_n)
    e_median_diff = zeros(ep_n)
    e_max_dif = zeros(ep_n)
    e_dev_mean = zeros(ep_n)

    @inbounds @simd for epoch_idx in 1:ep_n
        e_mean[epoch_idx] = @views mean(eeg.eeg_signals[:, :, epoch_idx])
        e_median[epoch_idx] = @views median(eeg.eeg_signals[:, :, epoch_idx])
        e_std[epoch_idx] = @views std(eeg.eeg_signals[:, :, epoch_idx])
        e_var[epoch_idx] = @views var(eeg.eeg_signals[:, :, epoch_idx])
        e_kurt[epoch_idx] = @views kurtosis(eeg.eeg_signals[:, :, epoch_idx])
        e_skew[epoch_idx] = @views skewness(eeg.eeg_signals[:, :, epoch_idx])
        e_mean_diff = @views mean(diff(eeg.eeg_signals[:, :, epoch_idx], dims=2))
        e_median_diff = @views median(diff(eeg.eeg_signals[:, :, epoch_idx], dims=2))
        e_max_dif = @views maximum(eeg.eeg_signals[:, :, epoch_idx]) - minimum(eeg.eeg_signals[:, :, epoch_idx])
        e_dev_mean = @views abs(mean(eeg.eeg_signals[:, :, epoch_idx])) - mean(eeg.eeg_signals[:, :, epoch_idx])
    end

    return (e_mean=e_mean, e_median=e_median, e_std=e_std, e_var=e_var, e_kurt=e_kurt, e_skew=e_skew, e_mean_diff=e_mean_diff, e_median_diff=e_median_diff, e_max_dif=e_max_dif, e_dev_mean=e_dev_mean)
end

"""
    eeg_spectrogram(eeg; channel, norm, mt, st, demean)

Calculate spectrogram.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `method::Symbol=:standard`: method of calculating spectrogram:
    - `:standard`: standard
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `pad::Int64=0`: number of zeros to add
- `frq_lim::Tuple{Real, Real}=(0, 0)`: frequency limits
- `frq_n::Int64=0`: number of frequencies, default is length(frq_lim[1]:frq_lim[2])
- `norm::Bool=true`: normalize powers to dB
- `demean::Bool=true`: demean signal prior to analysis
- `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n) for frq === :log or ncyc = linspace(ncyc[1], ncyc[2], frq_n) for frq === :lin
- `wt<:CWT=wavelet(Morlet(π), β=2)`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

Named tuple containing:
- `s_pow::Array{Float64, 3}`
- `s_frq::Vector{Float64}`
- `s_t::Vector{Float64}`
"""
function eeg_spectrogram(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), pad::Int64=0, frq_lim::Tuple{Real, Real}=(0, 0), frq_n::Int64=0, method::Symbol=:standard, norm::Bool=true, demean::Bool=true, frq::Symbol=:log, gw::Real=5, ncyc::Union{Int64, Tuple{Int64, Int64}}=6, wt::T=wavelet(Morlet(π), β=2)) where {T <: CWT}

    _check_var(method, [:standard, :stft, :mt, :mw, :gh, :cwt], "method")
    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    # get frequency range
    fs = eeg_sr(eeg)
    frq_lim == (0, 0) && (frq_lim = (0, div(fs, 2)))
    frq_n == 0 && (frq_n = length(frq_lim[1]:frq_lim[2]))

    if method === :standard
        p_tmp, s_frq, _ = @views s_spectrogram(eeg.eeg_signals[1, :, 1], fs=fs, norm=norm, mt=false, st=false, demean=demean)
    elseif method === :mt
        p_tmp, s_frq, _ = @views s_spectrogram(eeg.eeg_signals[1, :, 1], fs=fs, norm=norm, mt=true, st=false, demean=demean)
    elseif method === :stft
        p_tmp, s_frq, _ = @views s_spectrogram(eeg.eeg_signals[1, :, 1], fs=fs, norm=norm, mt=false, st=true, demean=demean)
    elseif method === :mw
    _, p_tmp, _, s_frq = @views s_wspectrogram(eeg.eeg_signals[1, :, 1], pad=pad, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc, demean=demean)
    elseif method === :gh
        p_tmp, s_frq = @views s_ghspectrogram(eeg.eeg_signals[1, :, 1], fs=fs, frq_lim=frq_lim, frq_n=frq_n, norm=norm, frq=frq, demean=demean, gw=gw)
    elseif method === :cwt
        p_tmp, s_frq = @views s_cwtspectrogram(eeg.eeg_signals[1, :, 1], wt=wt, fs=fs, frq_lim=frq_lim, norm=norm, demean=demean)
    end

    s_t = linspace(0, (eeg_epoch_len(eeg) / fs), size(p_tmp, 2))
    s_pow = zeros(size(p_tmp, 1), size(p_tmp, 2), ch_n, ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ch_n, 1))

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            if method === :standard
                s_pow[:, :, channel_idx, epoch_idx], _, _ = @views s_spectrogram(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], fs=fs, norm=norm, mt=false, st=false, demean=demean)
            elseif method === :mt
                s_pow[:, :, channel_idx, epoch_idx], _, _ = @views s_spectrogram(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], fs=fs, norm=norm, mt=true, st=false, demean=demean)
            elseif method === :stft
                s_pow[:, :, channel_idx, epoch_idx], _, _ = @views s_spectrogram(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], fs=fs, norm=norm, mt=false, st=true, demean=demean)
            elseif method === :mw
                _, s_pow[:, :, channel_idx, epoch_idx], _, _ = @views s_wspectrogram(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], pad=pad, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc, demean=demean)
            elseif method === :gh
                s_pow[:, :, channel_idx, epoch_idx], _, _ = @views s_ghspectrogram(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], fs=fs, frq_lim=frq_lim, frq_n=frq_n, norm=norm, frq=frq, demean=demean, gw=gw)
            elseif method === :cwt
                s_pow[:, :, channel_idx, epoch_idx], _ = @views s_cwtspectrogram(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], wt=wt, fs=fs, frq_lim=frq_lim, norm=norm, demean=demean)
            end

            # update progress bar
            progress_bar == true && next!(p)
        end
    end

    s_frq = round.(s_frq, digits=2)
    s_t = round.(s_t, digits=2)
    s_t .+= eeg.eeg_epoch_time[1]

    return (s_pow=s_pow, s_frq=s_frq, s_t=s_t)
end

"""
    eeg_spectrum(eeg; channel, pad, h)

Calculate FFT/Hilbert transformation components, amplitudes, powers and phases.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `pad::Int64=0`: number of zeros to add signal for FFT
- `h::Bool=false`: use Hilbert transform for calculations instead of FFT
- `norm::Bool=false`: normalize do dB

# Returns

Named tuple containing:
- `c::Array{ComplexF64, 3}`: Fourier or Hilbert components
- `amp::Array{Float64, 3}`: amplitudes
- `pow::Array{Float64, 3}`: powers
- `pha::Array{Float64, 3}: phase angles
"""
function eeg_spectrum(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), pad::Int64=0, h::Bool=false, norm::Bool=false)

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    fft_size = eeg_epoch_len(eeg) + pad

    s_c = zeros(ComplexF64, ch_n, fft_size, ep_n)
    s_pha = zeros(ch_n, fft_size, ep_n)
    if h == true
        s_amp = zeros(ch_n, fft_size, ep_n)
        s_pow = zeros(ch_n, fft_size, ep_n)
    else
        s_amp = zeros(ch_n, fft_size ÷ 2, ep_n)
        s_pow = zeros(ch_n, fft_size ÷ 2, ep_n)
    end        

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            if h == true
                s_c[channel_idx, :, epoch_idx], s_amp[channel_idx, :, epoch_idx], s_pow[channel_idx, :, epoch_idx], s_pha[channel_idx, :, epoch_idx] = @views s_hspectrum(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], pad=pad, norm=norm)
            else
                s_c[channel_idx, :, epoch_idx], s_amp[channel_idx, :, epoch_idx], s_pow[channel_idx, :, epoch_idx], s_pha[channel_idx, :, epoch_idx] = @views s_spectrum(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], pad=pad, norm=norm)
            end
        end
    end

    return (c=s_c, amp=s_amp, pow=s_pow, pha=s_pha)
end

"""
    eeg_s2t(eeg; t)

Convert time in samples to seconds.

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

Convert time in seconds to samples.

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
    eeg_channel_stats(eeg)

Calculate channels statistics per epoch.

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
function eeg_channel_stats(eeg::NeuroAnalyzer.EEG)

    ch_n = eeg_channel_n(eeg)
    ep_n = eeg_epoch_n(eeg)
    
    c_mean = zeros(ch_n, ep_n)
    c_median = zeros(ch_n, ep_n)
    c_std = zeros(ch_n, ep_n)
    c_var = zeros(ch_n, ep_n)
    c_kurt = zeros(ch_n, ep_n)
    c_skew = zeros(ch_n, ep_n)
    c_mean_diff = zeros(ch_n, ep_n)
    c_median_diff = zeros(ch_n, ep_n)
    c_max_dif = zeros(ch_n, ep_n)
    c_dev_mean = zeros(ch_n, ep_n)

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            c_mean[channel_idx, epoch_idx] = @views mean(eeg.eeg_signals[channel_idx, :, epoch_idx])
            c_median[channel_idx, epoch_idx] = @views median(eeg.eeg_signals[channel_idx, :, epoch_idx])
            c_std[channel_idx, epoch_idx] = @views std(eeg.eeg_signals[channel_idx, :, epoch_idx])
            c_var[channel_idx, epoch_idx] = @views var(eeg.eeg_signals[channel_idx, :, epoch_idx])
            c_kurt[channel_idx, epoch_idx] = @views kurtosis(eeg.eeg_signals[channel_idx, :, epoch_idx])
            c_skew[channel_idx, epoch_idx] = @views skewness(eeg.eeg_signals[channel_idx, :, epoch_idx])
            c_mean_diff[channel_idx, epoch_idx] = @views mean(diff(eeg.eeg_signals[channel_idx, :, epoch_idx]))
            c_median_diff[channel_idx, epoch_idx] = @views median(diff(eeg.eeg_signals[channel_idx, :, epoch_idx]))
            c_max_dif[channel_idx, epoch_idx] = @views maximum(eeg.eeg_signals[channel_idx, :, epoch_idx]) - minimum(eeg.eeg_signals[channel_idx, :, epoch_idx])
            c_dev_mean[channel_idx, epoch_idx] = @views abs(mean(eeg.eeg_signals[channel_idx, :, epoch_idx])) - mean(eeg.eeg_signals[channel_idx, :, epoch_idx])
        end
    end

    return (c_mean=c_mean, c_median=c_median, c_std=c_std, c_var=c_var, c_kurt=c_kurt, c_skew=c_skew, c_mean_diff=c_mean_diff, c_median_diff=c_median_diff, c_max_dif=c_max_dif, c_dev_mean=c_dev_mean)
end

"""
    eeg_snr(eeg; channel, type)

Calculate SNR of `eeg` channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `type::Symbol=:rms`: SNR type:
    - `:mean`: mean-based
    - `:rms`: RMS-based

# Returns

Named tuple containing:
- `snr::Matrix(Float64)`: SNR for each channel over frequencies 1:Nyquist
- `hz::Vector(Float64)`: frequencies
"""
function eeg_snr(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), type::Symbol=:rms)

    _check_var(type, [:mean, :rms], "type")
    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    ep_n == 1 && throw(ArgumentError("EEG must contain ≥ 2 epochs."))

    hz, _ = s_freqs(eeg.eeg_epoch_time)
    amp = zeros(ch_n, length(hz), ep_n)
    snr = zeros(ch_n, length(hz))

    # create spectrum for each channel
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            _, amp[channel_idx, :, epoch_idx], _, _ = @views s_spectrum(eeg.eeg_signals[channel[channel_idx], :, epoch_idx])
        end
    end

    # calculate SNR for each channel spectrum
    @inbounds @simd for hz_idx in 1:length(hz)
        Threads.@threads for channel_idx in 1:ch_n
            if type === :mean
                snr[channel_idx, hz_idx] = @views s_snr(amp[channel_idx, hz_idx, :])
            else
                snr[channel_idx, hz_idx] = @views s_snr2(amp[channel_idx, hz_idx, :])
            end
        end
    end

    return (snr=snr, hz=hz)
end

"""
    eeg_standardize(eeg; channel)

Standardize channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels

# Returns

- `eeg_new::NeuroAnalyzer.EEG`: standardized EEG
- `scaler::Matrix{Float64}`: standardizing matrix
"""
function eeg_standardize(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg))
    
    _check_channels(eeg, channel)
    ep_n = eeg_epoch_n(eeg)

    scaler = Vector{Any}()

    eeg_new = deepcopy(eeg)
    @inbounds @simd for epoch_idx in 1:ep_n
        @views push!(scaler, StatsBase.fit(ZScoreTransform, eeg.eeg_signals[channel, :, epoch_idx], dims=2)) 
        @views eeg_new.eeg_signals[channel,:, epoch_idx] = StatsBase.transform(scaler[epoch_idx], eeg.eeg_signals[channel, :, epoch_idx])
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_standardize(EEG)")

    return eeg_new, scaler
end

"""
    eeg_standardize!(eeg; channel)

Standardize channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels

# Returns

- `scaler::Matrix{Float64}`: standardizing matrix
"""
function eeg_standardize!(eeg::NeuroAnalyzer.EEG)

    eeg_tmp, scaler = eeg_standardize(eeg, channel=channel)
    eeg.eeg_signals = eeg_tmp.eeg_signals
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg_reset_components!(eeg)

    return scaler
end

"""
    eeg_fconv(eeg; channel, kernel, norm)

Perform convolution in the frequency domain.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`: kernel for convolution
- `norm::Bool=true`: normalize kernel to keep the post-convolution results in the same scale as the original data
- `pad::Int64=0`: number of zeros to add signal for FFT

# Returns

- `s_convoluted::Array{Float64, 3}`: convoluted signal
"""
function eeg_fconv(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), kernel::Union{Vector{<:Real}, Vector{ComplexF64}}, norm::Bool=true, pad::Int64=0)

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    s_convoluted = zeros(ch_n, eeg_epoch_len(eeg), ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ch_n, 1))

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            s_convoluted[channel_idx, :, epoch_idx] = @views s_fconv(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], kernel=kernel, norm=norm, pad=pad)
            
            # update progress bar
            progress_bar == true && next!(p)
        end
    end

    return s_convoluted
end

"""
    eeg_tconv(eeg; channel, kernel)

Perform convolution in the time domain.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`: kernel used for convolution

# Returns

- `s_convoluted::Array{Float64, 3}`: convoluted signal
"""
function eeg_tconv(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), kernel::Union{Vector{<:Real}, Vector{ComplexF64}})

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    s_convoluted = zeros(ch_n, eeg_epoch_len(eeg), ep_n)

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            s_convoluted[channel_idx, :, epoch_idx] = s_tconv(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], kernel=kernel)
        end
    end

    return s_convoluted
end

"""
    eeg_dft(eeg; channel)

Return FFT and DFT sample frequencies for a DFT.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `pad::Int64=0`: number of zeros to add signal for FFT

# Returns

Named tuple containing:
- `sfft::Array{ComplexF64, 3}`: FFT
- `sf::Vector{Float64}`: sample frequencies
"""
function eeg_dft(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), pad::Int64=0)

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    fs = eeg_sr(eeg)
    sfft = zeros(ComplexF64, ch_n, eeg_epoch_len(eeg), ep_n)
    sf = nothing

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            sfft[channel_idx, :, epoch_idx], sf = @views s_dft(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], fs=fs, pad=pad)
        end
    end

    return (sfft=sfft, sf=sf)
end

"""
    eeg_msci95(eeg; channel, n, method)

Calculate mean, standard deviation and 95% confidence interval for EEG channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `n::Int64=3`: number of bootstraps
- `method::Symbol=:normal`: use normal (`:normal`) method or `n`-times bootstrapping (`:boot`)

# Returns

Named tuple containing:
- `s_m::Matrix{Float64}`: mean
- `s_s::Matrix{Float64}`: standard deviation
- `s_u::Matrix{Float64}`: upper 95% CI
- `s_l::Matrix{Float64}`: lower 95% CI
"""
function eeg_msci95(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), n::Int64=3, method::Symbol=:normal)

    _check_var(method, [:normal, :boot], "method")

    _check_channels(eeg, channel)
    epoch_len = eeg_epoch_len(eeg)
    ep_n = eeg_epoch_n(eeg)

    s_m = zeros(ep_n, epoch_len)
    s_s = zeros(ep_n, epoch_len)
    s_u = zeros(ep_n, epoch_len)
    s_l = zeros(ep_n, epoch_len)

    Threads.@threads for epoch_idx in 1:ep_n
        s_m[epoch_idx, :], s_s[epoch_idx, :], s_u[epoch_idx, :], s_l[epoch_idx, :] = @views s_msci95(eeg.eeg_signals[channel, :, epoch_idx], n=n, method=method)
    end

    return (mean=s_m, sd=s_s, upper=s_u, lower=s_l)
end

"""
    eeg_mean(eeg1, eeg2; channel1, channel2, epoch1, epoch2)

Calculates mean, standard deviation and 95% confidence interval for two EEG channels.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2:NeuroAnalyzer.EEG`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs

# Returns

Named tuple containing:
- `s_m::Matrix{Float64}`: mean by epochs
- `s_s::Matrix{Float64}`: standard deviation by epochs
- `s_u::Matrix{Float64}`: upper 95% CI bound by epochs
- `s_l::Matrix{Float64}`: lower 95% CI bound by epochs
"""
function eeg_mean(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg2.eeg_header[:signal_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2)))

    # check channels
    _check_channels(eeg1, channel1)
    _check_channels(eeg2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    # check epochs
    _check_epochs(eeg1, epoch1)
    _check_epochs(eeg2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    eeg_epoch_len(eeg1) == eeg_epoch_len(eeg2) || throw(ArgumentError("eeg1 and eeg2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    epoch_len = eeg_epoch_len(eeg1)

    s_m = zeros(ep_n, epoch_len)
    s_s = zeros(ep_n, epoch_len)
    s_u = zeros(ep_n, epoch_len)
    s_l = zeros(ep_n, epoch_len)

    Threads.@threads for epoch_idx in 1:ep_n
        s1_mean = @views mean(eeg1.eeg_signals[channel1, :, epoch1[epoch_idx]], dims=1)
        s2_mean = @views mean(eeg2.eeg_signals[channel2, :, epoch2[epoch_idx]], dims=1)
        s_m[epoch_idx, :] = s1_mean - s2_mean
        s1_sd = @views std(eeg1.eeg_signals[channel1, :, epoch1[epoch_idx]], dims=1) / sqrt(size(eeg1.eeg_signals[channel1, :, epoch1[epoch_idx]], 2))
        s2_sd = @views std(eeg2.eeg_signals[channel2, :, epoch2[epoch_idx]], dims=1) / sqrt(size(eeg2.eeg_signals[channel2, :, epoch2[epoch_idx]], 2))
        s_s[epoch_idx, :] = sqrt.(s1_sd.^2 .+ s2_sd.^2)
        s_u[epoch_idx, :] = @. s_m[epoch_idx, :] + 1.96 * s_s[epoch_idx, :]
        s_l[epoch_idx, :] = @. s_m[epoch_idx, :] - 1.96 * s_s[epoch_idx, :]
    end

    return (s_m=s_m, s_s=s_s, s_u=s_u, s_l=s_l)
end

"""
    eeg_difference(eeg1, eeg2; channel1, channel2, epoch1, epoch2, n, method)

Calculates mean difference and 95% confidence interval for two EEG channels.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2:NeuroAnalyzer.EEG`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs
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
function eeg_difference(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2)), n::Int64=3, method::Symbol=:absdiff)

    _check_channels(eeg1, channel1)
    _check_channels(eeg2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(eeg1, epoch1)
    _check_epochs(eeg2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    eeg_epoch_len(eeg1) == eeg_epoch_len(eeg2) || throw(ArgumentError("eeg1 and eeg2 epoch lengths must be equal."))

    ep_n = length(epoch1)

    s_stat = zeros(ep_n, length(channel1) * n)
    s_stat_single = zeros(ep_n)
    p = zeros(ep_n)

    Threads.@threads for epoch_idx in 1:ep_n
        s_stat[epoch_idx, :], s_stat_single[epoch_idx], p[epoch_idx] = @views s2_difference(eeg1.eeg_signals[channel1, :, epoch1[epoch_idx]], eeg2.eeg_signals[channel2, :, epoch2[epoch_idx]], n=n, method=method)
    end

    return (s_stat=s_stat, statsitic_single=s_stat_single, p=p)
end

"""
   eeg_acov(eeg; channel, lag=1, demean=false, norm=false)

Calculate autocovariance.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `lag::Int64`: lags range is `-lag:lag`
- `demean::Bool`: demean eeg prior to analysis
- `norm::Bool`: normalize autocovariance

# Returns

Named tuple containing:
- `acov::Matrix{Float64}`
- `lags::Vector{Float64}`: lags in ms
"""
function eeg_acov(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), lag::Int64=1, demean::Bool=false, norm::Bool=false)

    lag < 1 && throw(ArgumentError("lag must be ≥ 1."))

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    acov = zeros(ch_n, length(-lag:lag), ep_n)
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            acov[channel_idx, :, epoch_idx], _ = @views s_acov(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], lag=lag, demean=demean, norm=norm)
        end
    end

    lags = 1/eeg_sr(eeg) .* collect(-lag:lag) .* 1000

    return (acov=acov, acov_lags=lags)
end

"""
    eeg_tenv(eeg; channel, d)

Calculate temporal envelope (amplitude).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `t_env::Array{Float64, 3}`: temporal envelope
- `s_t::Vector{Float64}`: signal time
"""
function eeg_tenv(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), d::Int64=32)
    
    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    t_env = zeros(ch_n, eeg_epoch_len(eeg), ep_n)
    s_t = eeg.eeg_epoch_time

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            s = @view eeg.eeg_signals[channel[channel_idx], :, epoch_idx]
            # find peaks
            p_idx = s_findpeaks(s, d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(s))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(s_t[p_idx], s[p_idx])
                try
                    t_env[channel_idx, :, epoch_idx] = model(s_t)
                catch
                    @error "CubicSpline error, using Loess."
                    model = loess(s_t[p_idx], s[p_idx], span=0.5)
                    t_env[channel_idx, :, epoch_idx] = Loess.predict(model, s_t)
                end
            else
                _info("Less than 5 peaks detected, using Loess.")
                model = loess(s_t[p_idx], s[p_idx], span=0.5)
                t_env[channel_idx, :, epoch_idx] = Loess.predict(model, s_t)
            end
            t_env[channel_idx, 1, epoch_idx] = t_env[channel_idx, 2, epoch_idx]
        end
    end
    
    return (t_env=t_env, s_t=s_t)
end

"""
    eeg_tenv_mean(eeg; channel, dims, d)

Calculate temporal envelope (amplitude): mean and 95% CI.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `t_env_m::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: mean
- `t_env_u::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI upper bound
- `t_env_l::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function eeg_tenv_mean(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), dims::Int64, d::Int64=32)
    
    if dims == 1
        eeg_channel_n(eeg) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        eeg_epoch_n(eeg) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        eeg_channel_n(eeg) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        eeg_epoch_n(eeg) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_a, s_t = eeg_tenv(eeg, channel=channel, d=d)
    ch_n = size(s_a, 1)
    ep_n = size(s_a, 3)

    if dims == 1
        # mean over channels

        t_env_m = zeros(length(s_t), ep_n)
        t_env_u = zeros(length(s_t), ep_n)
        t_env_l = zeros(length(s_t), ep_n)

        @inbounds @simd for epoch_idx in 1:ep_n
            t_env_m[:, epoch_idx] = mean(s_a[:, :, epoch_idx], dims=1)
            s = std(t_env_m[:, epoch_idx]) / sqrt(length(t_env_m[:, epoch_idx]))
            t_env_u[:, epoch_idx] = @. t_env_m[:, epoch_idx] + 1.96 * s
            t_env_l[:, epoch_idx] = @. t_env_m[:, epoch_idx] - 1.96 * s
        end
    elseif dims == 2
        # mean over epochs

        t_env_m = zeros(length(s_t), ch_n)
        t_env_u = zeros(length(s_t), ch_n)
        t_env_l = zeros(length(s_t), ch_n)

        @inbounds @simd for channel_idx in 1:ch_n
            t_env_m[:, channel_idx] = mean(s_a[channel_idx, :, :], dims=2)
            s = std(t_env_m[:, channel_idx]) / sqrt(length(t_env_m[:, channel_idx]))
            t_env_u[:, channel_idx] = @. t_env_m[:, channel_idx] + 1.96 * s
            t_env_l[:, channel_idx] = @. t_env_m[:, channel_idx] - 1.96 * s
        end
    else
        # mean over channels and epochs

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
    eeg_tenv_median(eeg; channel, dims, d)

Calculate temporal envelope (amplitude): median and 95% CI.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `dims::Int64`: median over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `t_env_m::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: median
- `t_env_u::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI upper bound
- `t_env_l::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function eeg_tenv_median(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), dims::Int64, d::Int64=32)
    
    if dims == 1
        eeg_channel_n(eeg) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        eeg_epoch_n(eeg) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        eeg_channel_n(eeg) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        eeg_epoch_n(eeg) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_a, s_t = eeg_tenv(eeg, channel=channel, d=d)
    ch_n = size(s_a, 1)
    ep_n = size(s_a, 3)

    if dims == 1
        # median over channels

        t_env_m = zeros(length(s_t), ep_n)
        t_env_u = zeros(length(s_t), ep_n)
        t_env_l = zeros(length(s_t), ep_n)

        @inbounds @simd for epoch_idx in 1:ep_n
            t_env_m[:, epoch_idx] = median(s_a[:, :, epoch_idx], dims=1)
            t_idx = s_findpeaks(t_env_m[:, epoch_idx], d=d)
            pushfirst!(t_idx, 1)
            push!(t_idx, length(t_env_m[:, epoch_idx]))
            if length(t_idx) > 4
                model = CubicSpline(s_t[t_idx], t_env_m[t_idx])
                try
                    t_env_m[:, epoch_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(t_env_m[:, epoch_idx]) / sqrt(length(t_env_m[:, epoch_idx]))
            t_env_u[:, epoch_idx] = @. t_env_m[:, epoch_idx] + 1.96 * s
            t_env_l[:, epoch_idx] = @. t_env_m[:, epoch_idx] - 1.96 * s
        end
    elseif dims == 2
        # median over epochs

        t_env_m = zeros(length(s_t), ch_n)
        t_env_u = zeros(length(s_t), ch_n)
        t_env_l = zeros(length(s_t), ch_n)

        @inbounds @simd for channel_idx in 1:ch_n
            t_env_m[:, idx] = median(s_a[channel_idx, :, :], dims=2)
            t_idx = s_findpeaks(t_env_m[:, channel_idx], d=d)
            pushfirst!(t_idx, 1)
            push!(t_idx, length(t_env_m[:, channel_idx]))
            if length(t_idx) > 4
                model = CubicSpline(s_t[t_idx], t_env_m[t_idx])
                try
                    t_env_m[:, channel_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(t_env_m[:, channel_idx]) / sqrt(length(t_env_m[:, channel_idx]))
            t_env_u[:, channel_idx] = @. t_env_m[:, channel_idx] + 1.96 * s
            t_env_l[:, channel_idx] = @. t_env_m[:, channel_idx] - 1.96 * s
        end
    else
        # median over channels and epochs

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
    eeg_penv(eeg; channel, d)

Calculate power (in dB) envelope.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns

Named tuple containing:
- `p_env::Array{Float64, 3}`: power spectrum envelope
- `p_env_frq::Vector{Float64}`: frequencies for each envelope
"""
function eeg_penv(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), d::Int64=8, mt::Bool=false, nt::Int64=8)
    
    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    fs = eeg_sr(eeg)

    psd_tmp, frq = s_psd(eeg.eeg_signals[1, :, 1], fs=fs, mt=mt, nt=nt)
    p_env = zeros(ch_n, length(psd_tmp), ep_n)

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            psd_pow, _ = s_psd(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], fs=fs, mt=mt, norm=true, nt=nt)
            # find peaks
            p_idx = s_findpeaks(psd_pow, d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(psd_pow))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(frq[p_idx], psd_pow[p_idx])
                try
                    p_env[channel_idx, :, epoch_idx] = model(frq)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
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
    eeg_penv_mean(eeg; channel, dims, d)

Calculate power (in dB) envelope: mean and 95% CI.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
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
function eeg_penv_mean(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), dims::Int64, d::Int64=8, mt::Bool=false)
    
    if dims == 1
        eeg_channel_n(eeg) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        eeg_epoch_n(eeg) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        eeg_channel_n(eeg) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        eeg_epoch_n(eeg) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_p, s_f = eeg_psd(eeg, channel=channel, norm=true, mt=mt)
    ch_n = size(s_p, 1)
    ep_n = size(s_p, 3)

    if dims == 1
        # mean over channels

        p_env_m = zeros(length(s_f), ep_n)
        p_env_u = zeros(length(s_f), ep_n)
        p_env_l = zeros(length(s_f), ep_n)

        @inbounds @simd for epoch_idx in 1:ep_n
            p_env_m[:, epoch_idx] = mean(s_p[:, :, epoch_idx], dims=1)
            # find peaks
            p_idx = s_findpeaks(p_env_m[:, epoch_idx], d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(p_env_m[:, epoch_idx]))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(s_f[p_idx], p_env_m[p_idx])
                try
                    p_env_m[:, epoch_idx] = model(s_f)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = std(p_env_m[:, epoch_idx]) / sqrt(length(p_env_m[:, epoch_idx]))
            p_env_u[:, epoch_idx] = @. p_env_m[:, epoch_idx] + 1.96 * s
            p_env_l[:, epoch_idx] = @. p_env_m[:, epoch_idx] - 1.96 * s
        end
    elseif dims == 2
        # mean over epochs

        p_env_m = zeros(length(s_f), ch_n)
        p_env_u = zeros(length(s_f), ch_n)
        p_env_l = zeros(length(s_f), ch_n)

        @inbounds @simd for channel_idx in 1:ch_n
            p_env_m[:, channel_idx] = mean(s_p[channel_idx, :, :], dims=2)
            # find peaks
            p_idx = s_findpeaks(p_env_m[:, channel_idx], d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(p_env_m[:, channel_idx]))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(s_f[p_idx], p_env_m[p_idx])
                try
                    p_env_m[:, channel_idx] = model(s_f)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = std(p_env_m[:, channel_idx]) / sqrt(length(p_env_m[:, channel_idx]))
            p_env_u[:, channel_idx] = @. p_env_m[:, channel_idx] + 1.96 * s
            p_env_l[:, channel_idx] = @. p_env_m[:, channel_idx] - 1.96 * s
        end
    else
        # mean over channels and epochs

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
    eeg_penv_median(eeg; channel, dims, d)

Calculate power (in dB) envelope: median and 95% CI.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
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
function eeg_penv_median(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), dims::Int64, d::Int64=8, mt::Bool=false)
    
    if dims == 1
        eeg_channel_n(eeg) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        eeg_epoch_n(eeg) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        eeg_channel_n(eeg) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        eeg_epoch_n(eeg) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_p, s_f = eeg_psd(eeg, channel=channel, norm=true, mt=mt)
    ch_n = size(s_p, 1)
    ep_n = size(s_p, 3)

    if dims == 1
        # median over channels

        p_env_m = zeros(length(s_f), ep_n)
        p_env_u = zeros(length(s_f), ep_n)
        p_env_l = zeros(length(s_f), ep_n)

        @inbounds @simd for epoch_idx in 1:ep_n
            p_env_m[:, epoch_idx] = median(s_p[:, :, epoch_idx], dims=1)
            # find peaks
            p_idx = s_findpeaks(p_env_m[:, epoch_idx], d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(p_env_m[:, epoch_idx]))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(s_f[p_idx], p_env_m[p_idx])
                try
                    p_env_m[:, epoch_idx] = model(s_f)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(p_env_m[:, epoch_idx]) / sqrt(length(p_env_m[:, epoch_idx]))
            p_env_u[:, epoch_idx] = @. p_env_m[:, epoch_idx] + 1.96 * s
            p_env_l[:, epoch_idx] = @. p_env_m[:, epoch_idx] - 1.96 * s
        end
    elseif dims == 2
        # median over epochs

        p_env_m = zeros(length(s_f), ch_n)
        p_env_u = zeros(length(s_f), ch_n)
        p_env_l = zeros(length(s_f), ch_n)

        @inbounds @simd for channel_idx in 1:ch_n
            p_env_m[:, channel_idx] = median(s_p[channel_idx, :, :], dims=2)
            # find peaks
            p_idx = s_findpeaks(p_env_m[:, channel_idx], d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(p_env_m[:, channel_idx]))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(s_f[p_idx], p_env_m[p_idx])
                try
                    p_env_m[:, channel_idx] = model(s_f)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(p_env_m[:, channel_idx]) / sqrt(length(p_env_m[:, channel_idx]))
            p_env_u[:, channel_idx] = @. p_env_m[:, channel_idx] + 1.96 * s
            p_env_l[:, channel_idx] = @. p_env_m[:, channel_idx] - 1.96 * s
        end
    else
        # median over channels and epochs
        
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
    eeg_senv(eeg; channel, d, mt, t)

Calculate spectral envelope.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
- `mt::Bool=false`: if true use multi-tapered spectrogram
- `t::Union{Real, Nothing}=nothing`: spectrogram threshold (maximize all powers > t)

# Returns

Named tuple containing:
- `s_env::Array{Float64, 3}`: spectral envelope
- `s_env_t::Vector{Float64}`: spectrogram time
"""
function eeg_senv(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), d::Int64=2, mt::Bool=false, t::Union{Real, Nothing}=nothing)
    
    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    fs = eeg_sr(eeg)
    s_tmp = @view eeg.eeg_signals[1, :, 1]
    interval = fs
    overlap = round(Int64, fs * 0.75)
    # for short signals always use multi-taper
    length(s_tmp) < 4 * fs && (mt = true)
    if mt == true
        spec_tmp = mt_spectrogram(s_tmp, fs=fs)
    else
        spec_tmp = spectrogram(s_tmp, interval, overlap, nfft=length(s_tmp), fs=fs, window=hanning)
    end
    sp_t = collect(spec_tmp.time)
    sp_t .+= eeg.eeg_epoch_time[1]

    s_env = zeros(ch_n, length(sp_t), ep_n)

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            # prepare spectrogram
            if mt == true
                spec = @views mt_spectrogram(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], fs=fs)
            else
                spec = @views spectrogram(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], interval, overlap, nfft=length(s_tmp), fs=fs, window=hanning)
            end
            s_frq = Vector(spec.freq)
            s_p = pow2db.(spec.power)

            # maximize all powers above threshold (t)
            if t !== nothing
                s_p[s_p .> t] .= 0
                reverse!(s_p)
                reverse!(s_frq)
            end
            
            f_idx = zeros(length(spec.time))
            m = maximum(s_p, dims=1)
            for idx2 in eachindex(m)
                f_idx[idx2] = s_frq[vsearch(m[idx2], s_p[:, idx2])]
            end
            # find peaks
            p_idx = s_findpeaks(f_idx, d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(spec.time))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(sp_t[p_idx], f_idx[p_idx])
                try
                    s_env[channel_idx, :, epoch_idx] = model(sp_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
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
    eeg_senv_mean(eeg; channel, dims, d, mt, t)

Calculate spectral envelope: mean and 95% CI.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
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
function eeg_senv_mean(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), dims::Int64, d::Int64=2, mt::Bool=false, t::Union{Real, Nothing}=nothing)

    if dims == 1
        eeg_channel_n(eeg) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        eeg_epoch_n(eeg) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        eeg_channel_n(eeg) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        eeg_epoch_n(eeg) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_p, s_t = eeg_senv(eeg, channel=channel, d=d, mt=mt, t=t)
    ch_n = size(s_p, 1)
    ep_n = size(s_p, 3)

    if dims == 1
        # mean over channels

        s_env_m = zeros(length(s_t), ep_n)
        s_env_u = zeros(length(s_t), ep_n)
        s_env_l = zeros(length(s_t), ep_n)

        @inbounds @simd for epoch_idx in 1:ep_n
            s_env_m[:, epoch_idx] = mean(s_p[:, :, epoch_idx], dims=1)
            # find peaks
            s_idx = s_findpeaks(s_env_m[:, epoch_idx], d=d)
            # add first time-point
            pushfirst!(s_idx, 1)
            # interpolate peaks using cubic spline or loess
            push!(s_idx, length(s_env_m[:, epoch_idx]))
            if length(s_idx) > 4
                model = CubicSpline(s_t[s_idx], s_env_m[s_idx])
                try
                    s_env_m[:, epoch_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = std(s_env_m[:, epoch_idx]) / sqrt(length(s_env_m[:, epoch_idx]))
            s_env_u[:, epoch_idx] = @. s_env_m[:, epoch_idx] + 1.96 * s
            s_env_l[:, epoch_idx] = @. s_env_m[:, epoch_idx] - 1.96 * s
        end
    elseif dims == 2
        # mean over epochs

        s_env_m = zeros(length(s_t), ch_n)
        s_env_u = zeros(length(s_t), ch_n)
        s_env_l = zeros(length(s_t), ch_n)

        @inbounds @simd for channel_idx in 1:ch_n
            s_env_m[:, channel_idx] = mean(s_p[channel_idx, :, :], dims=2)
            # find peaks
            s_idx = s_findpeaks(s_env_m[:, channel_idx], d=d)
            # add first time-point
            pushfirst!(s_idx, 1)
            # add last time-point
            push!(s_idx, length(s_env_m[:, channel_idx]))
            # interpolate peaks using cubic spline or loess
            if length(s_idx) > 4
                model = CubicSpline(s_t[s_idx], s_env_m[s_idx])
                try
                    s_env_m[:, channel_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = std(s_env_m[:, channel_idx]) / sqrt(length(s_env_m[:, channel_idx]))
            s_env_u[:, channel_idx] = @. s_env_m[:, channel_idx] + 1.96 * s
            s_env_l[:, channel_idx] = @. s_env_m[:, channel_idx] - 1.96 * s
        end
    else
        # mean over channels and epochs

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
    eeg_senv_median(eeg; channel, dims, d, mt)

Calculate spectral envelope: median and 95% CI.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `dims::Int64`: median over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
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
function eeg_senv_median(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), dims::Int64, d::Int64=2, mt::Bool=false, t::Union{Real, Nothing}=nothing)
    
    if dims == 1
        eeg_channel_n(eeg) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        eeg_epoch_n(eeg) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        eeg_channel_n(eeg) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        eeg_epoch_n(eeg) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_p, s_t = eeg_senv(eeg, channel=channel, d=d, mt=mt, t=t)
    ch_n = size(s_p, 1)
    ep_n = size(s_p, 3)

    if dims == 1
        # median over channels

        s_env_m = zeros(length(s_t), ep_n)
        s_env_u = zeros(length(s_t), ep_n)
        s_env_l = zeros(length(s_t), ep_n)

        @inbounds @simd for epoch_idx in 1:ep_n
            s_env_m[:, epoch_idx] = median(s_p[:, :, epoch_idx], dims=1)
            # find peaks
            s_idx = s_findpeaks(s_env_m[:, epoch_idx], d=d)
            # add first time-point
            pushfirst!(s_idx, 1)
            # add last time-point
            push!(s_idx, length(s_env_m[:, epoch_idx]))
            # interpolate peaks using cubic spline or loess
            if length(s_idx) > 4
                model = CubicSpline(s_t[s_idx], s_env_m[s_idx])
                try
                    s_env_m[:, epoch_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(s_env_m[:, epoch_idx]) / sqrt(length(s_env_m[:, epoch_idx]))
            s_env_u[:, epoch_idx] = @. s_env_m[:, epoch_idx] + 1.96 * s
            s_env_l[:, epoch_idx] = @. s_env_m[:, epoch_idx] - 1.96 * s
        end
    elseif dims == 2
        # median over epochs

        s_env_m = zeros(length(s_t), ch_n)
        s_env_u = zeros(length(s_t), ch_n)
        s_env_l = zeros(length(s_t), ch_n)

        @inbounds @simd for channel_idx in 1:ch_n
            s_env_m[:, channel_idx] = median(s_p[channel_idx, :, :], dims=2)
            # find peaks
            s_idx = s_findpeaks(s_env_m[:, channel_idx], d=d)
            # add first time-point
            pushfirst!(s_idx, 1)
            # add last time-point
            push!(s_idx, length(s_env_m[:, channel_idx]))
            # interpolate peaks using cubic spline or loess
            if length(s_idx) > 4
                model = CubicSpline(s_t[s_idx], s_env_m[s_idx])
                try
                    s_env_m[:, channel_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(s_env_m[:, channel_idx]) / sqrt(length(s_env_m[:, channel_idx]))
            s_env_u[:, channel_idx] = @. s_env_m[:, channel_idx] + 1.96 * s
            s_env_l[:, channel_idx] = @. s_env_m[:, channel_idx] - 1.96 * s
        end
    else
        # median over channels and epochs

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

Calculate ISPC (Inter-Site-Phase Clustering).

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2::NeuroAnalyzer.EEG`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs

# Returns

Named tuple containing:
- `ispc::Array{Float64, 2}`: ISPC value
- `ispc_angle::Array{Float64, 2}`: ISPC angle
- `signal_diff::Array{Float64, 3}`: signal difference (signal2 - signal1)
- `phase_diff::Array{Float64, 3}`: phase difference (signal2 - signal1)
- `s1_phase::Array{Float64, 3}`: signal 1 phase
- `s2_phase::Array{Float64, 3}`: signal 2 phase
"""
function eeg_ispc(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2)))

    _check_channels(eeg1, channel1)
    _check_channels(eeg2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(eeg1, epoch1)
    _check_epochs(eeg2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    eeg_epoch_len(eeg1) == eeg_epoch_len(eeg2) || throw(ArgumentError("eeg1 and eeg2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)

    ispc = zeros(ch_n, ep_n)
    ispc_angle = zeros(ch_n, ep_n)
    signal_diff = zeros(ch_n, eeg_epoch_len(eeg1), ep_n)
    phase_diff = zeros(ch_n, eeg_epoch_len(eeg1), ep_n)
    s1_phase = zeros(ch_n, eeg_epoch_len(eeg1), ep_n)
    s2_phase = zeros(ch_n, eeg_epoch_len(eeg1), ep_n)

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            ispc[channel_idx, epoch_idx], ispc_angle[channel_idx, epoch_idx], signal_diff[channel_idx, :, epoch_idx], phase_diff[channel_idx, :, epoch_idx], s1_phase[channel_idx, :, epoch_idx], s2_phase[channel_idx, :, epoch_idx] = @views s2_ispc(eeg1.eeg_signals[channel1[channel_idx], :, epoch1[epoch_idx]], eeg2.eeg_signals[channel2[channel_idx], :, epoch2[epoch_idx]])
        end
    end

    return (ispc=ispc, ispc_angle=ispc_angle, signal_diff=signal_diff, phase_diff=phase_diff, s1_phase=s1_phase, s2_phase=s2_phase)
end

"""
    eeg_itpc(eeg; channel, t, w)

Calculate ITPC (Inter-Trial-Phase Clustering) at sample number `t` over epochs/trials.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `t::Int64`: time point (sample number) at which ITPC is calculated
- `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation

# Returns

Named tuple containing:
- `itpc::Vector{Float64}`: ITPC or wITPC value
- `itpcz::Vector{Float64}`: Rayleigh's ITPC Z value
- `itpc_angle::Vector{Float64}`: ITPC angle
- `itpc_phases::Array{Float64, 2}`: phase difference (channel2 - channel1)
"""
function eeg_itpc(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), t::Int64, w::Union{Vector{<:Real}, Nothing}=nothing)

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)
    t < 1 && throw(ArgumentError("t must be ≥ 1."))
    t > eeg_epoch_len(eeg) && throw(ArgumentError("t must be ≤ $(eeg_epoch_len(eeg))."))
    ep_n < 2 && throw(ArgumentError("EEG must contain ≥ 2 epochs."))

    itpc = zeros(ch_n)
    itpcz = zeros(ch_n)
    itpc_angle = zeros(ch_n)
    itpc_phases = zeros(ch_n, ep_n)

    Threads.@threads for channel_idx in 1:ch_n
        @inbounds itpc[channel_idx], itpcz[channel_idx], itpc_angle[channel_idx], itpc_phases[channel_idx, :] = @views s_itpc(reshape(eeg.eeg_signals[channel[channel_idx], :, :], 1, :, ep_n), t=t, w=w)
    end
    return (itpc=itpc, itpcz=itpcz, itpc_angle=itpc_angle, itpc_phases=itpc_phases)
end

"""
    eeg_pli(eeg1, eeg2; channel1, channel2, epoch1, epoch2)

Calculate PLI (Phase Lag Index) between `eeg1` and `eeg2`.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2::NeuroAnalyzer.EEG`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs

# Returns

Named tuple containing:
- `pli::Array{Float64, 2}`: PLI value
- `signal_diff::Array{Float64, 3}`: signal difference (signal2 - signal1)
- `phase_diff::Array{Float64, 3}`: phase difference (signal2 - signal1)
- `s1_phase::Array{Float64, 3}`: signal 1 phase
- `s2_phase::Array{Float64, 3}`: signal 2 phase
"""
function eeg_pli(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2)))

    _check_channels(eeg1, channel1)
    _check_channels(eeg2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(eeg1, epoch1)
    _check_epochs(eeg2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    eeg_epoch_len(eeg1) == eeg_epoch_len(eeg2) || throw(ArgumentError("eeg1 and eeg2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)

    pli = zeros(ch_n, ep_n)
    signal_diff = zeros(ch_n, eeg_epoch_len(eeg1), ep_n)
    phase_diff = zeros(ch_n, eeg_epoch_len(eeg1), ep_n)
    s1_phase = zeros(ch_n, eeg_epoch_len(eeg1), ep_n)
    s2_phase = zeros(ch_n, eeg_epoch_len(eeg1), ep_n)

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            pli[channel_idx, epoch_idx], signal_diff[channel_idx, :, epoch_idx], phase_diff[channel_idx, :, epoch_idx], s1_phase[channel_idx, :, epoch_idx], s2_phase[channel_idx, :, epoch_idx] = @views s2_pli(eeg1.eeg_signals[channel1[channel_idx], :, epoch1[epoch_idx]], eeg2.eeg_signals[channel2[channel_idx], :, epoch2[epoch_idx]])
        end
    end

    return (pli=pli, signal_diff=signal_diff, phase_dif=phase_diff, s1_phase=s1_phase, s2_phase=s2_phase)
end

"""
    eeg_pli(eeg; channel)

Calculate PLIs (Phase Lag Index).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels

# Returns

- `pli_m::Array{Float64, 3}`: PLI value matrices over epochs
"""
function eeg_pli(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg))

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    pli_m = zeros(ch_n, ch_n, ep_n)

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx1 in 1:ch_n
            for channel_idx2 in 1:channel_idx1
                pli_m[channel_idx1, channel_idx2, epoch_idx], _, _, _, _ = @views s2_pli(eeg.eeg_signals[channel[channel_idx1], :, epoch_idx], eeg.eeg_signals[channel[channel_idx2], :, epoch_idx])
            end
        end
        Threads.@threads for channel_idx1 in 1:(ch_n - 1)
            for channel_idx2 in (channel_idx1 + 1):ch_n
                pli_m[channel_idx1, channel_idx2, epoch_idx] = @views pli_m[channel_idx2, channel_idx1, epoch_idx]
            end
        end
    end

    return pli_m
end

"""
    eeg_ispc(eeg; channel)

Calculate ISPCs (Inter-Site-Phase Clustering).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels

# Returns

- `ispc_m::Array{Float64, 3}`: ISPC value matrices over epochs
"""
function eeg_ispc(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg))

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    ispc_m = zeros(ch_n, ch_n, ep_n)

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx1 in 1:ch_n
            for channel_idx2 in 1:channel_idx1
                ispc_m[channel_idx1, channel_idx2, epoch_idx], _, _, _, _, _ = @views s2_ispc(eeg.eeg_signals[channel[channel_idx1], :, epoch_idx], eeg.eeg_signals[channel[channel_idx2], :, epoch_idx])
            end
        end

        Threads.@threads for channel_idx1 in 1:(ch_n - 1)
            for channel_idx2 in (channel_idx1 + 1):ch_n
                ispc_m[channel_idx1, channel_idx2, epoch_idx] = @views ispc_m[channel_idx2, channel_idx1, epoch_idx]
            end
        end
    end

    return ispc_m
end

"""
    eeg_ec(eeg1, eeg2; type, channel1, channel2, epoch1, epoch2)

Calculate envelope correlation.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2::NeuroAnalyzer.EEG`
- `type::Symbol=:amp`: envelope type:
    - `:amp`: amplitude
    - `:pow`: power
    - `:spec`: spectrogram
    - `:hamp`: Hilbert spectrum amplitude
- `channel1::Int64`
- `channel2::Int64`
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs

# Returns

Named tuple containing:
- `ec::Vector{Float64}`: power correlation value
- `ec_p::Vector{Float64}`: power correlation p-value
"""
function eeg_ec(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; type::Symbol=:amp, channel1::Int64, channel2::Int64, epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2)))

    _check_var(type, [:amp, :pow, :spec, :hamp], "type")

    _check_channels(eeg1, channel1)
    _check_channels(eeg2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))

    _check_epochs(eeg1, epoch1)
    _check_epochs(eeg2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    eeg_epoch_len(eeg1) == eeg_epoch_len(eeg2) || throw(ArgumentError("eeg1 and eeg2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    
    ec_r = zeros(ep_n)
    ec_p = zeros(ep_n)

    # calculate envelopes
    if type === :amp
        s1, _ = eeg_tenv(eeg1, channel=channel1)
        s2, _ = eeg_tenv(eeg2, channel=channel2)
    elseif type === :pow
        s1, _ = eeg_penv(eeg1, channel=channel1)
        s2, _ = eeg_penv(eeg2, channel=channel2)
    elseif type === :spec
        s1, _ = eeg_senv(eeg1, channel=channel1)
        s2, _ = eeg_senv(eeg2, channel=channel2)
    elseif type === :hamp
        s1, _ = eeg_henv(eeg1, channel=channel1)
        s2, _ = eeg_henv(eeg2, channel=channel2)
    end
    s1 = s1[:, :, epoch1]
    s2 = s2[:, :, epoch2]
    
    # compare envelopes per epochs
    Threads.@threads for epoch_idx in 1:ep_n
        ec = CorrelationTest(vec(s1[:, :, epoch_idx]), vec(s2[:, :, epoch_idx]))
        @inbounds ec_r[epoch_idx] = ec.r
        @inbounds ec_p[epoch_idx] = pvalue(ec)
    end

    return (ec=ec_r, ec_p=ec_p)
end

"""
    eeg_ged(eeg1, eeg2; channel1, channel2, epoch1, epoch2)

Perform generalized eigendecomposition.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`: signal data to be analyzed
- `eeg2::NeuroAnalyzer.EEG`: original signal data
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs

# Returns

- `sged::Array{Float64, 3}`
- `ress::Matrix{Float64}`
- `ress_normalized::Matrix{Float64}`
"""
function eeg_ged(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2)))

    _check_channels(eeg1, channel1)
    _check_channels(eeg2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(eeg1, epoch1)
    _check_epochs(eeg2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    eeg_epoch_len(eeg1) == eeg_epoch_len(eeg2) || throw(ArgumentError("eeg1 and eeg2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)

    sged = zeros(ch_n, eeg_epoch_len(eeg1), ep_n)
    ress = zeros(ch_n, ep_n)
    ress_normalized = zeros(ch_n, ep_n)

    Threads.@threads for epoch_idx in 1:ep_n
        sged[:, :, epoch_idx], ress[:, epoch_idx], ress_normalized[:, epoch_idx] = @views s2_ged(eeg1.eeg_signals[channel1, :, epoch1[epoch_idx]], eeg2.eeg_signals[channel2, :, epoch2[epoch_idx]])
    end

    return (sged=sged, ress=ress, ress_normalized=ress_normalized)
end

"""
    eeg_frqinst(eeg; channel)

Calculate instantaneous frequency.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels

# Returns

- `frqinst::Array{Float64, 3}`
"""
function eeg_frqinst(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg))

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    frqinst = zeros(ch_n, eeg_epoch_len(eeg), ep_n)
    fs = eeg_sr(eeg)

    _info("eeg_frqinst() uses Hilbert transform, the signal should be narrowband for best results.")

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ch_n, 1))
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            frqinst[channel_idx, :, epoch_idx] = @views s_frqinst(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], fs=fs)
        end
        # update progress bar
        progress_bar == true && next!(p)
    end
    return frqinst
end

"""
    eeg_itpc_s(eeg; <keyword arguments>)

Calculate spectrogram of ITPC (Inter-Trial-Phase Clustering).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64`
- `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
- `frq_n::Int64`: number of frequencies
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation

# Returns

Named tuple containing:
- `itpc_s::Array{Float64, 3}`: spectrogram of ITPC values
- `itpc_z_s::Array{Float64, 3}`: spectrogram ITPCz values
- `itpc_frq::Vector{Float64}`: frequencies list
"""
function eeg_itpc_s(eeg::NeuroAnalyzer.EEG; channel::Int64, frq_lim::Tuple{Real, Real}, frq_n::Int64, frq::Symbol=:log, w::Union{Vector{<:Real}, Nothing}=nothing)

    _check_var(frq, [:log, :lin], "frq")
    frq_lim = tuple_order(frq_lim)
    frq_lim[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
    frq_lim[2] > eeg_sr(eeg) ÷ 2 && throw(ArgumentError("Upper frequency bound must be ≤ $(eeg_sr(eeg) ÷ 2)."))
    frq_n < 2 && throw(ArgumentError("frq_n must be ≥ 2."))
    if frq === :log
        frq_lim[1] == 0 && (frq_lim = (0.01, frq_lim[2]))
        frq_lim = (frq_lim[1], frq_lim[2])
        frq_list = logspace(log10(frq_lim[1]), log10(frq_lim[2]), frq_n)
    else
        frq_list = linspace(frq_lim[1], frq_lim[2], frq_n)
    end

    _check_channels(eeg, channel)
    ep_n = eeg_epoch_n(eeg)
    epoch_len = eeg_epoch_len(eeg)
    ep_n < 2 && throw(ArgumentError("eeg must contain ≥ 2 epochs."))

    itpc_s = zeros(frq_n, epoch_len)
    itpc_z_s = zeros(frq_n, epoch_len)

    # initialize progress bar
    progress_bar == true && (p = Progress(frq_n, 1))

    Threads.@threads for frq_idx in 1:frq_n
        # create Morlet wavelet
        kernel = generate_morlet(eeg_sr(eeg), frq_list[frq_idx], 1, ncyc=10)
        half_kernel = floor(Int64, length(kernel) / 2) + 1
        s_conv = zeros(Float64, 1, epoch_len, ep_n)
        # convolute with Morlet wavelet
        @inbounds @simd for epoch_idx in 1:ep_n
            s_conv[1, :, epoch_idx] = @views DSP.conv(eeg.eeg_signals[channel, :, epoch_idx], kernel)[(half_kernel - 1):(end - half_kernel)]
        end
        # calculate ITPC of the convoluted signals
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
    eeg_tkeo(eeg; channel)

Calculate Teager-Kaiser energy-tracking operator: y(t) = x(t)^2 - x(t-1) × x(t+1)

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels

# Returns

- `tkeo::Array{Float64, 3}`
"""
function eeg_tkeo(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg))

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    tkeo = zeros(ch_n, eeg_epoch_len(eeg), ep_n)
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            tkeo[channel_idx, :, epoch_idx] = @views s_tkeo(eeg.eeg_signals[channel[channel_idx], :, epoch_idx])
        end
    end

    return tkeo
end

"""
    eeg_mwpsd(eeg; channel, pad, norm, frq_lim, frq_n, frq, ncyc)

Calculate power spectrum using Morlet wavelet convolution.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `pad::Int64`: pad the `signal` with `pad` zeros
- `norm::Bool`=true: normalize powers to dB
- `frq_lim::Tuple{Real, Real}=(0, 0)`: frequency bounds for the spectrogram
- `frq_n::Int64=10`: number of frequencies
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n) for frq === :log or ncyc = linspace(ncyc[1], ncyc[2], frq_n) for frq === :lin

# Returns

Named tuple containing:
- `w_pow::Array{Float64, 4}`
- `w_frq::Matrix{Float64}`
"""
function eeg_mwpsd(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), pad::Int64=0, norm::Bool=true, frq_lim::Tuple{Real, Real}=(0, 0), frq_n::Int64=0, frq::Symbol=:lin, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    fs = eeg_sr(eeg)
    p_tmp, w_frq = @views s_mwpsd(eeg.eeg_signals[1, :, 1], fs=fs, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)
    w_pow = zeros(ch_n, length(p_tmp), ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ch_n, 1))

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            w_pow[channel_idx, :, epoch_idx], _ = @views s_mwpsd(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], pad=pad, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)

            # update progress bar
            progress_bar == true && next!(p)
        end
    end

    return (w_pow=w_pow, w_frq=w_frq)
end

"""
    eeg_fcoherence(eeg1, eeg2; channel1, channel2, epoch1, epoch2, frq_lim)

Calculate coherence (mean over frequencies) and MSC (magnitude-squared coherence).

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2::NeuroAnalyzer.EEG`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs
- `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: return coherence only for the given frequency range

# Returns

Named tuple containing:
- `c::Array{Float64, 3}`: coherence
- `msc::Array{Float64, 3}`: MSC
- `f::Vector{Float64}`: frequencies
"""
function eeg_fcoherence(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Union{Int64, Vector{Int64}, AbstractRange}=0, channel2::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch1::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch2::Union{Int64, Vector{Int64}, AbstractRange}=0, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)

    _check_channels(eeg1, channel1)
    _check_channels(eeg2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(eeg1, epoch1)
    _check_epochs(eeg2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    eeg_epoch_len(eeg1) == eeg_epoch_len(eeg2) || throw(ArgumentError("eeg1 and eeg2 epoch lengths must be equal."))

    eeg_sr(eeg1) == eeg_sr(eeg2) || throw(ArgumentError("EEG1 and EEG2 must have the same sampling rate."))

    c_tmp, _, f = @views s2_fcoherence(eeg1.eeg_signals[1, :, 1], eeg1.eeg_signals[1, :, 1], fs=eeg_sr(eeg1), frq_lim=frq_lim)
    c = zeros(length(channel1), length(c_tmp), length(epoch1))
    msc = zeros(length(channel1), length(c_tmp), length(epoch1))
    f = zeros(length(channel1), length(c_tmp), length(epoch1))
    @inbounds @simd for epoch_idx in eachindex(epoch1)
        Threads.@threads for channel_idx in eachindex(channel1)
            c[channel_idx, :, epoch_idx], msc[channel_idx, :, epoch_idx], _ = @views s2_fcoherence(eeg1.eeg_signals[channel1[channel_idx], :, epoch1[epoch_idx]], eeg2.eeg_signals[channel2[channel_idx], :, epoch2[epoch_idx]], fs=eeg_sr(eeg1), frq_lim=frq_lim)
        end
    end

    return (c=c, msc=msc, f=f)
end

"""
    eeg_vartest(eeg; channel)

Calculate variance F-test.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels

# Returns

Named tuple containing:
- `f::Array{Float64, 3}`
- `p::Array{Float64, 3}`
"""
function eeg_vartest(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg))

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    f = zeros(ch_n, ch_n, ep_n)
    p = zeros(ch_n, ch_n, ep_n)
    @inbounds @simd for epoch_idx in 1:ep_n
       Threads.@threads for channel_idx1 in 1:ch_n
            # create half of the matrix
            for channel_idx2 in 1:channel_idx1
                ftest = @views VarianceFTest(eeg.eeg_signals[channel[channel_idx1], :, epoch_idx], eeg.eeg_signals[channel[channel_idx2], :, epoch_idx])
                f[channel_idx1, channel_idx2, epoch_idx] = ftest.F
                p[channel_idx1, channel_idx2, epoch_idx] = pvalue(ftest)
            end
        end
        # copy to the other half
        Threads.@threads for channel_idx1 in 1:(ch_n - 1)
            for channel_idx2 in (channel_idx1 + 1):ch_n
                f[channel_idx1, channel_idx2, epoch_idx] = @views f[channel_idx2, channel_idx1, epoch_idx]
                p[channel_idx1, channel_idx2, epoch_idx] = @views p[channel_idx2, channel_idx1, epoch_idx]
            end
        end
    end

    return (f=f, p=p)
end

"""
    eeg_vartest(eeg1, eeg2; channel1, channel2, epoch1, epoch2)

Calculate variance F-test.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2::NeuroAnalyzer.EEG`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs

# Returns

Named tuple containing:
- `f::Array{Float64, 3}`
- `p::Array{Float64, 3}`
"""
function eeg_vartest(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2)))

    _check_channels(eeg1, channel1)
    _check_channels(eeg2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(eeg1, epoch1)
    _check_epochs(eeg2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    eeg_epoch_len(eeg1) == eeg_epoch_len(eeg2) || throw(ArgumentError("eeg1 and eeg2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)

    f = zeros(ch_n, ch_n, ep_n)
    p = zeros(ch_n, ch_n, ep_n)

    @inbounds @simd for epoch_idx in 1:ep_n
       Threads.@threads for channel_idx1 in 1:ch_n
            for channel_idx2 in 1:ch_n
                ftest = @views VarianceFTest(eeg1.eeg_signals[channel1[channel_idx1], :, epoch1[epoch_idx]], eeg2.eeg_signals[channel2[channel_idx2], :, epoch2[epoch_idx]])
                f[channel_idx1, channel_idx2, epoch_idx] = ftest.F
                p[channel_idx1, channel_idx2, epoch_idx] = pvalue(ftest)
            end
        end
    end

    return (f=f, p=p)
end

"""
    eeg_band_mpower(eeg; channel, f, mt)

Calculate mean and maximum band power and its frequency.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

Named tuple containing:
- `mbp::Matrix{Float64}`: mean band power [μV^2/Hz] per channel per epoch
- `maxfrq::Matrix{Float64}`: frequency of maximum band power [Hz] per channel per epoch
- `maxbp::Matrix{Float64}`: power at maximum band frequency [μV^2/Hz] per channel per epoch
"""
function eeg_band_mpower(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), f::Tuple{Real, Real}, mt::Bool=false)

    fs = eeg_sr(eeg)
    f = tuple_order(f)
    f[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
    f[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be < $(fs / 2)."))

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    mbp = zeros(ch_n, ep_n)
    maxfrq = zeros(ch_n, ep_n)
    maxbp = zeros(ch_n, ep_n)

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            mbp[channel_idx, epoch_idx], maxfrq[channel_idx, epoch_idx], maxbp[channel_idx, epoch_idx] = @views s_band_mpower(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], fs=fs, f=f, mt=mt)
        end
    end

    return (mbp=mbp, maxfrq=maxfrq, maxbp=maxbp)
end

"""
    eeg_rel_psd(eeg; channel, norm, mt, f)

Calculate relative power spectrum density.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `f::Union{Tuple{Real, Real}, Nothing}=nothing`: calculate power relative to frequency range or total power

# Returns

Named tuple containing:
- `psd_pow::Array{Float64, 3}`:powers
- `psd_frq::Array{Float64, 3}`: frequencies
"""
function eeg_rel_psd(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), norm::Bool=false, mt::Bool=false, f::Union{Tuple{Real, Real}, Nothing}=nothing)

    fs = eeg_sr(eeg)
    if f !== nothing
        f = tuple_order(f)
        f[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
        f[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be < $(fs / 2)."))
    end

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    psd_tmp, psd_frq = @views s_rel_psd(eeg.eeg_signals[1, :, 1], fs=fs, norm=norm, mt=mt, f=f)
    psd_pow = zeros(ch_n, length(psd_tmp), ep_n)
    @inbounds @simd for epoch_idx in 1:ep_n
        for channel_idx in 1:ch_n
            psd_pow[channel_idx, :, epoch_idx], _ = @views s_rel_psd(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], fs=fs, norm=norm, mt=mt, f=f)
        end
    end

    return (psd_pow=psd_pow, psd_frq=psd_frq)
end

"""
    eeg_fbsplit(eeg; channel, order, window)

Split EEG signal into frequency bands.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `order::Int64=8`: number of taps for FIR band-pass filter
- `window::Union{Nothing, AbstractVector, Int64}=nothing`: window for `:fir` filter; default is Hamming window, number of taps is calculated using fred harris' rule-of-thumb

# Returns

Named tuple containing:
- `band_names::Vector{Symbol}`
- `band_frq::Vector{Tuple{Real, Real}}`
- `signal_split::Array{Float64, 4}`
"""
function eeg_fbsplit(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), order::Int64=8, window::Union{Nothing, AbstractVector, Int64}=nothing)
    
    band = [:delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher]

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    fs = eeg_sr(eeg)
    signal_split = zeros(length(band), ch_n, eeg_epoch_len(eeg), ep_n)
    band_frq = Vector{Tuple{Real, Real}}()

    @inbounds for band_idx in eachindex(band)
        band_f = eeg_band(eeg, band=band[band_idx])
        push!(band_frq, band_f)
        flt = s_filter_create(fs=fs, fprototype=:fir, ftype=:bp, cutoff=band_f, order=order, window=window, n=eeg_epoch_len(eeg))
        @inbounds @simd for epoch_idx in 1:ep_n
            Threads.@threads for channel_idx in 1:ch_n
                signal_split[band_idx, channel_idx, :, epoch_idx] = @views s_filter_apply(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], flt=flt)
            end
        end
    end

    return (band_names=band, band_frq=band_frq, signal_split=signal_split)
end

"""
    eeg_chdiff(eeg1, eeg2; channel1, channel2, epoch1, epoch2)

Subtract channels.

# Arguments

- `eeg1::NeuroAnalyzer.EEG`
- `eeg2::NeuroAnalyzer.EEG`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs

# Returns

- `ch_diff::Matrix{Float64}`
"""
function eeg_chdiff(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2)))

    _check_channels(eeg1, channel1)
    _check_channels(eeg2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(eeg1, epoch1)
    _check_epochs(eeg2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    eeg_epoch_len(eeg1) == eeg_epoch_len(eeg2) || throw(ArgumentError("eeg1 and eeg2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)

    ch_diff = zeros(ch_n, eeg_epoch_len(eeg1), ep_n)
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            ch_diff[channel_idx, :, epoch_idx] = @views eeg1.eeg_signals[channel1[channel_idx], :, epoch1[epoch_idx]] .- eeg2.eeg_signals[channel2[channel_idx], :, epoch2[epoch_idx]]
        end
    end

    return ch_diff
end

"""
    eeg_cps(eeg; channel, norm)

Calculate cross power spectrum.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `norm::Bool=true`: normalize do dB

# Returns

Named tuple containing:
- `cps_pw::Array{Float64, 4}`: cross power spectrum power
- `cps_ph::Array{Float64, 4}`: cross power spectrum phase (in radians)
- `cps_fq::Vector{Float64}`: cross power spectrum frequencies
"""
function eeg_cps(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), norm::Bool=true)

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    fs = eeg_sr(eeg)
    
    cps_pw_tmp, cps_ph_tmp, cps_fq = @views s2_cps(eeg.eeg_signals[1, :, 1], eeg.eeg_signals[1, :, 1], fs=fs)
    cps_pw = zeros(ch_n, ch_n, length(cps_pw_tmp), ep_n)
    cps_ph = zeros(ch_n, ch_n, length(cps_ph_tmp), ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ch_n, 1))

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx1 in 1:ch_n
           for channel_idx2 in 1:channel_idx1
                cps_pw[channel_idx1, channel_idx2, :, epoch_idx], cps_ph[channel_idx1, channel_idx2, :, epoch_idx], _ = @views s2_cps(eeg.eeg_signals[channel[channel_idx1], :, epoch_idx], eeg.eeg_signals[channel[channel_idx2], :, epoch_idx], fs=fs, norm=norm)
            end

        # update progress bar
        progress_bar == true && next!(p)
        end
    end

    @inbounds @simd for time_idx in 1:size(cps_pw, 3)
        Threads.@threads for epoch_idx in 1:ep_n
            for channel_idx1 in 1:(ch_n - 1)
                for channel_idx2 in (channel_idx1 + 1):ch_n
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

Calculate cross power spectrum.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type]))`: index of channels, default is all EEG channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2))`: default use all epochs
- `norm::Bool=true`: normalize do dB

# Returns

Named tuple containing:
- `cps_pw::Vector{Float64}`: cross power spectrum power
- `cps_ph::Vector{Float64}`: cross power spectrum phase (in radians)
- `cps_fq::Vector{Float64}`: cross power spectrum frequencies
"""
function eeg_cps(eeg1::NeuroAnalyzer.EEG, eeg2::NeuroAnalyzer.EEG; channel1::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg1, type=Symbol(eeg1.eeg_header[:signal_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg2, type=Symbol(eeg2.eeg_header[:signal_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_epoch_n(eeg2)), norm::Bool=true)

    eeg_sr(eeg1) == eeg_sr(eeg2) || throw(ArgumentError("EEG1 and EEG2 must have the same sampling rate."))

    _check_channels(eeg1, channel1)
    _check_channels(eeg2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(eeg1, epoch1)
    _check_epochs(eeg2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    eeg_epoch_len(eeg1) == eeg_epoch_len(eeg2) || throw(ArgumentError("eeg1 and eeg2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)
    fs = eeg_sr(eeg1)

    cps_pw, cps_ph, cps_fq = @views s2_cps(eeg1.eeg_signals[1, :, 1], eeg2.eeg_signals[1, :, 1], fs=fs, norm=norm)

    cps_pw = zeros(ch_n, length(cps_pw), ep_n)
    cps_ph = zeros(ch_n, length(cps_ph), ep_n)
    cps_fq = zeros(ch_n, length(cps_fq), ep_n)
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            cps_pw[channel_idx, :, epoch_idx], cps_ph[channel_idx, :, epoch_idx], cps_fq[channel_idx, :, epoch_idx] = @views s2_cps(eeg1.eeg_signals[channel1[channel_idx], :, epoch1[epoch_idx]], eeg2.eeg_signals[channel2[channel_idx], :, epoch2[epoch_idx]], fs=fs, norm=norm)
        end
    end

    return (cps_pw=cps_pw, cps_ph=cps_ph, cps_fq=cps_fq)
end

"""
    eeg_phdiff(eeg; channel, pad, h)

Calculate phase difference between EEG channels and mean phase of reference `channel`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of reference channels, default is all EEG/MEG channels except the analyzed one
- `avg::Symbol=:phase`: method of averaging:
    - `:phase`: phase is calculated for each reference channel separately and then averaged
    - `:signal`: signals are averaged prior to phase calculation
- `pad::Int64=0`: pad signals with 0s
- `h::Bool=false`: use FFT or Hilbert transformation

# Returns
 
- `ph_diff::Array{Float64, 3}`
"""
function eeg_phdiff(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), avg::Symbol=:phase, pad::Int64=0, h::Bool=false)

    avg in [:phase, :signal] || throw(ArgumentError("avg must be :phase or :signal."))

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    ph_diff = zeros(ch_n, eeg_epoch_len(eeg), ep_n)
    if avg === :phase
        @inbounds @simd for epoch_idx in 1:ep_n
            Threads.@threads for channel_idx in 1:ch_n
                ref_channels = setdiff(channel, channel_idx)
                ph_ref = zeros(length(ref_channels), eeg_epoch_len(eeg))
                for ref_idx in eachindex(ref_channels)
                    if h == true
                        _, _, _, ph = @views s_hspectrum(eeg.eeg_signals[ref_channels[ref_idx], :, epoch_idx], pad=pad)
                    else
                        _, _, _, ph = @views s_spectrum(eeg.eeg_signals[ref_channels[ref_idx], :, epoch_idx], pad=pad)
                    end
                    ph_ref[ref_idx, :] = ph
                end
                ph_ref = vec(mean(ph_ref, dims=1))
                if h == true
                    _, _, _, ph = @views s_hspectrum(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], pad=pad)
                else
                    _, _, _, ph = @views s_spectrum(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], pad=pad)
                end
                ph_diff[channel_idx, :, epoch_idx] = ph - ph_ref
            end
        end
    else
        @inbounds @simd for epoch_idx in 1:ep_n
            Threads.@threads for channel_idx in 1:ch_n
                ref_channels = setdiff(channel, channel_idx)
                signal_m = @views vec(mean(eeg.eeg_signals[ref_channels, :, epoch_idx], dims=1))
                ph_diff[channel_idx, :, epoch_idx] = @views s_phdiff(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], signal_m)
            end
        end
    end

    return ph_diff
end

"""
    eeg_ampdiff(eeg; channel)

Calculate amplitude difference between each channel and mean amplitude.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of reference channels, default is all EEG/MEG channels except the analyzed one

# Returns
 
- `amp_diff::Array{Float64, 3}`
"""
function eeg_ampdiff(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg))

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    amp_diff = zeros(ch_n, eeg_epoch_len(eeg), ep_n)
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            ref_channels = setdiff(channel, channel_idx)
            amp_ref = @views vec(mean(eeg.eeg_signals[ref_channels, :, epoch_idx], dims=1))
            amp_diff[channel_idx, :, epoch_idx] = @views eeg.eeg_signals[channel[channel_idx], :, epoch_idx] - amp_ref
        end
    end

    return amp_diff
end

"""
    eeg_dwt(eeg; channel, wt, type, l)

Perform discrete wavelet transformation (DWT).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: transformation type: 
    - `:sdwt`: Stationary Wavelet Transforms
    - `:acdwt`: Autocorrelation Wavelet Transforms
- `l::Int64=0`: number of levels, default is maximum number of levels available or total transformation

# Returns
 
- `dwt_c::Array{Float64, 4}`: DWT coefficients cAl, cD1, ..., cDl (by rows)
"""
function eeg_dwt(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), wt::T, type::Symbol, l::Int64=0) where {T <: DiscreteWavelet}

    if l == 0
        l = maxtransformlevels(eeg.eeg_signals[1, :, 1])
        _info("Calculating DWT using maximum level: $l.")
    end

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    dwt_c = zeros(ch_n, (l + 1), eeg_epoch_len(eeg), ep_n)
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            dwt_c[channel_idx, :, :, epoch_idx] = @views s_dwt(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], wt=wt, type=type, l=l)
        end
    end

    return dwt_c
end

"""
    eeg_cwt(eeg; channel, wt)

Perform continuous wavelet transformation (CWT).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns
 
- `cwt_c::Array{Float64, 4}`: CWT coefficients (by rows)
"""
function eeg_cwt(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), wt::T) where {T <: CWT}

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    l = size(ContinuousWavelets.cwt(eeg.eeg_signals[1, :, 1], wt), 2)
    cwt_c = zeros(ch_n, l, eeg_epoch_len(eeg), ep_n)
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            cwt_c[channel_idx, :, :, epoch_idx] = @views s_cwt(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], wt=wt)
        end
    end

    return cwt_c
end

"""
    eeg_psdslope(eeg; channel, f, norm, mt)

Calculate PSD linear fit and slope.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `f::Tuple{Real, Real}=(0, eeg_sr(eeg)/2)`: calculate slope of the total power (default) or frequency range f[1] to f[2]
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns

Named tuple containing:
- `lf::Array{Float64, 3}`: linear fit for each channel and epoch
- `psd_slope::Array{Float64, 2}`: slopes of each linear fit
- `frq::Vector{Float64}`: range of frequencies for the linear fits
"""
function eeg_psdslope(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), f::Tuple{Real, Real}=(0, eeg_sr(eeg)/2), norm::Bool=false, mt::Bool=false, nt::Int64=8)

    fs = eeg_sr(eeg)
    f = tuple_order(f)
    f[1] < 0 && throw(ArgumentError("Lower frequency bound must be be ≥ 0."))
    f[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be be < $(fs / 2)."))

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    _, frq = s_psd(eeg.eeg_signals[1, :, 1], fs=fs, norm=norm, mt=mt, nt=nt)
    f1_idx = vsearch(f[1], frq)
    f2_idx = vsearch(f[2], frq)
    lf = zeros(ch_n, length(frq[f1_idx:f2_idx]), ep_n)
    psd_slope = zeros(ch_n, ep_n)

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            pow, _ = s_psd(eeg.eeg_signals[channel[channel_idx], :, epoch_idx], fs=fs, norm=norm, mt=mt, nt=nt)
            _, _, _, _, _, _, lf[channel_idx, :, epoch_idx] = @views linreg(frq[f1_idx:f2_idx], pow[f1_idx:f2_idx])
            psd_slope[channel_idx, epoch_idx] = lf[channel_idx, 2, epoch_idx] - lf[channel_idx, 1, epoch_idx]
        end
    end

    return (lf=lf, psd_slope=psd_slope, frq=frq[f1_idx:f2_idx])
end

"""
    eeg_henv(eeg; channel, d)

Calculate Hilbert spectrum amplitude envelope.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `h_env::Array{Float64, 3}`: Hilbert spectrum amplitude envelope
- `s_t::Vector{Float64}`: signal time
"""
function eeg_henv(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), d::Int64=32)

    _check_channels(eeg, channel)
    _, signal, _, _ = @views eeg_spectrum(eeg_keep_channel(eeg, channel=channel), h=true)

    ch_n = size(signal, 1)
    ep_n = size(signal, 3)
    h_env = similar(signal)
    s_t = eeg.eeg_epoch_time

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            s = @view signal[channel_idx, :, epoch_idx]
            # find peaks
            p_idx = s_findpeaks(s, d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(s))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(s_t[p_idx], s[p_idx])
                try
                    h_env[channel_idx, :, epoch_idx] = model(s_t)
                catch
                    @error "CubicSpline error, using Loess."
                    model = loess(s_t[p_idx], s[p_idx], span=0.5)
                    h_env[channel_idx, :, epoch_idx] = Loess.predict(model, s_t)
                end
            else
                _info("Less than 5 peaks detected, using Loess.")
                model = loess(s_t[p_idx], s[p_idx], span=0.5)
                h_env[channel_idx, :, epoch_idx] = Loess.predict(model, s_t)
            end
            h_env[channel_idx, 1, epoch_idx] = h_env[channel_idx, 2, epoch_idx]
        end
    end
    
    return (h_env=h_env, s_t=s_t)
end

"""
    eeg_henv_mean(eeg; channel, dims, d)

Calculate Hilbert spectrum amplitude envelope: mean and 95% CI.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `h_env_m::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: mean
- `h_env_u::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI upper bound
- `h_env_l::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function eeg_henv_mean(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), dims::Int64, d::Int64=32)
    
    if dims == 1
        eeg_channel_n(eeg) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        eeg_epoch_n(eeg) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        eeg_channel_n(eeg) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        eeg_epoch_n(eeg) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_a, s_t = eeg_henv(eeg, channel=channel, d=d)
    ch_n = size(s_a, 1)
    ep_n = size(s_a, 3)

    if dims == 1
        # mean over channels

        h_env_m = zeros(length(s_t), ep_n)
        h_env_u = zeros(length(s_t), ep_n)
        h_env_l = zeros(length(s_t), ep_n)

        @inbounds @simd for epoch_idx in 1:ep_n
            h_env_m[:, epoch_idx] = mean(s_a[:, :, epoch_idx], dims=1)
            s = std(h_env_m[:, epoch_idx]) / sqrt(length(h_env_m[:, epoch_idx]))
            h_env_u[:, epoch_idx] = @. h_env_m[:, epoch_idx] + 1.96 * s
            h_env_l[:, epoch_idx] = @. h_env_m[:, epoch_idx] - 1.96 * s
        end
    elseif dims == 2
        # mean over epochs

        h_env_m = zeros(length(s_t), ch_n)
        h_env_u = zeros(length(s_t), ch_n)
        h_env_l = zeros(length(s_t), ch_n)

        @inbounds @simd for channel_idx in 1:ch_n
            h_env_m[:, channel_idx] = mean(s_a[channel_idx, :, :], dims=2)
            s = std(h_env_m[:, channel_idx]) / sqrt(length(h_env_m[:, channel_idx]))
            h_env_u[:, channel_idx] = @. h_env_m[:, channel_idx] + 1.96 * s
            h_env_l[:, channel_idx] = @. h_env_m[:, channel_idx] - 1.96 * s
        end
    else
        # mean over channels and epochs

        h_env_m, h_env_u, h_env_l, _ = eeg_henv_mean(eeg, dims=1, d=d)
        h_env_m = mean(h_env_m, dims=2)
        h_env_u = mean(h_env_u, dims=2)
        h_env_l = mean(h_env_l, dims=2)
        h_env_m = reshape(h_env_m, size(h_env_m, 1))
        h_env_u = reshape(h_env_u, size(h_env_u, 1))
        h_env_l = reshape(h_env_l, size(h_env_l, 1))
    end

    return (h_env_m=h_env_m, h_env_u=h_env_u, h_env_l=h_env_l, s_t=s_t)
end

"""
    eeg_henv_median(eeg; channel, dims, d)

Calculate Hilbert spectrum amplitude envelope of `eeg`: median and 95% CI.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `dims::Int64`: median over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `h_env_m::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: median
- `h_env_u::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI upper bound
- `h_env_l::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function eeg_henv_median(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), dims::Int64, d::Int64=32)
    
    if dims == 1
        eeg_channel_n(eeg) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        eeg_epoch_n(eeg) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        eeg_channel_n(eeg) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        eeg_epoch_n(eeg) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_a, s_t = eeg_henv(eeg, channel=channel, d=d)
    ch_n = size(s_a, 1)
    ep_n = size(s_a, 3)

    if dims == 1
        # median over channels

        h_env_m = zeros(length(s_t), ep_n)
        h_env_u = zeros(length(s_t), ep_n)
        h_env_l = zeros(length(s_t), ep_n)

        @inbounds @simd for epoch_idx in 1:ep_n
            h_env_m[:, epoch_idx] = median(s_a[:, :, epoch_idx], dims=1)
            t_idx = s_findpeaks(h_env_m[:, epoch_idx], d=d)
            pushfirst!(t_idx, 1)
            push!(t_idx, length(h_env_m[:, epoch_idx]))
            if length(t_idx) > 4
                model = CubicSpline(s_t[t_idx], h_env_m[t_idx])
                try
                    h_env_m[:, epoch_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(h_env_m[:, epoch_idx]) / sqrt(length(h_env_m[:, epoch_idx]))
            h_env_u[:, epoch_idx] = @. h_env_m[:, epoch_idx] + 1.96 * s
            h_env_l[:, epoch_idx] = @. h_env_m[:, epoch_idx] - 1.96 * s
        end
    elseif dims == 2
        # median over epochs

        h_env_m = zeros(length(s_t), ch_n)
        h_env_u = zeros(length(s_t), ch_n)
        h_env_l = zeros(length(s_t), ch_n)

        @inbounds @simd for channel_idx in 1:ch_n
            h_env_m[:, channel_idx] = median(s_a[channel_idx, :, :], dims=2)
            t_idx = s_findpeaks(h_env_m[:, channel_idx], d=d)
            pushfirst!(t_idx, 1)
            push!(t_idx, length(h_env_m[:, channel_idx]))
            if length(t_idx) > 4
                model = CubicSpline(s_t[t_idx], h_env_m[t_idx])
                try
                    h_env_m[:, channel_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(h_env_m[:, channel_idx]) / sqrt(length(h_env_m[:, channel_idx]))
            h_env_u[:, channel_idx] = @. h_env_m[:, channel_idx] + 1.96 * s
            h_env_l[:, channel_idx] = @. h_env_m[:, channel_idx] - 1.96 * s
        end
    else
        # median over channels and epochs

        h_env_m, h_env_u, h_env_l, _ = eeg_henv_median(eeg, dims=1, d=d)
        h_env_m = median(h_env_m, dims=2)
        h_env_u = median(h_env_u, dims=2)
        h_env_l = median(h_env_l, dims=2)
        h_env_m = reshape(h_env_m, size(h_env_m, 1))
        h_env_u = reshape(h_env_u, size(h_env_u, 1))
        h_env_l = reshape(h_env_l, size(h_env_l, 1))
    end

    return (h_env_m=h_env_m, h_env_u=h_env_u, h_env_l=h_env_l, s_t=s_t)
end

"""
    eeg_apply(eeg; channel, f)

Apply custom function.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all EEG channels
- `f::String`: function to be applied, e.g. `f="mean(eeg, dims=3)"; EEG signal is given using variable `eeg` here.

# Returns

- `out::Array{Float64, 3}`
"""
function eeg_apply(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), f::String)

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    f_tmp = replace(f, "eeg" => "$(eeg.eeg_signals[1, :, 1])")
    out_tmp = eval(Meta.parse(f_tmp))
    out = zeros(eltype(out_tmp), ch_n, length(out_tmp), ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ch_n * ep_n, 1))
    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            f_tmp = replace(f, "eeg" => "$(eeg.eeg_signals[channel[channel_idx], :, epoch_idx])")
            try
                out[channel_idx, :, epoch_idx] = eval(Meta.parse(f_tmp))

            catch
                @error "Formula is incorrect."
            end
            # update progress bar
            progress_bar == true && next!(p)
        end
    end
    return out
end

"""
    eeg_channels_cluster(eeg, cluster)

Return channels belonging to a `cluster` of channels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `cluster::Symbol`: available clusters are:
    - `:f1`: left frontal (F1, F3, F5, F7, AF3, AF7)
    - `:f2`: right frontal (F2, F4, F6, F8, AF4, AF8)
    - `:t1`: left temporal (C3, C5, T7, FC3, FC5, FT7)
    - `:t2`: right temporal (C4, C6, T8, FC4, FC6, FT8)
    - `:c1`: anterior central (Cz, C1, C2, FC1, FC2, FCz)
    - `:c2`: posterior central (Pz, P1, P2, CP1, CP2, CPz)
    - `:p1`: left parietal (P3, P5, P7, CP3, CP5, TP7)
    - `:p2`: right parietal (P4, P6, P8, CP4, CP6, TP8)
    - `:o`: occipital (Oz, O1, O2, POz, PO3, PO4)

# Returns

- `channels::Vector{Int64}`: list of channel numbers belonging to a given cluster of channels
"""
function eeg_channel_cluster(eeg::NeuroAnalyzer.EEG; cluster::Symbol)

    length(eeg_labels(eeg)) == 0 && throw(ArgumentError("EEG does not contain channel labels."))

    _check_var(cluster, [:f1, :f2, :t1, :t2, :c1, :c2, :p1, :p2, :o], "cluster")
    labels = lowercase.(eeg_labels(eeg))
    channels = Int64[]

    cluster === :f1 && (cluster = ["fp1", "f1", "f3", "f5", "f7", "f9", "af3", "af7"])
    cluster === :f2 && (cluster = ["fp2", "f2", "f4", "f6", "f8", "f10", "af4", "af8"])
    cluster === :t1 && (cluster = ["c3", "c5", "t7", "t9", "fc3", "fc5", "ft7", "ft9"])
    cluster === :t2 && (cluster = ["c4", "c6", "t8", "t10", "fc4", "fc6", "ft8", "ft10"])
    cluster === :c1 && (cluster = ["cz", "c1", "c2", "fc1", "fc2", "fcz"])
    cluster === :c2 && (cluster = ["pz", "p1", "p2", "cp1", "cp2", "cpz"])
    cluster === :p1 && (cluster = ["p3", "p5", "p7", "p9", "cp3", "cp5", "tp7", "tp9"])
    cluster === :p2 && (cluster = ["p4", "p6", "p8", "p10", "cp4", "cp6", "tp8", "tp10"])
    cluster === :o && (cluster = ["o1", "o2", "poz", "po3", "po4", "po7", "po8", "po9", "po10"])

    for idx in cluster
        idx in labels && push!(channels, eeg_get_channel(eeg, channel=idx))
    end

    return channels
end

"""
    eeg_erp_peaks(eeg)

Detect a pair of positive and negative peaks of ERP.

# Arguments

- `eeg::NeuroAnalyzer.EEG`:

# Returns
 
- `p::Array{Int64, 2}`: peaks: channels × positive peak position, negative peak position
"""
function eeg_erp_peaks(eeg::NeuroAnalyzer.EEG)

    channels = eeg_signal_channels(eeg)
    erp = eeg_erp(eeg).eeg_signals[channels, :]

    ch_n = size(erp, 1)
    p = zeros(Int64, ch_n, 2)
    @inbounds @simd for channel_idx in 1:ch_n
        pp_pos = @views maximum(erp[channel_idx, :])
        pp_neg = @views minimum(erp[channel_idx, :])
        p[channel_idx, :] = @views [vsearch(pp_pos, erp[channel_idx, :]), vsearch(pp_neg, erp[channel_idx, :])]
    end

    return p
end

"""
    eeg_bands_dwt(eeg; channel, wt, type, n)

Split EEG channel into bands using discrete wavelet transformation (DWT).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64}`: channel number
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: transformation type: 
    - `:sdwt`: Stationary Wavelet Transforms
    - `:acdwt`: Autocorrelation Wavelet Transforms
- `n::Int64=0`: number of bands, default is maximum number of bands available or total transformation

# Returns
 
- `bands::Array{Float64, 4}`: bands from lowest to highest frequency (by rows)
"""
function eeg_bands_dwt(eeg::NeuroAnalyzer.EEG; channel::Int64, wt::T, type::Symbol, n::Int64=0) where {T <: DiscreteWavelet}

    n -= 1
    if n == 0
        n = maxtransformlevels(eeg.eeg_signals[1, :, 1])
        _info("Calculating DWT using maximum level: $n.")
    end
    n < 2 && throw(ArgumentError("n must be ≥ 2."))

    _check_channels(eeg, channel)
    ep_n = eeg_epoch_n(eeg)

    dwt_c = zeros((n + 1), eeg_epoch_len(eeg), ep_n)
    Threads.@threads for epoch_idx in 1:ep_n
        @inbounds dwt_c[:, :, epoch_idx] = @views s_dwt(eeg.eeg_signals[channel, :, epoch_idx], wt=wt, type=type, l=n)
    end
    
    bands = similar(dwt_c)
    bands[1, :, :] = dwt_c[1, :, :]
    bands[2:end, :, :] = dwt_c[end:-1:2, :, :]

    return bands
end