"""
    signal_derivative(signal)

Returns the derivative of the `signal` vector with length same as the signal.

# Arguments

- `signal::Vector{Float64}` - the signal vector to analyze
"""
signal_derivative(signal::Vector{Float64}) = vcat(diff(signal), diff(signal)[end])

"""
    signal_derivative(signal)

Returns the derivative of each the `signal` matrix channels with length same as the signal.

# Arguments
- `signal::Matrix{Float64}` - the signal matrix to analyze (rows: channels, columns: time)
"""
function signal_derivative(signal::Matrix{Float64})
    channels_no = size(signal, 1)
    signal_der = zeros(size(signal))

    for idx in 1:channels_no
        signal_der[idx, :] = signal_derivative(signal[idx, :])
    end

    return signal_der
end

"""
    signal_total_power(signal)

Calculates total power for the `signal` vector.

# Arguments
- `signal::Vector{Float64}` - the signal vector to analyze
"""
function signal_total_power(signal::Vector{Float64}, fs)
    psd = welch_pgram(signal, 4*fs, fs=fs)
    # dx: frequency resolution
    dx = psd.freq[2] - psd.freq[1]
    stp = simpson(psd.power, dx=dx)

    return stp
end

"""
    signal_total_power(signal)

Calculates total power for each the `signal` matrix channels.

# Arguments

- `signal::Matrix{Float64}` - the signal matrix to analyze (rows: channels, columns: time)
"""
function signal_total_power(signal::Matrix{Float64}, fs)
    channels_no = size(signal, 1)
    stp = zeros(size(signal, 1))

    for idx in 1:channels_no
        stp[idx] = signal_total_power(signal[idx, :])
    end

    return stp
end

"""
    signal_band_power(signal, fs, f1, f2)

Calculates absolute band power between frequencies `f1` and `f2` for the `signal` vector.

# Arguments

- `signal::Vector{Float64}` - the signal vector to analyze
- `fs::Float64` - Sampling rate of the signal
- `f1::Float64` - Lower frequency bound
- `f2::Float64` - Upper frequency bound
"""
function signal_band_power(signal::Vector{Float64}, fs::Float64, f1::Float64, f2::Float64)
    psd = welch_pgram(signal, 4*fs, fs=fs)
    frq_idx = [vsearch(Vector(psd.freq), f1), vsearch(Vector(psd.freq), f2)]
    # dx: frequency resolution
    dx = psd.freq[2] - psd.freq[1]
    sbp = simpson(psd.power[frq_idx[1]:frq_idx[2]], psd.freq[frq_idx[1]:frq_idx[2]], dx=dx)

    return sbp
end

"""
    signal_band_power(signal, fs, f1, f2)

Calculates absolute band power between frequencies `f1` and `f2` for each the `signal` matrix channels.

# Arguments

- `signal::Matrix{Float64}` - the signal matrix to analyze
- `fs::Float64` - Sampling rate of the signal
- `f1::Float64` - Lower frequency bound
- `f2::Float64` - Upper frequency bound
"""
function signal_band_power(signal::Matrix{Float64}, fs, f1, f2)
    channels_no = size(signal, 1)
    sbp = zeros(size(signal, 1))

    for idx in 1:channels_no
        sbp[idx] = signal_band_power(signal[idx, :], fs, f1, f2)
    end

    return sbp
end

"""
    signal_make_spectrum(signal, fs)

Returns FFT and DFT sample frequencies for a DFT for the `signal` vector.

# Arguments

- `signal::Vector{Float64}` - the signal vector to analyze
- `fs::Float64` - Sampling rate of the signal
"""
function signal_make_spectrum(signal::Vector{Float64}, fs)
    signal_fft = fft(signal)
    # number of samples
    n = length(signal)
    # time between samples
    d = 1 / fs
    signal_sf = fftfreq(n, d)

    return signal_fft, signal_sf
end

"""
    signal_make_spectrum(signal, fs)

Returns FFT and DFT sample frequencies for a DFT for each the `signal` matrix channels.

# Arguments

- `signal::Matrix{Float64}` - the signal matrix to analyze
- `fs::Float64` - Sampling rate of the signal
"""
function signal_make_spectrum(signal::Matrix{Float64}, fs)
    channels_no = size(signal, 1)
    signal_fft = zeros(ComplexF64, size(signal))
    signal_sf = zeros(size(signal))

    for idx in 1:channels_no
        signal_fft[idx, :], signal_sf[idx, :] = signal_make_spectrum(signal[idx, :], fs)
    end

    return signal_fft, signal_sf
end

"""
    signal_detrend(signal; type=:linear)

Removes linear trend from the `signal` vector.

# Arguments

- `signal::Vector{Float64}` - the signal vector to analyze
- `type::Symbol[:linear, :constant]`, optional
    - `linear` - the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `constant` - the mean of `signal` is subtracted
"""
function signal_detrend(signal::Vector{Float64}; type=:linear)
    type in [:linear, :constant] || throw(ArgumentError("""Trend type must be ":linear" or ":constant"."""))

    if type == :constant
        signal_det = demean(signal)
    else
        A = ones(length(signal))
        coef = A \ signal
        signal_det = @. signal - dot(A, coef)
    end

    return signal_det
end

"""
    signal_detrend(signal; type=:linear)

Removes linear trend for each the `signal` matrix channels.

# Arguments

- `signal::Matrix{Float64}` the signal matrix to analyze
- `type::Symbol[:linear, :constant]`, optional
    - `linear` - the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `constant` - the mean of `signal` is subtracted
"""
function signal_detrend(signal::Matrix{Float64}; type=:linear)
    channels_no = size(signal, 1)
    signal_det = zeros(size(signal))

    for idx in 1:channels_no
        signal_det[idx, :] = signal_detrend(signal[idx, :], type=type)
    end

    return signal_det
end

"""
    signal_ci95(signal; n=3, method=:normal)

Calculates mean, std and 95% confidence interval for each the `signal` matrix channels.

# Arguments

- `signal::Matrix{Float64}` - the signal matrix to analyze
- `n::Int` - number of bootstraps
- `method::Symbol[:normal, :boot]` - use normal method or `n`-times boostrapping
"""
function signal_ci95(signal::Matrix{Float64}; n=3, method=:normal)
    method in [:normal, :boot] || throw(ArgumentError("""Method must be ":normal" or ":boot"."""))

    if method === :normal
        signal_mean = mean(signal, dims=1)'
        signal_sd = std(signal, dims=1)' / sqrt(size(signal, 1))
        upper_bound = signal_mean + 1.96 * signal_sd
        lower_bound = signal_mean - 1.96 * signal_sd
    else
        signal_tmp1 = zeros(size(signal, 1) * n, size(signal, 2))
        Threads.@threads for idx1 in 1:size(signal, 1) * n
            signal_tmp2 = zeros(size(signal))
            sample_idx = rand(1:size(signal, 1), size(signal, 1))
            for idx2 in 1:size(signal, 1)
                signal_tmp2[idx2, :] = signal[sample_idx[idx2], :]'
            end
            signal_tmp1[idx1, :] = mean(signal_tmp2, dims=1)
        end
        signal_mean = mean(signal_tmp1, dims=1)'
        signal_sd = std(signal_tmp1, dims=1)' / sqrt(size(signal_tmp1, 1))
        signal_sorted = sort(signal_tmp1, dims=1)
        lower_bound = signal_sorted[round(Int, 0.025 * size(signal_tmp1, 1)), :]
        upper_bound = signal_sorted[round(Int, 0.975 * size(signal_tmp1, 1)), :]
    end

    return Vector(signal_mean[:, 1]), Vector(signal_sd[:, 1]), Vector(upper_bound[:, 1]), Vector(lower_bound[:, 1])
end

"""
    signal_mean(signal1, signal2)

Calculates mean and 95% confidence interval for 2 signals.

# Arguments

- `signal1::Matrix{Float64}` - the signal 1 matrix to analyze
- `signal2:Matrix{Float64}` - the signal 2 matrix to analyze
"""
function signal_mean(signal1::Matrix{Float64}, signal2::Matrix{Float64})
    size(signal1) != size(signal2) && throw(ArgumentError("Both matrices must be of the same as size."))

    signal1_mean = mean(signal1, dims=1)'
    signal2_mean = mean(signal2, dims=1)'
    signals_mean = signal1_mean - signal2_mean
    signal1_sd = std(signal1, dims=1) / sqrt(size(signal1, 1))
    signal2_sd = std(signal2, dims=1) / sqrt(size(signal2, 1))
    signals_mean_sd = sqrt.(signal1_sd.^2 .+ signal2_sd.^2)'

    return Vector(signals_mean[:, 1]), Vector(signals_mean_sd[:, 1]), Vector((signals_mean + 1.96 * signals_mean_sd)[:, 1]), Vector((signals_mean - 1.96 * signals_mean_sd)[:, 1])
end

"""
    signal_difference(signal1::Matrix, signal2::Matrix; n=3, method=:absdiff)

Calculates mean difference and 95% confidence interval for 2 signals.

# Arguments

- `signal1::Matrix` - the signal 1 matrix to analyze
- `signal2:Matrix` - the signal 2 matrix to analyze
- `n::Int` - number of bootstraps.
- `method::Symbol[:absdiff, :diff2int]`
    - `:absdiff` - maximum difference
    - `:diff2int` - integrated area of the squared difference
"""
function signal_difference(signal1::Matrix, signal2::Matrix; n=3, method=:absdiff)
    size(signal1) != size(signal2) && throw(ArgumentError("Both matrices must be of the same as size."))
    method in [:absdiff, :diff2int] || throw(ArgumentError("""Method must be ":absdiff" or ":diff2int"."""))

    signal1_mean = mean(signal1, dims=1)'
    signal2_mean = mean(signal2, dims=1)'

    if method === :absdiff
        # statistic: maximum difference
        signals_diff = signal1_mean - signal2_mean
        signals_statistic_single = maximum(abs.(signals_diff))
    else
        # statistic: integrated area of the squared difference
        signals_diff_squared = (signal1_mean - signal2_mean).^2
        signals_statistic_single = simpson(signals_diff_squared)
    end

    signals = [signal1; signal2]
    signals_statistic = zeros(size(signal1, 1) * n)

    Threads.@threads for idx1 in 1:(size(signal1, 1) * n)
        signals_tmp1 = zeros(size(signal1, 1), size(signal1, 2))
        sample_idx = rand(1:size(signals, 1), size(signals, 1))
        # sample_idx = sample_idx[1:1000]
        for idx2 in 1:size(signal1, 1)
            signals_tmp1[idx2, :] = signals[sample_idx[idx2], :]'
        end
        signal1_mean = mean(signals_tmp1, dims=1)
        signals_tmp1 = zeros(size(signal1, 1), size(signal1, 2))
        sample_idx = rand(1:size(signals, 1), size(signals, 1))
        # sample_idx = sample_idx[1:1000]
        for idx2 in 1:size(signal1, 1)
            signals_tmp1[idx2, :] = signals[sample_idx[idx2], :]'
        end
        signal2_mean = mean(signals_tmp1, dims=1)
        if method === :absdiff
            # statistic: maximum difference
            signals_diff = signal1_mean - signal2_mean
            signals_statistic[idx1] = maximum(abs.(signals_diff))
        else
            # statistic: integrated area of the squared difference
            signals_diff_squared = (signal1_mean - signal2_mean).^2
            signals_statistic[idx1] = simpson(signals_diff_squared)
        end
    end

    p = length(signals_statistic[signals_statistic .> signals_statistic_single]) / size(signal1, 1) * n

    return signals_statistic, signals_statistic_single, p
end

"""
   signal_autocov(signal; lag=1, remove_dc=false, normalize=false)

Calculates autocovariance of the `signal` vector.

# Arguments

- `signal::Vector{Float64}` - the signal vector to analyze
- `lag::Int` - lags range is `-lag:lag`
- `remove_dc::Bool[true, false]` - demean `signal` prior to calculations
- `normalize::Bool[true, false]` - normalize autocovariance
"""
function signal_autocov(signal::Vector{Float64}; lag=1, remove_dc=false, normalize=false)
    signal_lags = collect(-lag:lag)

    if remove_dc == true
        signal_demeaned = signal .- mean(signal)
    else
        signal_demeaned = signal
    end

    signal_ac = zeros(length(signal_lags))

    for idx in 1:length(signal_lags)
        if signal_lags[idx] == 0
            # no lag
            signal_lagged = signal_demeaned
            signals_mul = signal_demeaned .* signal_lagged
        elseif signal_lags[idx] > 0
            # positive lag
            signal_lagged = signal_demeaned[1:(end - signal_lags[idx])]
            signals_mul = signal_demeaned[(1 + signal_lags[idx]):end] .* signal_lagged
        elseif signal_lags[idx] < 0
            # negative lag
            signal_lagged = signal_demeaned[(1 + abs(signal_lags[idx])):end]
            signals_mul = signal_demeaned[1:(end - abs(signal_lags[idx]))] .* signal_lagged
        end
        signals_sum = sum(signals_mul)
        if normalize == true
            signal_ac[idx] = signals_sum / length(signal)
        else
            signal_ac[idx] = signals_sum
        end
    end

    return signal_ac, signal_lags
end

"""
   signal_autocov(signal; lag=1, remove_dc=false, normalize=false)

Calculates autocovariance of each the `signal` matrix channels.

# Arguments

- `signal::Matrix{Float64}` - the signal vector to analyze
- `lag::Int` - lags range is `-lag:lag`
- `remove_dc::Bool` - demean signal prior to analysis
- `normalize::Bool` - normalize autocovariance
"""
function signal_autocov(signal::Matrix{Float64}; lag=1, remove_dc=false, normalize=false)
    signal_lags = collect(-lag:lag)
    channels_no = size(signal, 1)
    signal_ac = zeros(channels_no, length(signal_lags))

    for idx in 1:channels_no
        signal_ac[idx, :], _ = signal_autocov(signal[idx, :],
                                              lag=lag,
                                              remove_dc=remove_dc,
                                              normalize=normalize)
    end

    return signal_ac, signal_lags
end

"""
   signal_crosscov(signal1, signal2; lag=1, remove_dc=false, normalize=false)

Calculates cross-covariance between `signal1` and `signal2` vectors.

# Arguments

- `signal1::Vector{Float64}` - the signal 1 vector to analyze
- `signal2::Vector{Float64}` - the signal 2 vector to analyze
- `lag::Int` - lags range is `-lag:lag`
- `remove_dc::Bool` - demean signal prior to analysis
- `normalize::Bool` - normalize cross-covariance
"""
function signal_crosscov(signal1::Vector{Float64}, signal2::Vector{Float64}; lag=1, remove_dc=false, normalize=false)
    length(signal1) != length(signal2) && throw(ArgumentError("Both vectors must be of the same as length."))

    lags = collect(-lag:lag)

    if remove_dc == true
        signal_demeaned1 = signal1 .- mean(signal1)
        signal_demeaned2 = signal2 .- mean(signal2)
    else
        signal_demeaned1 = signal1
        signal_demeaned2 = signal2
    end

    ac = zeros(length(lags))

    for idx in 1:length(lags)
        if lags[idx] == 0
            # no lag
            signal_lagged = signal_demeaned2
            signals_mul = signal_demeaned1 .* signal_lagged
        elseif lags[idx] > 0
            # positive lag
            signal_lagged = signal_demeaned2[1:(end - lags[idx])]
            signals_mul = signal_demeaned1[(1 + lags[idx]):end] .* signal_lagged
        elseif lags[idx] < 0
            # negative lag
            signal_lagged = signal_demeaned2[(1 + abs(lags[idx])):end]
            signals_mul = signal_demeaned1[1:(end - abs(lags[idx]))] .* signal_lagged
        end
        signals_sum = sum(signals_mul)
        if normalize == true
            ac[idx] = signals_sum / length(signal1)
        else
            ac[idx] = signals_sum
        end
    end

    return ac, lags
end

"""
   signal_crosscov(signal1, signal2; lag=1, remove_dc=false, normalize=false)

Calculates cross-covariance between same channels in `signal1` and `signal2` matrices.

# Arguments

- `signal1::Matrix{Float64}` - the signal 1 matrix to analyze
- `signal2::Matrix{Float64}` - the signal 2 matrix to analyze
- `lag::Int` - lags range is `-lag:lag`
- `remove_dc::Bool` - demean signal prior to analysis
- `normalize::Bool` - normalize cross-covariance
"""
function signal_crosscov(signal1::Matrix{Float64}, signal2::Matrix{Float64}; lag=1, remove_dc=false, normalize=false)
    size(signal1) != size(signal2) && throw(ArgumentError("Both matrices must be of the same as size."))

    signal_lags = collect(-lag:lag)
    channels_no = size(signal1, 1)
    signal_ac = zeros(channels_no, length(signal_lags))

    for idx in 1:channels_no
        signal_ac[idx, :], _ = signal_crosscov(signal1[idx, :],
                                               signal2[idx, :],
                                               lag=lag,
                                               remove_dc=remove_dc,
                                               normalize=normalize)
    end

    return signal_ac, signal_lags
end

"""
   signal_crosscov(signal; lag=1, remove_dc=false, normalize=false)

Calculates cross-covariance for all channels in the `signal` matrix. Return matrix of cross-covariances:
signal_1_channel_1 vs signal_2_channel_1, signal_1_channel_1 vs signal_2_channel_2, signal_1_channel_1 vs signal_2_channel_3, ..., signal_1_channel_n vs signal_2_channel_n.

# Arguments

- `signal::Matrix{Float64}` - the signal matrix to analyze
- `lag::Int` - lags range is `-lag:lag`
- `remove_dc::Bool` - demean `signal` prior to analysis
- `normalize::Bool` - normalize cross-covariance
"""
function signal_crosscov(signal::Matrix{Float64}; lag=1, remove_dc=false, normalize=false)
    signal_lags = collect(-lag:lag)
    channels_no = size(signal, 1)
    signal_ac_packed = Array{Vector{Float64}}(undef, channels_no, channels_no)
    signal_ac = zeros(channels_no^2, channels_no)

    for idx1 in 1:channels_no
        for idx2 in 1:channels_no
            signal_ac_packed[idx1, idx2], _ = signal_crosscov(signal[idx1, :],
                                                              signal[idx2, :],
                                                              lag=lag,
                                                              remove_dc=remove_dc,
                                                              normalize=normalize)
        end
    end

    for idx in 1:channels_no^2
        signal_ac[idx, :] = signal_ac_packed[idx]
    end

    return reverse(signal_ac), signal_lags
end

"""
    signal_spectrum(signal; pad=0)

Calculates FFT, amplitudes, powers and phases of the `signal` vector.

# Arguments

- `signal::Vector{Float64}` - the signal vector to analyze
- `pad::Int` - pad the `signal` with `pad` zeros
"""
function signal_spectrum(signal::Vector{Float64}; pad=0)
    pad < 0 && throw(ArgumentError("""Value of "pad" cannot be negative."""))

    if pad == 0
        signal_fft = fft(signal)
    else
        signal_fft = fft0(signal, pad)
    end

    # normalize
    signal_fft ./= length(signal)

    # amplitudes
    signal_amplitudes = @. 2 * abs(signal_fft)

    # power
    signal_powers = signal_amplitudes.^2

    # phases
    signal_phases = atan.(imag(signal_fft), real(signal_fft))

    return signal_fft, signal_amplitudes, signal_powers, signal_phases
end

"""
    signal_spectrum(signal; pad=0)

Calculates FFT, amplitudes, powers and phases for each channel of the `signal` matrix.

# Arguments

- `signal::Vector{Float64}` - the signal matrix to analyze
- `pad::Int` - pad the `signal` with `pad` zeros
"""
function signal_spectrum(signal::Matrix{Float64}; pad=0)
    pad < 0 && throw(ArgumentError("""Value of "pad" cannot be negative."""))

    channels_no = size(signal, 1)

    signal_fft = zeros(ComplexF64, size(signal))
    signal_amplitudes = zeros(size(signal))
    signal_powers = zeros(size(signal))
    signal_phases = zeros(size(signal))

    for idx in 1:channels_no
        signal_fft[idx, :], signal_amplitudes[idx, :], signal_powers[idx, :], signal_phases[idx, :] = 
        signal_spectrum(signal[idx, :], pad=pad)
    end

    return signal_fft, signal_amplitudes, signal_powers, signal_phases
end

"""
    signal_epoch(signal; epoch_no, epoch_len, average=true)

Splits `signal` vector into epochs.

# Arguments

- `signal::Vector{Float64}` - the signal vector to analyze
- `epoch_no::Int` - number of epochs
- `epoch_len::Int` - epoch length in samples
- `average::Bool` - average all epochs, returns one averaged epoch; if false than returns array of epochs, each row is one epoch
"""
function signal_epoch(signal::Vector{Float64}; epoch_no=nothing, epoch_len=nothing, average=true)
    (epoch_len === nothing && epoch_no === nothing) && throw(ArgumentError("Either number of epochs or epoch length must be set."))
    (epoch_len != nothing && epoch_no != nothing) && throw(ArgumentError("Both number of epochs and epoch length cannot be set."))

    if epoch_no === nothing
        epoch_no = length(signal) ÷ epoch_len
    else
        epoch_len = length(signal) ÷ epoch_no
    end

    epochs = zeros(epoch_no, epoch_len)

    idx1 = 1
    for idx2 in 1:epoch_len:(epoch_no * epoch_len - 1)
        epochs[idx1, :] = signal[idx2:(idx2 + epoch_len - 1)]
        idx1 += 1
    end

    if average == true
        epochs = Vector(mean(epochs, dims=1)[1, :])
    end

    return epochs
end

"""
    signal_epoch(signal, n; epoch_no=nothing, epoch_len=nothing, average=true)

Splits `signal` matrix into epochs.

# Arguments

- `signal::Matrix{Float64}` - the signal matrix to analyze
- `epoch_no::Int` - number of epochs
- `epoch_len::Int` - epoch length in samples
- `average::Bool` - average all epochs, returns one averaged epoch; if false than returns array of epochs, each row is one epoch
"""
function signal_epoch(signal::Matrix; epoch_no=nothing, epoch_len=nothing, average=true)
    (epoch_len === nothing && epoch_no === nothing) && throw(ArgumentError("Either number of epochs or epoch length must be set."))
    (epoch_len != nothing && epoch_no != nothing) && throw(ArgumentError("Both number of epochs and epoch length cannot be set."))

    channels_no = size(signal, 1)

    if epoch_no === nothing
        epoch_no = size(signal, 2) ÷ epoch_len
    else
        epoch_len = size(signal, 2) ÷ epoch_no
    end

    epochs = zeros(channels_no, epoch_len, epoch_no)

    idx1 = 1
    for idx2 in 1:epoch_len:(epoch_no * epoch_len - 1)
        epochs[:, :, idx1] = signal[:, idx2:(idx2 + epoch_len - 1)]
        idx1 += 1
    end

    if average == true
        epochs = mean(epochs, dims=3)[:, :]
    end

    return epochs
end

"""
    signal_filter_butter(signal; filter_type, cutoff, fs, poles=8)

Filters `signal` vector using Butterworth filter.

# Arguments

- `signal::Vector{Float64}` - the signal vector to analyze
- `filter_type::Symbol[:lp, :hp, :bp, :bs]` - filter type
- `cutoff::Float64` - filter cutoff in Hz (tuple or vector for `:bp` and `:bs`)
- `fs::Float64` - sampling rate
- `poles::Int` - filter pole
"""
function signal_filter_butter(signal::Vector{Float64}; filter_type, cutoff, fs, poles=8)
    filter_type in [:lp, :hp, :bp, :bs] || throw(ArgumentError("""Filter type must be ":bp", ":hp", ":bp" or ":bs"."""))

    if filter_type == :lp
        responsetype = Lowpass(cutoff; fs=fs)
        prototype = Butterworth(poles)
    elseif filter_type == :hp
        responsetype = Highpass(cutoff; fs=fs)
        prototype = Butterworth(poles)
    elseif filter_type == :bp
        length(cutoff) < 2 && throw(ArgumentError("For band-pass filter two frequencies must be given."))
        responsetype = Bandpass(cutoff[1], cutoff[2]; fs=fs)
        prototype = Butterworth(poles)
    elseif filter_type == :bs
        length(cutoff) < 2 && throw(ArgumentError("For band-stop filter two frequencies must be given."))
        responsetype = Bandstop(cutoff[1], cutoff[2]; fs=fs)
        prototype = Butterworth(poles)
    end

    filter = digitalfilter(responsetype, prototype)
    signal_filtered = filt(filter, signal)

    return signal_filtered
end

"""
    signal_filter_butter(signal; filter_type, cutoff, fs, poles=8)

Filters `signal` matrix using Butterworth filter.

# Arguments

- `signal::Matrix{Float64}` - the signal matrix to analyze
- `filter_type::Symbol[:lp, :hp, :bp, :bs]` - filter type
- `cutoff::Float64` - filter cutoff in Hz (tuple or vector for `:bp` and `:bs`)
- `fs::Float64` - sampling rate
- `poles::Int` - filter pole
"""
function signal_filter_butter(signal::Matrix{Float64}; filter_type, cutoff, fs, poles=8)
    filter_type in [:lp, :hp, :bp, :bs] || throw(ArgumentError("""Filter type must be ":bp", ":hp", ":bp" or ":bs"."""))

    channels_no = size(signal, 1)
    signal_filtered = zeros(size(signal))

    for idx in 1:channels_no
        signal_filtered[idx, :] = signal_filter_butter(signal[idx, :], filter_type=filter_type, cutoff=cutoff, fs=fs, poles=poles)
    end

    return signal_filtered
end

"""
    signal_plot(t, signal; offset=1, labels=[], normalize=false, xlabel="Time [s]", ylabel="Amplitude [μV]", yamp=nothing, figure::String="")

Plots `signal` against time vector `t`.

# Arguments

- `t::Vector{Float64}` - the time vector
- `signal::Vector{Float64}` - the signal vector
- `offset::Float64` - displayed segment offset in samples
- `labels::Vector{String}` - channel labels vector
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis lable
- `yamp::Float64` - y-axis limits (-yamp:yamp)
- `figure::String` - name of the output figure file
"""
function signal_plot(t, signal::Vector{Float64}; offset=1, labels=[], xlabel="Time [s]", ylabel="Amplitude [μV]", yamp=nothing, figure::String="")

    if typeof(t) == UnitRange{Int64}
        t = float(collect(t))
    end

    if yamp === nothing
        yamp = maximum(signal)
        yamp = ceil(Int64, yamp)
    end

    p = plot(t, signal[ofset:(offset + length(t))], xlabel=xlabel, ylabel=ylabel, legend=false, t=:line, c=:black, ylims=(-yamp, yamp))

    # TO DO: catching error while saving
    figure !== "" && (savefig(p, figure))

    return p, signal
end

"""
    signal_plot(t, signal; offset=1, labels=[], normalize=false, xlabel="Time [s]", ylabel="Channels", figure::String="")

Plots `signal` matrix against time vector `t`.

# Arguments

- `t::Vector{Float64}` - the time vector
- `signal::Matrix{Float64}` - the signal matrix
- `offset::Float64` - displayed segment offset in samples
- `labels::Vector{String}` - channel labels vector
- `normalize::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis label
- `figure::String` - name of the output figure file
"""
function signal_plot(t, signal::Matrix{Float64}; offset=1, labels=[], normalize=false, xlabel="Time [s]", ylabel="Channels", figure::String="")
    
    if typeof(t) == UnitRange{Int64}
        t = float(collect(t))
    end
    
    channels_no = size(signal, 1)

    # reverse so 1st channel is on top
    signal = reverse(signal, dims = 1)

    if normalize == true
        # normalize and shift so all channels are visible
        variances = var(signal, dims=2)
        mean_variance = mean(variances)
        for idx in 1:channels_no
            signal[idx, :] = (signal[idx, :] .- mean(signal[idx, :])) ./ mean_variance .+ (idx - 1)
        end
    end

    # plot channels
    p = plot(xlabel=xlabel, ylabel=ylabel, ylim=(-0.5, channels_no-0.5))
    for idx in 1:channels_no
        p = plot!(t, signal[idx, offset:(offset + length(t))], legend=false, t=:line, c=:black)
    end
    p = plot!(p, yticks = (channels_no-1:-1:0, labels))

    plot(p)

    # TO DO: catching error while saving
    figure !== "" && (savefig(p, figure))

    return p, signal
end

"""
    signal_drop_channel(signal, channels)

Removes `channels` from the `signal` matrix.

# Arguments

- `signal::Matrix{Float64}` - the signal matrix
- `channels::Float64` - channels to be removed, vector of numbers or range
"""
function signal_drop_channel(signal::Matrix, channels)
    if typeof(channels) == UnitRange{Int64}
        channels = collect(channels)
    end

    channels = sort!(channels, rev=true)
    signal = signal[setdiff(1:end, (channels)), :]

    return signal
end

"""
    signal_reference_channel(signal, reference)

Re-references channels of the `signal` matrix to specific signal channel.

# Arguments

- `signal::Matrix{Float64}` - the signal matrix
- `reference::Float64` - index of channels used as reference; if multiple channels are specified, their average is used as the reference
"""
function signal_reference_channel(signal::Matrix, reference_idx)
    if typeof(reference_idx) == UnitRange{Int64}
        reference_idx = collect(reference_idx)
    end

    channels_no = size(signal, 1)
    signal_referenced = zeros(size(signal))

    channels_list = collect(1:channels_no)
    for idx in 1:length(reference_idx)
        if (reference_idx[idx] in channels_list) == false
            throw(ArgumentError("Reference channel index does not match signal channels."))
        end
    end

    reference_channel = vec(mean(signal[reference_idx, :], dims=1))

    for idx in 1:channels_no
        signal_referenced[idx, :] = signal[idx, :] .- reference_channel
    end

    return signal_referenced
end

"""
    signal_reference_channel(signal)

Re-references channels of the `signal` matrix to common average reference.

# Arguments

- `signal::Matrix{Float64}` - the signal matrix
"""
function signal_reference_car(signal::Matrix{Float64})
    channels_no = size(signal, 1)
    signal_referenced = zeros(size(signal))
    reference_channel = vec(mean(signal, dims=1))

    for idx in 1:channels_no
        signal_referenced[idx, :] = signal[idx, :] .- reference_channel
    end

    return signal_referenced
end

"""
    signal_taper(signal, taper)

Taper channels of the `signal` matrix with `taper`.

# Arguments

- `signal::Matrix{Float64}` - the signal matrix
- `taper::Vector`
"""
function signal_taper(signal::Matrix, taper::Vector)
    channels_no = size(signal, 1)

    length(taper) == size(signal, 2) || throw(ArgumentError("Taper length and signal length must be equal."))

    for idx in 1:channels_no
        signal_tapered[idx, :] = signal[idx, :] .* taper
    end

    return signal_tapered
end

"""
    signal_taper(signal, taper)

Taper channels of the `signal` vector with `taper`.

# Arguments

- `signal::Matrix{Float64}` - the signal matrix
- `taper::Vector`
"""
function signal_taper(signal::Vector, taper::Vector)
    length(taper) == length(signal) || throw(ArgumentError("Taper length and signal length must be equal."))
    signal_tapered = signal .* taper

    return signal_tapered
end

"""
    signal_demean(signal)

Removes mean value (DC offset) from the `signal` vector.

# Arguments

- `signal::Vector{Float64}` the signal matrix to analyze
"""
signal_demean(signal::Vector{Float64}) = signal .- mean(signal)

"""
    signal_demean(signal)

Removes mean value (DC offset) for each the `signal` matrix channels.

# Arguments

- `signal::Matrix{Float64}` the signal matrix to analyze
"""
function signal_demean(signal::Matrix{Float64})
    channels_no = size(signal, 1)
    signal_demeaned = zeros(size(signal))

    for idx in 1:channels_no
        signal_demeaned[idx, :] = signal_demean(signal[idx, :])
    end

    return signal_demeaned
end

"""
    signal_normalize_mean(signal)

Normalize (scales around the mean) `signal` vector.
"""
signal_normalize_mean(signal::Vector{Float64}) = (signal .- mean(signal)) ./ std(signal)

"""
    signal_normalize_mean(signal)

Normalize (scales around the mean) `signal` matrix.
"""
function signal_normalize_mean(signal::Matrix{Float64})
    channels_no = size(signal, 1)
    signal_normalized = zeros(size(signal))

    for idx in 1:channels_no
        signal_normalized[idx, :] = signal_normalize_mean(signal[idx, :])
    end

    return signal_normalized
end

"""
    signal_normalize_minmax(signal)

Normalize (to 0…1) `signal` vector.
"""
signal_normalize_minmax(signal::Vector{Float64}) = (signal .- minimum(signal)) ./ (maximum(signal) - minimum(signal))

"""
    signal_normalize_minmax(signal)

Normalize (to 0…1) each the `signal` matrix channel.
"""
function signal_normalize_minmax(signal::Matrix{Float64})
    channels_no = size(signal, 1)
    signal_normalized = zeros(size(signal))

    for idx in 1:channels_no
        signal_normalized[idx, :] = signal_normalize_minmax(signal[idx, :])
    end

    return signal_normalized
end

"""
   signal_cov(signal1, signal2; normalize=true)

Calculates covariance between `signal1` and `signal2` vectors.

# Arguments

- `signal1::Vector{Float64}` - the signal 1 vector to analyze
- `signal2::Vector{Float64}` - the signal 2 vector to analyze
- `normalize::Bool` - normalize covariance
"""
function signal_cov(signal1::Vector{Float64}, signal2::Vector{Float64}; normalize=true)
    length(signal1) != length(signal2) && throw(ArgumentError("Both vectors must be of the same as length."))

    result = cov(signal1 * signal2')

    # divide so that components are centered at (0, 0)
    normalize == true && (result = result ./ length(signal1))

    return result
end

"""
   signal_cov(signal; normalize=true)

Calculates covariance between all channels of the `signal` matrix.

# Arguments

- `signal::Matrix{Float64}` - the signal matrix to analyze
- `normalize::Bool` - normalize covariance
"""
function signal_cov(signal::Matrix{Float64}; normalize=true)
    # channels-vs-channels
    result = cov(signal')

    # divide so that components are centered at (0, 0)
    normalize == true && (result = result ./ size(signal, 2))

    return result
end

"""
   signal_cor(signal)

Calculates correlation coefficients between all channels of the `signal` matrix.

# Arguments

- `signal::Matrix{Float64}` - the signal matrix to analyze
"""
function signal_cor(signal::Matrix{Float64})
    # channels-vs-channels
    result = cor(signal')

    return result
end

"""
    signal_add_noise(signal)

Add random noise to the `signal` vector.

# Arguments

- `signal::Vector{Float64}` - the signal vector
"""
signal_add_noise(x::Vector{Float64}) = x .+ mean(x) .* rand(length(x))

"""
    signal_add_noise(signal)

Add random noise to each the `signal` matrix channel.

# Arguments

- `signal::Matrix{Float64}` - the signal matrix
"""
function signal_add_noise(signal::Matrix{Float64})
    channels_no = size(signal, 1)
    signal_noise = zeros(size(signal))

    for idx in 1:channels_no
        signal_noise[idx, :] = signal_add_noise(signal[idx, :])
    end

    return signal_noise
end

"""
    signal_upsample(signal; t, new_sr)

Upsamples the`signal` vector to `new_sr` sampling frequency.

# Arguments

- `signal::Vector{Float64}` - the signal vector
- `t::AbstractRange` - the time range
- `new_sr::Int64` - new sampling rate
"""
function signal_upsample(signal::Vector{Float64}; t::AbstractRange, new_sr::Int64)
    # sampling interval
    dt = t[2] - t[1]
    # sampling rate
    sr = 1 / dt
    new_sr < sr && throw(ArgumentError("New sampling rate mu be larger than signal sampling rate."))
    new_sr == sr && return(signal)

    # interpolate
    signal_interpolation = CubicSplineInterpolation(t, signal)
    t = t[1]:1/new_sr:t[end]
    signal_upsampled = signal_interpolation(t)

    return signal_upsampled
end

"""
    signal_upsample(signal; t, new_sr)

Upsamples all channels of the`signal` matrix to `new_sr` sampling frequency.

# Arguments

- `signal::Matrix{Float64}` - the signal vector
- `t::AbstractRange` - the time range
- `new_sr::Int64` - new sampling rate
"""
function signal_upsample(signal::Matrix{Float64}; t::AbstractRange, new_sr::Int64)
    channels_no = size(signal, 1)

    for idx in 1:channels_no
        signal_upsampled[idx, :] = signal_upsample(signal[idx, :], t=t, new_sr=new_sr)
    end

    return signal_upsampled
end