"""
    signal_derivative(signal)

Returns the derivative of the `signal` with length same as the signal.

# Arguments

- `signal::Union{Vector{Int64}, Vector{Float64}}`

# Returns

`signal_derivative::Union{Vector{Int64}, Vector{Float64}}`

"""
signal_derivative(signal::Union{Vector{Int64}, Vector{Float64}}) = vcat(diff(signal), diff(signal)[end])

"""
    signal_derivative(signal)

Returns the derivative of each the `signal` channels with length same as the signal.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `signal_derivative::Array{Float64, 3}`
"""
function signal_derivative(signal::Array{Float64, 3})
    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)
    signal_der = zeros(size(signal))

    Threads.@threads for epoch in 1:epochs_no
        for idx in 1:channels_no
            signal_der[idx, :, epoch] = signal_derivative(signal[idx, :, epoch])
        end
    end

    return signal_der
end

"""
    signal_total_power(signal; fs)

Calculates total power for the `signal`.

# Arguments

- `signal::Vector{Float64}`
- `fs::Int64` - sampling rate

# Returns

- `stp::Float64`
"""
function signal_total_power(signal::Vector{Float64}; fs::Int64)
    fs < 1 && throw(ArgumentError("Sampling rate must be ≥ 1 Hz."))

    psd = welch_pgram(signal, 4*fs, fs=fs)
    # dx: frequency resolution
    dx = psd.freq[2] - psd.freq[1]
    stp = simpson(psd.power, dx=dx)

    return stp
end

"""
    signal_total_power(signal; fs)

Calculates total power for each the `signal` channels.

# Arguments

- `signal::Array{Float64, 3}`
- `fs::Int64` - sampling rate

# Returns

- `stp::Matrix{Float64}`
"""
function signal_total_power(signal::Array{Float64, 3}; fs::Int64)
    fs < 1 && throw(ArgumentError("Sampling rate must be ≥ 1 Hz."))

    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)
    stp = zeros(channels_no, epochs_no)

    Threads.@threads for epoch in 1:epochs_no
        for idx in 1:channels_no
            stp[idx, epoch] = signal_total_power(signal[idx, :, epoch], fs=fs)
        end
    end

    return stp
end

"""
    signal_band_power(signal; fs, f1, f2)

Calculates absolute band power between frequencies `f1` and `f2` for the `signal`.

# Arguments

- `signal::Vector{Float64}`
- `fs::Int64` - sampling rate of the signal
- `f1::Union{Int64, Float64}` - lower frequency bound
- `f2::Union{Int64, Float64}` - upper frequency bound

# Returns

- `sbp::Float64`
"""
function signal_band_power(signal::Vector{Float64}; fs::Int64, f1::Union{Int64, Float64}, f2::Union{Int64, Float64})
    fs < 1 && throw(ArgumentError("Sampling rate must be ≥ 1 Hz."))
    f1 > f2 && throw(ArgumentError("Lower frequency bound must be lower than higher frequency bound."))
    f1 < 0 && throw(ArgumentError("Lower frequency bound must be be ≥ 1 Hz."))
    f2 > fs / 2 && throw(ArgumentError("Upper frequency bound must be be < Nyquist frequency ($(fs / 2) Hz)."))

    psd = welch_pgram(signal, 4*fs, fs=fs)
    psd_freq = Vector(psd.freq)
    frq_idx = [vsearch(psd_freq, f1), vsearch(psd_freq, f2)]
    # dx: frequency resolution
    dx = psd_freq[2] - psd_freq[1]
    sbp = simpson(psd.power[frq_idx[1]:frq_idx[2]], psd_freq[frq_idx[1]:frq_idx[2]], dx=dx)

    return sbp
end

"""
    signal_band_power(signal, fs, f1, f2)

Calculates absolute band power between frequencies `f1` and `f2` for each the `signal` channels.

# Arguments

- `signal::Array{Float64, 3}`
- `fs::Int64` - sampling rate
- `f1::Float64` - lower frequency bound
- `f2::Float64` - upper frequency bound

# Returns

- `sbp::Matrix{Float64}`
"""
function signal_band_power(signal::Array{Float64, 3}; fs::Int64, f1::Union{Int64, Float64}, f2::Union{Int64, Float64})
    fs < 1 && throw(ArgumentError("Sampling rate must be ≥ 1 Hz."))
    f1 > f2 && throw(ArgumentError("Lower frequency bound must be lower than higher frequency bound."))
    f1 < 0 && throw(ArgumentError("Lower frequency bound must be be ≥ 1 Hz."))
    f2 > fs / 2 && throw(ArgumentError("Upper frequency bound must be be < Nyquist frequency ($(fs / 2) Hz)."))

    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)
    sbp = zeros(channels_no, epochs_no)

    Threads.@threads for epoch in 1:epochs_no
        for idx in 1:channels_no
            sbp[idx, epoch] = signal_band_power(signal[idx, :, epoch], fs=fs, f1=f1, f2=f2)
        end
    end

    return sbp
end

"""
    signal_make_spectrum(signal, fs)

Returns FFT and DFT sample frequencies for a DFT for the `signal`.

# Arguments

- `signal::Vector{Float64}`
- `fs::Int64` - sampling rate

# Returns

- `signal_fft::Vector{ComplexF64}`
- `signal_sf::Vector{Float64}`
"""
function signal_make_spectrum(signal::Vector{Float64}, fs::Int64)
    fs < 1 && throw(ArgumentError("Sampling rate must be ≥ 1 Hz."))

    signal_fft = fft(signal)
    # number of samples
    n = length(signal)
    # time between samples
    d = 1 / fs
    signal_sf = fftfreq(n, d)

    return signal_fft, Vector(signal_sf)
end

"""
    signal_make_spectrum(signal, fs)

Returns FFT and DFT sample frequencies for a DFT for each the `signal` channels.

# Arguments

- `signal::Array{Float64, 3}`
- `fs::Int64` - sampling rate

# Returns

- `signal_fft::Array{ComplexF64, 3}`
- `signal_sf::Array{Float64, 3}`
"""
function signal_make_spectrum(signal::Array{Float64, 3}, fs::Int64)
    fs < 1 && throw(ArgumentError("Sampling rate must be ≥ 1 Hz."))

    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)
    signal_fft = zeros(ComplexF64, size(signal))
    signal_sf = zeros(size(signal))

    Threads.@threads for epoch in 1:epochs_no
        for idx in 1:channels_no
            signal_fft[idx, :, epoch], signal_sf[idx, :, epoch] = signal_make_spectrum(signal[idx, :, epoch], fs)
        end
    end

    return signal_fft, signal_sf
end

"""
    signal_detrend(signal; type=:linear)

Removes linear trend from the `signal`.

# Arguments

- `signal::Vector{Float64}`
- `type::Symbol[:linear, :constant]`, optional
    - `linear` - the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `constant` - the mean of `signal` is subtracted

# Returns

- `signal_detrended::Vector{Float64}`
"""
function signal_detrend(signal::Vector{Float64}; type::Symbol=:linear)
    type in [:linear, :constant] || throw(ArgumentError("Trend type must be :linear or :constant."))

    if type === :constant
        signal_det = signal_demean(signal)
    else
        A = ones(length(signal))
        coef = A \ signal
        signal_det = @. signal - dot(A, coef)
    end

    return signal_det
end

"""
    signal_detrend(signal; type=:linear)

Removes linear trend for each the `signal` channels.

# Arguments

- `signal::Array{Float64, 3}`
- `type::Symbol[:linear, :constant]`, optional
    - `linear` - the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `constant` - the mean of `signal` is subtracted

# Returns

- `signal_detrended::Array{Float64, 3}`
"""
function signal_detrend(signal::Array{Float64, 3}; type::Symbol=:linear)
    type in [:linear, :constant] || throw(ArgumentError("Trend type must be :linear or :constant."))

    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)
    signal_det = zeros(size(signal))
    
    Threads.@threads for epoch in 1:epochs_no
        for idx in 1:channels_no
            signal_det[idx, :, epoch] = signal_detrend(signal[idx, :, epoch], type=type)
        end
    end

    return signal_det
end

"""
    signal_ci95(signal; n=3, method=:normal)

Calculates mean, std and 95% confidence interval for `signal`.

# Arguments

- `signal::Vector{Float64}`
- `n::Int` - number of bootstraps
- `method::Symbol[:normal, :boot]` - use normal method or `n`-times boostrapping

# Returns

- `signal_m::Float64`
- `signal_s::Float64`
- `signal_u::Float64`
- `signal_l::Float64`
"""
function signal_ci95(signal::Vector{Float64}; n::Int=3, method::Symbol=:normal)
    method === :normal || throw(ArgumentError("For vector signal method must be :normal."))
    n < 1 && throw(ArgumentError("n must be ≥ 1."))

    signal_m = mean(signal)
    signal_s = std(signal) / sqrt(length(signal))
    signal_u = signal_m + 1.96 * signal_s
    signal_l = signal_m - 1.96 * signal_s

    return signal_m, signal_s, signal_u, signal_l
end

"""
    signal_ci95(signal::Matrix{Float64}; n::Int=3, method::Symbol=:normal)

Calculates mean, std and 95% confidence interval for each the `signal` channels.

# Arguments

- `signal::Matrix{Float64}`
- `n::Int` - number of bootstraps
- `method::Symbol[:normal, :boot]` - use normal method or `n`-times boostrapping

# Returns

- `signal_m::Vector{Float64}`
- `signal_s::Vector{Float64}`
- `signal_u::Vector{Float64}`
- `signal_l::Vector{Float64}`
"""
function signal_ci95(signal::Matrix{Float64}; n::Int=3, method::Symbol=:normal)
    method in [:normal, :boot] || throw(ArgumentError("Method must be :normal or :boot."))

    if method === :normal
        signal_m = mean(signal, dims=1)'
        signal_s = std(signal, dims=1)' / sqrt(size(signal, 1))
        signal_u = signal_m + 1.96 * signal_s
        signal_l = signal_m - 1.96 * signal_s
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
        signal_m = mean(signal_tmp1, dims=1)'
        signal_s = std(signal_tmp1, dims=1)' / sqrt(size(signal_tmp1, 1))
        signal_sorted = sort(signal_tmp1, dims=1)
        signal_l = signal_sorted[round(Int, 0.025 * size(signal_tmp1, 1)), :]
        signal_u = signal_sorted[round(Int, 0.975 * size(signal_tmp1, 1)), :]
    end

    return vec(signal_m[:, 1]), vec(signal_s[:, 1]), vec(signal_u[:, 1]), vec(signal_l[:, 1])
end

"""
    signal_ci95(signal; n=3, method=:normal)

Calculates mean, std and 95% confidence interval for each the `signal` channels.

# Arguments

- `signal::Array{Float64, 3}`
- `n::Int` - number of bootstraps
- `method::Symbol[:normal, :boot]` - use normal method or `n`-times boostrapping

# Returns

- `signal_m::Matrix{Float64}`
- `signal_s::Matrix{Float64}`
- `signal_u::Matrix{Float64}`
- `signal_l::Matrix{Float64}`
"""
function signal_ci95(signal::Array{Float64, 3}; n::Int=3, method::Symbol=:normal)
    method in [:normal, :boot] || throw(ArgumentError("Method must be :normal or :boot."))
    n < 1 && throw(ArgumentError("n must be ≥ 1."))

    signal_m = zeros(size(signal, 3), size(signal, 2))
    signal_s = zeros(size(signal, 3), size(signal, 2))
    signal_u = zeros(size(signal, 3), size(signal, 2))
    signal_l = zeros(size(signal, 3), size(signal, 2))

    epochs_no = size(signal, 3)

    Threads.@threads for epoch in 1:epochs_no
        signal_m[epoch, :], signal_s[epoch, :], signal_u[epoch, :], signal_l[epoch, :] = signal_ci95(signal[:, :, epoch])
    end

    return signal_m, signal_s, signal_u, signal_l
end

"""
    signal_mean(signal1, signal2)

Calculates mean and 95% confidence interval for 2 signals.

# Arguments

- `signal1::Vector{Float64}`
- `signal2:Vector{Float64}`

# Returns

- `signal_mean::Float64`
- `signal_sd::Float64`
- `signal_u::Float64`
- `signal_l::Float64`
"""
function signal_mean(signal1::Vector{Float64}, signal2::Vector{Float64})
    length(signal1) != length(signal2) && throw(ArgumentError("Both signals must be of the same as size."))

    signals_m = zeros(length(signal1))
    signals_s = zeros(length(signal1))
    signals_u = zeros(length(signal1))
    signals_l = zeros(length(signal1))

    signal1_mean = mean(signal1)
    signal2_mean = mean(signal2)
    signals_m = signal1_mean - signal2_mean
    signal1_sd = std(signal1) / sqrt(length(signal1))
    signal2_sd = std(signal2) / sqrt(length(signal2))
    signals_s = sqrt(signal1_sd^2 + signal2_sd^2)
    signals_u = signals_m + 1.96 * signals_s
    signals_l = signals_m - 1.96 * signals_s

    return signals_m, signals_s, signals_u, signals_l
end

"""
    signal_mean(signal1, signal2)

Calculates mean and 95% confidence interval for 2 signals.

# Arguments

- `signal1::Array{Float64, 3}`
- `signal2:Array{Float64, 3}`

# Returns

- `signal_mean::Matrix{Float64}`
- `signal_sd::Matrix{Float64}`
- `signal_u::Matrix{Float64}`
- `signal_l::Matrix{Float64}`
"""
function signal_mean(signal1::Array{Float64, 3}, signal2::Array{Float64, 3})
    size(signal1) != size(signal2) && throw(ArgumentError("Both signals must be of the same as size."))

    signals_m = zeros(size(signal1, 3), size(signal1, 2))
    signals_s = zeros(size(signal1, 3), size(signal1, 2))
    signals_u = zeros(size(signal1, 3), size(signal1, 2))
    signals_l = zeros(size(signal1, 3), size(signal1, 2))
    epochs_no = size(signal1, 3)

    Threads.@threads for epoch in 1:epochs_no
        signal1_mean = mean(signal1[:, :, epoch], dims=1)
        signal2_mean = mean(signal2[:, :, epoch], dims=1)
        signals_m[epoch, :] = signal1_mean - signal2_mean
        signal1_sd = std(signal1[:, :, epoch], dims=1) / sqrt(size(signal1[:, :, epoch], 2))
        signal2_sd = std(signal2[:, :, epoch], dims=1) / sqrt(size(signal2[:, :, epoch], 2))
        signals_s[epoch, :] = sqrt.(signal1_sd.^2 .+ signal2_sd.^2)
        signals_u[epoch, :] = @. signals_m[epoch, :] + 1.96 * signals_s[epoch, :]
        signals_l[epoch, :] = @. signals_m[epoch, :] - 1.96 * signals_s[epoch, :]
    end

    return signals_m, signals_s, signals_u, signals_l
end

"""
    signal_difference(signal1, signal2; n=3, method=:absdiff)

Calculates mean difference and 95% confidence interval for 2 signals.

# Arguments

- `signal1::Matrix{Float64}`
- `signal2:Matrix{Float64}`
- `n::Int` - number of bootstraps.
- `method::Symbol[:absdiff, :diff2int]`
    - `:absdiff` - maximum difference
    - `:diff2int` - integrated area of the squared difference

# Returns

- `signals_statistic::Vector{Float64}`
- `signals_statistic_single::Float64`
- `p::Float64`
"""
function signal_difference(signal1::Matrix{Float64}, signal2::Matrix{Float64}; n=3, method::Symbol=:absdiff)
    size(signal1) != size(signal2) && throw(ArgumentError("Both signals must be of the same size."))
    method in [:absdiff, :diff2int] || throw(ArgumentError("Method must be :absdiff or :diff2int."))

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
    p > 1 && (p = 1.0)

    return signals_statistic, signals_statistic_single, p
end

"""
    signal_difference(signal1, signal2; n=3, method=:absdiff)

Calculates mean difference and 95% confidence interval for 2 signals.

# Arguments

- `signal1::Array{Float64, 3}`
- `signal2:Array{Float64, 3}`
- `n::Int` - number of bootstraps.
- `method::Symbol[:absdiff, :diff2int]`
    - `:absdiff` - maximum difference
    - `:diff2int` - integrated area of the squared difference

# Returns

- `signals_statistic::Matrix{Float64}`
- `signals_statistic_single::Vector{Float64}`
- `p::Vector{Float64}`
"""
function signal_difference(signal1::Array{Float64, 3}, signal2::Array{Float64, 3}; n=3, method::Symbol=:absdiff)
    size(signal1) != size(signal2) && throw(ArgumentError("Both signals must be of the same size."))
    method in [:absdiff, :diff2int] || throw(ArgumentError("Method must be :absdiff or :diff2int."))

    epochs_no = size(signal1, 3)
    signals_statistic = zeros(epochs_no, size(signal1, 1) * n)
    signals_statistic_single = zeros(epochs_no)
    p = zeros(epochs_no)

    Threads.@threads for epoch in 1:epochs_no
        signals_statistic[epoch, :], signals_statistic_single[epoch], p[epoch] = signal_difference(signal1[:, :, epoch], signal2[:, :, epoch])
    end

    return signals_statistic, signals_statistic_single, p
end

"""
   signal_autocov(signal; lag=1, demean=false, normalize=false)

Calculates autocovariance of the `signal`.

# Arguments

- `signal::Vector{Float64}`
- `lag::Int64` - lags range is `-lag:lag`
- `demean::Bool` - demean `signal` prior to calculations
- `normalize::Bool` - normalize autocovariance

# Returns

- `acov::Vector{Float64}`
- `lags::Vector{Int64}`
"""
function signal_autocov(signal::Vector{Float64}; lag::Int64=1, demean::Bool=false, normalize::Bool=false)
    lag < 1 && throw(ArgumentError("Lag must be ≥ 1."))

    lags = collect(-lag:lag)

    if demean == true
        signal_demeaned = signal .- mean(signal)
    else
        signal_demeaned = signal
    end

    acov = zeros(length(lags))

    for idx in 1:length(lags)
        if lags[idx] == 0
            # no lag
            signal_lagged = signal_demeaned
            signals_mul = signal_demeaned .* signal_lagged
        elseif lags[idx] > 0
            # positive lag
            signal_lagged = signal_demeaned[1:(end - lags[idx])]
            signals_mul = signal_demeaned[(1 + lags[idx]):end] .* signal_lagged
        elseif lags[idx] < 0
            # negative lag
            signal_lagged = signal_demeaned[(1 + abs(lags[idx])):end]
            signals_mul = signal_demeaned[1:(end - abs(lags[idx]))] .* signal_lagged
        end
        signals_sum = sum(signals_mul)
        if normalize == true
            acov[idx] = signals_sum / length(signal)
        else
            acov[idx] = signals_sum
        end
    end

    return acov, lags
end

"""
   signal_autocov(signal; lag=1, demean=false, normalize=false)

Calculates autocovariance of each the `signal` channels.

# Arguments

- `signal::Array{Float64, 3}`
- `lag::Int64` - lags range is `-lag:lag`
- `demean::Bool` - demean signal prior to analysis
- `normalize::Bool` - normalize autocovariance

# Returns

- `acov::Matrix{Float64}`
- `lags::Vector{Int64}`
"""
function signal_autocov(signal::Array{Float64, 3}; lag::Int64=1, demean::Bool=false, normalize::Bool=false)
    lag < 1 && throw(ArgumentError("Lag must be ≥ 1."))

    lags = collect(-lag:lag)
    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)
    acov = zeros(channels_no, length(lags), epochs_no)

    Threads.@threads for epoch in 1:epochs_no
        for idx in 1:channels_no
            acov[idx, :, epoch], lags = signal_autocov(signal[idx, :, epoch],
                                                 lag=lag,
                                                 demean=demean,
                                                 normalize=normalize)
        end
    end

    return acov, lags
end

"""
   signal_crosscov(signal1, signal2; lag=1, demean=false, normalize=false)

Calculates cross-covariance between `signal1` and `signal2`.

# Arguments

- `signal1::Vector{Float64}`
- `signal2::Vector{Float64}`
- `lag::Int64` - lags range is `-lag:lag`
- `demean::Bool` - demean signal prior to analysis
- `normalize::Bool` - normalize cross-covariance

# Returns

- `ccov::Vector{Float64}`
- `lags::Vector{Int64}`
"""
function signal_crosscov(signal1::Vector{Float64}, signal2::Vector{Float64}; lag::Int64=1, demean::Bool=false, normalize::Bool=false)
    length(signal1) != length(signal2) && throw(ArgumentError("Both vectors must be of the same as length."))
    lag < 1 && throw(ArgumentError("Lag must be ≥ 1."))

    lags = collect(-lag:lag)

    if demean == true
        signal_demeaned1 = signal1 .- mean(signal1)
        signal_demeaned2 = signal2 .- mean(signal2)
    else
        signal_demeaned1 = signal1
        signal_demeaned2 = signal2
    end

    ccov = zeros(length(lags))

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
            ccov[idx] = signals_sum / length(signal1)
        else
            ccov[idx] = signals_sum
        end
    end

    return ccov, lags
end

"""
   signal_crosscov(signal; lag=1, demean=false, normalize=false)

Calculates cross-covariance between all channels in the `signal`.

# Arguments

- `signal::Matrix{Float64}` - the signal
- `lag::Int64` - lags range is `-lag:lag`
- `demean::Bool` - demean `signal` prior to analysis
- `normalize::Bool` - normalize cross-covariance

# Returns

- `ccov::Array{Float64, 3}`
- `lags::Vector{Int64}`
"""
function signal_crosscov(signal::Array{Float64, 3}; lag::Int64=1, demean::Bool=false, normalize::Bool=false)
    lag < 1 && throw(ArgumentError("Lag must be ≥ 1 Hz."))

    lags = collect(-lag:lag)
    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)
    ccov = zeros(channels_no^2, length(lags), epochs_no)

    Threads.@threads for epoch in 1:epochs_no
        ccov_packed = Array{Vector{Float64}}(undef, channels_no, channels_no)
        for idx1 in 1:channels_no
            for idx2 in 1:channels_no
                ccov_packed[idx1, idx2], lags = signal_crosscov(signal[idx1, :, epoch],
                                                          signal[idx2, :, epoch],
                                                          lag=lag,
                                                          demean=demean,
                                                          normalize=normalize)
            end
        end
        for idx in 1:channels_no^2
            ccov[idx, :, epoch] = ccov_packed[idx]
        end
    end

    return ccov, lags
end

"""
   signal_crosscov(signal1, signal2; lag=1, demean=false, normalize=false)

Calculates cross-covariance between same channels in `signal1` and `signal2`.

# Arguments

- `signal1::Array{Float64, 3}`
- `signal2::Array{Float64, 3}`
- `lag::Int64` - lags range is `-lag:lag`
- `demean::Bool` - demean signal prior to analysis
- `normalize::Bool` - normalize cross-covariance

# Returns

- `ccov::Array{Float64, 3}`
- `lags::Vector{Int64}`
"""
function signal_crosscov(signal1::Array{Float64, 3}, signal2::Array{Float64, 3}; lag::Int64=1, demean::Bool=false, normalize::Bool=false)
    size(signal1) != size(signal2) && throw(ArgumentError("Both arrays must be of the same as size."))
    lag < 1 && throw(ArgumentError("Lag must be ≥ 1 Hz."))

    lags = collect(-lag:lag)
    channels_no = size(signal1, 1)
    epochs_no = size(signal1, 3)
    ccov = zeros(channels_no, length(lags), epochs_no)

    Threads.@threads for epoch in 1:epochs_no
        for idx in 1:channels_no
            ccov[idx, :, epoch], lags = signal_crosscov(signal1[idx, :, epoch],
                                                  signal2[idx, :, epoch],
                                                  lag=lag,
                                                  demean=demean,
                                                  normalize=normalize)
        end
    end

    return ccov, lags
end

"""
    signal_spectrum(signal; pad=0)

Calculates FFT, amplitudes, powers and phases of the `signal`.

# Arguments

- `signal::Vector{Float64}`
- `pad::Int64` - pad the `signal` with `pad` zeros

# Returns

- `signal_fft::Vector(ComplexF64}`
- `signal_amplitudes::Vector{Float64}`
- `signal_powers::Vector{Float64}`
- `signal_phases::Vector{Float64}
"""
function signal_spectrum(signal::Vector{Float64}; pad::Int64=0)
    pad < 0 && throw(ArgumentError("Pad cannot be negative."))

    signal_fft = fft0(signal, pad)

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

- `signal::Array{Float64, 3}` - the signal
- `pad::Int64` - pad the `signal` with `pad` zeros

# Returns

- `signal_fft::Array{ComplexF64, 3}`
- `signal_amplitudes::Array{Float64, 3}`
- `signal_powers::Array{Float64, 3}`
- `signal_phases::Array{Float64, 3}
"""
function signal_spectrum(signal::Array{Float64, 3}; pad::Int64=0)
    pad < 0 && throw(ArgumentError("Pad cannot be negative."))

    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)

    signal_fft = zeros(ComplexF64, channels_no, size(signal, 2) + pad, epochs_no)
    signal_amplitudes = zeros(channels_no, size(signal, 2) + pad, epochs_no)
    signal_powers = zeros(channels_no, size(signal, 2) + pad, epochs_no)
    signal_phases = zeros(channels_no, size(signal, 2) + pad, epochs_no)

    Threads.@threads for epoch in 1:epochs_no
        for idx in 1:channels_no
            signal_fft[idx, :, epoch], signal_amplitudes[idx, :, epoch], signal_powers[idx, :, epoch], signal_phases[idx, :, epoch] = signal_spectrum(signal[idx, :, epoch], pad=pad)
        end
    end

    return signal_fft, signal_amplitudes, signal_powers, signal_phases
end

"""
    signal_epochs(signal; epochs_no, epochs_len, average=true)

Splits `signal` into epochs.

# Arguments

- `signal::Vector{Float64}`
- `epochs_no::Union{Int64, Nothing}` - number of epochs
- `epochs_len::Union{Int64, Nothing}` - epoch length in samples
- `average::Bool` - average all epochs, returns one averaged epoch; if false than returns array of epochs, each row is one epoch

# Returns

- `epochs::Matrix{Float64}`
"""
function signal_epochs(signal::Vector{Float64}; epochs_no::Union{Int64, Nothing}=nothing, epochs_len::Union{Int64, Nothing}=nothing, average::Bool=false)
    (epochs_len === nothing && epochs_no === nothing) && throw(ArgumentError("Either number of epochs or epoch length must be set."))
    (epochs_len !== nothing && epochs_no !== nothing) && throw(ArgumentError("Both number of epochs and epoch length cannot be set."))
    (epochs_len != nothing && epochs_len < 1) && throw(ArgumentError("Epoch length rate must be ≥ 1."))
    (epochs_no != nothing && epochs_no < 1) && throw(ArgumentError("Epoch length rate must be ≥ 1."))

    if epochs_no === nothing
        epochs_no = length(signal) ÷ epochs_len
    end
    if epochs_len === nothing
        epochs_len = length(signal) ÷ epochs_no
    end

    epochs = zeros(epochs_no, epochs_len)

    signal = signal[1:epochs_len * epochs_no]
    epochs = reshape(signal, epochs_no, epochs_len)'

    if average == true
        epochs = vec(mean(epochs, dims=1)[1, :])
    end

    return epochs
end

"""
    signal_epochs(signal; epochs_no=nothing, epochs_len=nothing, average=true)

Splits `signal` into epochs.

# Arguments

- `signal::Array{Float64, 3}`
- `epochs_no::Union{Int64, Nothing}` - number of epochs
- `epochs_len::Union{Int64, Nothing}` - epoch length in samples
- `average::Bool` - average all epochs, returns one averaged epoch; if false than returns array of epochs, each row is one epoch

# Returns

- `epochs::Array{Float64, 3}`
"""
function signal_epochs(signal::Matrix{Float64}; epochs_no::Union{Int64, Nothing}=nothing, epochs_len::Union{Int64, Nothing}=nothing, average::Bool=false)
    (epochs_len === nothing && epochs_no === nothing) && throw(ArgumentError("Either number of epochs or epoch length must be set."))
    (epochs_len !== nothing && epochs_no !== nothing) && throw(ArgumentError("Both number of epochs and epoch length cannot be set."))
    (epochs_len != nothing && epochs_len < 1) && throw(ArgumentError("Epoch length rate must be ≥ 1."))
    (epochs_no != nothing && epochs_no < 1) && throw(ArgumentError("Epoch length rate must be ≥ 1."))

    channels_no = size(signal, 1)

    if epochs_no === nothing
        epochs_no = size(signal, 2) ÷ epochs_len
    else
        epochs_len = size(signal, 2) ÷ epochs_no
    end

    epochs = zeros(channels_no, epochs_len, epochs_no)

    # signal = signal[1:epochs_len * epochs_no]
    # epochs = reshape(signal, epochs_no, epochs_len)'

    idx1 = 1
    for idx2 in 1:epochs_len:(epochs_no * epochs_len - 1)
        epochs[:, :, idx1] = signal[:, idx2:(idx2 + epochs_len - 1), 1]
        idx1 += 1
    end

    if average == true
        epochs = mean(epochs, dims=3)[:, :]
    end

    return epochs
end

"""
    signal_drop_channel(signal, channels)

Removes `channels` from the `signal`.

# Arguments

- `signal::Matrix{Float64}`
- `channels::Union{Int64, Vector{Int64}, AbstractRange}` - channels to be removed, vector of numbers or range

# Returns

- `signal_modified::Matrix{Float64}`
"""
function signal_drop_channel(signal::Matrix{Float64}, channels::Union{Int64, Vector{Int64}, AbstractRange})
    if typeof(channels) <: AbstractRange
        channels = collect(channels)
    end
    for idx in length(channels):-1:1
        channels[idx] > size(signal, 2) && throw(ArgumentError("Channels index does not match signal."))
    end

    length(channels) > 1 && (channels = sort!(channels, rev=true))
    signal_new = signal[setdiff(1:end, (channels)), :]

    return signal_new
end

"""
    signal_drop_channel(signal, channels)

Removes `channels` from the `signal`.

# Arguments

- `signal::Array{Float64, 3}`
- `channels::Union{Int64, Vector{Int64}, AbstractRange}` - channels to be removed, vector of numbers or range

# Returns

- `signal_modified::Matrix{Float64}`
"""
function signal_drop_channel(signal::Array{Float64, 3}, channels::Union{Int64, Vector{Int64}, AbstractRange})
    if typeof(channels) <: AbstractRange
        channels = collect(channels)
    end
    for idx in length(channels):-1:1
        channels[idx] > size(signal, 2) && throw(ArgumentError("Channels index does not match signal."))
    end

    length(channels) > 1 && (channels = sort!(channels, rev=true))
    signal_new = signal[setdiff(1:end, (channels)), :, :]

    return signal_new
end
"""
    signal_reference_channel(signal, reference)

Re-references channels of the `signal` to specific signal channel.

# Arguments

- `signal::Matrix{Float64}`
- `reference::Union{Int64, Vector{Int64}, AbstractRange}}` - index of channels used as reference; if multiple channels are specified, their average is used as the reference

# Returns

- `signal_referenced::Matrix{Float64}`
"""
function signal_reference_channel(signal::Matrix{Float64}, reference_idx::Union{Int64, Vector{Int64}, AbstractRange})
    if typeof(reference_idx) <: AbstractRange
        reference_idx = collect(reference_idx)
    end

    channels_no = size(signal, 1)
    signal_ref = zeros(size(signal))

    channels_list = collect(1:channels_no)
    for idx in 1:length(reference_idx)
        if (reference_idx[idx] in channels_list) == false
            throw(ArgumentError("Reference channel index does not match signal."))
        end
    end

    if length(reference_idx) == 1
        reference_channel = mean(signal[reference_idx, :], dims=2)
    else
        reference_channel = vec(mean(signal[reference_idx, :], dims=1))
    end
    for idx in 1:channels_no
        signal_ref[idx, :] = signal[idx, :] .- reference_channel
    end
    length(reference_idx) == 1 && (signal_ref[reference_idx, :] = reference_channel)

    return signal_ref
end

"""
    signal_reference_channel(signal, reference)

Re-references channels of the `signal` to specific signal channel.

# Arguments

- `signal::Array{Float64, 3}`
- `reference::Union{Int64, Vector{Int64}, AbstractRange}}` - index of channels used as reference; if multiple channels are specified, their average is used as the reference

# Returns

- `signal_referenced::Matrix{Float64}`
"""
function signal_reference_channel(signal::Array{Float64, 3}, reference_idx::Union{Int64, Vector{Int64}, AbstractRange})
    if typeof(reference_idx) <: AbstractRange
        reference_idx = collect(reference_idx)
    end

    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)
    signal_ref = zeros(size(signal))

    channels_list = collect(1:channels_no)
    for idx in 1:length(reference_idx)
        if (reference_idx[idx] in channels_list) == false
            throw(ArgumentError("Reference channel index does not match signal."))
        end
    end
    Threads.@threads for epoch in 1:epochs_no
        if length(reference_idx) == 1
            reference_channel = mean(signal[reference_idx, :, epoch], dims=2)
        else
            reference_channel = vec(mean(signal[reference_idx, :, epoch], dims=1))
        end
        for idx in 1:channels_no
            signal_ref[idx, :, epoch] = signal[idx, :, epoch] .- reference_channel
        end
        length(reference_idx) == 1 && (signal_ref[reference_idx, :, epoch] = reference_channel)
    end

    return signal_ref
end

"""
    signal_reference_car(signal)

Re-references channels of the `signal` to common average reference.

# Arguments

- `signal::Matrix{Float64}`

# Returns

- `signal_referenced::Matrix{Float64}`
"""
function signal_reference_car(signal::Matrix{Float64})
    channels_no = size(signal, 1)
    reference_channel = vec(mean(signal, dims=1))
    signal_ref = zeros(size(signal))

    for idx in 1:channels_no
        signal_ref[idx, :] = signal[idx, :] .- reference_channel
    end

    return signal_ref
end

"""
    signal_reference_car(signal)

Re-references channels of the `signal` to common average reference.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `signal_referenced::Array{Float64, 3}`
"""
function signal_reference_car(signal::Array{Float64, 3})
    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)
    signal_ref = zeros(size(signal))

    Threads.@threads for epoch in 1:epochs_no
        reference_channel = vec(mean(signal[:, :, epoch], dims=1))
        for idx in 1:channels_no
            signal_ref[idx, :, epoch] = signal[idx, :, epoch] .- reference_channel
        end
    end

    return signal_ref
end

"""
    signal_taper(signal, taper)

Taper the `signal` with `taper`.

# Arguments

- `signal::Vector{Float64}`
- `taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `signal_tapered::Vector{Union{Float64, ComplexF64}}`
"""
function signal_taper(signal::Vector{Float64}, taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})
    length(taper) == length(signal) || throw(ArgumentError("Taper length and signal length must be equal."))
    signal_tap = signal .* taper

    return signal_tap
end

"""
    signal_taper(signal, taper)

Tapers channels of the `signal` with `taper`.

# Arguments

- `signal::Array{Float64, 3}`
- `taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `signal_tapered::Union{Array{Float64, 3}, Array{ComplexF64, 3}}`
"""
function signal_taper(signal::Array{Float64, 3}, taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})
    length(taper) == size(signal, 2) || throw(ArgumentError("Taper length and signal length must be equal."))

    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)

    signal_tap = zeros(eltype(taper), size(signal))

    Threads.@threads for epoch in 1:epochs_no
        for idx in 1:channels_no
            signal_tap[idx, :, epoch] = signal_taper(signal[idx, :, epoch], taper)
        end
    end

    return signal_tap
end

"""
    signal_demean(signal)

Removes mean value (DC offset) from the `signal`.

# Arguments

- `signal::Vector{Float64}`

# Returns

- `signal_demeaned::Vector{Float64}`
"""
function signal_demean(signal::Vector{Float64})
    signal_dem = signal .- mean(signal)
    return signal_dem
end

"""
    signal_demean(signal)

Removes mean value (DC offset) for each the `signal` channels.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `signal_demeaned::Array{Float64, 3}`
"""
function signal_demean(signal::Array{Float64, 3})
    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)
    signal_dem = zeros(size(signal))

    Threads.@threads for epoch in 1:epochs_no
        for idx in 1:channels_no
            signal_dem[idx, :, epoch] = signal_demean(signal[idx, :, epoch])
        end
    end

    return signal_dem
end

"""
    signal_normalize_zscore(signal)

Normalize (by z-score) `signal`.

# Arguments

- `signal::Vector{Float64}`

# Returns

- `signal_normalized::Vector{Float64}`
"""
function signal_normalize_zscore(signal::Vector{Float64})
    signal_norm = (signal .- mean(signal)) ./ std(signal)
    return signal_norm
end

"""
    signal_normalize_zscore(signal)

Normalize (by z-score) each the `signal` channel.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `signal_normalized::Array{Float64, 3}`
"""
function signal_normalize_zscore(signal::Array{Float64, 3})
    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)
    signal_norm = zeros(size(signal))

    Threads.@threads for epoch in 1:epochs_no
        for idx in 1:channels_no
            signal_norm[idx, :, epoch] = signal_normalize_zscore(signal[idx, :, epoch])
        end
    end

    return signal_norm
end

"""
    signal_normalize_minmax(signal)

Normalize (to 0…1) `signal`.

# Arguments

- `signal::Vector{Float64}`

# Returns

- `signal_normalized::Vector{Float64}`
"""
function signal_normalize_minmax(signal::Vector{Float64})
    signal_norm = (signal .- minimum(signal)) ./ (maximum(signal) - minimum(signal))
    return signal_norm
end

"""
    signal_normalize_minmax(signal)

Normalize (to 0…1) each the `signal` channel.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `signal_normalized::Array{Float64, 3}`
"""
function signal_normalize_minmax(signal::Array{Float64, 3})
    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)
    signal_norm = zeros(size(signal))

    Threads.@threads for epoch in 1:epochs_no
        for idx in 1:channels_no
            signal_norm[idx, :, epoch] = signal_normalize_minmax(signal[idx, :, epoch])
        end
    end

    return signal_norm
end

"""
   signal_cov(signal1, signal2; normalize=true)

Calculates covariance between `signal1` and `signal2`.

# Arguments

- `signal1::Vector{Float64}`
- `signal2::Vector{Float64}`
- `normalize::Bool` - normalize covariance

# Returns

- `cov_mat::Matrix{Float64}`
"""
function signal_cov(signal1::Vector{Float64}, signal2::Vector{Float64}; normalize::Bool=true)
    length(signal1) != length(signal2) && throw(ArgumentError("Both vectors must be of the same as length."))

    cov_mat = cov(signal1 * signal2')

    # divide so that components are centered at (0, 0)
    normalize == true && (cov_mat = cov_mat ./ length(signal1))

    return cov_mat
end

"""
   signal_cov(signal; normalize=true)

Calculates covariance between all channels of the `signal`.

# Arguments

- `signal::Array{Float64, 3}`
- `normalize::Bool` - normalize covariance

# Returns

- `cov_mat::Array{Float64, 3}`
"""
function signal_cov(signal::Array{Float64, 3}; normalize::Bool=true)
    epochs_no = size(signal, 3)
    cov_mat = zeros(size(signal, 1), size(signal, 1), epochs_no)

    Threads.@threads for epoch in 1:epochs_no
        # channels-vs-channels
        cov_mat[:, :, epoch] = cov(signal[:, :, epoch]')
    end

    # divide so that components are centered at (0, 0)
    normalize == true && (cov_mat = cov_mat ./ size(cov_mat, 2))

    return cov_mat
end

"""
   signal_cor(signal)

Calculates correlation coefficients between all channels of the `signal`.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `cor_mat::Array{Float64, 3}`
"""
function signal_cor(signal::Array{Float64, 3})
    epochs_no = size(signal, 3)
    cor_mat = zeros(size(signal, 1), size(signal, 1), epochs_no)

    Threads.@threads for epoch in 1:epochs_no
        # channels-vs-channels
        cor_mat[:, :, epoch] = cor(signal[:, :, epoch]')
    end

    return cor_mat
end

"""
    signal_add_noise(signal)

Adds random noise to the `signal`.

# Arguments

- `signal::Vector{Float64}`

# Returns

- `signal_noisy::Vector{Float64}`
"""
function signal_add_noise(signal::Vector{Float64})
    signal_noise = signal .+ mean(signal) .* rand(length(signal))

    return signal_noise
end

"""
    signal_add_noise(signal)

Add random noise to the `signal` channels.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `signal_noisy::Array{Float64, 3}`
"""
function signal_add_noise(signal::Array{Float64, 3})
    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)
    signal_noise = zeros(size(signal))

    Threads.@threads for epoch in 1:epochs_no
        for idx in 1:channels_no
            signal_noise[idx, :, epoch] = signal_add_noise(signal[idx, :, epoch])
        end
    end

    return signal_noise
end

"""
    signal_upsample(signal; t, new_sr)

Upsamples `signal` to `new_sr` sampling frequency.

# Arguments

- `signal::Vector{Float64}`
- `t::AbstractRange`
- `new_sr::Int64` - new sampling rate
# Returns

- `signal_upsampled::Vector{Float64}`
- `t_upsampled::AbstractRange`
"""
function signal_upsample(signal::Vector{Float64}; t::AbstractRange, new_sr::Int64)
    new_sr < 1 && throw(ArgumentError("New sampling rate must be positive."))

    # sampling interval
    dt = t[2] - t[1]
    # sampling rate
    sr = 1 / dt
    new_sr < sr && throw(ArgumentError("New sampling rate must be larger than signal sampling rate."))
    new_sr == sr && return(signal)

    # interpolate
    signal_interpolation = CubicSplineInterpolation(t, signal)
    t_upsampled = t[1]:1/new_sr:t[end]
    signal_upsampled = signal_interpolation(t_upsampled)

    return signal_upsampled, t_upsampled
end

"""
    signal_upsample(signal; t, new_sr)

Upsamples all channels of `signal` to `new_sr` sampling frequency.

# Arguments

- `signal::Array{Float64, 3}`
- `t::AbstractRange`
- `new_sr::Int64` - new sampling rate

# Returns

- `signal_upsampled::Array{Float64, 3}`
- `t_upsampled::AbstractRange`
"""
function signal_upsample(signal::Array{Float64, 3}; t::AbstractRange, new_sr::Int64)
    new_sr < 1 && throw(ArgumentError("New sampling rate must be positive."))

    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)

    signal_upsampled_length = length(signal_upsample(signal[1, :, 1], t=t, new_sr=new_sr)[1])
    signal_upsampled = zeros(channels_no, signal_upsampled_length, epochs_no) 

    t_upsampled = nothing
    Threads.@threads for epoch in 1:epochs_no
        for idx in 1:channels_no
            signal_upsampled[idx, :, epoch], t_upsampled = signal_upsample(signal[idx, :, epoch], t=t, new_sr=new_sr)
        end
    end

    return signal_upsampled, t_upsampled
end

"""
    signal_tconv(signal, kernel)

Performs convolution in the time domain between `signal` and `kernel`.

# Arguments

- `signal::Vector{Float64}`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `signal_convoluted::Union{Vector{Float64}, Vector{ComplexF64}}`
"""
function signal_tconv(signal::Vector{Float64}, kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})
    signal_tconv = conv(signal, kernel)
    half_kernel = floor(Int, length(kernel) / 2)
    signal_tconv = signal[(half_kernel + 1):(end - half_kernel)]

    return signal_tconv
end

"""
    signal_tconv(signal, kernel)

Performs convolution in the time domain between `signal` and `kernel`.

# Arguments

- `signal::Array{Float64, 3}`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `signal_convoluted::Union{Array{Float64, 3}, Array{ComplexF64, 3}}`
"""
function signal_tconv(signal::Array{Float64, 3}, kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})
    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)

    if typeof(kernel) == Vector{ComplexF64}
        signal_conv = zeros(ComplexF64, size(signal))
    else
        signal_conv = zeros(channels_no, length(signal_tconv(signal[1, :, 1], kernel)), epochs_no)
    end

    Threads.@threads for epoch in 1:epochs_no
        for idx in 1:channels_no
            signal_conv[idx, :, epoch] = signal_tconv(signal[idx, :, epoch], kernel)
        end
    end

    return signal_conv
end

"""
    signal_filter(signal; fprototype, ftype, cutoff, fs, order=8, window=hanning(64))

Filters `signal` using zero phase distortion filter.

# Arguments

- `signal::Vector{Float64}`
- `fprototype::Symbol[:butterworth, :fir]
- `ftype::Symbol[:lp, :hp, :bp, :bs]` - filter type
- `cutoff::Union{Int64, Float64, Vector{Int64}, Vector{Float64}}` - filter cutoff in Hz (vector for `:bp` and `:bs`)
- `fs::Int64` - sampling rate
- `order::Int64` - filter order
- `rp::Float64` - dB ripple in the passband
- `rs::Float64` - dB attentuation in the stopband
- `dir:Symbol[:oenpass, :onepass_reverse, :twopass]` - filter direction
- `window::Union{Vector{Float64}, Nothing} - window, required for FIR filter

# Returns

- `signal_filtered::Vector{Float64}`
"""
function signal_filter(signal::Vector{Float64}; fprototype::Symbol, ftype::Symbol, cutoff::Union{Int64, Float64, Vector{Int64}, Vector{Float64}}, fs::Int64, order::Int64, rp::Union{Int64, Float64, Nothing}=nothing, rs::Union{Int64, Float64, Nothing}=nothing, dir::Symbol=:twopass, window::Union{Vector{Float64}, Nothing}=nothing)
    (fprototype === :fir && (window === nothing || length(window) > length(signal))) && throw(ArgumentError("For FIR filter window must be longer that signal."))
    ftype in [:lp, :hp, :bp, :bs] || throw(ArgumentError("Filter type must be :bp, :hp, :bp or :bs."))
    dir in [:onepass, :onepass_reverse, :twopass] || throw(ArgumentError("Filter direction must be :onepass, :onepass_reverse or :twopass."))
    fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir] || throw(ArgumentError("Filter prototype must be :butterworth, :chebyshev1, :chebyshev2, :elliptic or :fir."))
    mod(order, 2) != 0 && throw(ArgumentError("Filter order must be even."))
    order < 0 && throw(ArgumentError("Order must be positive."))

    if ftype === :lp
        length(cutoff) > 1 && throw(ArgumentError("For low-pass filter one frequency must be given."))
        responsetype = Lowpass(cutoff; fs=fs)
    elseif ftype === :hp
        length(cutoff) > 1 && throw(ArgumentError("For high-pass filter one frequency must be given."))
        responsetype = Highpass(cutoff; fs=fs)
    elseif ftype === :bp
        responsetype = Bandpass(cutoff[1], cutoff[2]; fs=fs)
    elseif ftype === :bs
        (length(cutoff) < 2 || length(cutoff) > 2) && throw(ArgumentError("For band-stop filter two frequencies must be given."))
        responsetype = Bandstop(cutoff[1], cutoff[2]; fs=fs)
    end

    fprototype === :butterworth && (prototype = Butterworth(order))
    if fprototype === :fir
        window === nothing && throw(ArgumentError("For FIR filter window must be given."))
        if ftype === :hp || ftype === :bp || ftype === :bs
            mod(length(window), 2) == 0 && (window = vcat(window[1:((length(window) ÷ 2) - 1)], window[((length(window) ÷ 2) + 1):end]))
        end
        prototype = FIRWindow(window)
    end
    if fprototype === :chebyshev1
        rs === nothing && throw(ArgumentError("For Chebyshev1 filter rs must be given."))
        prototype = Chebyshev1(order, rs)
    end
    if fprototype === :chebyshev2
        rp === nothing && throw(ArgumentError("For Chebyshev2 filter rp must be given."))
        prototype = Chebyshev2(order, rp)
    end
    if fprototype === :elliptic
        rs === nothing && throw(ArgumentError("For Elliptic filter rs must be given."))
        rp === nothing && throw(ArgumentError("For Elliptic filter rp must be given."))
        rp > rs && throw(ArgumentError("For Elliptic filter rp must be less than rs."))
        prototype = Elliptic(order, rp, rs)
    end

    eeg_filter = digitalfilter(responsetype, prototype)

    dir === :twopass && (signal_filtered = filtfilt(eeg_filter, signal))
    dir === :onepass && (signal_filtered = filt(eeg_filter, signal))
    dir === :onepass_reverse && (signal_filtered = filt(eeg_filter, reverse(signal)))

    return signal_filtered
end

"""
    signal_filter(signal; fprototype, ftype, cutoff, fs, order=8, window=hanning(64))

Filters `signal` using zero phase distortion filter.

# Arguments

- `signal::Array{Float64, 3}`
- `fprototype::Symbol[:butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir]
- `ftype::Symbol[:lp, :hp, :bp, :bs]` - filter type
- `cutoff::Union{Int64, Float64, Vector{Int64}, Vector{Float64}}` - filter cutoff in Hz (vector for `:bp` and `:bs`)
- `fs::Int64` - sampling rate
- `order::Int64` - filter order
- `rp::Union{Int64, Float64, Nothing}` - dB ripple in the passband
- `rs::Union{Int64, Float64, Nothing}` - dB attentuation in the stopband
- `dir:Symbol[:oenpass, :onepass_reverse, :twopass]` - filter direction
- `window::Union{Vector{Float64}, Nothing}` - window, required for FIR filter

# Returns

- `signal_filtered::Array{Float64, 3}`
"""
function signal_filter(signal::Array{Float64, 3}; fprototype::Symbol, ftype::Symbol, cutoff::Union{Int64, Float64, Vector{Int64}, Vector{Float64}}, fs::Int64, order::Int64, rp::Union{Int64, Float64, Nothing}=nothing, rs::Union{Int64, Float64, Nothing}=nothing, dir::Symbol=:twopass, window::Union{Vector{Float64}, Nothing}=nothing)
    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)
    signal_filtered = zeros(size(signal))

    Threads.@threads for epoch in 1:epochs_no
        for idx in 1:channels_no
            signal_filtered[idx, :, epoch] = signal_filter(signal[idx, :, epoch],
                                                           fprototype=fprototype,
                                                           ftype=ftype,
                                                           cutoff=cutoff,
                                                           fs=fs,
                                                           order=order,
                                                           rp=rp,
                                                           rs=rs,
                                                           dir=dir,
                                                           window=window)
        end
    end

    return signal_filtered
end

"""
    signal_downsample(signal; t, new_sr)

Downsamples the`signal` to `new_sr` sampling frequency.

# Arguments

- `signal::Vector{Float64}`
- `t::AbstractRange`
- `new_sr::Int64` - new sampling rate

# Returns

- `signal_downsampled::Vector{Float64}`
- `t_downsampled::Vector{Float64}`
"""
function signal_downsample(signal::Vector{Float64}; t::AbstractRange, new_sr::Int64)
    new_sr < 1 && throw(ArgumentError("New sampling rate must be positive."))

    # sampling interval
    dt = t[2] - t[1]
    # sampling rate
    sr = 1 / dt
    new_sr > sr && throw(ArgumentError("New sampling rate must be lower than signal sampling rate."))
    new_sr == sr && return(signal)
    sr_ratio = new_sr / sr
    # downsample
    signal_downsampled = resample(signal, sr_ratio)
    t = collect(t)
    t_downsampled = t[1]:1/new_sr:t[end]

    return signal_downsampled, t_downsampled
end

"""
    signal_downsample(signal; t, new_sr)

Downsamples all channels of the`signal` to `new_sr` sampling frequency.

# Arguments

- `signal::Array{Float64, 3}`
- `new_sr::Int64` - new sampling rate
- `t::AbstractRange`

# Returns

- `signal_downsampled::Array{Float64, 3}`
- `t_downsampled::AbstractRange`
"""
function signal_downsample(signal::Array{Float64, 3}; t::AbstractRange, new_sr::Int64)
    new_sr < 1 && throw(ArgumentError("New sampling rate must be positive."))
    
    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)

    signal_downsampled_length = length(signal_downsample(signal[1, :, 1], t=t, new_sr=new_sr)[1])
    signal_downsampled = zeros(channels_no, signal_downsampled_length, epochs_no) 

    t_downsampled = nothing
    Threads.@threads for epoch in 1:epochs_no
        for idx in 1:channels_no
            signal_downsampled[idx, :, epoch], t_downsampled = signal_downsample(signal[idx, :, epoch], t=t, new_sr=new_sr)
        end
    end

    return signal_downsampled, t_downsampled
end

"""
    signal_psd(signal; fs, normalize=false)

Calculates power spectrum density of the `signal`.

# Arguments
- `signal::Vector{Float64}`
- `fs::Int64` - sampling rate
- `normalize::Bool` - normalize do dB

# Returns

- `psd_pow::Vector{Float64}`
- `psd_frq::Vector{Float64}`
"""
function signal_psd(signal::Vector{Float64}; fs::Int64, normalize::Bool=false)
    fs < 1 && throw(ArgumentError("Sampling rate must be positive."))
    
    psd = welch_pgram(signal, 4*fs, fs=fs)
    psd_pow = power(psd)
    psd_frq = freq(psd)
    normalize == true && (psd_pow = pow2db.(psd_pow))

    return psd_pow, Vector(psd_frq)
end

"""
    signal_psd(signal; fs, normalize=false)

Calculates power spectrum density for each the `signal` channels.

# Arguments

- `signal::Matrix{Float64}`
- `fs::Int64` - sampling rate
- `normalize::Bool` - normalize do dB

# Returns

- `psd_pow::Matrix{Float64}`
- `psd_frq::Matrix{Float64}`
"""
function signal_psd(signal::Matrix{Float64}; fs::Int64, normalize::Bool=false)
    fs < 1 && throw(ArgumentError("Sampling rate must be positive."))

    channels_no = size(signal, 1)
    psd_length, _ = signal_psd(signal[1, :], fs=fs, normalize=normalize)
    psd_pow = zeros(channels_no, length(psd_length))
    psd_frq = zeros(channels_no, length(psd_length))
    Threads.@threads for idx in 1:channels_no
        psd_pow[idx, :], psd_frq[idx, :] = signal_psd(signal[idx, :],
                                                      fs=fs,
                                                      normalize=normalize)
    end

    return psd_pow, psd_frq
end

"""
    signal_psd(signal; fs, normalize=false)

Calculates power spectrum density for each the `signal` channels.

# Arguments

- `signal::Array{Float64, 3}`
- `fs::Int64` sampling rate
- `normalize::Bool` - normalize do dB

# Returns

- `psd_pow::Array{Float64, 3}`
- `psd_frq::Array{Float64, 3}`
"""
function signal_psd(signal::Array{Float64, 3}; fs::Int64, normalize::Bool=false)
    fs < 1 && throw(ArgumentError("Sampling rate must be positive."))

    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)
    psd_length, _ = signal_psd(signal[1, :, 1], fs=fs, normalize=normalize)
    psd_pow = zeros(channels_no, length(psd_length), epochs_no)
    psd_frq = zeros(channels_no, length(psd_length), epochs_no)
    Threads.@threads for epoch in 1:epochs_no
        for idx in 1:channels_no
            psd_pow[idx, :, epoch], psd_frq[idx, :, epoch] = signal_psd(signal[idx, :, epoch],
                                                                        fs=fs,
                                                                        normalize=normalize)
        end
    end

    return psd_pow, psd_frq
end

"""
    signal_stationarity_hilbert(signal::Vector{Float64})

Calculates phase stationarity using Hilbert transformation.

# Arguments

- `signal::Vector{Float64}`

# Returns

- `phase_stationarity::Vector{Float64}`

"""
function signal_stationarity_hilbert(signal::Vector{Float64})
    
    phase_stationarity = diff(DSP.unwrap(angle.(hilbert(signal))))
    
    return phase_stationarity
end

"""
    signal_stationarity_hilbert(signal)

Calculates phase stationarity using Hilbert transformation.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `phase_stationarity::Array{Float64, 3}`
"""
function signal_stationarity_hilbert(signal::Array{Float64, 3})
    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)
    phase_stationarity = zeros(channels_no, size(signal, 2) - 1, epochs_no)
    Threads.@threads for epoch in 1:epochs_no
        for idx in 1:channels_no
            phase_stationarity[idx, :, epoch] = signal_stationarity_hilbert(signal[idx, :, epoch])
        end
    end

    return phase_stationarity
end

"""
    signal_stationarity_mean(signal::Vector{Float64})

Calculates mean stationarity.

# Arguments

- `signal::Vector{Float64}`
- `window::Int64` - time window in samples

# Returns

- `mean_stationarity::Vector{Float64}`

"""
function signal_stationarity_mean(signal::Vector{Float64}; window::Int64)
    signal = signal[1:(window * floor(Int64, length(signal) / window))]
    signal = reshape(signal, Int(length(signal) / window), window)
    mean_stationarity = mean(signal, dims=1)

    return mean_stationarity
end


"""
    signal_stationarity_var(signal::Vector{Float64})

Calculates variance stationarity.

# Arguments

- `signal::Vector{Float64}`
- `window::Int64` - time window in samples

# Returns

- `var_stationarity::Vector{Float64}`

"""
function signal_stationarity_var(signal::Vector{Float64}; window::Int64)
    signal = signal[1:(window * floor(Int64, length(signal) / window))]
    signal = reshape(signal, Int(length(signal) / window), window)
    var_stationarity = var(signal, dims=1)

    return var_stationarity
end

"""
    signal_stationarity(signal; window=10, method=:euclid)

Calculates stationarity.

# Arguments

- `signal:Array{Float64, 3}`
- `window::Int64` - time window in samples
- `method::Symbol[:mean, :var, :euclid, :hilbert]

# Returns

- `stationarity::Union{Matrix{Float64}, Array{Float64, 3}}`

"""
function signal_stationarity(signal::Array{Float64, 3}; window::Int64=10, method::Symbol=:hilbert)
    method in [:mean, :var, :euclid, :hilbert] || throw(ArgumentError("Method must be must be :mean, :var, :euclid or :hilbert."))
    (typeof(window) == Int64 && window < 1) && throw(ArgumentError("Time window must be ≥ 1."))
    (typeof(window) == Int64 && window > size(signal, 2)) && throw(ArgumentError("Time window must be ≤ epoch duration in samples."))
    
    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)

    if method === :mean
        stationarity = zeros(channels_no, window, epochs_no)
        Threads.@threads for epoch in 1:epochs_no
            for idx = 1:channels_no
                stationarity[idx, :, epoch] = signal_stationarity_mean(signal[idx, :, epoch], window=window)
            end
        end
    end

    if method === :var
        stationarity = zeros(channels_no, window, epochs_no)
        Threads.@threads for epoch in 1:epochs_no
            for idx = 1:channels_no
                stationarity[idx, :, epoch] = signal_stationarity_var(signal[idx, :, epoch], window=window)
            end
        end
    end

    if method === :hilbert
        stationarity = zeros(channels_no, size(signal, 2) - 1, epochs_no)
        Threads.@threads for epoch in 1:epochs_no
            for idx = 1:channels_no
                stationarity[idx, :, epoch] = signal_stationarity_hilbert(signal[idx, :, epoch])
            end
        end
    end

    if method === :euclid
        # number of time windows per epoch
        window_n = size(signal, 2)
        cov_mat = zeros(channels_no, channels_no, window_n, epochs_no)
        stationarity = zeros(1 + length(2:window:window_n), epochs_no)

        Threads.@threads for epoch in 1:epochs_no
            for idx = 1:window_n
                cov_mat[:, :, idx, epoch] = signal_cov(signal[:, idx, epoch], signal[:, idx, epoch])
            end
        end

        Threads.@threads for epoch in 1:epochs_no
            phase_idx = 1
            for idx = 2:window:window_n
                stationarity[phase_idx, epoch] = euclidean(cov_mat[:, :, idx - 1, epoch],
                                                           cov_mat[:, :, idx, epoch])
                phase_idx += 1
            end
        end
    end

    size(stationarity, 3) == 1 && reshape(stationarity, size(stationarity, 1), size(stationarity, 2))

    return stationarity
end

"""
    signal_trim(signal::Vector{Float64}; trim_len::Int64)

Removes `trim_len` samples from the beginning (`from` = :start, default) or end (`from` = :end) of the `signal`.

# Arguments

- `signal::Vector{Float64}; trim_len::Int64`
- `offset::Int64` - offset from which trimming starts, only works for `from` = :start
- `from::Symbol[:start, :end]

# Returns

- `signal_trimmed::Vector{Float64}`

"""
function signal_trim(signal::Vector{Float64}; trim_len::Int64, offset::Int64=0, from::Symbol=:start)
    from in [:start, :end] || throw(ArgumentError("Argument from must be :start or :end."))
    trim_len < 0 && throw(ArgumentError("Trim length must be ≥ 1."))
    trim_len >= length(signal) && throw(ArgumentError("Trim length must be less than signal length."))
    offset < 0 && throw(ArgumentError("Offset must be ≥ 1."))
    offset >= length(signal) - 1 && throw(ArgumentError("Offset must be less than signal length."))
    (from ===:start && 1 + offset + trim_len > length(signal)) && throw(ArgumentError("Offset + trim length must be less than signal length."))
    
    from === :start && (signal_trimmed = vcat(signal[1:offset], signal[(1 + offset + trim_len):end]))
    from === :end && (signal_trimmed = signal[1:(end - trim_len)])
    
    return signal_trimmed::Vector{Float64}
end


"""
    signal_trim(signal::Array{Float64, 3}; trim_len, offset=0, from=:start)

Removes `trim_len` samples from the beginning (`from` = :start, default) or end (`from` = :end) of the `signal`.

# Arguments

- `signal::Array{Float64, 3}`
- `trim_len::Int64` - number of samples to remove
- `offset::Int64` - offset from which trimming starts, only works for `from` = :start
- `from::Symbol[:start, :end]`

# Returns

- `signal_trimmed::Array{Float64, 3}`

"""
function signal_trim(signal::Array{Float64, 3}; trim_len::Int64, offset::Int64=0, from::Symbol=:start)
    from in [:start, :end] || throw(ArgumentError("Argument from must be :start or :end."))
    trim_len < 0 && throw(ArgumentError("Trim length must be ≥ 1."))
    trim_len >= size(signal, 2) && throw(ArgumentError("Trim length must be less than signal length."))
    offset < 0 && throw(ArgumentError("Offset must be ≥ 1."))
    offset >= size(signal, 2) - 1 && throw(ArgumentError("Offset must be less than signal length."))
    (from ===:start && 1 + offset + trim_len > size(signal, 2)) && throw(ArgumentError("Offset + trim length must be less than signal length."))
    
    channels_no = size(signal, 1)
    epochs_no = size(signal, 3)

    signal_trimmed = zeros(channels_no, (size(signal, 2) - trim_len), epochs_no)

    Threads.@threads for epoch in 1:epochs_no
        for idx = 1:channels_no
            signal_trimmed[idx, :, epoch] = signal_trim(signal[idx, :, epoch], trim_len=trim_len, from=from)
        end
    end
    
    return signal_trimmed
end