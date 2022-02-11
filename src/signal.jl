"""
    signal_derivative(signal)

Returns the derivative of the `signal` with length same as the signal.

# Arguments

- `signal::Union{Vector{Int64}, Vector{Float64}}`

# Returns

- `s_derivative::Union{Vector{Int64}, Vector{Float64}}`

"""
signal_derivative(signal::Union{Vector{Int64}, Vector{Float64}}) = vcat(diff(signal), diff(signal)[end])

"""
    signal_derivative(signal)

Returns the derivative of each the `signal` channels with length same as the signal.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `s_derivative::Array{Float64, 3}`
"""
function signal_derivative(signal::Array{Float64, 3})
    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)
    s_der = zeros(size(signal))

    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
            s_der[idx, :, epoch] = signal_derivative(signal[idx, :, epoch])
        end
    end

    return s_der
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

    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)
    stp = zeros(channels_n, epochs_n)

    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
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

    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)
    sbp = zeros(channels_n, epochs_n)

    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
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

- `s_fft::Vector{ComplexF64}`
- `s_sf::Vector{Float64}`
"""
function signal_make_spectrum(signal::Vector{Float64}, fs::Int64)
    fs < 1 && throw(ArgumentError("Sampling rate must be ≥ 1 Hz."))

    s_fft = fft(signal)
    # number of samples
    n = length(signal)
    # time between samples
    d = 1 / fs
    s_sf = fftfreq(n, d)

    return s_fft, Vector(s_sf)
end

"""
    signal_make_spectrum(signal, fs)

Returns FFT and DFT sample frequencies for a DFT for each the `signal` channels.

# Arguments

- `signal::Array{Float64, 3}`
- `fs::Int64` - sampling rate

# Returns

- `s_fft::Array{ComplexF64, 3}`
- `s_sf::Array{Float64, 3}`
"""
function signal_make_spectrum(signal::Array{Float64, 3}, fs::Int64)
    fs < 1 && throw(ArgumentError("Sampling rate must be ≥ 1 Hz."))

    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)
    s_fft = zeros(ComplexF64, size(signal))
    s_sf = zeros(size(signal))

    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
            s_fft[idx, :, epoch], s_sf[idx, :, epoch] = signal_make_spectrum(signal[idx, :, epoch], fs)
        end
    end

    return s_fft, s_sf
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

- `s_detrended::Vector{Float64}`
"""
function signal_detrend(signal::Vector{Float64}; type::Symbol=:linear)
    type in [:linear, :constant] || throw(ArgumentError("Trend type must be :linear or :constant."))

    if type === :constant
        s_det = signal_demean(signal)
    else
        A = ones(length(signal))
        coef = A \ signal
        s_det = @. signal - dot(A, coef)
    end

    return s_det
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

- `s_detrended::Array{Float64, 3}`
"""
function signal_detrend(signal::Array{Float64, 3}; type::Symbol=:linear)
    type in [:linear, :constant] || throw(ArgumentError("Trend type must be :linear or :constant."))

    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)
    s_det = zeros(size(signal))
    
    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
            s_det[idx, :, epoch] = signal_detrend(signal[idx, :, epoch], type=type)
        end
    end

    return s_det
end

"""
    signal_ci95(signal; n::Int64=3, method=:normal)

Calculates mean, std and 95% confidence interval for `signal`.

# Arguments

- `signal::Vector{Float64}`
- `n::Int64` - number of bootstraps
- `method::Symbol[:normal, :boot]` - use normal method or `n`-times boostrapping

# Returns

- `s_m::Float64`
- `s_s::Float64`
- `s_u::Float64`
- `s_l::Float64`
"""
function signal_ci95(signal::Vector{Float64}; n::Int=3, method::Symbol=:normal)
    method === :normal || throw(ArgumentError("For vector signal method must be :normal."))
    n < 1 && throw(ArgumentError("n must be ≥ 1."))

    s_m = mean(signal)
    s_s = std(signal) / sqrt(length(signal))
    s_u = s_m + 1.96 * s_s
    s_l = s_m - 1.96 * s_s

    return s_m, s_s, s_u, s_l
end

"""
    signal_ci95(signal::Matrix{Float64}; n::Int=3, method::Symbol=:normal)

Calculates mean, std and 95% confidence interval for each the `signal` channels.

# Arguments

- `signal::Matrix{Float64}`
- `n::Int64` - number of bootstraps
- `method::Symbol[:normal, :boot]` - use normal method or `n`-times boostrapping

# Returns

- `s_m::Vector{Float64}`
- `s_s::Vector{Float64}`
- `s_u::Vector{Float64}`
- `s_l::Vector{Float64}`
"""
function signal_ci95(signal::Matrix{Float64}; n::Int=3, method::Symbol=:normal)
    method in [:normal, :boot] || throw(ArgumentError("Method must be :normal or :boot."))

    if method === :normal
        s_m = mean(signal, dims=1)'
        s_s = std(signal, dims=1)' / sqrt(size(signal, 1))
        s_u = s_m + 1.96 * s_s
        s_l = s_m - 1.96 * s_s
    else
        s_tmp1 = zeros(size(signal, 1) * n, size(signal, 2))
        Threads.@threads for idx1 in 1:size(signal, 1) * n
            s_tmp2 = zeros(size(signal))
            sample_idx = rand(1:size(signal, 1), size(signal, 1))
            for idx2 in 1:size(signal, 1)
                s_tmp2[idx2, :] = signal[sample_idx[idx2], :]'
            end
            s_tmp1[idx1, :] = mean(s_tmp2, dims=1)
        end
        s_m = mean(s_tmp1, dims=1)'
        s_s = std(s_tmp1, dims=1)' / sqrt(size(s_tmp1, 1))
        s_sorted = sort(s_tmp1, dims=1)
        s_l = s_sorted[round(Int, 0.025 * size(s_tmp1, 1)), :]
        s_u = s_sorted[round(Int, 0.975 * size(s_tmp1, 1)), :]
    end

    return vec(s_m[:, 1]), vec(s_s[:, 1]), vec(s_u[:, 1]), vec(s_l[:, 1])
end

"""
    signal_ci95(signal; n::Int64=3, method=:normal)

Calculates mean, std and 95% confidence interval for each the `signal` channels.

# Arguments

- `signal::Array{Float64, 3}`
- `n::Int64` - number of bootstraps
- `method::Symbol[:normal, :boot]` - use normal method or `n`-times boostrapping

# Returns

- `s_m::Matrix{Float64}`
- `s_s::Matrix{Float64}`
- `s_u::Matrix{Float64}`
- `s_l::Matrix{Float64}`
"""
function signal_ci95(signal::Array{Float64, 3}; n::Int=3, method::Symbol=:normal)
    method in [:normal, :boot] || throw(ArgumentError("Method must be :normal or :boot."))
    n < 1 && throw(ArgumentError("n must be ≥ 1."))

    s_m = zeros(size(signal, 3), size(signal, 2))
    s_s = zeros(size(signal, 3), size(signal, 2))
    s_u = zeros(size(signal, 3), size(signal, 2))
    s_l = zeros(size(signal, 3), size(signal, 2))

    epochs_n = size(signal, 3)

    Threads.@threads for epoch in 1:epochs_n
        s_m[epoch, :], s_s[epoch, :], s_u[epoch, :], s_l[epoch, :] = signal_ci95(signal[:, :, epoch])
    end

    return s_m, s_s, s_u, s_l
end

"""
    signal_mean(signal1, signal2)

Calculates mean and 95% confidence interval for 2 signals.

# Arguments

- `signal1::Vector{Float64}`
- `signal2:Vector{Float64}`

# Returns

- `s_mean::Float64`
- `s_sd::Float64`
- `s_u::Float64`
- `s_l::Float64`
"""
function signal_mean(signal1::Vector{Float64}, signal2::Vector{Float64})
    length(signal1) != length(signal2) && throw(ArgumentError("Both signals must be of the same as size."))

    s_m = zeros(length(signal1))
    s_s = zeros(length(signal1))
    s_u = zeros(length(signal1))
    s_l = zeros(length(signal1))

    s1_mean = mean(signal1)
    s2_mean = mean(signal2)
    s_m = s1_mean - s2_mean
    s1_sd = std(signal1) / sqrt(length(signal1))
    s2_sd = std(signal2) / sqrt(length(signal2))
    s_s = sqrt(s1_sd^2 + s2_sd^2)
    s_u = s_m + 1.96 * s_s
    s_l = s_m - 1.96 * s_s

    return s_m, s_s, s_u, s_l
end

"""
    signal_mean(signal1, signal2)

Calculates mean and 95% confidence interval for 2 signals.

# Arguments

- `signal1::Array{Float64, 3}`
- `signal2:Array{Float64, 3}`

# Returns

- `s_mean::Matrix{Float64}`
- `s_sd::Matrix{Float64}`
- `s_u::Matrix{Float64}`
- `s_l::Matrix{Float64}`
"""
function signal_mean(signal1::Array{Float64, 3}, signal2::Array{Float64, 3})
    size(signal1) != size(signal2) && throw(ArgumentError("Both signals must be of the same as size."))

    s_m = zeros(size(signal1, 3), size(signal1, 2))
    s_s = zeros(size(signal1, 3), size(signal1, 2))
    s_u = zeros(size(signal1, 3), size(signal1, 2))
    s_l = zeros(size(signal1, 3), size(signal1, 2))
    epochs_n = size(signal1, 3)

    Threads.@threads for epoch in 1:epochs_n
        s1_mean = mean(signal1[:, :, epoch], dims=1)
        s2_mean = mean(signal2[:, :, epoch], dims=1)
        s_m[epoch, :] = s1_mean - s2_mean
        s1_sd = std(signal1[:, :, epoch], dims=1) / sqrt(size(signal1[:, :, epoch], 2))
        s2_sd = std(signal2[:, :, epoch], dims=1) / sqrt(size(signal2[:, :, epoch], 2))
        s_s[epoch, :] = sqrt.(s1_sd.^2 .+ s2_sd.^2)
        s_u[epoch, :] = @. s_m[epoch, :] + 1.96 * s_s[epoch, :]
        s_l[epoch, :] = @. s_m[epoch, :] - 1.96 * s_s[epoch, :]
    end

    return s_m, s_s, s_u, s_l
end

"""
    signal_difference(signal1, signal2; n::Int64=3, method=:absdiff)

Calculates mean difference and 95% confidence interval for 2 signals.

# Arguments

- `signal1::Matrix{Float64}`
- `signal2:Matrix{Float64}`
- `n::Int64` - number of bootstraps
- `method::Symbol[:absdiff, :diff2int]`
    - `:absdiff` - maximum difference
    - `:diff2int` - integrated area of the squared difference

# Returns

- `s_statistic::Vector{Float64}`
- `s_statistic_single::Float64`
- `p::Float64`
"""
function signal_difference(signal1::Matrix{Float64}, signal2::Matrix{Float64}; n::Int64=3, method::Symbol=:absdiff)
    size(signal1) != size(signal2) && throw(ArgumentError("Both signals must be of the same size."))
    method in [:absdiff, :diff2int] || throw(ArgumentError("Method must be :absdiff or :diff2int."))

    s1_mean = mean(signal1, dims=1)'
    s2_mean = mean(signal2, dims=1)'

    if method === :absdiff
        # statistic: maximum difference
        s_diff = s1_mean - s2_mean
        s_statistic_single = maximum(abs.(s_diff))
    else
        # statistic: integrated area of the squared difference
        s_diff_squared = (s1_mean - s2_mean).^2
        s_statistic_single = simpson(s_diff_squared)
    end

    signals = [signal1; signal2]
    s_statistic = zeros(size(signal1, 1) * n)

    Threads.@threads for idx1 in 1:(size(signal1, 1) * n)
        s_tmp1 = zeros(size(signal1, 1), size(signal1, 2))
        sample_idx = rand(1:size(signals, 1), size(signals, 1))
        # sample_idx = sample_idx[1:1000]
        for idx2 in 1:size(signal1, 1)
            s_tmp1[idx2, :] = signals[sample_idx[idx2], :]'
        end
        s1_mean = mean(s_tmp1, dims=1)
        s_tmp1 = zeros(size(signal1, 1), size(signal1, 2))
        sample_idx = rand(1:size(signals, 1), size(signals, 1))
        for idx2 in 1:size(signal1, 1)
            s_tmp1[idx2, :] = signals[sample_idx[idx2], :]'
        end
        s2_mean = mean(s_tmp1, dims=1)
        if method === :absdiff
            # statistic: maximum difference
            s_diff = s1_mean - s2_mean
            s_statistic[idx1] = maximum(abs.(s_diff))
        else
            # statistic: integrated area of the squared difference
            s_diff_squared = (s1_mean - s2_mean).^2
            s_statistic[idx1] = simpson(s_diff_squared)
        end
    end

    p = length(s_statistic[s_statistic .> s_statistic_single]) / size(signal1, 1) * n
    p > 1 && (p = 1.0)

    return s_statistic, s_statistic_single, p
end

"""
    signal_difference(signal1, signal2; n::Int64=3, method=:absdiff)

Calculates mean difference and 95% confidence interval for 2 signals.

# Arguments

- `signal1::Array{Float64, 3}`
- `signal2:Array{Float64, 3}`
- `n::Int64` - number of bootstraps
- `method::Symbol[:absdiff, :diff2int]`
    - `:absdiff` - maximum difference
    - `:diff2int` - integrated area of the squared difference

# Returns

- `s_statistic::Matrix{Float64}`
- `s_statistic_single::Vector{Float64}`
- `p::Vector{Float64}`
"""
function signal_difference(signal1::Array{Float64, 3}, signal2::Array{Float64, 3}; n::Int64=3, method::Symbol=:absdiff)
    size(signal1) != size(signal2) && throw(ArgumentError("Both signals must be of the same size."))
    method in [:absdiff, :diff2int] || throw(ArgumentError("Method must be :absdiff or :diff2int."))

    epochs_n = size(signal1, 3)
    s_statistic = zeros(epochs_n, size(signal1, 1) * n)
    s_statistic_single = zeros(epochs_n)
    p = zeros(epochs_n)

    Threads.@threads for epoch in 1:epochs_n
        s_statistic[epoch, :], s_statistic_single[epoch], p[epoch] = signal_difference(signal1[:, :, epoch], signal2[:, :, epoch])
    end

    return s_statistic, s_statistic_single, p
end

"""
   signal_autocov(signal; lag=1, demean=false, norm=false)

Calculates autocovariance of the `signal`.

# Arguments

- `signal::Vector{Float64}`
- `lag::Int64` - lags range is `-lag:lag`
- `demean::Bool` - demean `signal` prior to calculations
- `norm::Bool` - normalize autocovariance

# Returns

- `acov::Vector{Float64}`
- `lags::Vector{Int64}`
"""
function signal_autocov(signal::Vector{Float64}; lag::Int64=1, demean::Bool=false, norm::Bool=false)
    lag < 1 && throw(ArgumentError("Lag must be ≥ 1."))

    lags = collect(-lag:lag)

    if demean == true
        s_demeaned = signal_demean(signal)
    else
        s_demeaned = signal
    end

    acov = zeros(length(lags))

    for idx in 1:length(lags)
        if lags[idx] == 0
            # no lag
            s_lagged = s_demeaned
            s_mul = s_demeaned .* s_lagged
        elseif lags[idx] > 0
            # positive lag
            s_lagged = s_demeaned[1:(end - lags[idx])]
            s_mul = s_demeaned[(1 + lags[idx]):end] .* s_lagged
        elseif lags[idx] < 0
            # negative lag
            s_lagged = s_demeaned[(1 + abs(lags[idx])):end]
            s_mul = s_demeaned[1:(end - abs(lags[idx]))] .* s_lagged
        end
        s_sum = sum(s_mul)
        if norm == true
            acov[idx] = s_sum / length(signal)
        else
            acov[idx] = s_sum
        end
    end

    return acov, lags
end

"""
   signal_autocov(signal; lag=1, demean=false, norm=false)

Calculates autocovariance of each the `signal` channels.

# Arguments

- `signal::Array{Float64, 3}`
- `lag::Int64` - lags range is `-lag:lag`
- `demean::Bool` - demean signal prior to analysis
- `norm::Bool` - normalize autocovariance

# Returns

- `acov::Matrix{Float64}`
- `lags::Vector{Int64}`
"""
function signal_autocov(signal::Array{Float64, 3}; lag::Int64=1, demean::Bool=false, norm::Bool=false)
    lag < 1 && throw(ArgumentError("Lag must be ≥ 1."))

    lags = collect(-lag:lag)
    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)
    acov = zeros(channels_n, length(lags), epochs_n)

    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
            acov[idx, :, epoch], lags = signal_autocov(signal[idx, :, epoch],
                                                 lag=lag,
                                                 demean=demean,
                                                 norm=norm)
        end
    end

    return acov, lags
end

"""
   signal_crosscov(signal1, signal2; lag=1, demean=false, norm=false)

Calculates cross-covariance between `signal1` and `signal2`.

# Arguments

- `signal1::Vector{Float64}`
- `signal2::Vector{Float64}`
- `lag::Int64` - lags range is `-lag:lag`
- `demean::Bool` - demean signal prior to analysis
- `norm::Bool` - normalize cross-covariance

# Returns

- `ccov::Vector{Float64}`
- `lags::Vector{Int64}`
"""
function signal_crosscov(signal1::Vector{Float64}, signal2::Vector{Float64}; lag::Int64=1, demean::Bool=false, norm::Bool=false)
    length(signal1) != length(signal2) && throw(ArgumentError("Both signals must be of the same as length."))
    lag < 1 && throw(ArgumentError("Lag must be ≥ 1."))

    lags = collect(-lag:lag)

    if demean == true
        s_demeaned1 = signal_demean(signal1)
        s_demeaned2 = signal_demean(signal2)
    else
        s_demeaned1 = signal1
        s_demeaned2 = signal2
    end

    ccov = zeros(length(lags))

    for idx in 1:length(lags)
        if lags[idx] == 0
            # no lag
            s_lagged = s_demeaned2
            s_mul = s_demeaned1 .* s_lagged
        elseif lags[idx] > 0
            # positive lag
            s_lagged = s_demeaned2[1:(end - lags[idx])]
            s_mul = s_demeaned1[(1 + lags[idx]):end] .* s_lagged
        elseif lags[idx] < 0
            # negative lag
            s_lagged = s_demeaned2[(1 + abs(lags[idx])):end]
            s_mul = s_demeaned1[1:(end - abs(lags[idx]))] .* s_lagged
        end
        s_sum = sum(s_mul)
        if norm == true
            ccov[idx] = s_sum / length(signal1)
        else
            ccov[idx] = s_sum
        end
    end

    return ccov, lags
end

"""
   signal_crosscov(signal; lag=1, demean=false, norm=false)

Calculates cross-covariance between all channels in the `signal`.

# Arguments

- `signal::Matrix{Float64}` - the signal
- `lag::Int64` - lags range is `-lag:lag`
- `demean::Bool` - demean `signal` prior to analysis
- `norm::Bool` - normalize cross-covariance

# Returns

- `ccov::Array{Float64, 3}`
- `lags::Vector{Int64}`
"""
function signal_crosscov(signal::Array{Float64, 3}; lag::Int64=1, demean::Bool=false, norm::Bool=false)
    lag < 1 && throw(ArgumentError("Lag must be ≥ 1 Hz."))

    lags = collect(-lag:lag)
    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)
    ccov = zeros(channels_n^2, length(lags), epochs_n)

    Threads.@threads for epoch in 1:epochs_n
        ccov_packed = Array{Vector{Float64}}(undef, channels_n, channels_n)
        for idx1 in 1:channels_n
            for idx2 in 1:channels_n
                ccov_packed[idx1, idx2], lags = signal_crosscov(signal[idx1, :, epoch],
                                                          signal[idx2, :, epoch],
                                                          lag=lag,
                                                          demean=demean,
                                                          norm=norm)
            end
        end
        for idx in 1:channels_n^2
            ccov[idx, :, epoch] = ccov_packed[idx]
        end
    end

    return ccov, lags
end

"""
   signal_crosscov(signal1, signal2; lag=1, demean=false, norm=false)

Calculates cross-covariance between same channels in `signal1` and `signal2`.

# Arguments

- `signal1::Array{Float64, 3}`
- `signal2::Array{Float64, 3}`
- `lag::Int64` - lags range is `-lag:lag`
- `demean::Bool` - demean signal prior to analysis
- `norm::Bool` - normalize cross-covariance

# Returns

- `ccov::Array{Float64, 3}`
- `lags::Vector{Int64}`
"""
function signal_crosscov(signal1::Array{Float64, 3}, signal2::Array{Float64, 3}; lag::Int64=1, demean::Bool=false, norm::Bool=false)
    size(signal1) != size(signal2) && throw(ArgumentError("Both arrays must be of the same as size."))
    lag < 1 && throw(ArgumentError("Lag must be ≥ 1 Hz."))

    lags = collect(-lag:lag)
    channels_n = size(signal1, 1)
    epochs_n = size(signal1, 3)
    ccov = zeros(channels_n, length(lags), epochs_n)

    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
            ccov[idx, :, epoch], lags = signal_crosscov(signal1[idx, :, epoch],
                                                  signal2[idx, :, epoch],
                                                  lag=lag,
                                                  demean=demean,
                                                  norm=norm)
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

- `fft::Vector(ComplexF64}`
- `amplitudes::Vector{Float64}`
- `powers::Vector{Float64}`
- `phases::Vector{Float64}
"""
function signal_spectrum(signal::Vector{Float64}; pad::Int64=0)
    pad < 0 && throw(ArgumentError("Pad cannot be negative."))

    s_fft = fft0(signal, pad)

    # normalize
    s_fft ./= length(signal)

    # amplitudes
    s_amplitudes = @. 2 * abs(s_fft)

    # power
    s_powers = s_amplitudes.^2

    # phases
    s_phases = angle.(s_fft)

    return s_fft, s_amplitudes, s_powers, s_phases
end

"""
    signal_spectrum(signal; pad=0)

Calculates FFT, amplitudes, powers and phases for each channel of the `signal` matrix.

# Arguments

- `signal::Array{Float64, 3}` - the signal
- `pad::Int64` - pad the `signal` with `pad` zeros

# Returns

- `fft::Array{ComplexF64, 3}`
- `amplitudes::Array{Float64, 3}`
- `powers::Array{Float64, 3}`
- `phases::Array{Float64, 3}
"""
function signal_spectrum(signal::Array{Float64, 3}; pad::Int64=0)
    pad < 0 && throw(ArgumentError("Pad cannot be negative."))

    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)

    s_fft = zeros(ComplexF64, channels_n, size(signal, 2) + pad, epochs_n)
    s_amplitudes = zeros(channels_n, size(signal, 2) + pad, epochs_n)
    s_powers = zeros(channels_n, size(signal, 2) + pad, epochs_n)
    s_phases = zeros(channels_n, size(signal, 2) + pad, epochs_n)

    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
            s_fft[idx, :, epoch], s_amplitudes[idx, :, epoch], s_powers[idx, :, epoch], s_phases[idx, :, epoch] = signal_spectrum(signal[idx, :, epoch], pad=pad)
        end
    end

    return s_fft, s_amplitudes, s_powers, s_phases
end

"""
    signal_epochs(signal; epochs_n, epochs_len, average=true)

Splits `signal` into epochs.

# Arguments

- `signal::Vector{Float64}`
- `epochs_n::Union{Int64, Nothing}` - number of epochs
- `epochs_len::Union{Int64, Nothing}` - epoch length in samples
- `average::Bool` - average all epochs, returns one averaged epoch; if false than returns array of epochs, each row is one epoch

# Returns

- `epochs::Matrix{Float64}`
"""
function signal_epochs(signal::Vector{Float64}; epochs_n::Union{Int64, Nothing}=nothing, epochs_len::Union{Int64, Nothing}=nothing, average::Bool=false)
    (epochs_len === nothing && epochs_n === nothing) && throw(ArgumentError("Either number of epochs or epoch length must be set."))
    (epochs_len !== nothing && epochs_n !== nothing) && throw(ArgumentError("Both number of epochs and epoch length cannot be set."))
    (epochs_len != nothing && epochs_len < 1) && throw(ArgumentError("Epoch length rate must be ≥ 1."))
    (epochs_n != nothing && epochs_n < 1) && throw(ArgumentError("Epoch length rate must be ≥ 1."))

    if epochs_n === nothing
        epochs_n = length(signal) ÷ epochs_len
    end
    if epochs_len === nothing
        epochs_len = length(signal) ÷ epochs_n
    end

    epochs = zeros(epochs_n, epochs_len)

    signal = signal[1:epochs_len * epochs_n]
    epochs = reshape(signal, epochs_n, epochs_len)'

    if average == true
        epochs = vec(mean(epochs, dims=1)[1, :])
    end

    return epochs
end

"""
    signal_epochs(signal; epochs_n=nothing, epochs_len=nothing, average=true)

Splits `signal` into epochs.

# Arguments

- `signal::Array{Float64, 3}`
- `epochs_n::Union{Int64, Nothing}` - number of epochs
- `epochs_len::Union{Int64, Nothing}` - epoch length in samples
- `average::Bool` - average all epochs, returns one averaged epoch; if false than returns array of epochs, each row is one epoch

# Returns

- `epochs::Array{Float64, 3}`
"""
function signal_epochs(signal::Matrix{Float64}; epochs_n::Union{Int64, Nothing}=nothing, epochs_len::Union{Int64, Nothing}=nothing, average::Bool=false)
    (epochs_len === nothing && epochs_n === nothing) && throw(ArgumentError("Either number of epochs or epoch length must be set."))
    (epochs_len !== nothing && epochs_n !== nothing) && throw(ArgumentError("Both number of epochs and epoch length cannot be set."))
    (epochs_len != nothing && epochs_len < 1) && throw(ArgumentError("Epoch length rate must be ≥ 1."))
    (epochs_n != nothing && epochs_n < 1) && throw(ArgumentError("Epoch length rate must be ≥ 1."))

    channels_n = size(signal, 1)

    if epochs_n === nothing
        epochs_n = size(signal, 2) ÷ epochs_len
    else
        epochs_len = size(signal, 2) ÷ epochs_n
    end

    epochs = zeros(channels_n, epochs_len, epochs_n)

    # signal = signal[1:epochs_len * epochs_n]
    # epochs = reshape(signal, epochs_n, epochs_len)'

    idx1 = 1
    for idx2 in 1:epochs_len:(epochs_n * epochs_len - 1)
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

- `s_modified::Matrix{Float64}`
"""
function signal_drop_channel(signal::Matrix{Float64}, channels::Union{Int64, Vector{Int64}, AbstractRange})
    if typeof(channels) <: AbstractRange
        channels = collect(channels)
    end
    for idx in length(channels):-1:1
        channels[idx] > size(signal, 2) && throw(ArgumentError("Channels index does not match signal."))
    end

    length(channels) > 1 && (channels = sort!(channels, rev=true))
    s_new = signal[setdiff(1:end, (channels)), :]

    return s_new
end

"""
    signal_drop_channel(signal, channels)

Removes `channels` from the `signal`.

# Arguments

- `signal::Array{Float64, 3}`
- `channels::Union{Int64, Vector{Int64}, AbstractRange}` - channels to be removed, vector of numbers or range

# Returns

- `s_modified::Matrix{Float64}`
"""
function signal_drop_channel(signal::Array{Float64, 3}, channels::Union{Int64, Vector{Int64}, AbstractRange})
    if typeof(channels) <: AbstractRange
        channels = collect(channels)
    end
    for idx in length(channels):-1:1
        channels[idx] > size(signal, 2) && throw(ArgumentError("Channels index does not match signal."))
    end

    length(channels) > 1 && (channels = sort!(channels, rev=true))
    s_new = signal[setdiff(1:end, (channels)), :, :]

    return s_new
end
"""
    signal_reference_channel(signal, reference)

Re-references channels of the `signal` to specific signal channel.

# Arguments

- `signal::Matrix{Float64}`
- `reference::Union{Int64, Vector{Int64}, AbstractRange}}` - index of channels used as reference; if multiple channels are specified, their average is used as the reference

# Returns

- `s_referenced::Matrix{Float64}`
"""
function signal_reference_channel(signal::Matrix{Float64}, reference_idx::Union{Int64, Vector{Int64}, AbstractRange})
    if typeof(reference_idx) <: AbstractRange
        reference_idx = collect(reference_idx)
    end

    channels_n = size(signal, 1)
    s_ref = zeros(size(signal))

    channels_list = collect(1:channels_n)
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
    for idx in 1:channels_n
        s_ref[idx, :] = signal[idx, :] .- reference_channel
    end
    length(reference_idx) == 1 && (s_ref[reference_idx, :] = reference_channel)

    return s_ref
end

"""
    signal_reference_channel(signal, reference)

Re-references channels of the `signal` to specific signal channel.

# Arguments

- `signal::Array{Float64, 3}`
- `reference::Union{Int64, Vector{Int64}, AbstractRange}}` - index of channels used as reference; if multiple channels are specified, their average is used as the reference

# Returns

- `s_referenced::Matrix{Float64}`
"""
function signal_reference_channel(signal::Array{Float64, 3}, reference_idx::Union{Int64, Vector{Int64}, AbstractRange})
    if typeof(reference_idx) <: AbstractRange
        reference_idx = collect(reference_idx)
    end

    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)
    s_ref = zeros(size(signal))

    channels_list = collect(1:channels_n)
    for idx in 1:length(reference_idx)
        if (reference_idx[idx] in channels_list) == false
            throw(ArgumentError("Reference channel index does not match signal."))
        end
    end
    Threads.@threads for epoch in 1:epochs_n
        if length(reference_idx) == 1
            reference_channel = mean(signal[reference_idx, :, epoch], dims=2)
        else
            reference_channel = vec(mean(signal[reference_idx, :, epoch], dims=1))
        end
        for idx in 1:channels_n
            s_ref[idx, :, epoch] = signal[idx, :, epoch] .- reference_channel
        end
        length(reference_idx) == 1 && (s_ref[reference_idx, :, epoch] = reference_channel)
    end

    return s_ref
end

"""
    signal_reference_car(signal)

Re-references channels of the `signal` to common average reference.

# Arguments

- `signal::Matrix{Float64}`

# Returns

- `s_referenced::Matrix{Float64}`
"""
function signal_reference_car(signal::Matrix{Float64})
    channels_n = size(signal, 1)
    reference_channel = vec(mean(signal, dims=1))
    s_ref = zeros(size(signal))

    for idx in 1:channels_n
        s_ref[idx, :] = signal[idx, :] .- reference_channel
    end

    return s_ref
end

"""
    signal_reference_car(signal)

Re-references channels of the `signal` to common average reference.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `s_referenced::Array{Float64, 3}`
"""
function signal_reference_car(signal::Array{Float64, 3})
    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)
    s_ref = zeros(size(signal))

    Threads.@threads for epoch in 1:epochs_n
        reference_channel = vec(mean(signal[:, :, epoch], dims=1))
        for idx in 1:channels_n
            s_ref[idx, :, epoch] = signal[idx, :, epoch] .- reference_channel
        end
    end

    return s_ref
end

"""
    signal_taper(signal, taper)

Taper the `signal` with `taper`.

# Arguments

- `signal::Vector{Float64}`
- `taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `s_tapered::Vector{Union{Float64, ComplexF64}}`
"""
function signal_taper(signal::Vector{Float64}, taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})
    length(taper) == length(signal) || throw(ArgumentError("Taper length and signal length must be equal."))
    s_tap = signal .* taper

    return s_tap
end

"""
    signal_taper(signal, taper)

Tapers channels of the `signal` with `taper`.

# Arguments

- `signal::Array{Float64, 3}`
- `taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `s_tapered::Union{Array{Float64, 3}, Array{ComplexF64, 3}}`
"""
function signal_taper(signal::Array{Float64, 3}, taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})
    length(taper) == size(signal, 2) || throw(ArgumentError("Taper length and signal length must be equal."))

    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)

    s_tap = zeros(eltype(taper), size(signal))

    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
            s_tap[idx, :, epoch] = signal_taper(signal[idx, :, epoch], taper)
        end
    end

    return s_tap
end

"""
    signal_demean(signal)

Removes mean value (DC offset) from the `signal`.

# Arguments

- `signal::Vector{Float64}`

# Returns

- `s_demeaned::Vector{Float64}`
"""
function signal_demean(signal::Vector{Float64})
    s_dem = signal .- mean(signal)

    return s_dem
end

"""
    signal_demean(signal)

Removes mean value (DC offset) for each the `signal` channels.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `s_demeaned::Array{Float64, 3}`
"""
function signal_demean(signal::Array{Float64, 3})
    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)
    s_dem = zeros(size(signal))

    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
            s_dem[idx, :, epoch] = signal_demean(signal[idx, :, epoch])
        end
    end

    return s_dem
end

"""
    signal_normalize_zscore(signal)

Normalize (by z-score) `signal`.

# Arguments

- `signal::Vector{Float64}`

# Returns

- `s_normalized::Vector{Float64}`
"""
function signal_normalize_zscore(signal::Vector{Float64})
    s_norm = (signal .- mean(signal)) ./ std(signal)
    return s_norm
end

"""
    signal_normalize_zscore(signal)

Normalize (by z-score) each the `signal` channel.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `s_normalized::Array{Float64, 3}`
"""
function signal_normalize_zscore(signal::Array{Float64, 3})
    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)
    s_norm = zeros(size(signal))

    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
            s_norm[idx, :, epoch] = signal_normalize_zscore(signal[idx, :, epoch])
        end
    end

    return s_norm
end

"""
    signal_normalize_minmax(signal)

Normalize (to 0…1) `signal`.

# Arguments

- `signal::Vector{Float64}`

# Returns

- `s_normalized::Vector{Float64}`
"""
function signal_normalize_minmax(signal::Vector{Float64})
    s_norm = (signal .- minimum(signal)) ./ (maximum(signal) - minimum(signal))
    return s_norm
end

"""
    signal_normalize_minmax(signal)

Normalize (to 0…1) each the `signal` channel.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `s_normalized::Array{Float64, 3}`
"""
function signal_normalize_minmax(signal::Array{Float64, 3})
    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)
    s_norm = zeros(size(signal))

    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
            s_norm[idx, :, epoch] = signal_normalize_minmax(signal[idx, :, epoch])
        end
    end

    return s_norm
end

"""
   signal_cov(signal1, signal2; norm=true)

Calculates covariance between `signal1` and `signal2`.

# Arguments

- `signal1::Vector{Float64}`
- `signal2::Vector{Float64}`
- `norm::Bool` - normalize covariance

# Returns

- `cov_mat::Matrix{Float64}`
"""
function signal_cov(signal1::Vector{Float64}, signal2::Vector{Float64}; norm::Bool=false)
    length(signal1) == length(signal2) || throw(ArgumentError("Both signals must be of the same length."))

    # channels-vs-channels
    cov_mat = cov(signal1 * signal2')

    # normalize
    norm == true && (cov_mat = cov_mat ./ (size(cov_mat, 2) - 1))

    return cov_mat
end

"""
   signal_cov(signal; norm=true)

Calculates covariance between all channels of the `signal`.

# Arguments

- `signal::Matrix{Float64}`
- `norm::Bool` - normalize covariance

# Returns

- `cov_mat::Matrix{Float64}`
"""
function signal_cov(signal::Matrix{Float64}; norm::Bool=false)
    signal = signal'
    
    # channels-vs-channels
    cov_mat = cov(signal)

    # normalize
    norm == true && (cov_mat = cov_mat ./ (size(cov_mat, 2) - 1))

    return cov_mat
end

"""
   signal_cov(signal; norm=true)

Calculates covariance between all channels of the `signal`.

# Arguments

- `signal::Array{Float64, 3}`
- `norm::Bool` - normalize covariance

# Returns

- `cov_mat::Array{Float64, 3}`
"""
function signal_cov(signal::Array{Float64, 3}; norm::Bool=false)
    epochs_n = size(signal, 3)
    cov_mat = zeros(size(signal, 1), size(signal, 1), epochs_n)

    Threads.@threads for epoch in 1:epochs_n
        # channels-vs-channels
        cov_mat[:, :, epoch] = cov(signal[:, :, epoch]')
    end

    # normalize
    norm == true && (cov_mat = cov_mat ./ (size(cov_mat, 2) - 1))

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
    epochs_n = size(signal, 3)
    cor_mat = zeros(size(signal, 1), size(signal, 1), epochs_n)

    Threads.@threads for epoch in 1:epochs_n
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

- `s_noisy::Vector{Float64}`
"""
function signal_add_noise(signal::Vector{Float64})
    s_noise = signal .+ rand(length(signal))

    return s_noise
end

"""
    signal_add_noise(signal)

Add random noise to the `signal` channels.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `s_noisy::Array{Float64, 3}`
"""
function signal_add_noise(signal::Array{Float64, 3})
    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)
    s_noise = zeros(size(signal))

    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
            s_noise[idx, :, epoch] = signal_add_noise(signal[idx, :, epoch])
        end
    end

    return s_noise
end

"""
    signal_upsample(signal; t, new_sr)

Upsamples `signal` to `new_sr` sampling frequency.

# Arguments

- `signal::Vector{Float64}`
- `t::AbstractRange`
- `new_sr::Int64` - new sampling rate
# Returns

- `s_upsampled::Vector{Float64}`
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
    s_upsampled = signal_interpolation(t_upsampled)

    return s_upsampled, t_upsampled
end

"""
    signal_upsample(signal; t, new_sr)

Upsamples all channels of `signal` to `new_sr` sampling frequency.

# Arguments

- `signal::Array{Float64, 3}`
- `t::AbstractRange`
- `new_sr::Int64` - new sampling rate

# Returns

- `s_upsampled::Array{Float64, 3}`
- `t_upsampled::AbstractRange`
"""
function signal_upsample(signal::Array{Float64, 3}; t::AbstractRange, new_sr::Int64)
    new_sr < 1 && throw(ArgumentError("New sampling rate must be positive."))

    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)

    s_upsampled_len = length(signal_upsample(signal[1, :, 1], t=t, new_sr=new_sr)[1])
    s_upsampled = zeros(channels_n, s_upsampled_len, epochs_n) 

    t_upsampled = nothing
    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
            s_upsampled[idx, :, epoch], t_upsampled = signal_upsample(signal[idx, :, epoch], t=t, new_sr=new_sr)
        end
    end

    return s_upsampled, t_upsampled
end

"""
    signal_tconv(signal, kernel)

Performs convolution in the time domain between `signal` and `kernel`.

# Arguments

- `signal::Vector{Float64}`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `s_conv::Union{Vector{Float64}, Vector{ComplexF64}}`
"""
function signal_tconv(signal::Vector{Float64}, kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})
    s_conv = conv(signal, kernel)
    half_kernel = floor(Int, length(kernel) / 2)

    # remove in- and out- edges
    if mod(length(kernel), 2) == 0 
        s_conv = s_conv[half_kernel:(end - half_kernel)]
    else
        s_conv = s_conv[half_kernel:(end - half_kernel - 1)]
    end

    return s_conv
end

"""
    signal_tconv(signal, kernel)

Performs convolution in the time domain between `signal` and `kernel`.

# Arguments

- `signal::Array{Float64, 3}`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `s_conv::Union{Array{Float64, 3}, Array{ComplexF64, 3}}`
"""
function signal_tconv(signal::Array{Float64, 3}, kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})
    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)

    if typeof(kernel) == Vector{ComplexF64}
        s_conv = zeros(ComplexF64, size(signal))
    else
        s_conv = zeros(size(signal))
    end

    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
            s_conv[idx, :, epoch] = signal_tconv(signal[idx, :, epoch], kernel)
        end
    end

    return s_conv
end

"""
    signal_filter(signal; fprototype, ftype, cutoff, fs, order, rp, rs, dir=:twopass, d=1, window)

Filters `signal` using zero phase distortion filter.

# Arguments

- `signal::Vector{Float64}`
- `fprototype::Symbol[:mavg, :mmed, :poly, :butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir]` - filter prototype:
    - `:mavg` - moving average (with threshold and/or weight window)
    - `:mmed` - moving median (with threshold and/or weight window)
    - `:poly` - polynomial of `order` order
- `ftype::Symbol[:lp, :hp, :bp, :bs]` - filter type
- `cutoff::Union{Int64, Float64, Vector{Int64}, Vector{Float64}, Tuple, Nothing}` - filter cutoff in Hz (vector for `:bp` and `:bs`)
- `fs::Union{Int64, Nothing}` - sampling rate
- `order::Union{Int64, Nothing}` - filter order
- `rp::Union{Float64, Nothing}` - dB ripple in the passband
- `rs::Union{Float64, Nothing}` - dB attentuation in the stopband
- `dir:Symbol[:onepass, :onepass_reverse, :twopass]` - filter direction
- `d::Int64` - window length for mean average and median average filter
- `t::Union{Int64, Float64}` - threshold for :mavg and :mmed filters; threshold = threshold * std(signal) + mean(signal) for :mavg or threshold = threshold * std(signal) + median(signal) for :mmed filter
- `window::Union{Vector{Float64}, Nothing} - window, required for FIR filter

# Returns

- `s_filtered::Vector{Float64}`
"""
function signal_filter(signal::Vector{Float64}; fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Int64, Float64, Vector{Int64}, Vector{Float64}, Tuple, Nothing}=nothing, fs::Union{Int64, Nothing}=nothing, order::Union{Int64, Nothing}=nothing, rp::Union{Int64, Float64, Nothing}=nothing, rs::Union{Int64, Float64, Nothing}=nothing, dir::Symbol=:twopass, d::Int64=1, t::Union{Int64, Float64}=0, window::Union{Vector{Float64}, Nothing}=nothing)
    fprototype in [:mavg, :mmed, :poly, :butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir] || throw(ArgumentError("Filter prototype must be :mavg, :mmed,:butterworth, :chebyshev1, :chebyshev2, :elliptic or :fir."))
    (fprototype === :fir && (window === nothing || length(window) < length(signal))) && throw(ArgumentError("For FIR filter window must be shorter than signal."))
    (fprototype !== :mavg && fprototype !== :mmed && fprototype !== :poly) && (ftype in [:lp, :hp, :bp, :bs] || throw(ArgumentError("Filter type must be :bp, :hp, :bp or :bs.")))
    (fprototype !== :mavg && fprototype !== :mmed) && (fs === nothing && throw(ArgumentError("Sampling frequency must be given.")))
    dir in [:onepass, :onepass_reverse, :twopass] || throw(ArgumentError("Filter direction must be :onepass, :onepass_reverse or :twopass."))
    (order !== nothing && fprototype !== :poly) && (mod(order, 2) != 0 && throw(ArgumentError("Filter order must be even.")))
    (order === 0 && fprototype === :poly) && throw(ArgumentError("Filter order must be ≥ 1."))
    order !== nothing && (order < 0 && throw(ArgumentError("Filter order must be positive.")))
    d > length(signal) && throw(ArgumentError("Value of d cannot be higher than signal length."))

    if fprototype === :mavg
        if window === nothing
            s_filtered = signal
            for idx in d:-1:1
                if t > 0
                    if signal[idx] > t * std(signal) + mean(signal)
                        s_filtered[idx] = mean(signal[idx:idx+1])
                    end
                else
                    s_filtered[idx] = mean(signal[idx:idx+1])
                end
            end
            for idx in (1 + d):(length(signal) - d)
                if t > 0
                    if signal[idx] > t * std(signal) + mean(signal)
                        s_filtered[idx] = mean(signal[(idx - d):(idx + d)])
                    end
                else
                    s_filtered[idx] = mean(signal[(idx - d):(idx + d)])
                end
            end
            for idx in (length(signal) - d + 1):length(signal)
                if t > 0
                    if signal[idx] > t * std(signal) + mean(signal)
                        s_filtered[idx] = mean(signal[idx-1:idx])
                    end
                else
                    s_filtered[idx] = mean(signal[idx-1:idx])
                end
            end
        else
            s_filtered = signal_tconv(signal, window)
        end

        return s_filtered
    end

    if fprototype === :mmed
        s_filtered = signal
        for idx in d:-1:1
            if t > 0
                if signal[idx] > t * std(signal) + median(signal)
                    s_filtered[idx] = median(signal[idx:idx+1])
                end
            else
                s_filtered[idx] = median(signal[idx:idx+1])
            end
        end
        for idx in (1 + d):(length(signal) - d)
            if t > 0
                if signal[idx] > t * std(signal) + median(signal)
                    s_filtered[idx] = median(signal[(idx - d):(idx + d)])
                end
            else
                s_filtered[idx] = median(signal[(idx - d):(idx + d)])
            end
        end
        for idx in (length(signal) - d + 1):length(signal)
            if t > 0
                if signal[idx] > t * std(signal) + median(signal)
                    s_filtered[idx] = median(signal[idx-1:idx])
                end
            else
                s_filtered[idx] = median(signal[idx-1:idx])
            end
        end

        return s_filtered
    end

    if fprototype === :poly
        t = collect(0:1/fs:(length(signal) - 1) / fs)        
        p = Polynomials.fit(t, signal, order)
        s_filtered = zeros(length(signal))
        for idx in 1:length(signal)
            s_filtered[idx] = p(t[idx])
        end

        return s_filtered
    end

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

    dir === :twopass && (s_filtered = filtfilt(eeg_filter, signal))
    dir === :onepass && (s_filtered = filt(eeg_filter, signal))
    dir === :onepass_reverse && (s_filtered = filt(eeg_filter, reverse(signal)))

    return s_filtered
end

"""
    signal_filter(signal; fprototype, ftype, cutoff, fs, order, rp, rs, dir=:twopass, d=1, window)

Filters `signal` using zero phase distortion filter.

# Arguments

- `signal::Array{Float64, 3}`
- `fprototype::Symbol[:mavg, :mmed, :poly, :butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir]` - filter prototype:
    - `:mavg` - moving average (with threshold and/or weight window)
    - `:mmed` - moving median (with threshold and/or weight window)
    - `:poly` - polynomial of `order` order
- `ftype::Symbol[:lp, :hp, :bp, :bs]` - filter type
- `cutoff::Union{Int64, Float64, Vector{Int64}, Vector{Float64}, Tuple, Nothing}` - filter cutoff in Hz (vector for `:bp` and `:bs`)
- `fs::Union{Int64, Nothing}` - sampling rate
- `order::Union{Int64, Nothing}` - filter order
- `rp::Union{Float64, Nothing}` - dB ripple in the passband
- `rs::Union{Float64, Nothing}` - dB attentuation in the stopband
- `dir:Symbol[:onepass, :onepass_reverse, :twopass]` - filter direction
- `d::Int64` - window length for mean average and median average filter
- `t::Union{Int64, Float64}` - threshold for :mavg and :mmed filters; threshold = threshold * std(signal) + mean(signal) for :mavg or threshold = threshold * std(signal) + median(signal) for :mmed filter
- `window::Union{Vector{Float64}, Nothing} - window, required for :fir and :mavg filters

# Returns

- `s_filtered::Array{Float64, 3}`
"""
function signal_filter(signal::Array{Float64, 3}; fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Int64, Float64, Vector{Int64}, Vector{Float64}, Tuple, Nothing}=nothing, fs::Union{Int64, Nothing}=nothing, order::Union{Int64, Nothing}=nothing, rp::Union{Int64, Float64, Nothing}=nothing, rs::Union{Int64, Float64, Nothing}=nothing, dir::Symbol=:twopass, d::Int64=1, t::Union{Int64, Float64}=0, window::Union{Vector{Float64}, Nothing}=nothing)
    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)
    s_filtered = zeros(size(signal))

    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
            s_filtered[idx, :, epoch] = signal_filter(signal[idx, :, epoch],
                                                           fprototype=fprototype,
                                                           ftype=ftype,
                                                           cutoff=cutoff,
                                                           fs=fs,
                                                           order=order,
                                                           rp=rp,
                                                           rs=rs,
                                                           dir=dir,
                                                           d=d,
                                                           t=t,
                                                           window=window)
        end
    end

    return s_filtered
end

"""
    signal_downsample(signal; t, new_sr)

Downsamples the`signal` to `new_sr` sampling frequency.

# Arguments

- `signal::Vector{Float64}`
- `t::AbstractRange`
- `new_sr::Int64` - new sampling rate

# Returns

- `s_downsampled::Vector{Float64}`
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
    s_downsampled = resample(signal, sr_ratio)
    t = collect(t)
    t_downsampled = t[1]:1/new_sr:t[end]

    return s_downsampled, t_downsampled
end

"""
    signal_downsample(signal; t, new_sr)

Downsamples all channels of the`signal` to `new_sr` sampling frequency.

# Arguments

- `signal::Array{Float64, 3}`
- `new_sr::Int64` - new sampling rate
- `t::AbstractRange`

# Returns

- `s_downsampled::Array{Float64, 3}`
- `t_downsampled::AbstractRange`
"""
function signal_downsample(signal::Array{Float64, 3}; t::AbstractRange, new_sr::Int64)
    new_sr < 1 && throw(ArgumentError("New sampling rate must be positive."))
    
    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)

    s_downsampled_len = length(signal_downsample(signal[1, :, 1], t=t, new_sr=new_sr)[1])
    s_downsampled = zeros(channels_n, s_downsampled_len, epochs_n) 

    t_downsampled = nothing
    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
            s_downsampled[idx, :, epoch], t_downsampled = signal_downsample(signal[idx, :, epoch], t=t, new_sr=new_sr)
        end
    end

    return s_downsampled, t_downsampled
end

"""
    signal_psd(signal; fs, norm=false)

Calculates power spectrum density of the `signal`.

# Arguments
- `signal::Vector{Float64}`
- `fs::Int64` - sampling rate
- `norm::Bool` - normalize do dB

# Returns

- `psd_pow::Vector{Float64}`
- `psd_frq::Vector{Float64}`
"""
function signal_psd(signal::Vector{Float64}; fs::Int64, norm::Bool=false)
    fs < 1 && throw(ArgumentError("Sampling rate must be positive."))
    
    psd = welch_pgram(signal, 4*fs, fs=fs)
    psd_pow = power(psd)
    psd_frq = freq(psd)
    norm == true && (psd_pow = pow2db.(psd_pow))

    return psd_pow, Vector(psd_frq)
end

"""
    signal_psd(signal; fs, norm=false)

Calculates power spectrum density for each the `signal` channels.

# Arguments

- `signal::Matrix{Float64}`
- `fs::Int64` - sampling rate
- `norm::Bool` - normalize do dB

# Returns

- `psd_pow::Matrix{Float64}`
- `psd_frq::Matrix{Float64}`
"""
function signal_psd(signal::Matrix{Float64}; fs::Int64, norm::Bool=false)
    fs < 1 && throw(ArgumentError("Sampling rate must be positive."))

    channels_n = size(signal, 1)
    psd_len, _ = signal_psd(signal[1, :], fs=fs, norm=norm)
    psd_pow = zeros(channels_n, length(psd_len))
    psd_frq = zeros(channels_n, length(psd_len))
    Threads.@threads for idx in 1:channels_n
        psd_pow[idx, :], psd_frq[idx, :] = signal_psd(signal[idx, :],
                                                      fs=fs,
                                                      norm=norm)
    end

    return psd_pow, psd_frq
end

"""
    signal_psd(signal; fs, norm=false)

Calculates power spectrum density for each the `signal` channels.

# Arguments

- `signal::Array{Float64, 3}`
- `fs::Int64` sampling rate
- `norm::Bool` - normalize do dB

# Returns

- `psd_pow::Array{Float64, 3}`
- `psd_frq::Array{Float64, 3}`
"""
function signal_psd(signal::Array{Float64, 3}; fs::Int64, norm::Bool=false)
    fs < 1 && throw(ArgumentError("Sampling rate must be positive."))

    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)
    psd_len, _ = signal_psd(signal[1, :, 1], fs=fs, norm=norm)
    psd_pow = zeros(channels_n, length(psd_len), epochs_n)
    psd_frq = zeros(channels_n, length(psd_len), epochs_n)
    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
            psd_pow[idx, :, epoch], psd_frq[idx, :, epoch] = signal_psd(signal[idx, :, epoch],
                                                                        fs=fs,
                                                                        norm=norm)
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
    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)
    phase_stationarity = zeros(channels_n, size(signal, 2) - 1, epochs_n)
    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
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
    
    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)

    if method === :mean
        stationarity = zeros(channels_n, window, epochs_n)
        Threads.@threads for epoch in 1:epochs_n
            for idx = 1:channels_n
                stationarity[idx, :, epoch] = signal_stationarity_mean(signal[idx, :, epoch], window=window)
            end
        end
    end

    if method === :var
        stationarity = zeros(channels_n, window, epochs_n)
        Threads.@threads for epoch in 1:epochs_n
            for idx = 1:channels_n
                stationarity[idx, :, epoch] = signal_stationarity_var(signal[idx, :, epoch], window=window)
            end
        end
    end

    if method === :hilbert
        stationarity = zeros(channels_n, size(signal, 2) - 1, epochs_n)
        Threads.@threads for epoch in 1:epochs_n
            for idx = 1:channels_n
                stationarity[idx, :, epoch] = signal_stationarity_hilbert(signal[idx, :, epoch])
            end
        end
    end

    if method === :euclid
        # number of time windows per epoch
        window_n = size(signal, 2)
        cov_mat = zeros(channels_n, channels_n, window_n, epochs_n)
        stationarity = zeros(1 + length(2:window:window_n), epochs_n)

        Threads.@threads for epoch in 1:epochs_n
            for idx = 1:window_n
                cov_mat[:, :, idx, epoch] = signal_cov(signal[:, idx, epoch], signal[:, idx, epoch])
            end
        end

        Threads.@threads for epoch in 1:epochs_n
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
    signal_trim(signal; trim_len::Int64)

Removes `trim_len` samples from the beginning (`from` = :start, default) or end (`from` = :end) of the `signal`.

# Arguments

- `signal::Vector{Float64}; trim_len::Int64`
- `offset::Int64` - offset from which trimming starts, only works for `from` = :start
- `from::Symbol[:start, :end]

# Returns

- `s_trimmed::Vector{Float64}`

"""
function signal_trim(signal::Vector{Float64}; trim_len::Int64, offset::Int64=0, from::Symbol=:start)
    from in [:start, :end] || throw(ArgumentError("Argument from must be :start or :end."))
    trim_len < 0 && throw(ArgumentError("Trim length must be ≥ 1."))
    trim_len >= length(signal) && throw(ArgumentError("Trim length must be less than signal length."))
    offset < 0 && throw(ArgumentError("Offset must be ≥ 1."))
    offset >= length(signal) - 1 && throw(ArgumentError("Offset must be less than signal length."))
    (from ===:start && 1 + offset + trim_len > length(signal)) && throw(ArgumentError("Offset + trim length must be less than signal length."))
    
    from === :start && (s_trimmed = vcat(signal[1:offset], signal[(1 + offset + trim_len):end]))
    from === :end && (s_trimmed = signal[1:(end - trim_len)])
    
    return s_trimmed::Vector{Float64}
end


"""
    signal_trim(signal; trim_len, offset=0, from=:start)

Removes `trim_len` samples from the beginning (`from` = :start, default) or end (`from` = :end) of the `signal`.

# Arguments

- `signal::Array{Float64, 3}`
- `trim_len::Int64` - number of samples to remove
- `offset::Int64` - offset from which trimming starts, only works for `from` = :start
- `from::Symbol[:start, :end]`

# Returns

- `s_trimmed::Array{Float64, 3}`

"""
function signal_trim(signal::Array{Float64, 3}; trim_len::Int64, offset::Int64=0, from::Symbol=:start)
    from in [:start, :end] || throw(ArgumentError("Argument from must be :start or :end."))
    trim_len < 0 && throw(ArgumentError("Trim length must be ≥ 1."))
    trim_len >= size(signal, 2) && throw(ArgumentError("Trim length must be less than signal length."))
    offset < 0 && throw(ArgumentError("Offset must be ≥ 1."))
    offset >= size(signal, 2) - 1 && throw(ArgumentError("Offset must be less than signal length."))
    (from ===:start && 1 + offset + trim_len > size(signal, 2)) && throw(ArgumentError("Offset + trim length must be less than signal length."))
    
    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)

    s_trimmed = zeros(channels_n, (size(signal, 2) - trim_len), epochs_n)

    Threads.@threads for epoch in 1:epochs_n
        for idx = 1:channels_n
            s_trimmed[idx, :, epoch] = signal_trim(signal[idx, :, epoch], trim_len=trim_len, from=from)
        end
    end
    
    return s_trimmed
end

"""
    signal_mi(signal1::Vector{Float64}, signal2::Vector{Float64})

Calculates mutual information between `signal1` and `signal2`.

# Arguments

- `signal1::Vector{Float64}`
- `signal2::Vector{Float64}`

# Returns

- `mi::Float64`

"""
function signal_mi(signal1::Vector{Float64}, signal2::Vector{Float64})
    mi = get_mutual_information(signal1, signal2)
    
    return mi
end

"""
    signal_mi(signal)

Calculates mutual information between each the `signal` channels.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `mi::Array{Float64, 3}`
"""
function signal_mi(signal::Array{Float64, 3})
    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)
    mi = zeros(channels_n, channels_n, epochs_n)

    Threads.@threads for epoch in 1:epochs_n
        for idx1 in 1:channels_n
            for idx2 in 1:channels_n
                mi[idx1, idx2, epoch] = signal_mi(signal[idx1, :, epoch], signal[idx2, :, epoch])
            end
        end
    end

    return mi
end

"""
    signal_mi(signal1, signal2)

Calculates mutual information between each the `signal1` and `signal2` channels.

# Arguments

- `signal1::Array{Float64, 3}`
- `signal2::Array{Float64, 3}`

# Returns

- `mi::Array{Float64, 3}`
"""
function signal_mi(signal1::Array{Float64, 3}, signal2::Array{Float64, 3})
    size(signal1, 1) == size(signal2, 1) || throw(ArgumentError("Both signals must have the same number of channels."))
    size(signal1, 3) == size(signal2, 3) || throw(ArgumentError("Both signals must have the same number of epochs."))
    channels_n = size(signal1, 1)
    epochs_n = size(signal1, 3)
    mi = zeros(channels_n, channels_n, epochs_n)

    Threads.@threads for epoch in 1:epochs_n
        for idx1 in 1:channels_n
            for idx2 in 1:channels_n
                mi[idx1, idx2, epoch] = signal_mi(signal1[idx1, :, epoch], signal2[idx2, :, epoch])
            end
        end
    end

    return mi
end

"""
    signal_entropy(signal::Vector{Float64})

Calculates entropy of `signal`.

# Arguments

- `signal::Vector{Float64}`

# Returns

- `ent::Float64`

"""
function signal_entropy(signal::Vector{Float64})
    n = length(signal)
    maxmin_range = maximum(signal) - minimum(signal)
    fd_bins = ceil(Int64, maxmin_range/(2.0 * iqr(signal) * n^(-1/3))) # Freedman-Diaconis

    # recompute entropy with optimal bins for comparison
    h = StatsKit.fit(Histogram, signal, nbins=fd_bins)
    hdat1 = h.weights ./ sum(h.weights)

    # convert histograms to probability values
    ent = -sum(hdat1 .* log2.(hdat1 .+ eps()))    

    return ent
end


"""
    signal_entropy(signal)

Calculates mutual information between each the `signal1` and `signal2` channels.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `s_entropy::Matrix{Float64}`
"""
function signal_entropy(signal::Array{Float64, 3})
    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)
    s_entropy = zeros(channels_n, epochs_n)

    Threads.@threads for epoch in 1:epochs_n
        for idx1 in 1:channels_n
                s_entropy[idx1, epoch] = signal_entropy(signal[idx1, :, epoch])
        end
    end

    return s_entropy
end

"""
    signal_average(signal1, signal2)

Averages `signal1` and `signal2`.

# Arguments

- `signal1::Vector{Float64}`
- `signal2::Vector{Float64}`

# Returns

- `s_averaged::Vector{Float64}`

"""
function signal_average(signal1::Vector{Float64}, signal2::Vector{Float64})
    length(signal1) == length(signal2) || throw(ArgumentError("Both signals must have the same length."))

    s_averaged = zeros(length(signal1))
    for idx in 1:length(signal1)
        s_averaged[idx] = mean([signal1[idx], signal2[idx]])
    end

    return s_averaged
end

"""
    signal_average(signal1, signal2)

Averages `signal1` and `signal2`.

# Arguments

- `signal1::Array{Float64, 3}`
- `signal2::Array{Float64, 3}`

# Returns

- `s_averaged::Array{Float64, 3}`
"""
function signal_average(signal1::Array{Float64, 3}, signal2::Array{Float64, 3})
    size(signal1) == size(signal2) || throw(ArgumentError("Both signals must have the same size."))
    channels_n = size(signal1, 1)
    epochs_n = size(signal1, 3)
    s_averaged = zeros(size(signal1))

    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
            s_averaged[idx, :, epoch] = signal_average(signal1[idx, :, epoch], signal2[idx, :, epoch])
        end
    end

    return s_averaged
end

"""
    signal_coherence(signal1, signal2)

Calculates coherence between `signal1` and `signal2`.

# Arguments

- `signal1::Vector{Float64}`
- `signal2::Vector{Float64}`

# Returns

- `coherence::Vector{ComplexF64`

"""
function signal_coherence(signal1::Vector{Float64}, signal2::Vector{Float64})
    length(signal1) == length(signal2) || throw(ArgumentError("Both signals must have the same length."))

    s1_fft = fft(signal1) ./ length(signal1)
    s2_fft = fft(signal2) ./ length(signal2)

    coherence = (abs.((s1_fft) .* conj.(s2_fft)).^2) ./ (s1_fft .* s2_fft)

    return coherence
end

"""
    signal_coherence(signal1, signal2)

Calculates coherence between `signal1` and `signal2`.

# Arguments

- `signal1::Array{Float64, 3}`
- `signal2::Array{Float64, 3}`

# Returns

- `coherence::Array{ComplexF64, 3}`
"""
function signal_coherence(signal1::Array{Float64, 3}, signal2::Array{Float64, 3})
    size(signal1) == size(signal2) || throw(ArgumentError("Both signals must have the same size."))
    channels_n = size(signal1, 1)
    epochs_n = size(signal1, 3)
    coherence = zeros(ComplexF64, size(signal1))

    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
            coherence[idx, :, epoch] = signal_coherence(signal1[idx, :, epoch], signal2[idx, :, epoch])
        end
    end

    return coherence
end

"""
    signal_pca(signal, n)

Calculates `n` first PCs for `signal`.

# Arguments

- `signal::Array{Float64, 3}`
- `n::Int64` - number of PCs

# Returns

- `pc::Array{Float64, 3}:` - PC(1)..PC(n) × epoch
- `pc_var::Matrix{Float64}` - PC_VAR(1)..PC_VAR(n) × epoch
"""
function signal_pca(signal::Array{Float64, 3}; n::Int64)
    n < 0 && throw(ArgumentError("Number of PCs must be ≥ 1."))
    n > size(signal, 1) && throw(ArgumentError("Number of PCs cannot be higher than signal rows."))

    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)
    pc = zeros(n, size(signal, 2), epochs_n)
    pc_var = zeros(n, epochs_n)

    Threads.@threads for epoch in 1:epochs_n
        m_cov = signal_cov(signal[:, :, epoch])
        eig_val, eig_vec = eigen(m_cov)
        eig_val_idx = sortperm(eig_val, rev=true)
        eig_val = eig_val[eig_val_idx]
        eig_vec = matrix_sort(eig_vec, eig_val_idx)
        eig_val = 100 .* eig_val / sum(eig_val) # convert to %

        for idx in 1:n
            pc_var[idx, epoch] = eig_val[idx]
            pc[idx, :, epoch] = (eig_vec[:, idx] .* signal[:, :, epoch])[idx, :]
        end
    end

    return pc, pc_var
end

"""
    signal_fconv(signal, kernel)

Performs convolution in the frequency domain between `signal` and `kernel`.

# Arguments

- `signal::Vector{Float64}`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `s_conv::Vector{ComplexF64}`
"""
function signal_fconv(signal::Vector{Float64}, kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})
    n_signal = length(signal)
    n_kernel = length(kernel)
    n_conv = n_signal + n_kernel - 1
    half_kernel = floor(Int64, n_kernel / 2)
    s_fft = fft0(signal, n_conv)
    kernel_fft = fft0(kernel, n_conv)
    s_conv = ifft(s_fft .* kernel_fft)
    
    # remove in- and out- edges
    if mod(n_kernel, 2) == 0 
        s_conv = s_conv[half_kernel:(end - half_kernel)]
    else
        s_conv = s_conv[half_kernel:(end - half_kernel - 1)]
    end

    return s_conv
end

"""
    signal_fconv(signal, kernel)

Performs convolution in the frequency domain between `signal` and `kernel`.

# Arguments

- `signal::Array{Float64, 3}`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `s_conv::Array{ComplexF64, 3}`
"""
function signal_fconv(signal::Array{Float64, 3}, kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})
    channels_n = size(signal, 1)
    epochs_n = size(signal, 3)

    s_conv = zeros(ComplexF64, size(signal))

    Threads.@threads for epoch in 1:epochs_n
        for idx in 1:channels_n
            s_conv[idx, :, epoch] = signal_fconv(signal[idx, :, epoch], kernel)
        end
    end

    return s_conv
end