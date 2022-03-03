"""
    signal_derivative(signal)

Returns the derivative of the `signal` with length same as the signal.

# Arguments

- `signal::AbstractArray`

# Returns

- `s_derivative::Union{Vector{Int64}, Vector{Float64}}`
"""
function signal_derivative(signal::AbstractArray)

    s_der = diff(signal)
    s_der = vcat(s_der, s_der[end])
    
    return s_der
end

"""
    signal_derivative(signal)

Returns the derivative of each the `signal` channels with length same as the signal.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `s_derivative::Array{Float64, 3}`
"""
function signal_derivative(signal::Array{Float64, 3})

    channel_n, _, epoch_n = size(signal)
    s_der = zeros(size(signal))
    
    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            s_der[idx, :, epoch] = signal_derivative(s)
        end
    end

    return s_der
end

"""
    signal_total_power(signal; fs)

Calculates total power for the `signal`.

# Arguments

- `signal::AbstractArray`
- `fs::Int64`: sampling rate

# Returns

- `stp::Float64`
"""
function signal_total_power(signal::AbstractArray; fs::Int64)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))

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
- `fs::Int64`: sampling rate

# Returns

- `stp::Matrix{Float64}`
"""
function signal_total_power(signal::Array{Float64, 3}; fs::Int64)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))

    channel_n, _, epoch_n = size(signal)
    stp = zeros(channel_n, epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            stp[idx, epoch] = signal_total_power(s, fs=fs)
        end
    end

    return stp
end

"""
    signal_band_power(signal; fs, f)

Calculates absolute band power between frequencies `f[1]` and `f[2]` for the `signal`.

# Arguments

- `signal::AbstractArray`
- `fs::Int64`: sampling rate of the signal
- `f::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}`: lower and upper frequency bound

# Returns

- `sbp::Float64`
"""
function signal_band_power(signal::AbstractArray; fs::Int64, f::Tuple)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    length(f) != 2 && throw(ArgumentError("f must contain two frequencies."))
    f = tuple_order(f)
    f[1] <= 0 && throw(ArgumentError("Lower frequency bound must be be > 0."))
    f[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be be < $(fs / 2)."))

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
    signal_band_power(signal; fs, f)

Calculates absolute band power between frequencies `f[1]` and `f[2]` for the `signal`.

# Arguments

- `signal::Array{Float64, 3}`
- `fs::Int64`: sampling rate
- `f::Tuple{Union{Int64, Float64}, Union{Int64, Float64}}`: lower and upper frequency bound

# Returns

- `sbp::Matrix{Float64}`
"""
function signal_band_power(signal::Array{Float64, 3}; fs::Int64, f::Tuple)
    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    length(f) != 2 && throw(ArgumentError("f must contain two frequencies."))
    f = tuple_order(f)
    f[1] <= 0 && throw(ArgumentError("Lower frequency bound must be be > 0."))
    f[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be be < $(fs / 2)."))

    channel_n, _, epoch_n = size(signal)
    sbp = zeros(channel_n, epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            sbp[idx, epoch] = signal_band_power(s, fs=fs, f=f)
        end
    end

    return sbp
end

"""
    signal_make_spectrum(signal; fs)

Returns FFT and DFT sample frequencies for a DFT for the `signal`.

# Arguments

- `signal::AbstractArray`
- `fs::Int64`: sampling rate

# Returns

- `s_fft::Vector{ComplexF64}`
- `s_sf::Vector{Float64}`
"""
function signal_make_spectrum(signal::AbstractArray; fs::Int64)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))

    s_fft = fft(signal)
    # number of samples
    n = length(signal)
    # time between samples
    d = 1 / fs
    s_sf = fftfreq(n, d)

    return s_fft, Vector(s_sf)
end

"""
    signal_make_spectrum(signal; fs)

Returns FFT and DFT sample frequencies for a DFT for each the `signal` channels.

# Arguments

- `signal::Array{Float64, 3}`
- `fs::Int64`: sampling rate

# Returns

- `s_fft::Array{ComplexF64, 3}`
- `s_sf::Array{Float64, 3}`
"""
function signal_make_spectrum(signal::Array{Float64, 3}; fs::Int64)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))

    channel_n, _, epoch_n = size(signal)
    s_fft = zeros(ComplexF64, size(signal))
    s_sf = zeros(size(signal))

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            s_fft[idx, :, epoch], s_sf[idx, :, epoch] = signal_make_spectrum(s, fs=fs)
        end
    end

    return s_fft, s_sf
end

"""
    signal_detrend(signal; type=:linear)

Removes linear trend from the `signal`.

# Arguments

- `signal::AbstractArray`
- `type::Symbol[:linear, :constant]`, optional
    - `linear`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `constant`: the mean of `signal` is subtracted

# Returns

- `s_detrended::Vector{Float64}`
"""
function signal_detrend(signal::AbstractArray; type::Symbol=:linear)

    type in [:linear, :constant] || throw(ArgumentError("type must be :linear or :constant."))

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
    - `linear`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `constant`: the mean of `signal` is subtracted

# Returns

- `s_detrended::Array{Float64, 3}`
"""
function signal_detrend(signal::Array{Float64, 3}; type::Symbol=:linear)

    type in [:linear, :constant] || throw(ArgumentError("type must be :linear or :constant."))

    channel_n, _, epoch_n = size(signal)
    s_det = zeros(size(signal))
    
    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            s_det[idx, :, epoch] = signal_detrend(s, type=type)
        end
    end

    return s_det
end

"""
    signal_ci95(signal; n=3, method=:normal)

Calculates mean, std and 95% confidence interval for `signal`.

# Arguments

- `signal::Vector{Float64}`
- `n::Int64`: number of bootstraps
- `method::Symbol[:normal, :boot]`: use normal method or `n`-times boostrapping

# Returns

- `s_m::Float64`
- `s_s::Float64`
- `s_u::Float64`
- `s_l::Float64`
"""
function signal_ci95(signal::Vector{Float64}; n::Int64=3, method::Symbol=:normal)

    method === :normal || throw(ArgumentError("method must be :normal."))
    n < 1 && throw(ArgumentError("n must be ≥ 1."))

    s_m = mean(signal)
    s_s = std(signal) / sqrt(length(signal))
    s_u = s_m + 1.96 * s_s
    s_l = s_m - 1.96 * s_s

    return s_m, s_s, s_u, s_l
end

"""
    signal_ci95(signal; n=3, method=:normal)

Calculates mean, std and 95% confidence interval for each the `signal` channels.

# Arguments

- `signal::AbstractArray`
- `n::Int64`: number of bootstraps
- `method::Symbol[:normal, :boot]`: use normal method or `n`-times boostrapping

# Returns

- `s_m::Vector{Float64}`
- `s_s::Vector{Float64}`
- `s_u::Vector{Float64}`
- `s_l::Vector{Float64}`
"""
function signal_ci95(signal::AbstractArray; n::Int64=3, method::Symbol=:normal)

    method in [:normal, :boot] || throw(ArgumentError("method must be :normal or :boot."))

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
    signal_ci95(signal; n::=3, method=:normal)

Calculates mean, std and 95% confidence interval for each the `signal` channels.

# Arguments

- `signal::Array{Float64, 3}`
- `n::Int64`: number of bootstraps
- `method::Symbol[:normal, :boot]`: use normal method or `n`-times boostrapping

# Returns

- `s_m::Matrix{Float64}`
- `s_s::Matrix{Float64}`
- `s_u::Matrix{Float64}`
- `s_l::Matrix{Float64}`
"""
function signal_ci95(signal::Array{Float64, 3}; n::Int64=3, method::Symbol=:normal)

    method in [:normal, :boot] || throw(ArgumentError("method must be :normal or :boot."))
    n < 1 && throw(ArgumentError("n must be ≥ 1."))

    s_m = zeros(size(signal, 3), size(signal, 2))
    s_s = zeros(size(signal, 3), size(signal, 2))
    s_u = zeros(size(signal, 3), size(signal, 2))
    s_l = zeros(size(signal, 3), size(signal, 2))

    _, _, epoch_n = size(signal)

    Threads.@threads for epoch in 1:epoch_n
        s = @view signal[:, :, epoch]
        s_m[epoch, :], s_s[epoch, :], s_u[epoch, :], s_l[epoch, :] = signal_ci95(s)
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

    _, s_len, epoch_n = size(signal1)

    s_m = zeros(epoch_n, s_len)
    s_s = zeros(epoch_n, s_len)
    s_u = zeros(epoch_n, s_len)
    s_l = zeros(epoch_n, s_len)

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

    return s_m, s_s, s_u, s_l
end

"""
    signal_difference(signal1, signal2; n=3, method=:absdiff)

Calculates mean difference and 95% confidence interval for 2 signals.

# Arguments

- `signal1::AbstractArray`
- `signal2::AbstractArray`
- `n::Int64`: number of bootstraps
- `method::Symbol[:absdiff, :diff2int]`
    - `:absdiff`: maximum difference
    - `:diff2int`: integrated area of the squared difference

# Returns

- `s_statistic::Vector{Float64}`
- `s_statistic_single::Float64`
- `p::Float64`
"""
function signal_difference(signal1::AbstractArray, signal2::AbstractArray; n::Int64=3, method::Symbol=:absdiff)

    size(signal1) != size(signal2) && throw(ArgumentError("Both signals must be of the same size."))
    method in [:absdiff, :diff2int] || throw(ArgumentError("method must be :absdiff or :diff2int."))

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
            s = @view signals[sample_idx[idx2], :]
            s_tmp1[idx2, :] = s'
        end
        s1_mean = mean(s_tmp1, dims=1)
        s_tmp1 = zeros(size(signal1, 1), size(signal1, 2))
        sample_idx = rand(1:size(signals, 1), size(signals, 1))
        for idx2 in 1:size(signal1, 1)
            s = @view signals[sample_idx[idx2], :]
            s_tmp1[idx2, :] = s'
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
    signal_difference(signal1, signal2; n=3, method=:absdiff)

Calculates mean difference and 95% confidence interval for 2 signals.

# Arguments

- `signal1::Array{Float64, 3}`
- `signal2:Array{Float64, 3}`
- `n::Int64`: number of bootstraps
- `method::Symbol[:absdiff, :diff2int]`
    - `:absdiff`: maximum difference
    - `:diff2int`: integrated area of the squared difference

# Returns

- `s_statistic::Matrix{Float64}`
- `s_statistic_single::Vector{Float64}`
- `p::Vector{Float64}`
"""
function signal_difference(signal1::Array{Float64, 3}, signal2::Array{Float64, 3}; n::Int64=3, method::Symbol=:absdiff)

    size(signal1) != size(signal2) && throw(ArgumentError("Both signals must be of the same size."))
    method in [:absdiff, :diff2int] || throw(ArgumentError("method must be :absdiff or :diff2int."))

    _, _, epoch_n = size(signal1)
    s_statistic = zeros(epoch_n, size(signal1, 1) * n)
    s_statistic_single = zeros(epoch_n)
    p = zeros(epoch_n)

    Threads.@threads for epoch in 1:epoch_n
        s1 = @view signal1[:, :, epoch]
        s2 = @view signal2[:, :, epoch]
        s_statistic[epoch, :], s_statistic_single[epoch], p[epoch] = signal_difference(s1, s2)
    end

    return s_statistic, s_statistic_single, p
end

"""
   signal_autocov(signal; lag=1, demean=false, norm=false)

Calculates autocovariance of the `signal`.

# Arguments

- `signal::AbstractArray`
- `lag::Int64`: lags range is `-lag:lag`
- `demean::Bool`: demean `signal` prior to calculations
- `norm::Bool`: normalize autocovariance

# Returns

- `acov::Vector{Float64}`
- `lags::Vector{Int64}`
"""
function signal_autocov(signal::AbstractArray; lag::Int64=1, demean::Bool=false, norm::Bool=false)

    lag < 1 && throw(ArgumentError("lag must be ≥ 1."))

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
- `lag::Int64`: lags range is `-lag:lag`
- `demean::Bool`: demean signal prior to analysis
- `norm::Bool`: normalize autocovariance

# Returns

- `acov::Matrix{Float64}`
- `lags::Vector{Int64}`
"""
function signal_autocov(signal::Array{Float64, 3}; lag::Int64=1, demean::Bool=false, norm::Bool=false)

    lag < 1 && throw(ArgumentError("lag must be ≥ 1."))

    lags = collect(-lag:lag)
    channel_n, _, epoch_n = size(signal)
    acov = zeros(channel_n, length(lags), epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            acov[idx, :, epoch], lags = signal_autocov(s,
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

- `signal1::AbstractArray`
- `signal2::AbstractArray`
- `lag::Int64`: lags range is `-lag:lag`
- `demean::Bool`: demean signal prior to analysis
- `norm::Bool`: normalize cross-covariance

# Returns

- `ccov::Vector{Float64}`
- `lags::Vector{Int64}`
"""
function signal_crosscov(signal1::AbstractArray, signal2::AbstractArray; lag::Int64=1, demean::Bool=false, norm::Bool=false)

    length(signal1) != length(signal2) && throw(ArgumentError("Both signals must be of the same as length."))
    lag < 1 && throw(ArgumentError("lag must be ≥ 1."))

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

- `signal::Array{Float64, 3}`: the signal
- `lag::Int64`: lags range is `-lag:lag`
- `demean::Bool`: demean `signal` prior to analysis
- `norm::Bool`: normalize cross-covariance

# Returns

- `ccov::Array{Float64, 3}`
- `lags::Vector{Int64}`
"""
function signal_crosscov(signal::Array{Float64, 3}; lag::Int64=1, demean::Bool=false, norm::Bool=false)

    lag < 1 && throw(ArgumentError("lag must be ≥ 1 Hz."))

    lags = collect(-lag:lag)
    channel_n, _, epoch_n = size(signal)
    ccov = zeros(channel_n^2, length(lags), epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        ccov_packed = Array{Vector{Float64}}(undef, channel_n, channel_n)
        Threads.@threads for idx1 in 1:channel_n
            for idx2 in 1:channel_n
                s1 = @view signal[idx1, :, epoch]
                s2 = @view signal[idx2, :, epoch]
                ccov_packed[idx1, idx2], lags = signal_crosscov(s1,
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

    return ccov, lags
end

"""
   signal_crosscov(signal1, signal2; lag=1, demean=false, norm=false)

Calculates cross-covariance between same channels in `signal1` and `signal2`.

# Arguments

- `signal1::Array{Float64, 3}`
- `signal2::Array{Float64, 3}`
- `lag::Int64`: lags range is `-lag:lag`
- `demean::Bool`: demean signal prior to analysis
- `norm::Bool`: normalize cross-covariance

# Returns

- `ccov::Array{Float64, 3}`
- `lags::Vector{Int64}`
"""
function signal_crosscov(signal1::Array{Float64, 3}, signal2::Array{Float64, 3}; lag::Int64=1, demean::Bool=false, norm::Bool=false)

    size(signal1) != size(signal2) && throw(ArgumentError("Both arrays must be of the same as size."))
    lag < 1 && throw(ArgumentError("lag must be ≥ 1 Hz."))

    lags = collect(-lag:lag)
    channel_n, _, epoch_n = size(signal1)
    ccov = zeros(channel_n, length(lags), epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s1 = @view signal1[idx, :, epoch]
            s2 = @view signal2[idx, :, epoch]
            ccov[idx, :, epoch], lags = signal_crosscov(s1,
                                                        s2,
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

- `signal::AbstractArray`
- `pad::Int64`: pad the `signal` with `pad` zeros

# Returns

- `fft::Vector(ComplexF64}`
- `amplitudes::Vector{Float64}`
- `powers::Vector{Float64}`
- `phases::Vector{Float64}
"""
function signal_spectrum(signal::AbstractArray; pad::Int64=0)

    pad < 0 && throw(ArgumentError("pad must be ≥ 0."))

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

- `signal::Array{Float64, 3}`: the signal
- `pad::Int64`: pad the `signal` with `pad` zeros

# Returns

- `fft::Array{ComplexF64, 3}`
- `amplitudes::Array{Float64, 3}`
- `powers::Array{Float64, 3}`
- `phases::Array{Float64, 3}
"""
function signal_spectrum(signal::Array{Float64, 3}; pad::Int64=0)

    pad < 0 && throw(ArgumentError("pad must be ≥ 0."))

    channel_n, _, epoch_n = size(signal)

    s_fft = zeros(ComplexF64, channel_n, size(signal, 2) + pad, epoch_n)
    s_amplitudes = zeros(channel_n, size(signal, 2) + pad, epoch_n)
    s_powers = zeros(channel_n, size(signal, 2) + pad, epoch_n)
    s_phases = zeros(channel_n, size(signal, 2) + pad, epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            s_fft[idx, :, epoch], s_amplitudes[idx, :, epoch], s_powers[idx, :, epoch], s_phases[idx, :, epoch] = signal_spectrum(s, pad=pad)
        end
    end

    return s_fft, s_amplitudes, s_powers, s_phases
end

"""
    signal_epochs(signal; epoch_n, epoch_len, average=true)

Splits `signal` into epochs.

# Arguments

- `signal::Vector{Float64}`
- `epoch_n::Union{Int64, Nothing}`: number of epochs
- `epoch_len::Union{Int64, Nothing}`: epoch length in samples
- `average::Bool`: average all epochs, returns one averaged epoch; if false than returns array of epochs, each row is one epoch

# Returns

- `epochs::Matrix{Float64}`
"""
function signal_epochs(signal::Vector{Float64}; epoch_n::Union{Int64, Nothing}=nothing, epoch_len::Union{Int64, Nothing}=nothing, average::Bool=false)

    (epoch_len === nothing && epoch_n === nothing) && throw(ArgumentError("Either epoch_n or epoch_len must be set."))
    (epoch_len !== nothing && epoch_n !== nothing) && throw(ArgumentError("Both epoch_n and epoch_len cannot be set."))
    (epoch_len != nothing && epoch_len < 1) && throw(ArgumentError("epoch_len must be ≥ 1."))
    (epoch_n != nothing && epoch_n < 1) && throw(ArgumentError("epoch_n must be ≥ 1."))

    if epoch_n === nothing
        epoch_n = length(signal) ÷ epoch_len
    end
    if epoch_len === nothing
        epoch_len = length(signal) ÷ epoch_n
    end

    epochs = zeros(epoch_n, epoch_len)

    signal = signal[1:epoch_len * epoch_n]
    epochs = reshape(signal, epoch_n, epoch_len)'

    if average == true
        epochs = vec(mean(epochs, dims=1)[1, :])
    end

    return epochs
end

"""
    signal_epochs(signal; epoch_n=nothing, epoch_len=nothing, average=true)

Splits `signal` into epochs.

# Arguments

- `signal::Array{Float64, 3}`
- `epoch_n::Union{Int64, Nothing}`: number of epochs
- `epoch_len::Union{Int64, Nothing}`: epoch length in samples
- `average::Bool`: average all epochs, returns one averaged epoch; if false than returns array of epochs, each row is one epoch

# Returns

- `epochs::Array{Float64, 3}`
"""
function signal_epochs(signal::Matrix{Float64}; epoch_n::Union{Int64, Nothing}=nothing, epoch_len::Union{Int64, Nothing}=nothing, average::Bool=false)

    (epoch_len === nothing && epoch_n === nothing) && throw(ArgumentError("Either epoch_n or epoch_len must be set."))
    (epoch_len !== nothing && epoch_n !== nothing) && throw(ArgumentError("Both epoch_n and epoch_len cannot be set."))
    (epoch_len != nothing && epoch_len < 1) && throw(ArgumentError("epoch_len must be ≥ 1."))
    (epoch_n != nothing && epoch_n < 1) && throw(ArgumentError("epoch_n must be ≥ 1."))

    channel_n, _ = size(signal)

    if epoch_n === nothing
        epoch_n = size(signal, 2) ÷ epoch_len
    else
        epoch_len = size(signal, 2) ÷ epoch_n
    end

    epochs = zeros(channel_n, epoch_len, epoch_n)

    # signal = signal[1:epoch_len * epoch_n]
    # epochs = reshape(signal, epoch_n, epoch_len)'

    idx1 = 1
    for idx2 in 1:epoch_len:(epoch_n * epoch_len - 1)
        epochs[:, :, idx1] = signal[:, idx2:(idx2 + epoch_len - 1), 1]
        idx1 += 1
    end

    if average == true
        epochs = mean(epochs, dims=3)[:, :]
    end

    return epochs
end

"""
    signal_delete_channel(signal; channel)

Removes `channel` from the `signal`.

# Arguments

- `signal::Matrix{Float64}`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel to be removed, vector of numbers or range

# Returns

- `signal::Matrix{Float64}`
"""
function signal_delete_channel(signal::Matrix{Float64}; channel::Union{Int64, Vector{Int64}, AbstractRange})

    if typeof(channel) <: AbstractRange
        channel = collect(channel)
    end
    for idx in length(channel):-1:1
        channel[idx] > size(signal, 2) && throw(ArgumentError("channel does not match signal channels."))
    end

    length(channel) > 1 && (channel = sort!(channel, rev=true))
    s_new = signal[setdiff(1:end, (channel)), :]

    return s_new
end

"""
    signal_delete_channel(signal; channel)

Removes `channel` from the `signal`.

# Arguments

- `signal::Array{Float64, 3}`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel to be removed, vector of numbers or range

# Returns

- `signal::Matrix{Float64}`
"""
function signal_delete_channel(signal::Array{Float64, 3}; channel::Union{Int64, Vector{Int64}, AbstractRange})

    if typeof(channel) <: AbstractRange
        channel = collect(channel)
    end
    for idx in length(channel):-1:1
        channel[idx] > size(signal, 2) && throw(ArgumentError("channel does not match signal channels."))
    end

    length(channel) > 1 && (channel = sort!(channel, rev=true))
    s_new = signal[setdiff(1:end, (channel)), :, :]

    return s_new
end
"""
    signal_reference_channel(signal, reference)

Re-references channels of the `signal` to specific signal channel.

# Arguments

- `signal::Matrix{Float64}`
- `reference::Union{Int64, Vector{Int64}, AbstractRange}}`: index of channels used as reference; if multiple channels are specified, their average is used as the reference

# Returns

- `s_referenced::Matrix{Float64}`
"""
function signal_reference_channel(signal::Matrix{Float64}; channel::Union{Int64, Vector{Int64}, AbstractRange})

    if typeof(channel) <: AbstractRange
        channel = collect(channel)
    end

    channel_n, _, = size(signal)
    s_ref = zeros(size(signal))

    channel_list = collect(1:channel_n)
    for idx in 1:length(channel)
        if (channel[idx] in channel_list) == false
            throw(ArgumentError("channel does not match signal channels."))
        end
    end

    if length(channel) == 1
        reference_channel = mean(signal[channel, :], dims=2)
    else
        reference_channel = vec(mean(signal[channel, :], dims=1))
    end
    for idx in 1:channel_n
        s_ref[idx, :] = signal[idx, :] .- reference_channel
    end
    length(channel) == 1 && (s_ref[channel, :] = reference_channel)

    return s_ref
end

"""
    signal_reference_channel(signal, channel)

Re-references channels of the `signal` to specific signal channel.

# Arguments

- `signal::Array{Float64, 3}`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}}`: index of channels used as reference; if multiple channels are specified, their average is used as the reference

# Returns

- `s_referenced::Matrix{Float64}`
"""
function signal_reference_channel(signal::Array{Float64, 3}; channel::Union{Int64, Vector{Int64}, AbstractRange})

    if typeof(channel) <: AbstractRange
        channel = collect(channel)
    end

    channel_n, _, epoch_n = size(signal)
    s_ref = zeros(size(signal))

    channel_list = collect(1:channel_n)
    for idx in 1:length(channel)
        if (channel[idx] in channel_list) == false
            throw(ArgumentError("channel does not match signal channels."))
        end
    end

    @inbounds @simd for epoch in 1:epoch_n
        s = @view signal[channel, :, epoch]
        if length(channel) == 1
            reference_channel = mean(s, dims=2)
        else
            reference_channel = vec(mean(s, dims=1))
        end
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            s_ref[idx, :, epoch] = s .- reference_channel
        end
        length(channel) == 1 && (s_ref[channel, :, epoch] = reference_channel)
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

    channel_n, _ = size(signal)
    reference_channel = vec(mean(signal, dims=1))
    s_ref = zeros(size(signal))

    @inbounds @simd for idx in 1:channel_n
        s = @view signal[idx, :]
        s_ref[idx, :] = s .- reference_channel
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

    channel_n, _, epoch_n = size(signal)
    s_ref = zeros(size(signal))

    @inbounds @simd for epoch in 1:epoch_n
        reference_channel = vec(mean(signal[:, :, epoch], dims=1))
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            s_ref[idx, :, epoch] = s .- reference_channel
        end
    end

    return s_ref
end

"""
    signal_taper(signal; taper)

Taper the `signal` with `taper`.

# Arguments

- `signal::AbstractArray`
- `taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `s_tapered::Vector{Union{Float64, ComplexF64}}`
"""
function signal_taper(signal::AbstractArray; taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

    length(taper) == length(signal) || throw(ArgumentError("Taper and signal lengths must be equal."))
    s_tap = signal .* taper

    return s_tap
end

"""
    signal_taper(signal; taper)

Tapers channels of the `signal` with `taper`.

# Arguments

- `signal::Array{Float64, 3}`
- `taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `s_tapered::Union{Array{Float64, 3}, Array{ComplexF64, 3}}`
"""
function signal_taper(signal::Array{Float64, 3}; taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

    length(taper) == size(signal, 2) || throw(ArgumentError("Taper and signal lengths must be equal."))

    channel_n, _, epoch_n = size(signal)

    s_tap = zeros(eltype(taper), size(signal))

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            s_tap[idx, :, epoch] = signal_taper(s, taper=taper)
        end
    end

    return s_tap
end

"""
    signal_demean(signal)

Removes mean value (DC offset) from the `signal`.

# Arguments

- `signal::AbstractArray`

# Returns

- `s_demeaned::Vector{Float64}`
"""
function signal_demean(signal::AbstractArray)

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

    channel_n, _, epoch_n = size(signal)
    s_dem = zeros(size(signal))

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            s_dem[idx, :, epoch] = signal_demean(s)
        end
    end

    return s_dem
end

"""
    signal_normalize_zscore(signal)

Normalize (by z-score) `signal`.

# Arguments

- `signal::AbstractArray`

# Returns

- `s_normalized::Vector{Float64}`
"""
function signal_normalize_zscore(signal::AbstractArray)

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

    channel_n, _, epoch_n = size(signal)
    s_norm = zeros(size(signal))

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            s_norm[idx, :, epoch] = signal_normalize_zscore(s)
        end
    end

    return s_norm
end

"""
    signal_normalize_minmax(signal)

Normalize `signal` in [-1, +1].

# Arguments

- `signal::AbstractArray`

# Returns

- `s_normalized::Vector{Float64}`
"""
function signal_normalize_minmax(signal::AbstractArray)

    s_norm = 2 .* (signal .- minimum(signal)) ./ (maximum(signal) - minimum(signal)) .- 1

    return s_norm
end

"""
    signal_normalize_minmax(signal)

Normalize each the `signal` channel in [-1, +1].

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `s_normalized::Array{Float64, 3}`
"""
function signal_normalize_minmax(signal::Array{Float64, 3})

    channel_n, _, epoch_n = size(signal)
    s_norm = zeros(size(signal))

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            s_norm[idx, :, epoch] = signal_normalize_minmax(s)
        end
    end

    return s_norm
end

"""
   signal_cov(signal1, signal2; norm=true)

Calculates covariance between `signal1` and `signal2`.

# Arguments

- `signal1::AbstractArray`
- `signal2::AbstractArray`
- `norm::Bool`: normalize covariance

# Returns

- `cov_mat::Matrix{Float64}`
"""
function signal_cov(signal1::AbstractArray, signal2::AbstractArray; norm::Bool=false)

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

- `signal::AbstractArray`
- `norm::Bool`: normalize covariance

# Returns

- `cov_mat::Matrix{Float64}`
"""
function signal_cov(signal::AbstractArray; norm::Bool=false)

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
- `norm::Bool`: normalize covariance

# Returns

- `cov_mat::Array{Float64, 3}`
"""
function signal_cov(signal::Array{Float64, 3}; norm::Bool=false)

    _, _, epoch_n = size(signal)
    cov_mat = zeros(size(signal, 1), size(signal, 1), epoch_n)

    Threads.@threads for epoch in 1:epoch_n
        # channels-vs-channels
        s = @view signal[:, :, epoch]
        cov_mat[:, :, epoch] = cov(s')
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

    _, _, epoch_n = size(signal)
    cor_mat = zeros(size(signal, 1), size(signal, 1), epoch_n)

    Threads.@threads for epoch in 1:epoch_n
        # channels-vs-channels
        s = @view signal[:, :, epoch]
        cor_mat[:, :, epoch] = cor(s')
    end

    return cor_mat
end

"""
    signal_add_noise(signal)

Adds random noise to the `signal`.

# Arguments

- `signal::AbstractArray`

# Returns

- `s_noisy::Vector{Float64}`
"""
function signal_add_noise(signal::AbstractArray)

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

    channel_n, _, epoch_n = size(signal)
    s_noise = zeros(size(signal))

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            s_noise[idx, :, epoch] = signal_add_noise(s)
        end
    end

    return s_noise
end

"""
    signal_upsample(signal; t, new_sr)

Upsamples `signal` to `new_sr` sampling frequency.

# Arguments

- `signal::AbstractArray`
- `t::AbstractRange`
- `new_sr::Int64`: new sampling rate
# Returns

- `s_upsampled::Vector{Float64}`
- `t_upsampled::AbstractRange`
"""
function signal_upsample(signal::AbstractArray; t::AbstractRange, new_sr::Int64)

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
- `new_sr::Int64`: new sampling rate

# Returns

- `s_upsampled::Array{Float64, 3}`
- `t_upsampled::AbstractRange`
"""
function signal_upsample(signal::Array{Float64, 3}; t::AbstractRange, new_sr::Int64)

    new_sr < 1 && throw(ArgumentError("New sampling rate must be positive."))

    channel_n, _, epoch_n = size(signal)

    s_upsampled_len = length(signal_upsample(signal[1, :, 1], t=t, new_sr=new_sr)[1])
    s_upsampled = zeros(channel_n, s_upsampled_len, epoch_n) 

    t_upsampled = nothing
    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            s_upsampled[idx, :, epoch], t_upsampled = signal_upsample(s, t=t, new_sr=new_sr)
        end
    end

    return s_upsampled, t_upsampled
end

"""
    signal_tconv(signal; kernel)

Performs convolution in the time domain between `signal` and `kernel`.

# Arguments

- `signal::AbstractArray`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `s_conv::Union{Vector{Float64}, Vector{ComplexF64}}`
"""
function signal_tconv(signal::AbstractArray; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

    signal = Vector(signal)
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
    signal_tconv(signal; kernel)

Performs convolution in the time domain between `signal` and `kernel`.

# Arguments

- `signal::Array{Float64, 3}`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `s_conv::Union{Array{Float64, 3}, Array{ComplexF64, 3}}`
"""
function signal_tconv(signal::Array{Float64, 3}; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

    channel_n, _, epoch_n = size(signal)

    if typeof(kernel) == Vector{ComplexF64}
        s_conv = zeros(ComplexF64, size(signal))
    else
        s_conv = zeros(size(signal))
    end

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            s_conv[idx, :, epoch] = signal_tconv(s, kernel=kernel)
        end
    end

    return s_conv
end

"""
    signal_filter(signal; fprototype, ftype=nothing, cutoff, fs, order, rp, rs, dir=:twopass, d=1, window=nothing)

Filters `signal`.

# Arguments

- `signal::AbstractArray`
- `fprototype::Symbol[:mavg, :mmed, :poly, :butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir]`: filter prototype:
    - `:mavg`: moving average (with threshold and/or weight window)
    - `:mmed`: moving median (with threshold and/or weight window)
    - `:poly`: polynomial of `order` order
- `ftype::Union{Symbol[:lp, :hp, :bp, :bs], Nothing}`: filter type
- `cutoff::Union{Int64, Float64, Tuple}`: filter cutoff in Hz (vector for `:bp` and `:bs`)
- `order::Int64`: filter order
- `rp::Union{Int64, Float64}`: dB ripple in the passband
- `rs::Union{Int64, Float64}`: dB attentuation in the stopband
- `dir:Symbol[:onepass, :onepass_reverse, :twopass]`: filter direction
- `d::Int64`: window length for mean average and median average filter
- `t::Union{Int64, Float64}`: threshold for :mavg and :mmed filters; threshold = threshold * std(signal) + mean(signal) for :mavg or threshold = threshold * std(signal) + median(signal) for :mmed filter
- `window::Union{Vector{Float64}, Nothing} - window, required for FIR filter

# Returns

- `s_filtered::Vector{Float64}`
"""
function signal_filter(signal::AbstractArray; fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Int64, Float64, Tuple}=0, fs::Int64=0, order::Int64=0, rp::Union{Int64, Float64}=-1, rs::Union{Int64, Float64}=-1, dir::Symbol=:twopass, d::Int64=1, t::Union{Int64, Float64}=0, window::Union{Vector{Float64}, Nothing}=nothing)

    fprototype in [:mavg, :mmed, :poly, :butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir] || throw(ArgumentError("fprototype must be :mavg, :mmed,:butterworth, :chebyshev1, :chebyshev2, :elliptic or :fir."))
    (fprototype === :fir && (window === nothing || length(window) > length(signal))) && throw(ArgumentError("For :fir filter window must be shorter than signal."))
    (fprototype !== :mavg && fprototype !== :mmed && fprototype !== :poly) && (ftype in [:lp, :hp, :bp, :bs] || throw(ArgumentError("ftype must be :bp, :hp, :bp or :bs.")))
    (fprototype !== :mavg && fprototype !== :mmed) && (fs < 1 && throw(ArgumentError("fs must be > 0.")))
    dir in [:onepass, :onepass_reverse, :twopass] || throw(ArgumentError("direction must be :onepass, :onepass_reverse or :twopass."))
    (order < 2 && fprototype !== :poly) && (mod(order, 2) != 0 && throw(ArgumentError("order must be even and ≥ 2.")))
    (order < 1 && (fprototype !== :mavg && fprototype !== :mmed)) && throw(ArgumentError("order must be > 0."))
    d > length(signal) && throw(ArgumentError("d must be ≤ signal length."))

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
            s_filtered = signal_tconv(signal, kernel=window)
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
        length(cutoff) != 1 && throw(ArgumentError("For :lp filter one frequency must be given."))
        responsetype = Lowpass(cutoff; fs=fs)
    elseif ftype === :hp
        length(cutoff) != 1 && throw(ArgumentError("For :hp filter one frequency must be given."))
        responsetype = Highpass(cutoff; fs=fs)
    elseif ftype === :bp
        length(cutoff) != 2 && throw(ArgumentError("For :bp filter two frequencies must be given."))
        responsetype = Bandpass(cutoff[1], cutoff[2]; fs=fs)
    elseif ftype === :bs
        length(cutoff) != 2 && throw(ArgumentError("For :bs filter two frequencies must be given."))
        responsetype = Bandstop(cutoff[1], cutoff[2]; fs=fs)
    end

    fprototype === :butterworth && (prototype = Butterworth(order))
    if fprototype === :fir
        window === nothing && throw(ArgumentError("For :fir filter window must be given."))
        if ftype === :hp || ftype === :bp || ftype === :bs
            mod(length(window), 2) == 0 && (window = vcat(window[1:((length(window) ÷ 2) - 1)], window[((length(window) ÷ 2) + 1):end]))
        end
        prototype = FIRWindow(window)
    end
    if fprototype === :chebyshev1
        (rs < 0 || rs > eeg_sr(eeg) / 2) && throw(ArgumentError("For :chebyshev1 filter rs must be > 0 and ≤ $(fs / 2)."))
        prototype = Chebyshev1(order, rs)
    end
    if fprototype === :chebyshev2
        (rp < 0 || rp > fs / 2) && throw(ArgumentError("For :chebyshev2 filter rp must be > 0 and ≤ $(fs / 2)."))
        prototype = Chebyshev2(order, rp)
    end
    if fprototype === :elliptic
        (rs < 0 || rs > fs / 2) && throw(ArgumentError("For :elliptic filter rs must be > 0 and ≤ $(fs / 2)."))
        (rp < 0 || rp > fs / 2) && throw(ArgumentError("For :elliptic filter rp must be > 0 and ≤ $(fs / 2)."))
        prototype = Elliptic(order, rp, rs)
    end

    eeg_filter = digitalfilter(responsetype, prototype)

    dir === :twopass && (s_filtered = filtfilt(eeg_filter, signal))
    dir === :onepass && (s_filtered = filt(eeg_filter, signal))
    dir === :onepass_reverse && (s_filtered = filt(eeg_filter, reverse(signal)))

    return s_filtered
end

"""
    signal_filter(signal; fprototype, ftype=nothing, cutoff, fs, order, rp, rs, dir=:twopass, d=1, window=nothing)

Filters `signal` using zero phase distortion filter.

# Arguments

- `signal::Array{Float64, 3}`
- `fprototype::Symbol[:mavg, :mmed, :poly, :butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir]`: filter prototype:
    - `:mavg`: moving average (with threshold and/or weight window)
    - `:mmed`: moving median (with threshold and/or weight window)
    - `:poly`: polynomial of `order` order
- `ftype::Union{Symbol[:lp, :hp, :bp, :bs], Nothing}`: filter type
- `cutoff::Union{Int64, Float64, Tuple}`: filter cutoff in Hz (vector for `:bp` and `:bs`)
- `order::Int64`: filter order
- `rp::Union{Int64, Float64}`: dB ripple in the passband
- `rs::Union{Int64, Float64}`: dB attentuation in the stopband
- `dir:Symbol[:onepass, :onepass_reverse, :twopass]`: filter direction
- `d::Int64`: window length for mean average and median average filter
- `t::Union{Int64, Float64}`: threshold for :mavg and :mmed filters; threshold = threshold * std(signal) + mean(signal) for :mavg or threshold = threshold * std(signal) + median(signal) for :mmed filter
- `window::Union{Vector{Float64}, Nothing} - window, required for FIR filter

# Returns

- `s_filtered::Array{Float64, 3}`
"""
function signal_filter(signal::Array{Float64, 3}; fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Int64, Float64, Tuple}=0, fs::Int64=0, order::Int64=0, rp::Union{Int64, Float64}=-1, rs::Union{Int64, Float64}=-1, dir::Symbol=:twopass, d::Int64=1, t::Union{Int64, Float64}=0, window::Union{Vector{Float64}, Nothing}=nothing)

    channel_n, _, epoch_n = size(signal)
    s_filtered = zeros(size(signal))

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            s_filtered[idx, :, epoch] = signal_filter(s,
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

- `signal::AbstractArray`
- `t::AbstractRange`
- `new_sr::Int64`: new sampling rate

# Returns

- `s_downsampled::Vector{Float64}`
- `t_downsampled::Vector{Float64}`
"""
function signal_downsample(signal::AbstractArray; t::AbstractRange, new_sr::Int64)

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
- `new_sr::Int64`: new sampling rate
- `t::AbstractRange`

# Returns

- `s_downsampled::Array{Float64, 3}`
- `t_downsampled::AbstractRange`
"""
function signal_downsample(signal::Array{Float64, 3}; t::AbstractRange, new_sr::Int64)

    new_sr < 1 && throw(ArgumentError("New sampling rate must be positive."))
    
    channel_n, _, epoch_n = size(signal)

    s_downsampled_len = length(signal_downsample(signal[1, :, 1], t=t, new_sr=new_sr)[1])
    s_downsampled = zeros(channel_n, s_downsampled_len, epoch_n) 

    t_downsampled = nothing
    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            s_downsampled[idx, :, epoch], t_downsampled = signal_downsample(s, t=t, new_sr=new_sr)
        end
    end

    return s_downsampled, t_downsampled
end

"""
    signal_psd(signal; fs, norm=false)

Calculates power spectrum density of the `signal`.

# Arguments
- `signal::AbstractArray`
- `fs::Int64`: sampling rate
- `norm::Bool`: normalize do dB

# Returns

- `psd_pow::Vector{Float64}`
- `psd_frq::Vector{Float64}`
"""
function signal_psd(signal::AbstractArray; fs::Int64, norm::Bool=false)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    
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
- `fs::Int64`: sampling rate
- `norm::Bool`: normalize do dB

# Returns

- `psd_pow::Matrix{Float64}`
- `psd_frq::Matrix{Float64}`
"""
function signal_psd(signal::Matrix{Float64}; fs::Int64, norm::Bool=false)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))

    channel_n, _, = size(signal)
    psd_len, _ = signal_psd(signal[1, :], fs=fs, norm=norm)
    psd_pow = zeros(channel_n, length(psd_len))
    psd_frq = zeros(channel_n, length(psd_len))
    Threads.@threads for idx in 1:channel_n
        s = @view signal[idx, :]
        psd_pow[idx, :], psd_frq[idx, :] = signal_psd(s,
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
- `norm::Bool`: normalize do dB

# Returns

- `psd_pow::Array{Float64, 3}`
- `psd_frq::Array{Float64, 3}`
"""
function signal_psd(signal::Array{Float64, 3}; fs::Int64, norm::Bool=false)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1"))

    channel_n, _, epoch_n = size(signal)
    psd_len, _ = signal_psd(signal[1, :, 1], fs=fs, norm=norm)
    psd_pow = zeros(channel_n, length(psd_len), epoch_n)
    psd_frq = zeros(channel_n, length(psd_len), epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = signal[idx, :, epoch]
            psd_pow[idx, :, epoch], psd_frq[idx, :, epoch] = signal_psd(s,
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

- `signal::AbstractArray`

# Returns

- `phase_stationarity::Vector{Float64}`
"""
function signal_stationarity_hilbert(signal::AbstractArray)
    
    phase_stationarity = diff(DSP.unwrap(angle.(hilbert(signal))))
    
    return phase_stationarity
end

"""
    signal_stationarity_mean(signal)

Calculates mean stationarity.

# Arguments

- `signal::AbstractArray`
- `window::Int64`: time window in samples

# Returns

- `mean_stationarity::Vector{Float64}`
"""
function signal_stationarity_mean(signal::AbstractArray; window::Int64)

    signal = signal[1:(window * floor(Int64, length(signal) / window))]
    signal = reshape(signal, Int(length(signal) / window), window)
    mean_stationarity = mean(signal, dims=1)

    return mean_stationarity
end

"""
    signal_stationarity_var(signal)

Calculates variance stationarity.

# Arguments

- `signal::AbstractArray`
- `window::Int64`: time window in samples

# Returns

- `var_stationarity::Vector{Float64}`
"""
function signal_stationarity_var(signal::AbstractArray; window::Int64)

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
- `window::Int64`: time window in samples
- `method::Symbol[:mean, :var, :euclid, :hilbert]

# Returns

- `stationarity::Union{Matrix{Float64}, Array{Float64, 3}}`
"""
function signal_stationarity(signal::Array{Float64, 3}; window::Int64=10, method::Symbol=:hilbert)

    method in [:mean, :var, :euclid, :hilbert] || throw(ArgumentError("Method must be must be :mean, :var, :euclid or :hilbert."))
    (typeof(window) == Int64 && window < 1) && throw(ArgumentError("window must be ≥ 1."))
    (typeof(window) == Int64 && window > size(signal, 2)) && throw(ArgumentError("window must be ≤ $(size(signal, 2))."))
    
    channel_n, _, epoch_n = size(signal)

    if method === :mean
        stationarity = zeros(channel_n, window, epoch_n)
        @inbounds @simd for epoch in 1:epoch_n
            Threads.@threads for idx in 1:channel_n
                s = @view signal[idx, :, epoch]
                stationarity[idx, :, epoch] = signal_stationarity_mean(s, window=window)
            end
        end
    end

    if method === :var
        stationarity = zeros(channel_n, window, epoch_n)
        @inbounds @simd for epoch in 1:epoch_n
            Threads.@threads for idx in 1:channel_n
                s = @view signal[idx, :, epoch]
                stationarity[idx, :, epoch] = signal_stationarity_var(s, window=window)
            end
        end
    end

    if method === :hilbert
        stationarity = zeros(channel_n, size(signal, 2) - 1, epoch_n)
        @inbounds @simd for epoch in 1:epoch_n
            Threads.@threads for idx in 1:channel_n
                s = @view signal[idx, :, epoch]
                stationarity[idx, :, epoch] = signal_stationarity_hilbert(s)
            end
        end
    end

    if method === :euclid
        # number of time windows per epoch
        window_n = size(signal, 2)
        cov_mat = zeros(channel_n, channel_n, window_n, epoch_n)
        stationarity = zeros(1 + length(2:window:window_n), epoch_n)

        @inbounds @simd for epoch in 1:epoch_n
            Threads.@threads for idx = 1:window_n
                s = @view signal[:, idx, epoch]
                cov_mat[:, :, idx, epoch] = signal_cov(s, s)
            end
        end

        @inbounds @simd for epoch in 1:epoch_n
            phase_idx = 1
            Threads.@threads for idx = 2:window:window_n
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
    signal_trim(signal; len::Int64)

Removes `len` samples from the beginning (`from` = :start, default) or end (`from` = :end) of the `signal`.

# Arguments

- `signal::AbstractArray`
- `len::Int64`: trimming length in samples
- `offset::Int64`: offset from which trimming starts, only works for `from` = :start
- `from::Symbol[:start, :end]

# Returns

- `s_trimmed::Vector{Float64}`
"""
function signal_trim(signal::AbstractArray; len::Int64, offset::Int64=1, from::Symbol=:start)

    from in [:start, :end] || throw(ArgumentError("from must be :start or :end."))
    len < 1 && throw(ArgumentError("len must be ≥ 1."))
    len >= length(signal) && throw(ArgumentError("len must be < $(length(signal))."))
    offset < 1 && throw(ArgumentError("offset must be ≥ 1."))
    offset >= length(signal) - 1 && throw(ArgumentError("offset must be < $(length(signal))."))
    (from ===:start && 1 + offset + len > length(signal)) && throw(ArgumentError("offset + len must be < $(length(signal))."))
    
    offset == 1 && (from === :start && (s_trimmed = signal[(offset + len):end]))
    offset > 1 && (from === :start && (s_trimmed = vcat(signal[1:offset], signal[(1 + offset + len):end])))
    from === :end && (s_trimmed = signal[1:(end - len)])
    
    return s_trimmed
end


"""
    signal_trim(signal; len, offset=0, from=:start)

Removes `len` samples from the beginning (`from` = :start, default) or end (`from` = :end) of the `signal`.

# Arguments

- `signal::Array{Float64, 3}`
- `len::Int64`: number of samples to remove
- `offset::Int64`: offset from which trimming starts, only works for `from` = :start
- `from::Symbol[:start, :end]`

# Returns

- `s_trimmed::Array{Float64, 3}`
"""
function signal_trim(signal::Array{Float64, 3}; len::Int64, offset::Int64=1, from::Symbol=:start)

    from in [:start, :end] || throw(ArgumentError("from must be :start or :end."))
    len < 0 && throw(ArgumentError("len must be ≥ 1."))
    len >= size(signal, 2) && throw(ArgumentError("len must be < $(size(signal, 2))."))
    offset < 0 && throw(ArgumentError("offset must be ≥ 1."))
    offset >= size(signal, 2) - 1 && throw(ArgumentError("Offset must be < $(size(signal, 2))."))
    (from ===:start && 1 + offset + len > size(signal, 2)) && throw(ArgumentError("offset + len must be < $(size(signal, 2))."))
    
    channel_n, _, epoch_n = size(signal)

    s_trimmed = zeros(channel_n, (size(signal, 2) - len), epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            s_trimmed[idx, :, epoch] = signal_trim(s, len=len, offset=offset, from=from)
        end
    end

    return s_trimmed
end

"""
    signal_mi(signal1::Vector{Float64}, signal2::Vector{Float64})

Calculates mutual information between `signal1` and `signal2`.

# Arguments

- `signal1::AbstractArray`
- `signal2::AbstractArray`

# Returns

- `mi::Float64`
"""
function signal_mi(signal1::AbstractArray, signal2::AbstractArray)

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

    channel_n, _, epoch_n = size(signal)
    mi = zeros(channel_n, channel_n, epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx1 in 1:channel_n
            for idx2 in 1:channel_n
                s = @view signal[idx2, :, epoch]
                mi[idx1, idx2, epoch] = signal_mi(s, s)
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
    channel_n, _, epoch_n = size(signal1)
    mi = zeros(channel_n, channel_n, epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx1 in 1:channel_n
            for idx2 in 1:channel_n
                s1 = signal1[idx2, :, epoch]
                s2 = signal2[idx2, :, epoch]
                mi[idx1, idx2, epoch] = signal_mi(s1, s2)
            end
        end
    end

    return mi
end

"""
    signal_entropy(signal)

Calculates entropy of `signal`.

# Arguments

- `signal::AbstractArray`

# Returns

- `ent::Float64`
"""
function signal_entropy(signal::AbstractArray)

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

    channel_n, _, epoch_n = size(signal)
    s_entropy = zeros(channel_n, epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            s_entropy[idx, epoch] = signal_entropy(s)
        end
    end

    return s_entropy
end

"""
    signal_average(signal1, signal2)

Averages `signal1` and `signal2`.

# Arguments

- `signal1::AbstractArray`
- `signal2::AbstractArray`

# Returns

- `s_averaged::Vector{Float64}`
"""
function signal_average(signal1::AbstractArray, signal2::AbstractArray)

    length(signal1) == length(signal2) || throw(ArgumentError("Both signals must have the same length."))

    s_averaged = mean(hcat(signal1, signal2), dims=2)

    return s_averaged
end

"""
    signal_average(signal)

Averages all channels of `signal`.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `s_averaged::Array{Float64, 3}`
"""
function signal_average(signal::Array{Float64, 3})

    s_averaged = mean(signal[:, :, :], dims=1)

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
    channel_n, _, epoch_n = size(signal1)
    s_averaged = zeros(size(signal1))

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s1 = @view signal1[idx, :, epoch]
            s2 = @view signal2[idx, :, epoch]
            s_averaged[idx, :, epoch] = signal_average(s1, s2)
        end
    end

    return s_averaged
end

"""
    signal_coherence(signal1, signal2)

Calculates coherence between `signal1` and `signal2`.

# Arguments

- `signal1::AbstractArray`
- `signal2::AbstractArray`

# Returns

- `coherence::Vector{ComplexF64}`
"""
function signal_coherence(signal1::AbstractArray, signal2::AbstractArray)

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

- `signal1::Matrix{Float64}`
- `signal2::Matrix{Float64}`

# Returns

- `coherence::Matrix{Float64}`
"""
function signal_coherence(signal1::Matrix{Float64}, signal2::Matrix{Float64})

    size(signal1) == size(signal2) || throw(ArgumentError("Both signals must have the same size."))
    channel_n, _, _ = size(signal1)
    coherence = zeros(ComplexF64, size(signal1))

    Threads.@threads for idx in 1:channel_n
        s1 = @view signal1[idx, :]
        s2 = @view signal2[idx, :]
        coherence[idx, :] = signal_coherence(s1, s2)
    end

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
    channel_n, _, epoch_n = size(signal1)
    coherence = zeros(ComplexF64, size(signal1))

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s1 = @view signal1[idx, :, epoch]
            s2 = @view signal2[idx, :, epoch]
            coherence[idx, :, epoch] = signal_coherence(s1, s2)
        end
    end

    return coherence
end

"""
    signal_pca(signal, n)

Calculates `n` first PCs for `signal`.

# Arguments

- `signal::Array{Float64, 3}`
- `n::Int64`: number of PCs

# Returns

- `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
- `pc_var::Matrix{Float64}`: PC_VAR(1)..PC_VAR(n) × epoch
"""
function signal_pca(signal::Array{Float64, 3}; n::Int64)

    n < 0 && throw(ArgumentError("n must be ≥ 1."))
    n > size(signal, 1) && throw(ArgumentError("Number of PCs must be ≤ $(size(signal, 1))."))

    channel_n, _, epoch_n = size(signal)
    pc = zeros(n, size(signal, 2), epoch_n)
    pc_var = zeros(n, epoch_n)
    pc_reconstructed = zeros(size(signal))
    M = []
    
    Threads.@threads for epoch in 1:epoch_n
        s = @view signal[:, :, epoch]
        # m_cov = signal_cov(s)
        # eig_val, eig_vec = eigen(m_cov)
        # eig_val_idx = sortperm(eig_val, rev=true)
        # eig_val = eig_val[eig_val_idx]
        # eig_vec = matrix_sort(eig_vec, eig_val_idx)
        # eig_val = 100 .* eig_val / sum(eig_val) # convert to %

        M = MultivariateStats.fit(PCA, s, maxoutdim=n)
        v = principalvars(M) ./ var(M) * 100

        for idx in 1:n
            pc_var[idx, epoch] = v[idx]
            # pc[idx, :, epoch] = (eig_vec[:, idx] .* s)[idx, :]
            pc[idx, :, epoch] = MultivariateStats.predict(M, s)[idx, :]
        end
    end

    return pc, pc_var, M
end

"""
    signal_fconv(signal; kernel)

Performs convolution in the frequency domain between `signal` and `kernel`.

# Arguments

- `signal::AbstractArray`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `s_conv::Vector{ComplexF64}`
"""
function signal_fconv(signal::AbstractArray; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

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
    signal_fconv(signal; kernel)

Performs convolution in the frequency domain between `signal` and `kernel`.

# Arguments

- `signal::Array{Float64, 3}`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `s_conv::Array{ComplexF64, 3}`
"""
function signal_fconv(signal::Array{Float64, 3}; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

    channel_n, _, epoch_n = size(signal)

    s_conv = zeros(ComplexF64, size(signal))

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            s_conv[idx, :, epoch] = signal_fconv(s, kernel=kernel)
        end
    end

    return s_conv
end

"""
    signal_ica(signal, n, tol=1.0e-6, iter=100, f=:tanh)

Calculates `n` first ICs for `signal`.

# Arguments

- `signal::Array{Float64, 3}`
- `n::Int64`: number of PCs
- `tol::Float64`: tolerance for ICA
- `iter::Int64`: maximum number of iterations
- `f::Symbol[:tanh, :gaus]`: neg-entropy functor

# Returns

- `ic::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
- `ic_mw::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
"""
function signal_ica(signal::Array{Float64, 3}; n::Int64, tol::Float64=1.0e-6, iter::Int64=100, f::Symbol=:tanh)

    f in [:tanh, :gaus] || throw(ArgumentError("f must be :tanh or :gaus."))
    n < 0 && throw(ArgumentError("n must be ≥ 1."))
    n > size(signal, 1) && throw(ArgumentError("Number of ICs must be ≤ $(size(signal, 1))."))
    channel_n, _, epoch_n = size(signal)
    ic = zeros(n, size(signal, 2), epoch_n)
    ic_mw = zeros(channel_n, n, epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        s = @view signal[:, :, epoch]

        f === :tanh && (M = MultivariateStats.fit(ICA, s, n, tol=tol, maxiter=iter, fun=MultivariateStats.Tanh(1.0)))
        f === :gaus && (M = MultivariateStats.fit(ICA, s, n, tol=tol, maxiter=iter, fun=MultivariateStats.Gaus()))

        n == size(signal, 1) && (mw = inv(M.W)')
        n < size(signal, 1) && (mw = pinv(M.W)')

        for idx in 1:n
            ic[idx, :, epoch] = MultivariateStats.predict(M, s)[idx, :]
        end

        ic_mw[:, :, epoch] = mw
    end

    return ic, ic_mw
end

"""
    signal_ica_reconstruct(signal, ic_activations, ic_mw, ic_v)

Reconstructs `signal` using removal of `ic_v` ICA components.

# Arguments

- `signal::Array{Float64, 3}`
- `ic_activation::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
- `ic_mw::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
- `ic_v::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove

# Returns

- `signal_reconstructed::Array{Float64, 3}`
"""
function signal_ica_reconstruct(signal::Array{Float64, 3}; ic_activations::Array{Float64, 3}, ic_mw::Array{Float64, 3}, ic_v::Union{Int64, Vector{Int64}, AbstractRange})

    typeof(ic_v) <: AbstractRange && (ic_v = collect(ic_v))
    if typeof(ic_v) == Vector{Int64}
        sort!(ic_v)
        for idx in 1:length(ic_v)
            (ic_v[idx] < 1 || ic_v[idx] > size(ic_mw, 2)) && throw(ArgumentError("ic_v must be ≥ 1 and ≤ $(size(ic_mw, 2))"))
        end
    else
        (ic_v < 1 || ic_v > size(ic_mw, 2)) && throw(ArgumentError("ic_v must be ≥ 1 and ≤ $(size(ic_mw, 2))"))
    end
    ic_removal = setdiff(1:size(ic_mw, 2), ic_v)

    s_reconstructed = zeros(size(signal))

    _, _, epoch_n = size(signal)

    @inbounds @simd for epoch in 1:epoch_n
        s_reconstructed[:, :, epoch] = ic_mw[:, ic_removal, epoch] * ic_activations[ic_removal, :, epoch]
    end

    return s_reconstructed
end

"""
    signal_epochs_var(signal)

Calculates variance for all `signal` epochs.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `var::Vector{Float64}`
"""
function signal_epochs_stats(signal::Array{Float64, 3})

    channel_n, _, epoch_n = size(signal)
    s_mean = zeros(epoch_n)
    s_median = zeros(epoch_n)
    s_sd = zeros(epoch_n)
    s_var = zeros(epoch_n)
    s_kurt = zeros(epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        s = @view signal[:, :, epoch]
        s_mean[epoch] = mean(s)
        s_median[epoch] = median(s)
        s_sd[epoch] = std(s)
        s_var[epoch] = var(s)
        s_kurt[epoch] = kurtosis(s)
    end

    return s_mean, s_median, s_sd, s_var, s_kurt
end

"""
    signal_spectrogram(signal; fs, norm=true, demean=true)

Calculates spectrogram of `signal`.

# Arguments

- `signal::AbstractArray`
- `fs::Int64`: sampling frequency
- `norm::Bool`: normalize powers to dB
- `demean::Bool`: demean signal prior to analysis

# Returns

- `spec.power::Matrix{Float64}`
- `spec.freq::Vector{Float64}`
- `spec.time::Vector{Float64}`
"""
function signal_spectrogram(signal::AbstractArray; fs::Int64, norm::Bool=true, demean::Bool=true)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1 Hz."))

    demean == true && (signal = signal_demean(signal))

    nfft = length(signal)
    interval = fs
    overlap = round(Int64, fs * 0.85)

    spec = spectrogram(signal, interval, overlap, nfft=nfft, fs=fs, window=hanning)
    s_t = collect(spec.time)
    s_frq = Vector(spec.freq)
    if norm == true
        s_pow = pow2db.(spec.power)
    else
        s_pow = Vector(spec.power)
    end

    return s_pow, s_frq, s_t
end


"""
    signal_spectrogram(signal; fs, norm=true, demean=true)

Calculates spectrogram of `signal`.

# Arguments

- `signal::Array{Float64, 3}`
- `fs::Int64`: sampling frequency
- `norm::Bool`: normalize powers to dB
- `demean::Bool`: demean signal prior to analysis

# Returns

- `spec.power::Array{Float64, 4}`
- `spec.freq::Matrix{Float64}`
- `spec.time::Matrix{Float64}`
"""
function signal_spectrogram(signal::Array{Float64, 3}; fs::Int64, norm::Bool=true, demean::Bool=true)

    channel_n, _, epoch_n = size(signal)
    p_tmp, f_tmp, t_tmp = signal_spectrogram(signal[1, :, 1], fs=fs, norm=norm, demean=demean)
    s_pow = zeros(size(p_tmp, 1), size(p_tmp, 2), channel_n, epoch_n)
    s_frq = zeros(length(f_tmp), epoch_n)
    s_t = zeros(length(t_tmp), epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            s_pow[:, :, idx, epoch], s_frq[:, epoch], s_t[:, epoch] = signal_spectrogram(s, fs=fs, norm=norm, demean=demean)
        end
    end

    return s_pow, s_frq, s_t
end

"""
    signal_band(fs, band)

Returns EEG band frequency limits.

# Arguments

- `fs::Int64`
- `band::Symbol`: EEG band name (:delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher)

# Returns

- `band_frequency::Tuple{Float64, Float64}`: lower and upper bounds
"""
function signal_band(fs::Union{Int64, Float64}, band::Symbol)

    fs <= 0 && throw(ArgumentError("fs must be > 0 Hz"))
    band in [:delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher] || throw(ArgumentError("band must be :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower or :gamma_higher."))

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
    
    band_frequency[1] > fs / 2 && (band_frequency = (fs / 2, band_frequency[2]))
    band_frequency[2] > fs / 2 && (band_frequency = (band_frequency[1], fs / 2))

    return band_frequency
end

"""
    signal_detect_epoch_flat(signal::Array{Float64, 3}, threshold=0.1)

Detect bad `signal` epochs based on: flat channel(s)

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch
"""
function signal_detect_epoch_flat(signal::Array{Float64, 3})
    
    channel_n, _, epoch_n = size(signal)

    bad_epochs_score = zeros(epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            idx=1
            epoch=1
            c_tmp_diff = abs.(diff(signal[idx, :, epoch]))
            # add tolerance around zero μV
            sum(c_tmp_diff) < eps() && (bad_epochs_score[epoch] += 1)
        end
    end

    bad_epochs_score = round.(bad_epochs_score ./ channel_n, digits=1)

    return bad_epochs_score
end

"""
    signal_detect_epoch_rmse(signal::Array{Float64, 3})

Detect bad `signal` epochs based on: RMSE vs average channel > 95%CI.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch
"""
function signal_detect_epoch_rmse(signal::Array{Float64, 3})
    
    channel_n, _, epoch_n = size(signal)

    bad_epochs_score = zeros(epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        ch_m = vec(median(signal[:, :, epoch], dims=1))
        rmse_ch = zeros(channel_n)
        Threads.@threads for idx in 1:channel_n
            rmse_ch[idx] = rmse(signal[idx, :, epoch], ch_m)
        end
        Threads.@threads for idx in 1:channel_n
            rmse_ch[idx] > HypothesisTests.confint(OneSampleTTest(rmse_ch))[2] && (bad_epochs_score[epoch] += 1)
        end
    end

    bad_epochs_score = round.(bad_epochs_score ./ channel_n, digits=1)

    return bad_epochs_score
end

"""
    signal_detect_epoch_rmsd(signal::Array{Float64, 3})

Detect bad `signal` epochs based on: RMSD vs average channel > 95%CI.
# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch
"""
function signal_detect_epoch_rmsd(signal::Array{Float64, 3})
    
    channel_n, _, epoch_n = size(signal)

    bad_epochs_score = zeros(epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        ch_m = median(signal[:, :, epoch], dims=1)
        rmsd_ch = zeros(channel_n)
        Threads.@threads for idx in 1:channel_n
            rmsd_ch[idx] = Distances.rmsd(signal[idx, :, epoch], ch_m)
        end
        Threads.@threads for idx in 1:channel_n
            rmsd_ch[idx] > HypothesisTests.confint(OneSampleTTest(rmsd_ch))[2] && (bad_epochs_score[epoch] += 1)
        end
    end

    bad_epochs_score = round.(bad_epochs_score ./ channel_n, digits=1)

    return bad_epochs_score
end

"""
    signal_detect_epoch_euclid(signal::Array{Float64, 3})

Detect bad `signal` epochs based on: Euclidean distance vs median channel > 95% CI.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch
"""
function signal_detect_epoch_euclid(signal::Array{Float64, 3})
    
    channel_n, _, epoch_n = size(signal)

    bad_epochs_score = zeros(epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        ch_m = median(signal[:, :, epoch], dims=1)
        ed_ch = zeros(channel_n)
        Threads.@threads for idx in 1:channel_n
            ed_ch[idx] = euclidean(signal[idx, :, epoch], ch_m)
        end
        Threads.@threads for idx in 1:channel_n
            ed_ch[idx] > HypothesisTests.confint(OneSampleTTest(ed_ch))[2] && (bad_epochs_score[epoch] += 1)
        end
    end

    bad_epochs_score = round.(bad_epochs_score ./ channel_n, digits=1)

    return bad_epochs_score
end

"""
    signal_detect_epoch_p2p(signal::Array{Float64, 3})

Detect bad `signal` epochs based on: p2p amplitude > upper 95% CI p2p amplitude.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch
"""
function signal_detect_epoch_p2p(signal::Array{Float64, 3})
    
    channel_n, _, epoch_n = size(signal)

    bad_epochs_score = zeros(epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        p2p = zeros(channel_n)
        Threads.@threads for idx in 1:channel_n
            p2p[idx] = maximum(signal[idx, :, epoch]) + abs(minimum(signal[idx, :, epoch]))
        end
        Threads.@threads for idx in 1:channel_n
            p2p[idx] > HypothesisTests.confint(OneSampleTTest(p2p))[2] && (bad_epochs_score[epoch] += 1)
        end
    end

    bad_epochs_score = round.(bad_epochs_score ./ channel_n, digits=1)

    return bad_epochs_score
end