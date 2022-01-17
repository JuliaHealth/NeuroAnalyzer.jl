"""
    signal_derivative(x)

Returns the derivative of the signal `x` with length same as the signal
"""
signal_derivative(x::Vector) = vcat(diff(x), diff(x)[end])

"""
    band_power(psd, f1, f2)

Calculates absolute band power between frequencies `f1` and `f2` for the `psd` power spectrum.
"""
function signal_band_power(psd, f1, f2)
    frq_idx = [vsearch(psd.freq, f1), vsearch(psd.freq, f2)]
    # dx: frequency resolution
    dx = psd.freq[2] - psd.freq[1]
    result = simpson(psd.power[frq_idx[1]:frq_idx[2]], start=frq_idx[1], stop=frq_idx[1], dx=dx)
    return result
end

"""
    make_spectrum(x, fs)

Returns FFT and DFT sample frequencies for a DFT of the `x` signal 
"""
function signal_make_spectrum(x::Vector, fs)
    hs = fft(x)
    n = length(x)               # number of samples
    d = 1/fs                    # time between samples
    fs = fftfreq(n, d)
    return hs, fs
end

"""
    signal_detrend(x, type=:linear)

Removes linear trend of signal `x`.

# Arguments
- `x::Vector` - the signal to detrend.
- `type::Symbol[:linear, :constant]`, optional.
    linear: the result of a linear least-squares fit to `y` is subtracted from `y`
    constant: the mean of `y` is subtracted.
"""
function signal_detrend(x::Vector; trend=:linear)
    trend in [:linear, :constant] || throw(ArgumentError("""Trend type must be ":linear" or ":constant"."""))
    if trend == :constant
        result = x .- mean(x)
    else
        A = ones(length(x))
        coef = A \ x
        result = @. x - dot(A, coef)
    end
    return result
end

"""
    signals_ci95(signals::Matrix, n=3; method=:normal)

Calculates mean, std and 95% confidence interval for the matrix of signals .

# Arguments

- `signals::Matrix` - the signal matrix to analyze (rows: trials, columns: time).
- `n::Int` - number of bootstraps.
- `method::Symbol[:normal, :boot]` - use normal method or `n`-times boostrapping.
"""
function signals_ci95(signals::Matrix, n=3; method=:normal)
    method in [:normal, :boot] || throw(ArgumentError("""Method must be ":normal" or ":boot"."""))
    if method === :normal
        signals_mean = mean(signals, dims=1)'
        signals_sd = std(signals, dims=1)' / sqrt(size(signals, 1))
        upper_bound = signals_mean + 1.96 * signals_sd
        lower_bound = signals_mean - 1.96 * signals_sd
    else
        signals_tmp1 = zeros(size(signals, 1) * n, size(signals, 2))
        Threads.@threads for idx1 in 1:size(signals, 1) * n
            signals_tmp2 = zeros(size(signals))
            sample_idx = rand(1:size(signals, 1), size(signals, 1))
            for idx2 in 1:size(signals, 1)
                signals_tmp2[idx2, :] = signals[sample_idx[idx2], :]'
            end
            signals_tmp1[idx1, :] = mean(signals_tmp2, dims=1)
        end
        signals_mean = mean(signals_tmp1, dims=1)'
        signals_sd = std(signals_tmp1, dims=1)' / sqrt(size(signals_tmp1, 1))
        signal_sorted = sort(signals_tmp1, dims=1)
        lower_bound = signal_sorted[round(Int, 0.025 * size(signals_tmp1, 1)), :]
        upper_bound = signal_sorted[round(Int, 0.975 * size(signals_tmp1, 1)), :]
    end
    return signals_mean, signals_sd, upper_bound, lower_bound
end

"""
    signals_mean(signals1::Matrix, signals2::Matrix)

Calculates mean difference and 95% confidence interval for 2 signals.

# Arguments

- `signal1::Matrix` - the signal matrix to analyze (rows: trials, columns: time).
- `signal2:Matrix` - the signal matrix to analyze (rows: trials, columns: time).
"""
function signals_mean(signals1::Matrix, signals2::Matrix)
    signals1_mean = mean(signals1, dims=1)'
    signals2_mean = mean(signals2, dims=1)'
    signals_mean = signals1_mean - signals2_mean
    signals1_sd = std(signals1, dims=1) / sqrt(size(signals1, 1))
    signals2_sd = std(signals2, dims=1) / sqrt(size(signals2, 1))
    signals_mean_sd = sqrt.(signals1_sd.^2 .+ signals2_sd.^2)'

    return signals_mean, signals_mean_sd, signals_mean + 1.96 * signals_mean_sd, signals_mean - 1.96 * signals_mean_sd
end

"""
    signals_difference(signals1::Matrix, signals2::Matrix, n=3; method=:absdiff)

Calculates mean difference and 95% confidence interval for 2 signals.

# Arguments

- `signal1::Matrix` - the signal matrix to analyze (rows: trials, columns: time).
- `signal2:Matrix` - the signal matrix to analyze (rows: trials, columns: time).
- `n::Int` - number of bootstraps.
- `method::Symbol[:absdiff, :diff2int]`
"""
function signals_difference(signals1::Matrix, signals2::Matrix, n::Int=3; method=:absdiff)
    method in [:absdiff, :diff2int] || throw(ArgumentError("""Method must be ":absdiff" or ":diff2int"."""))
    signals1_mean = mean(signals1, dims=1)'
    signals2_mean = mean(signals2, dims=1)'

    if method === :absdiff
        # statistic: maximum difference
        signals_diff = signals1_mean - signals2_mean
        signals_statistic_single = findmax(abs.(signals_diff))[1]
    else
        # statistic: integrated area of the squared difference
        signals_diff_squared = (signals1_mean - signals2_mean).^2
        signals_statistic_single = simpson(signals_diff_squared)
    end

    signals = [signals1; signals2]
    signals_statistic = zeros(size(signals1, 1) * n)

    Threads.@threads for idx1 in 1:(size(signals1, 1) * n)
        signals_tmp1 = zeros(size(signals1, 1), size(signals1, 2))
        sample_idx = rand(1:size(signals, 1), size(signals, 1))
        sample_idx = sample_idx[1:1000]
        for idx2 in 1:size(signals1, 1)
            signals_tmp1[idx2, :] = signals[sample_idx[idx2], :]'
        end
        signals1_mean = mean(signals_tmp1, dims=1)
        signals_tmp1 = zeros(size(signals1, 1), size(signals1, 2))
        sample_idx = rand(1:size(signals, 1), size(signals, 1))
        sample_idx = sample_idx[1:1000]
        for idx2 in 1:size(signals1, 1)
            signals_tmp1[idx2, :] = signals[sample_idx[idx2], :]'
        end
        signals2_mean = mean(signals_tmp1, dims=1)
        if method === :absdiff
            # statistic: maximum difference
            signals_diff = signals1_mean - signals2_mean
            signals_statistic[idx1] = findmax(abs.(signals_diff))[1]
        else
            # statistic: integrated area of the squared difference
            signals_diff_squared = (signals1_mean - signals2_mean).^2
            signals_statistic[idx1] = simpson(signals_diff_squared)
        end
    end

    p = length(signals_statistic[signals_statistic .> signals_statistic_single]) / size(signals1, 1) * n

    return signals_statistic, signals_statistic_single, p
end

"""
   signal_autocov(signal, lag=1, demean=true)

Calculates autocovariance for the `signal` vector for lags = -lag:lag.
"""
function signal_autocov(signal::Vector, lag=1, demean=false, normalize=false)
    lags = collect(-lag:lag)

    if demean == true
        signal_demeaned = signal .- mean(signal)
    else
        signal_demeaned = signal
    end

    ac = zeros(length(lags))

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
            ac[idx] = signals_sum / length(signal)
        else
            ac[idx] = signals_sum
        end
    end

    return ac, lags
end

"""
   signal_crosscov(signal1, signal2, lag=1, demean=true)

Calculates cross-covariance between `signal1` and `signal2` vectors for lags = -lag:lag.
"""
function signal_crosscov(signal1::Vector, signal2::Vector, lag=1, demean=false, normalize=false)
    lags = collect(-lag:lag)

    if demean == true
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
