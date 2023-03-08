export stationarity
export stationarity_hilbert
export stationarity_mean
export stationarity_var

"""
    stationarity_hilbert(signal)

Calculate phase stationarity using Hilbert transformation.

# Arguments

- `signal::AbstractVector`

# Returns

- `phase_stationarity::Vector{Float64}`
"""
function stationarity_hilbert(signal::AbstractVector)
    return diff(DSP.unwrap(angle.(hilbert(signal))))
end

"""
    stationarity_mean(signal; window)

Calculate mean stationarity. Signal is split into `window`-long windows and averaged across windows.

# Arguments

- `signal::AbstractVector`
- `window::Int64`: time window in samples

# Returns

- `mean_stationarity::Vector{Float64}`
"""
function stationarity_mean(signal::AbstractVector; window::Int64)
    window < 1 && throw(ArgumentError("window must be ≥ 1."))
    window > length(signal) && throw(ArgumentError("window must be ≤ $(length(signal))."))
    signal = signal[1:(window * floor(Int64, length(signal) / window))]
    signal = reshape(signal, Int(length(signal) / window), window)
    return mean(signal, dims=1)[:]
end

"""
    stationarity_var(signal; window)

Calculate variance stationarity. Signal is split into `window`-long windows and variance is calculated across windows.

# Arguments

- `signal::AbstractVector`
- `window::Int64`: time window in samples

# Returns

- `var_stationarity::Vector{Float64}`
"""
function stationarity_var(signal::AbstractVector; window::Int64)
    window < 1 && throw(ArgumentError("window must be ≥ 1."))
    window > length(signal) && throw(ArgumentError("window must be ≤ $(length(signal))."))
    signal = signal[1:(window * floor(Int64, length(signal) / window))]
    signal = reshape(signal, Int(length(signal) / window), window)
    return var(signal, dims=1)[:]
end

"""
    stationarity(obj; channel, window, method)

Calculate stationarity.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
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
function stationarity(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), window::Int64=10, method::Symbol=:hilbert)

    _check_var(method, [:mean, :var, :cov, :hilbert, :adf], "method")
    window < 1 && throw(ArgumentError("window must be ≥ 1."))
    window > epoch_len(obj) && throw(ArgumentError("window must be ≤ $(epoch_len(obj))."))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    if method === :mean
        s = zeros(ch_n, window, ep_n)
        @inbounds @simd for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
                s[ch_idx, :, ep_idx] = @views stationarity_mean(obj.data[channel[ch_idx], :, ep_idx], window=window)
            end
        end
        return s
    end

    if method === :var
        s = zeros(ch_n, window, ep_n)
        @inbounds @simd for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
                s[ch_idx, :, ep_idx] = @views stationarity_var(obj.data[channel[ch_idx], :, ep_idx], window=window)
            end
        end
        return s
    end

    if method === :hilbert
        s = zeros(ch_n, epoch_len(obj) - 1, ep_n)
        @inbounds @simd for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
                s[ch_idx, :, ep_idx] = @views stationarity_hilbert(obj.data[channel[ch_idx], :, ep_idx])
            end
        end
        return s
    end

    if method === :cov
        # number of time windows per epoch
        window_n = epoch_len(obj)
        cov_mat = zeros(ch_n, ch_n, window_n, ep_n)
        s = zeros(1 + length(2:window:window_n), ep_n)
        ch_n == 1 && throw(ArgumentError("For :cov method, number of channels must be ≥ 2."))

        # create covariance matrices per each window
        @inbounds @simd for ep_idx in 1:ep_n
            Threads.@threads for window_idx = 1:window_n
                cov_mat[:, :, window_idx, ep_idx] = @views covm(obj.data[channel, window_idx, ep_idx], obj.data[channel, window_idx, ep_idx])
            end
        end

        # calculate Euclidean distance between adjacent matrices
        @inbounds @simd for ep_idx in 1:ep_n
            w_idx = 1
            Threads.@threads for window_idx = 2:window:window_n
                s[w_idx, ep_idx] = @views euclidean(cov_mat[:, :, window_idx - 1, ep_idx], cov_mat[:, :, window_idx, ep_idx])
                w_idx += 1
            end
        end
        return s
    end

    if method === :adf
        s = zeros(ch_n, 2, ep_n)

        # initialize progress bar
        progress_bar == true && (pb = Progress(ep_n * ch_n, 1))

        # perform Augmented Dickey–Fuller test
        @inbounds @simd for ep_idx in 1:ep_n
            Threads.@threads for ch_idx = 1:ch_n
                adf = @views ADFTest(obj.data[channel[ch_idx], :, ep_idx], :constant, window)
                a = adf.stat
                p = pvalue(adf)
                p < eps() && (p = 0.0001)
                a = round(a, digits=2)
                p = round(p, digits=4)
                p == 0.0 && (p = 0.0001)
                s[ch_idx, :, ep_idx] = [a, p]

                # update progress bar
                progress_bar == true && next!(pb)
            end
        end
        return s
    end
end
