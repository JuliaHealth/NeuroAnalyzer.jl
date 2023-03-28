export stationarity
export stationarity_hilbert
export stationarity_mean
export stationarity_var

"""
    stationarity_hilbert(s)

Calculate phase stationarity using Hilbert transformation.

# Arguments

- `s::AbstractVector`

# Returns

- `stph::Vector{Float64}`
"""
function stationarity_hilbert(s::AbstractVector)

    stph = diff(DSP.unwrap(angle.(hilbert(s))))

    return stph

end

"""
    stationarity_mean(s; window)

Calculate mean stationarity. Signal is split into `window`-long windows and averaged across windows.

# Arguments

- `s::AbstractVector`
- `window::Int64`: time window in samples

# Returns

- `stm::Vector{Float64}`
"""
function stationarity_mean(s::AbstractVector; window::Int64)

    window < 1 && throw(ArgumentError("window must be ≥ 1."))
    window > length(s) && throw(ArgumentError("window must be ≤ $(length(s))."))

    s = s[1:(window * floor(Int64, length(s) / window))]
    s = reshape(s, Int(length(s) / window), window)

    stm = mean(s, dims=1)[:]

    return stm

end

"""
    stationarity_var(s; window)

Calculate variance stationarity. Signal is split into `window`-long windows and variance is calculated across windows.

# Arguments

- `s::AbstractVector`
- `window::Int64`: time window in samples

# Returns

- `stv::Vector{Float64}`
"""
function stationarity_var(s::AbstractVector; window::Int64)
    
    window < 1 && throw(ArgumentError("window must be ≥ 1."))
    window > length(s) && throw(ArgumentError("window must be ≤ $(length(s))."))

    s = s[1:(window * floor(Int64, length(s) / window))]
    s = reshape(s, Int(length(s) / window), window)

    stv = var(s, dims=1)[:]

    return stv 

end

"""
    stationarity(obj; ch, window, method)

Calculate stationarity.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
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
function stationarity(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), window::Int64=10, method::Symbol=:hilbert)

    _check_var(method, [:mean, :var, :cov, :hilbert, :adf], "method")
    window < 1 && throw(ArgumentError("window must be ≥ 1."))
    window > epoch_len(obj) && throw(ArgumentError("window must be ≤ $(epoch_len(obj))."))

    ch_n = channel_n(obj)
    ep_n = epoch_n(obj)

    if method === :mean

        s = zeros(ch_n, window, ep_n)
        
        @inbounds @simd for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
            end
        end
        return s
    end

    if method === :var

        s = zeros(ch_n, window, ep_n)

        @inbounds @simd for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
            end
        end
        return s
    end

    if method === :hilbert

        s = zeros(ch_n, epoch_len(obj) - 1, ep_n)

        @inbounds @simd for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
            end
        end

        return s

    end

    if method === :cov

        ch_n == 1 && throw(ArgumentError("For :cov method, number of channels must be ≥ 2."))
        
        # number of time windows per epoch
        window_n = epoch_len(obj)
        cov_mat = zeros(ch_n, ch_n, window_n, ep_n)
        s = zeros(1 + length(2:window:window_n), ep_n)

        # create covariance matrices per each window
        @inbounds @simd for ep_idx in 1:ep_n
            Threads.@threads for window_idx = 1:window_n
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
                adf = @views HypothesisTests.ADFTest(obj.data[ch_idx, :, ep_idx], :none, 1)
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
