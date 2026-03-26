export stationarity
export stationarity_hilbert
export stationarity_mean
export stationarity_var

"""
    stationarity_hilbert(s)

Calculate phase stationarity using Hilbert transformation for a 1-D signal vector.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `Vector{Float64}`
"""
function stationarity_hilbert(s::AbstractVector)::Vector{Float64}

    stph = diff(DSP.unwrap(DSP.angle.(hilbert(s))))

    return stph

end

"""
    stationarity_mean(s; <keyword arguments>)

Calculate mean stationarity for a 1-D signal vector.

Signal is split into `window`-long windows and averaged across windows.

# Arguments

- `s::AbstractVector`: signal vector
- `window::Int64`: time window in samples

# Returns

- `Vector{Float64}`
"""
function stationarity_mean(s::AbstractVector; window::Int64)::Vector{Float64}

    # validate
    window >= 1 || throw(ArgumentError("window must be ≥ 1."))
    window <= length(s) || throw(ArgumentError("window must be ≤ $(length(s))."))

    s = s[1:(window * floor(Int64, length(s) / window))]
    s = reshape(s, Int(length(s) / window), window)

    stm = mean(s, dims = 1)[:]

    return stm

end

"""
    stationarity_var(s; <keyword arguments>)

Calculate variance stationarity for a 1-D signal vector.

Signal is split into `window`-long windows and variance is calculated across windows.

# Arguments

- `s::AbstractVector`: signal vector
- `window::Int64`: time window in samples

# Returns

- `Vector{Float64}`
"""
function stationarity_var(s::AbstractVector; window::Int64)::Vector{Float64}

    # validate
    window >= 1 || throw(ArgumentError("window must be ≥ 1."))
    window <= length(s) || throw(ArgumentError("window must be ≤ $(length(s))."))

    s = s[1:(window * floor(Int64, length(s) / window))]
    s = reshape(s, Int(length(s) / window), window)

    stv = var(s, dims = 1)[:]

    return stv

end

"""
    stationarity(obj; <keyword arguments>)

Calculate stationarity for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}: list of channels
- `window::Int64=10`: time window in samples
- `method::Symbol=:euclid`: stationarity method:
    - `:mean`: mean across `window`-long windows
    - `:var`: variance across `window`-long windows
    - `:cov`: covariance stationarity based on Euclidean distance between covariance matrix of adjacent time windows
    - `:hilbert`: phase stationarity using Hilbert transformation
    - `:adf`: Augmented Dickey–Fuller test; returns ADF-test value and p-value (H0: signal is non-stationary; p-value < alpha means that signal is stationary)

# Returns

- `Union{Matrix{Float64}, Array{Float64, 3}}`
"""
function stationarity(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    window::Int64 = 10,
    method::Symbol = :hilbert
)::Union{Matrix{Float64}, Array{Float64, 3}}

    # validate
    _check_var(method, [:mean, :var, :cov, :hilbert, :adf], "method")
    window >= 1 || throw(ArgumentError("window must be ≥ 1."))
    window <= epoch_len(obj) || throw(ArgumentError("window must be ≤ $(epoch_len(obj))."))

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    ch_n = length(ch)
    ep_n = nepochs(obj)

    if method === :mean

        s = zeros(ch_n, window, ep_n)
        @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
            ch_idx, ep_idx = idx[1], idx[2]
            s[ch_idx, :, ep_idx] = stationarity_mean(@view(obj.data[ch[ch_idx], :, ep_idx]), window = window)
        end

        return s

    end

    if method === :var

        s = zeros(ch_n, window, ep_n)
        @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
            ch_idx, ep_idx = idx[1], idx[2]
            s[ch_idx, :, ep_idx] = stationarity_var(@view(obj.data[ch[ch_idx], :, ep_idx]), window = window)
        end

        return s

    end

    if method === :hilbert

        s = zeros(ch_n, epoch_len(obj) - 1, ep_n)
        @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
            ch_idx, ep_idx = idx[1], idx[2]
            s[ch_idx, :, ep_idx] = stationarity_hilbert(@view(obj.data[ch[ch_idx], :, ep_idx]))
        end

        return s

    end

    if method === :cov

        # validate
        ch_n >= 2 || throw(ArgumentError("For :cov method, number of channels must be ≥ 2."))

        # number of time windows per epoch
        window_n = epoch_len(obj)
        cov_mat = zeros(ch_n, ch_n, window_n, ep_n)
        s = zeros(1 + length(2:window:window_n), ep_n)

        # create covariance matrices per each window
        @inbounds Threads.@threads :static for idx in CartesianIndices((window_n, ep_n))
            window_idx, ep_idx = idx[1], idx[2]
            cov_mat[:, :, window_idx, ep_idx] = covm(
                @view(obj.data[ch, window_idx, ep_idx]),
                @view(obj.data[ch, window_idx, ep_idx])
            )
            end
        end

        # calculate Euclidean distance between adjacent matrices
        @inbounds for ep_idx in 1:ep_n
            w_idx = 1
            Threads.@threads :dynamic for window_idx in 2:window:window_n
                s[w_idx, ep_idx] = euclidean(
                    @view(cov_mat[:, :, window_idx - 1, ep_idx]),
                    @view(cov_mat[:, :, window_idx, ep_idx])
                )
                w_idx += 1
            end
        end

        return s

    end

    if method === :adf

        s = zeros(ch_n, 2, ep_n)

        # initialize progress bar
        progbar = Progress(ep_n * ch_n, dt = 1, barlen = 20, color = :white, enabled = progress_bar)

        # perform Augmented Dickey–Fuller test
        @inbounds for ep_idx in 1:ep_n
            Threads.@threads :dynamic for ch_idx in 1:ch_n
                adf = HypothesisTests.ADFTest(@view(obj.data[ch_idx, :, ep_idx]), :none, 1)
                a = adf.stat
                p = pvalue(adf)
                p < eps() && (p = 0.0001)
                a = round(a, digits = 2)
                p = round(p, digits = 4)
                p == 0.0 && (p = 0.0001)
                s[ch_idx, :, ep_idx] = [a, p]

                # update progress bar
                progress_bar && next!(progbar)
            end
        end

        return s

    end

end
