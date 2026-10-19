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

The signal is split into `window` consecutive, equal segments (trailing samples that do not fill a segment are dropped); the mean of each segment is returned.

# Arguments

- `s::AbstractVector`: signal vector
- `window::Int64`: number of segments

# Returns

- `Vector{Float64}`: segment means
"""
function stationarity_mean(s::AbstractVector; window::Int64)::Vector{Float64}
    # validate
    window >= 1 || throw(ArgumentError("window must be ≥ 1."))
    window <= length(s) || throw(ArgumentError("window must be ≤ $(length(s))."))

    s = s[1:(window * floor(Int64, length(s) / window))]
    s = reshape(s, Int(length(s) / window), window)

    stm = mean(s; dims = 1)[:]

    return stm
end

"""
    stationarity_var(s; <keyword arguments>)

Calculate variance stationarity for a 1-D signal vector.

The signal is split into `window` consecutive, equal segments (trailing samples that do not fill a segment are dropped); the variance of each segment is returned.

# Arguments

- `s::AbstractVector`: signal vector
- `window::Int64`: number of segments

# Returns

- `Vector{Float64}`: segment variances
"""
function stationarity_var(s::AbstractVector; window::Int64)::Vector{Float64}
    # validate
    window >= 1 || throw(ArgumentError("window must be ≥ 1."))
    window <= length(s) || throw(ArgumentError("window must be ≤ $(length(s))."))

    s = s[1:(window * floor(Int64, length(s) / window))]
    s = reshape(s, Int(length(s) / window), window)

    stv = var(s; dims = 1)[:]

    return stv
end

"""
    stationarity(obj; <keyword arguments>)

Calculate stationarity for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `window::Int64=10`: for `:mean` and `:var`: number of segments; for `:cov`: segment length in samples
- `method::Symbol=:hilbert`: stationarity method:
    - `:mean`: mean of each of `window` segments, shape (channels, window, epochs)
    - `:var`: variance of each of `window` segments, shape (channels, window, epochs)
    - `:cov`: covariance stationarity: Euclidean (Frobenius) distance between the channel covariance matrices of adjacent `window`-sample segments, shape (segments - 1, epochs)
    - `:hilbert`: phase stationarity using Hilbert transformation
    - `:adf`: Augmented Dickey–Fuller test; returns ADF-test value and p-value (H0: signal is non-stationary; p-value < alpha means that signal is stationary)

# Returns

- `Union{Matrix{Float64}, Array{Float64, 3}}`
"""
function stationarity(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    window::Int64 = 10,
    method::Symbol = :hilbert,
)::Union{Matrix{Float64}, Array{Float64, 3}}
    # validate
    _check_var(method, [:mean, :var, :cov, :hilbert, :adf], "method")
    window >= 1 || throw(ArgumentError("window must be ≥ 1."))
    window <= epoch_len(obj) || throw(ArgumentError("window must be ≤ $(epoch_len(obj))."))

    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)

    if method === :mean
        s = zeros(ch_n, window, ep_n)
        @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
            ch_idx, ep_idx = idx[1], idx[2]
            s[ch_idx, :, ep_idx] =
                stationarity_mean(@view(obj.data[ch[ch_idx], :, ep_idx]), window = window)
        end

        return s
    end

    if method === :var
        s = zeros(ch_n, window, ep_n)
        @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
            ch_idx, ep_idx = idx[1], idx[2]
            s[ch_idx, :, ep_idx] =
                stationarity_var(@view(obj.data[ch[ch_idx], :, ep_idx]), window = window)
        end

        return s
    end

    if method === :hilbert
        s = zeros(ch_n, epoch_len(obj) - 1, ep_n)
        @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
            ch_idx, ep_idx = idx[1], idx[2]
            s[ch_idx, :, ep_idx] =
                stationarity_hilbert(@view(obj.data[ch[ch_idx], :, ep_idx]))
        end

        return s
    end

    if method === :cov

        # validate
        ch_n >= 2 ||
            throw(ArgumentError("For :cov method, number of channels must be ≥ 2."))

        window >= 2 || throw(ArgumentError("For :cov method, window must be ≥ 2."))
        # number of complete windows per epoch
        window_n = epoch_len(obj) ÷ window
        window_n >= 2 || throw(ArgumentError("For :cov method, the epoch must contain ≥ 2 windows."))
        s = zeros(window_n - 1, ep_n)
        @inbounds Threads.@threads :dynamic for ep_idx = 1:ep_n
            # channel covariance matrix of each window
            cm = [
                cov(permutedims(@view(obj.data[ch, ((w - 1) * window + 1):(w * window), ep_idx]))) for
                w = 1:window_n
            ]
            # Euclidean distance between the covariance matrices of adjacent windows
            for w = 2:window_n
                s[w - 1, ep_idx] = euclidean(cm[w - 1], cm[w])
            end
        end
        return s
    end

    if method === :adf
        s = zeros(ch_n, 2, ep_n)

        # initialize progress bar
        progbar = Progress(
            ep_n * ch_n;
            dt = 1,
            barlen = 20,
            color = :white,
            enabled = progress_bar,
        )

        # perform Augmented Dickey–Fuller test
        @inbounds for ep_idx = 1:ep_n
            Threads.@threads :dynamic for ch_idx = 1:ch_n
                adf = HypothesisTests.ADFTest(@view(obj.data[ch[ch_idx], :, ep_idx]), :none, 1)
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
