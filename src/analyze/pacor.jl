export pacor

"""
    pacor(s; <keyword arguments>)

Computes the partial auto-correlation function (PACF) of a signal over symmetric lags −l … l using one of two methods:

- `:yw` (Yule-Walker) – fast, no matrix inversion, may need larger l
- `:reg` (regression) – via successive OLS models, more robust for small l

# Arguments

- `s::AbstractVector`: signal vector
- `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: symmetric lag range is `−l:l`
- `demean::Bool=true`: demean signal before computing PACF
- `method::Symbol=:yw`
- `:yw`: Yule-Walker equations
- `:reg`: successive regression models

# Returns

- `pac::Array{Float64, 3}`: shape `(1, 2l+1, 1)`

# Notes

If you get `ERROR: PosDefException: matrix is not positive definite; Cholesky factorization failed.`, try lowering `l` value or change method to `:yw`.
"""
function pacor(
    s::AbstractVector;
    l::Int64 = round(Int64, min(length(s) - 1, 10 * log10(length(s)))),
    demean::Bool = true,
    method::Symbol = :yw,
)::Array{Float64, 3}

    _check_var(method, [:reg, :yw], "method")

    # map short symbol names to the names expected by StatsBase.pacf()
    method === :reg && (method = :regression)
    method === :yw && (method = :yulewalker)

    s_tmp = demean ? remove_dc(s) : s

    # positive and negative lags computed separately (reverse trick for negative)
    pac = pacf(s_tmp, collect(0:l), method = method)
    pac_neg = pacf(reverse(s_tmp), collect(0:l), method = method)

    # concatenate to form the symmetric lag vector: [−l … −1, 0, +1 … +l].
    pac = vcat(reverse(pac_neg), pac[2:end])
    pac = round.(pac, digits = 3)

    return reshape(pac, 1, :, 1)

end

"""
    pacor(s; <keyword arguments>)

Calculate partial auto-correlation function (PACF) for each epoch of a matrix over symmetric lags −l … l using one of two methods:

- `:yw` (Yule-Walker) – fast, no matrix inversion, may need larger l
- `:reg` (regression) – via successive OLS models, more robust for small l

# Arguments

- `s::AbstractMatrix`: signal matrix (samples × epochs)
- `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: symmetric lag range is `−l:l`
- `demean::Bool=true`: demean signal before computing PACF
- `method::Symbol=:yw`: method of calculating auto-correlation:
- `:yw`: Yule-Walker equations
- `:reg`: successive regression models

# Returns

- `pac::Array{Float64, 3}`
"""
function pacor(
    s::AbstractMatrix;
    l::Int64 = round(Int64, min(size(s[:, 1], 1) - 1, 10 * log10(size(s[:, 1], 1)))),
    demean::Bool = true,
    method::Symbol = :yw,
)::Array{Float64, 3}

    # number of epochs
    ep_n = size(s, 2)

    # allocate output
    pac = zeros(1, length((-l):l), ep_n)

    @inbounds for ep_idx in 1:ep_n
        pac[1, :, ep_idx] = vec(pacor(@view(s[:, ep_idx]), l = l, demean = demean, method = method))
    end

    return pac

end

"""
    pacor(s; <keyword arguments>)

Calculate partial auto-correlation function (PACF) for a 3-D signal array over lags −l … l using one of two methods:

- `:yw` (Yule-Walker) – fast, no matrix inversion, may need larger l
- `:reg` (regression) – via successive OLS models, more robust for small l

# Arguments

- `s::AbstractArray`: signal array (channels, samples, epochs)
- `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: symmetric lag range is `−l:l`
- `demean::Bool=true`: demean signal before computing PACF
- `method::Symbol=:yw`: method of calculating auto-correlation:
- `:yw`: Yule-Walker equations
- `:reg`: successive regression models

# Returns

- `pac::Array{Float64, 3}`: shape `(channels, 2l+1, epochs)`
"""
function pacor(
    s::AbstractArray;
    l::Int64 = round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1)))),
    demean::Bool = true,
    method::Symbol = :yw,
)::Array{Float64, 3}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)
    !(size(s, 1) == 1) && throw(ArgumentError("s must have 1 channel."))

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    pac = zeros(ch_n, length((-l):l), ep_n)

    # calculate over channels and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        pac[ch_idx, :, ep_idx] = vec(pacor(@view(s[ch_idx, :, ep_idx]),
                                           l = l, demean = demean, method = method))
    end

    return pac

end

"""
    pacor(obj; <keyword arguments>)

Calculate partial auto-correlation function (PACF) over lags −l … l using one of two methods:

- `:yw` (Yule-Walker) – fast, no matrix inversion, may need larger l
- `:reg` (regression) – via successive OLS models, more robust for small l

For ERP objects, epoch 1 is the trial-averaged waveform and is prepended to the per-trial result.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `l::Real=1`: lag limit in seconds; lags vector is `−l:1/sr(obj):l`
- `demean::Bool=true`: demean signal before computing PACF
- `method::Symbol=:yw`: method of calculating auto-correlation:
- `:yw`: Yule-Walker equations
- `:reg`: successive regression models

# Returns

Named tuple:

- `pac::Array{Float64, 3}`: partial auto-correlations
- `l::Vector{Float64}`: lag values in seconds
"""
function pacor(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    l::Real = 1,
    demean::Bool = true,
    method::Symbol = :yw,
)::@NamedTuple{pac::Array{Float64, 3}, l::Vector{Float64}}

    !(!(method === :yw && l <= 1)) && throw(ArgumentError("For method=:yw, l must be > 1."))

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    # convert l from seconds to samples for the inner call
    l_samp = round(Int64, l * sr(obj))

    !(l_samp <= size(obj, 2)) && throw(ArgumentError("l must be <= $(size(obj, 2) / sr(obj)) s."))
    !(l_samp >= 0) && throw(ArgumentError("l must be >= 0."))

    if datatype(obj) == "erp"
        pac = pacor(
            @view(obj.data[ch, :, 2:end]),
            l = l_samp,
            demean = demean,
            method = method
        )
        pac = cat(mean(pac, dims = 3), pac, dims = 3)
    else
        pac = pacor(
            @view(obj.data[ch, :, :]),
            l = l_samp,
            demean = demean,
            method = method
        )
    end

    # build the symmetric lag vector in seconds
    lags = collect((-l_samp):l_samp) ./ sr(obj)

    return (pac = pac, l = lags)

end
