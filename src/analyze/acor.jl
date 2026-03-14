export acor

"""
    acor(s; <keyword arguments>)

Calculate auto-correlation.

# Arguments

- `s::AbstractVector`: signal vector
- `l::Int64=round(Int64, min(length(s) - 1, 10 * log10(length(s))))`: range of lags is `-l:l`
- `demean::Bool=true`: demean signal before computing auto-correlation
- `biased::Bool=true`: calculate biased or unbiased autocovariance
- `method::Symbol=:sum`: method of calculating auto-correlation:
    - `:sum`: `acf = sum(s[1:(end - l)] .* s[(1 + l):end]) ./ var(s)`
    - `:cor`: `acf = cov(s[1:(end - l)], s[(1 + l):end]) ./ var(s)`
    - `:stat`: use StatsBase `autocor()`, `biased` value is ignored

# Returns

- `ac::Vector{Float64}`: auto-correlation of length `2l + 1`
"""
function acor(
    s::AbstractVector;
    l::Int64 = round(Int64, min(length(s) - 1, 10 * log10(length(s)))),
    demean::Bool = true,
    biased::Bool = true,
    method::Symbol = :sum,
)::Vector{Float64}

    # reject any method symbol not in the supported set
    _check_var(method, [:sum, :cor, :stat], "method")

    if method === :sum

        # compute raw autocovariance via element-wise products, then divide by
        # the signal variance to obtain a correlation-scale (unitless) result
        ac = acov(
            s,
            l = l,
            demean = demean,
            biased = biased,
            method = :sum,
        )
        # normalize by variance to bring values onto [-1, 1] scale
        ac = round.(ac ./ Statistics.var(s), digits = 3)

    elseif method === :cor

        # similar to :sum but uses a covariance-based approach internally;
        # the final normalization step is the same
        ac = acov(
            s,
            l = l,
            demean = demean,
            biased = biased,
            method = :cov,
        )
        # normalize by variance to bring values onto [-1, 1] scale
        ac = round.(ac ./ Statistics.var(s), digits = 3)

    elseif method === :stat

        # delegate entirely to StatsBase.autocor, which handles normalization internally
        # the `biased` keyword is intentionally ignored here
        ac = StatsBase.autocor(
            s,
            0:l,
            demean = demean,
        )
        ac = round.(ac, digits = 3)
        ac = vcat(reverse(ac), ac[2:end])

    end

    return ac

end

"""
    acor(s; <keyword arguments>)

Calculate auto-correlation.

# Arguments

- `s::AbstractArray`: signal array (channels, samples, epochs)
- `l::Int64=round(Int64, min(size(s, 2) - 1, 10 * log10(size(s, 2))))`: range of lags is `-l:l`
- `demean::Bool=true`: demean signal before computing auto-correlation
- `biased::Bool=true`: calculate biased or unbiased autocovariance
- `method::Symbol=:sum`: method of calculating auto-correlation:
    - `:sum`: `acf = sum(s[1:(end - l)] .* s[(1 + l):end]) ./ var(s)`
    - `:cor`: `acf = cov(s[1:(end - l)], s[(1 + l):end]) ./ var(s)`
    - `:stat`: use StatsBase `autocor()`, `biased` value is ignored

# Returns

- `ac::Array{Float64, 3}`: auto-correlations, shape `(channels, 2l+1, epochs)`
"""
function acor(
    s::AbstractArray;
    l::Int64 = round(Int64, min(size(s, 2) - 1, 10 * log10(size(s, 2)))),
    demean::Bool = true,
    biased::Bool = true,
    method::Symbol = :sum,
)::Array{Float64, 3}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    ac = zeros(ch_n, length((-l):l), ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        ac[ch_idx, :, ep_idx] = acor(
            @view(s[ch_idx, :, ep_idx]),
            l = l,
            demean = demean,
            biased = biased,
            method = method,
        )
    end

    return ac

end

"""
    acor(obj; <keyword arguments>)

Calculate auto-correlation. For ERP return trial-averaged auto-correlation.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `l::Int64=round(Int64, min(size(obj.data, 2) - 1, 10 * log10(size(obj.data, 2))))`: range of lags is `-l:l`
- `demean::Bool=true`: demean signal before computing auto-correlation
- `biased::Bool=true`: calculate biased or unbiased autocovariance
- `method::Symbol=:sum`: method of calculating auto-correlation:
    - `:sum`: `acf = sum(s[1:(end - l)] .* s[(1 + l):end]) ./ var(s)`
    - `:cor`: `acf = cov(s[1:(end - l)], s[(1 + l):end]) ./ var(s)`
    - `:stat`: use StatsBase `autocor()`, `biased` value is ignored

# Returns

Named tuple:

- `ac::Array{Float64, 3}`: auto-correlations, shape `(channels, 2l+1, epochs)`
- `l::Vector{Float64}`: lags in seconds
"""
function acor(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    l::Int64=round(Int64, min(size(obj.data, 2) - 1, 10 * log10(size(obj.data, 2)))),
    demean::Bool = true,
    biased::Bool = true,
    method::Symbol = :sum,
)::@NamedTuple{ac::Array{Float64, 3}, l::Vector{Float64}}

    # validate lag bounds: must be non-negative and within the signal length
    !(l <= size(obj, 2)) && throw(ArgumentError("l must be ≤ $(size(obj, 2))."))
    !(l >= 0) && throw(ArgumentError("l must be ≥ 0."))

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    if datatype(obj) == "erp"

        # epoch 1 is the pre-computed average; epochs 2:end are individual trials
        # compute per-trial auto-correlations first, then prepend the mean across
        # trials as epoch 1 of the output (preserving the ERP convention)
        ac = acor(
            @view(obj.data[ch, :, 2:end]),
            l = l,
            demean = demean,
            biased = biased,
            method = method,
        )
        ac = cat(mean(ac, dims = 3), ac, dims = 3)

    else

        ac = acor(
            @view(obj.data[ch, :, :]),
            l = l,
            demean = demean,
            biased = biased,
            method = method,
        )

    end

    # convert integer lag indices (-l … l) to physical time in seconds using
    # the object's sampling rate
    return (ac = ac, l = collect((-l):l) .* 1 / sr(obj))

end
