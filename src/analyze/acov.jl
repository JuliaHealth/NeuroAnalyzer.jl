export acov

"""
    acov(s; <keyword arguments>)

Calculate auto-covariance.

# Arguments

- `s::AbstractVector`: signal vector
- `l::Int64=round(Int64, min(length(s) - 1, 10 * log10(length(s))))`: range of lags is `-l:l`
- `demean::Bool=true`: demean signal before computing auto-covariance
- `biased::Bool=true`: calculate biased or unbiased auto-covariance
- `method::Symbol=:sum`: method of calculating auto-covariance:
  - `:sum`: `acf = sum(s[1:(end - l)] .* s[(1 + l):end])`
  - `:cov`: `acf = cov(s[1:(end - l)], s[(1 + l):end])`
  - `:stat`: use StatsBase `autocov()`, `biased` value is ignored

# Returns

  - `ac::Vector{Float64}`: auto-covariance of length `2l + 1`
"""
function acov(
    s::AbstractVector;
    l::Int64 = round(Int64, min(length(s) - 1, 10 * log10(length(s)))),
    demean::Bool = true,
    biased::Bool = true,
    method::Symbol = :sum,
)::Vector{Float64}

    # reject any method symbol not in the supported set
    _check_var(method, [:sum, :cov, :stat], "method")

    # pre-allocate output for lags 0 … l (negative lags added later)
    ac = zeros(l + 1)

    if method === :sum

        # compute auto-covariance via dot products
        demean && (s = remove_dc(s))
        n = length(s)
        # hoist the biased branch: denominator formula is the only difference
        denom = biased ? (idx -> n) : (idx -> n - idx)
        @inbounds for idx in 0:l
            # dot avoids allocating the intermediate product array
            ac[idx + 1] = dot(@view(s[1:(end - idx)]),
                              @view(s[(1 + idx):end])) / denom(idx)
        end

    elseif method === :cov

        # similar to :sum but uses a covariance-based approach internally
        demean && (s = remove_dc(s))
        corrected = !biased
        @inbounds for idx in 0:l
            ac[idx + 1] = cov(
                            @view(s[1:(end - idx)]),
                            @view(s[(1 + idx):end]),
                            corrected = corrected,
                        )
        end

    elseif method === :stat

        # delegate entirely to StatsBase.autocor, which handles normalization internally
        # the `biased` keyword is intentionally ignored here
        ac = StatsBase.autocov(
                            s,
                            0:l,
                            demean = demean,
                        )

    end

    ac = round.(ac, digits = 3)
    ac = vcat(reverse(ac), ac[2:end])

    return ac

end

"""
    acov(s; <keyword arguments>)

Calculate auto-covariance.

# Arguments

- `s::AbstractArray`: signal array (channels × samples × epochs)
- `l::Int64=round(Int64, min(size(s, 2) - 1, 10 * log10(size(s, 2))))`: range of lags is `-l:l`
- `demean::Bool=true`: demean signal before computing auto-covariance
- `biased::Bool=true`: calculate biased or unbiased auto-covariance
- `method::Symbol=:sum`: method of calculating auto-covariance:
  - `:sum`: `acf = sum(s[1:(end - l)] .* s[1+l:end])`
  - `:cov`: `acf = cov(s[1:(end - l)], s[1+l:end])`
  - `:stat`: use StatsBase `autocor()`, `biased` value is ignored

# Returns

  - `ac::Array{Float64, 3}`: auto-covariances, shape `(channels, 2l+1, epochs)`
"""
function acov(
    s::AbstractArray;
    l::Int64 = round(Int64, min(size(s, 2) - 1, 10 * log10(size(s, 2)))),
    demean::Bool = true,
    biased::Bool = true,
    method::Symbol = :sum,
)::Array{Float64, 3}

    # validate that the input is a proper 3-D array (channels × samples × epochs)
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
        ac[ch_idx, :, ep_idx] = acov(
            @view(s[ch_idx, :, ep_idx]),
            l = l,
            demean = demean,
            biased = biased,
            method = method
        )
    end

    return ac

end

"""
    acov(obj; <keyword arguments>)

Calculate auto-covariance. For ERP return trial-averaged auto-covariance.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `l::Int64=round(Int64, min(size(obj.data, 2) - 1, 10 * log10(size(obj.data, 2))))`: range of lags is `-l:l`
- `demean::Bool=true`: demean signal before computing auto-covariance
- `biased::Bool=true`: calculate biased or unbiased auto-covariance
- `method::Symbol=:sum`: method of calculating auto-covariance:
  - `:sum`: `acf = sum(s[1:(end - l)] .* s[(1 + l):end])`
  - `:cov`: `acf = cov(s[1:(end - l)], s[(1 + l):end])`
  - `:stat`: use StatsBase `autocor()`, `biased` value is ignored

# Returns

Named tuple containing:

- `ac::Array{Float64, 3}`: auto-covariances of, shape `(channels, 2l+1, epochs)`
- `l::Vector{Float64}`: lags in seconds
"""
function acov(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    l::Int64=round(Int64, min(size(obj.data, 2) - 1, 10 * log10(size(obj.data, 2)))),
    demean::Bool = true,
    biased::Bool = true,
    method::Symbol = :sum,
)::@NamedTuple{ac::Array{Float64, 3}, l::Vector{Float64}}

    # validate lag bounds: must be non-negative and within the signal length
    @assert l <= size(obj, 2) "l must be ≤ $(size(obj, 2))."
    @assert l >= 0 "l must be ≥ 0."

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    if datatype(obj) == "erp"

        # epoch 1 is the pre-computed average; epochs 2:end are individual trials
        # compute per-trial auto-correlations first, then prepend the mean across
        # trials as epoch 1 of the output (preserving the ERP convention)
        ac = acov(
            @view(obj.data[ch, :, 2:end]),
            l = l,
            demean = demean,
            biased = biased,
            method = method,
        )
        ac = cat(mean(ac, dims = 3), ac, dims = 3)

    else
        ac = acov(
            @view(obj.data[ch, :, :]),
            l = l,
            demean = demean,
            biased = biased,
            method = method,
        )
    end

    return (ac = ac, l = collect((-l):l) .* 1 / sr(obj))

end
