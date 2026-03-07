export acov

"""
    acov(s; <keyword arguments>)

Calculate autocovariance.

# Arguments

  - `s::AbstractVector`
  - `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
  - `demean::Bool=true`: demean signal before computing autocovariance
  - `biased::Bool=true`: calculate biased or unbiased autocovariance
  - `method::Symbol=:sum`: method of calculating autocovariance:
      + `:sum`: `acf = Σ(s[1:end - l] .* s[1+l:end])`
      + `:cov`: `acf = cov(s[1:end - l], s[1+l:end])`
      + `:stat`: use StatsBase `autocov()`, `biased` value is ignored

# Returns

  - `ac::Array{Float64, 3}`
"""
function acov(
    s::AbstractVector;
    l::Int64 = round(Int64, min(length(s) - 1, 10 * log10(length(s)))),
    demean::Bool = true,
    biased::Bool = true,
    method::Symbol = :sum,
)::Vector{Float64}

    _check_var(method, [:sum, :cov, :stat], "method")

    ac = zeros(l + 1)

    if method === :sum

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

Calculate autocovariance.

# Arguments

  - `s::AbstractArray`
  - `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `0:l`
  - `demean::Bool=true`: demean signal before computing autocovariance
  - `biased::Bool=true`: calculate biased or unbiased autocovariance
  - `method::Symbol=:sum`: method of calculating autocovariance:
      + `:sum`: `acf = Σ(s[1:end - l] .* s[1+l:end])`
      + `:cor`: `acf = cov(s[1:end - l], s[1+l:end])`
      + `:stat`: use StatsBase `autocor()`, `biased` value is ignored

# Returns

  - `ac::Array{Float64, 3}`
"""
function acov(
    s::AbstractArray;
    l::Int64 = round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1)))),
    demean::Bool = true,
    biased::Bool = true,
    method::Symbol = :sum,
)::Array{Float64, 3}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    ac = zeros(ch_n, length((-l):l), ep_n)

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

Calculate autocovariance. For ERP return trial-averaged autocovariance.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `l::Int64=1`: lags range is `0:lag` [samples]
  - `demean::Bool=true`: demean signal before computing autocovariance
  - `biased::Bool=true`: calculate biased or unbiased autocovariance
  - `method::Symbol=:sum`: method of calculating autocovariance:
      + `:sum`: `acf = Σ(s[1:end - l] .* s[1+l:end])`
      + `:cor`: `acf = cov(s[1:end - l], s[1+l:end])`
      + `:stat`: use StatsBase `autocor()`, `biased` value is ignored

# Returns

Named tuple containing:

  - `ac::Array{Float64, 3}`
  - `l::Vector{Float64}`: lags [s]
"""
function acov(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex},
        l::Real = 1,
        demean::Bool = true,
        biased::Bool = true,
        method::Symbol = :sum,
    )::@NamedTuple{ac::Array{Float64, 3}, l::Vector{Float64}}

    @assert l <= size(obj, 2) "l must be ≤ $(size(obj, 2))."
    @assert l >= 0 "l must be ≥ 0."

    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    if datatype(obj) == "erp"
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
