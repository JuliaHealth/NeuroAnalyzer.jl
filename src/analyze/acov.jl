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
    - `:sum`: `acf = Σ(s[1:end - l] .* s[1+l:end])`
    - `:cov`: `acf = cov(s[1:end - l], s[1+l:end])`
    - `:stat`: use StatsBase `autocov()`, `biased` value is ignored

# Returns

- `ac::Array{Float64, 3}`
"""
function acov(s::AbstractVector; l::Int64=round(Int64, min(length(s) - 1, 10 * log10(length(s)))), demean::Bool=true, biased::Bool=true, method::Symbol=:sum)::Array{Float64, 3}

    _check_var(method, [:sum, :cov, :stat], "method")

    ac = zeros(l + 1)

    if demean
        s_tmp = delmean(s)
    else
        s_tmp = s
    end

    if method === :sum
        for idx in 0:l
            ac[idx + 1] = @views sum(s_tmp[1:(end - idx)] .* s_tmp[(1 + idx):end])
            if biased
                ac[idx + 1] /= length(s)
            else
                ac[idx + 1] /= (length(s) - idx)
            end
        end
    elseif method === :cov
        for idx in 0:l
            if !biased
                ac[idx + 1] = @views cov(s_tmp[1:(end - idx)], s_tmp[(1 + idx):end], corrected=true)
            else
                ac[idx + 1] = @views cov(s_tmp[1:(end - idx)], s_tmp[(1 + idx):end], corrected=false)
            end
        end
    elseif method === :stat
        ac = StatsBase.autocov(s, 0:l, demean=demean)
    end

    ac = round.(ac, digits=3)
    ac = vcat(reverse(ac), ac[2:end])

    return reshape(ac, 1, :, 1)

end

"""
    acov(s; <keyword arguments>)

Calculate autocovariance.

# Arguments

- `s::AbstractMatrix`
- `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing autocovariance
- `biased::Bool=true`: calculate biased or unbiased autocovariance
- `method::Symbol=:sum`: method of calculating autocovariance:
    - `:sum`: `acf = Σ(s[1:end - l] .* s[1+l:end])`
    - `:cor`: `acf = cov(s[1:end - l], s[1+l:end])`
    - `:stat`: use StatsBase `autocor()`, `biased` value is ignored

# Returns

- `ac::Array{Float64, 3}`
"""
function acov(s::AbstractMatrix; l::Int64=round(Int64, min(size(s[:, 1], 1) - 1, 10 * log10(size(s[:, 1], 1)))), demean::Bool=true, biased::Bool=true, method::Symbol=:sum)::Array{Float64, 3}

    ep_n = size(s, 2)

    ac = zeros(1, length(-l:l), ep_n)

    @inbounds for ep_idx in 1:ep_n
        ac[1, :, ep_idx] = @views reshape(acov(s[:, ep_idx], l=l, demean=demean, biased=biased, method=method), 1, :, ep_n)
    end

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
    - `:sum`: `acf = Σ(s[1:end - l] .* s[1+l:end])`
    - `:cor`: `acf = cov(s[1:end - l], s[1+l:end])`
    - `:stat`: use StatsBase `autocor()`, `biased` value is ignored

# Returns

- `ac::Array{Float64, 3}`
"""
function acov(s::AbstractArray; l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1)))), demean::Bool=true, biased::Bool=true, method::Symbol=:sum)::Array{Float64, 3}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    ac = zeros(ch_n, length(-l:l), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            ac[ch_idx, :, ep_idx] = @views acov(s[ch_idx, :, ep_idx], l=l, demean=demean, biased=biased, method=method)
        end
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
    - `:sum`: `acf = Σ(s[1:end - l] .* s[1+l:end])`
    - `:cor`: `acf = cov(s[1:end - l], s[1+l:end])`
    - `:stat`: use StatsBase `autocor()`, `biased` value is ignored

# Returns

Named tuple containing:
- `ac::Array{Float64, 3}`
- `l::Vector{Float64}`: lags [s]
"""
function acov(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, l::Real=1, demean::Bool=true, biased::Bool=true, method::Symbol=:sum)::@NamedTuple{ac::Array{Float64, 3}, l::Vector{Float64}}

    @assert l <= size(obj, 2) "l must be ≤ $(size(obj, 2))."
    @assert l >= 0 "l must be ≥ 0."

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")

    if datatype(obj) == "erp"
        ac = @views acov(obj.data[ch, :, 2:end], l=l, demean=demean, biased=biased, method=method)
        ac = cat(mean(ac, dims=3), ac, dims=3)
    else
        ac = @views acov(obj.data[ch, :, :], l=l, demean=demean, biased=biased, method=method)
    end

    return (ac=ac, l=collect(-l:l) .* 1/sr(obj))

end
