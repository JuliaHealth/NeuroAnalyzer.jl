export acor

"""
    acor(s; <keyword arguments>)

Calculate auto-correlation.

# Arguments

- `s::AbstractVector`
- `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing auto-correlation
- `biased::Bool=true`: calculate biased or unbiased autocovariance
- `method::Symbol=:sum`: method of calculating auto-correlation:
    - `:sum`: `acf = Σ(s[1:end - l] .* s[1+l:end]) ./ var(s)`
    - `:cor`: `acf = cor(s[1:end - l], s[1+l:end])`
    - `:stat`: use StatsBase `autocor()`, `biased` value is ignored

# Returns

- `ac::Array{Float64, 3}`
"""
function acor(s::AbstractVector; l::Int64=round(Int64, min(length(s) - 1, 10 * log10(length(s)))), demean::Bool=true, biased::Bool=true, method::Symbol=:sum)::Array{Float64, 3}

    _check_var(method, [:sum, :cor, :stat], "method")

    ac = zeros(l + 1)

    if method === :sum
        ac = acov(s, l=l, demean=demean, biased=biased, method=:sum)[1, :, 1]
        # normalize by the variance of s
        ac = round.(ac ./ StatsBase.var(s), digits=3)
    elseif method === :cor
        ac = acov(s, l=l, demean=demean, biased=biased, method=:cov)[1, :, 1]
        # normalize by the variance of s
        ac = round.(ac ./ StatsBase.var(s), digits=3)
    elseif method === :stat
        ac = StatsBase.autocor(s, 0:l, demean=demean)
        ac = round.(ac, digits=3)
        ac = vcat(reverse(ac), ac[2:end])
    end

    return reshape(ac, 1, :, 1)

end

"""
    acor(s; <keyword arguments>)

Calculate auto-correlation.

# Arguments

- `s::AbstractMatrix`
- `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing auto-correlation
- `biased::Bool=true`: calculate biased or unbiased autocovariance
- `method::Symbol=:sum`: method of calculating auto-correlation:
    - `:sum`: `acf = Σ(s[1:end - l] .* s[1+l:end]) ./ var(s)`
    - `:cor`: `acf = cor(s[1:end - l], s[1+l:end])`
    - `:stat`: use StatsBase `autocor()`, `biased` value is ignored

# Returns

- `ac::Array{Float64, 3}`
"""
function acor(s::AbstractMatrix; l::Int64=round(Int64, min(size(s[:, 1], 1) - 1, 10 * log10(size(s[:, 1], 1)))), demean::Bool=true, biased::Bool=true, method::Symbol=:sum)::Array{Float64, 3}

    ep_n = size(s, 2)
    ac = zeros(1, length(-l:l), ep_n)

    @inbounds for ep_idx in 1:ep_n
        ac[1, :, ep_idx] = @views acor(s[:, ep_idx], l=l, demean=demean, biased=biased, method=method)
    end

    return ac

end

"""
    acor(s; <keyword arguments>)

Calculate auto-correlation.

# Arguments

- `s::AbstractArray`
- `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing auto-correlation
- `biased::Bool=true`: calculate biased or unbiased autocovariance
- `method::Symbol=:sum`: method of calculating auto-correlation:
    - `:sum`: `acf = Σ(s[1:end - l] .* s[1+l:end]) ./ var(s)`
    - `:cor`: `acf = cor(s[1:end - l], s[1+l:end])`
    - `:stat`: use StatsBase `autocor()`, `biased` value is ignored

# Returns

- `ac::Array{Float64, 3}`
"""
function acor(s::AbstractArray; l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1)))), demean::Bool=true, biased::Bool=true, method::Symbol=:sum)::Array{Float64, 3}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    ac = zeros(ch_n, length(-l:l), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            ac[ch_idx, :, ep_idx] = @views acor(s[ch_idx, :, ep_idx], l=l, demean=demean, biased=biased, method=method)
        end
    end

    return ac

end

"""
    acor(obj; <keyword arguments>)

Calculate auto-correlation. For ERP return trial-averaged auto-correlation.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `l::Real=1`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing auto-correlation
- `biased::Bool=true`: calculate biased or unbiased autocovariance
- `method::Symbol=:sum`: method of calculating auto-correlation:
    - `:sum`: `acf = Σ(s[1:end - l] .* s[1+l:end]) ./ var(s)`
    - `:cor`: `acf = cor(s[1:end - l], s[1+l:end])`
    - `:stat`: use StatsBase `autocor()`, `biased` value is ignored

# Returns

Named tuple containing:
- `ac::Array{Float64, 3}`
- `l::Vector{Float64}`: lags [s]
"""
function acor(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, l::Real=1, demean::Bool=true, biased::Bool=true, method::Symbol=:sum)::@NamedTuple{ac::Array{Float64, 3}, l::Vector{Float64}}

    @assert l <= size(obj, 2) "l must be ≤ $(size(obj, 2))."
    @assert l >= 0 "l must be ≥ 0."

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")

    if datatype(obj) == "erp"
        ac = @views acor(obj.data[ch, :, 2:end], l=l, demean=demean, biased=biased, method=method)
        ac = cat(mean(ac, dims=3), ac, dims=3)
    else
        ac = @views acor(obj.data[ch, :, :], l=l, demean=demean, biased=biased, method=method)
    end

    return (ac=ac, l=collect(-l:l) .* 1/sr(obj))

end
