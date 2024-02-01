export acov

"""
   acov(s; l, demean, n, biased)

Calculate autocovariance.

# Arguments

- `s::AbstractVector`
- `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing autocovariance
- `biased::Bool=true`: calculate biased or unbiased autocovariance

# Returns

- `ac::Matrix{Float64}`
"""
function acov(s::AbstractVector; l::Int64=round(Int64, min(length(s) - 1, 10 * log10(length(s)))), demean::Bool=true, biased::Bool=true)

    ac = zeros(l + 1)

    ms = mean(s)
    if demean == true
        s_tmp = s .- ms
    else
        s_tmp = s
    end

    for idx in 0:l
        ac[idx + 1] = @views sum(s_tmp[1:(end - idx)] .* s_tmp[(1 + idx):end])
        if biased == true 
            ac[idx + 1] /= length(s)
        else
            ac[idx + 1] /= (length(s) - idx)
        end
    end
    
    ac = round.(ac, digits=8)
    ac = vcat(reverse(ac), ac[2:end])

    return reshape(ac, 1, :, 1)

end

"""
   acov(s; l, demean, n, biased)

Calculate autocovariance.

# Arguments

- `s::AbstractMatrix`
- `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing autocovariance
- `biased::Bool=true`: calculate biased or unbiased autocovariance

# Returns

- `ac::Matrix{Float64}`
"""
function acov(s::AbstractMatrix; l::Int64=round(Int64, min(size(s[:, 1], 1) - 1, 10 * log10(size(s[:, 1], 1)))), demean::Bool=true, biased::Bool=true)

    ep_n = size(s, 2)

    ac = zeros(1, length(-l:l), ep_n)

    @inbounds for ep_idx in 1:ep_n
        ac[1, :, ep_idx] = @views reshape(acov(s[:, ep_idx], l=l, demean=demean, biased=biased), 1, :, ep_n)
    end

    return ac

end

"""
   acov(s; l, demean, n, biased)

Calculate autocovariance.

# Arguments

- `s::AbstractArray`
- `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `0:l`
- `demean::Bool=true`: demean signal before computing autocovariance
- `biased::Bool=true`: calculate biased or unbiased autocovariance

# Returns

- `ac::Matrix{Float64}`
"""
function acov(s::AbstractArray; l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1)))), demean::Bool=true, biased::Bool=true)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    ac = zeros(ch_n, length(-l:l), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ac[ch_idx, :, ep_idx] = @views acov(s[ch_idx, :, ep_idx], l=l, demean=demean, biased=biased)
        end
    end

    return ac

end

"""
   acov(obj; ch, l, demean, n, biased)

Calculate autocovariance. For ERP return trial-averaged autocovariance.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `l::Int64=1`: lags range is `0:lag` [samples]
- `demean::Bool=true`: demean signal before computing autocovariance
- `biased::Bool=true`: calculate biased or unbiased autocovariance

# Returns

Named tuple containing:
- `ac::Array{Float64, 3}`
- `l::Vector{Float64}`: lags [s]
"""
function acov(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), l::Real=1, demean::Bool=true, biased::Bool=true)

    _check_channels(obj, ch)
    @assert l <= size(obj, 2) "l must be ≤ $(size(obj, 2))."
    @assert l >= 0 "l must be ≥ 0."

    if obj.header.recording[:data_type] == "erp"
        ac = @views acov(obj.data[ch, :, 2:end], l=l, demean=demean, biased=biased)
        ac = mean(ac, dims=3)
    else
        ac = @views acov(obj.data[ch, :, :], l=l, demean=demean, biased=biased)
    end

    return (ac=ac, l=round.(collect(-l:l) .* (1/sr(obj)), digits=5))

end
