export acor

"""
   acor(s; l, demean)

Calculate auto-correlation.

# Arguments

- `s::AbstractVector`
- `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing auto-correlation
- `biased::Bool=true`: calculate biased or unbiased autocovariance

# Returns

- `ac::Matrix{Float64}`
"""
function acor(s::AbstractVector; l::Int64=round(Int64, min(length(s) - 1, 10 * log10(length(s)))), demean::Bool=true, biased::Bool=true)

    if biased == true
        ac = autocor(s, 0:l, demean=demean)
    else
        ac = zeros(l + 1)
        ms = mean(s)
        if demean == true
            s_tmp = s .- ms
        else
            s_tmp = s
        end
        for idx in 0:l
            s1 = s_tmp[1:(end - idx)]
            s2 = s_tmp[(1 + idx):end]
            ac[idx + 1] = sum(s1 .* s2) / (length(s) - idx)
        end
        ac ./= std(s)^2
        ac = round.(ac, digits=8)
    end

    ac = vcat(reverse(ac), ac[2:end])

    return reshape(ac, 1, :, 1)

end

"""
   acor(s; l, demean)

Calculate auto-correlation.

# Arguments

- `s::AbstractMatrix`
- `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing auto-correlation
- `biased::Bool=true`: calculate biased or unbiased autocovariance

# Returns

- `ac::Matrix{Float64}`
"""
function acor(s::AbstractMatrix; l::Int64=round(Int64, min(size(s[:, 1], 1) - 1, 10 * log10(size(s[:, 1], 1)))), demean::Bool=true, biased::Bool=true)

    ep_n = size(s, 2)

    ac = zeros(1, length(-l:l), ep_n)

    @inbounds for ep_idx in 1:ep_n
        ac[1, :, ep_idx] = @views reshape(acor(s[:, ep_idx], l=l, demean=demean, biased=biased), 1, :, ep_n)
    end

    return ac

end

"""
   acor(s; l, demean)

Calculate auto-correlation.

# Arguments

- `s::AbstractArray`
- `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing auto-correlation
- `biased::Bool=true`: calculate biased or unbiased autocovariance

# Returns

- `ac::Matrix{Float64}`
"""
function acor(s::AbstractArray; l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1)))), demean::Bool=true, biased::Bool=true)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    ac = zeros(ch_n, length(-l:l), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ac[ch_idx, :, ep_idx] = @views acor(s[ch_idx, :, ep_idx], l=l, demean=demean, biased=biased)
        end
    end

    return ac

end

"""
   acor(obj; ch, lag, demean)

Calculate auto-correlation.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `l::Real=1`: lags range is `-l:l`
- `demean::Bool=true`: demean signal before computing auto-correlation
- `biased::Bool=true`: calculate biased or unbiased autocovariance

# Returns

Named tuple containing:
- `ac::Array{Float64, 3}`
- `l::Vector{Float64}`: lags [s]
"""
function acor(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), l::Real=1, demean::Bool=true, biased::Bool=true)

    _check_channels(obj, ch)
    @assert l <= size(obj, 2) "l must be ≤ $(size(obj, 2))."
    @assert l >= 0 "l must be ≥ 0."

    ac = @views acor(obj.data[ch, :, :], l=l, demean=demean, biased=biased)

    return (ac=ac, l=round.(collect(-l:l) .* (1/sr(obj)), digits=5))

end
