export acov

"""
   acov(s; l, demean)

Calculate auto-covariance.

# Arguments

- `s::AbstractVector`
- `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `0:l`
 - `demean::Bool=true`: demean signal before computing auto-covariance

# Returns

- `ac::Matrix{Float64}`
"""
function acov(s::AbstractVector; l::Int64=round(Int64, min(length(s) - 1, 10 * log10(length(s)))), demean::Bool=true)

    return reshape(autocov(s, 0:l, demean=demean), 1, :, 1)

end

"""
   acov(s; l, demean)

Calculate auto-covariance.

# Arguments

- `s::AbstractMatrix`
- `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `0:l`
 - `demean::Bool=true`: demean signal before computing auto-covariance

# Returns

- `ac::Matrix{Float64}`
"""
function acov(s::AbstractMatrix; l::Int64=round(Int64, min(size(s[:, 1], 1) - 1, 10 * log10(size(s[:, 1], 1)))), demean::Bool=true)

    ep_n = size(s, 2)

    ac = zeros(1, length(0:l), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        ac[1, :, ep_idx] = @views reshape(autocov(s[:, ep_idx], 0:l, demean=demean), 1, :, ep_n)
    end

    return ac

end

"""
   acov(s; l, demean)

Calculate auto-covariance.

# Arguments

- `s::AbstractArray`
- `l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1))))`: lags range is `0:l`
 - `demean::Bool=true`: demean signal before computing auto-covariance

# Returns

- `ac::Matrix{Float64}`
"""
function acov(s::AbstractArray; l::Int64=round(Int64, min(size(s[1, :, 1], 1) - 1, 10 * log10(size(s[1, :, 1], 1)))), demean::Bool=true)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    ac = zeros(ch_n, length(0:l), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ac[ch_idx, :, ep_idx] = @views autocov(s[ch_idx, :, ep_idx], 0:l, demean=demean)
        end
    end

    return ac

end

"""
   acov(obj; ch, lag, demean)

Calculate auto-covariance.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `lag::Real=1`: lags range is `0:lag` [s]
- `demean::Bool=true`: demean signal before computing auto-covariance

# Returns

Named tuple containing:
- `ac::Array{Float64, 3}`
- `lag::Vector{Float64}`: lags [s]
"""
function acov(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), lag::Real=1, demean::Bool=true)

    @assert lag <= obj.epoch_time[end] "lag must be ≤ $(obj.epoch_time[end])."
    _check_channels(obj, ch)
    @assert lag >= 0 "lag must be ≥ 0."

    l = vsearch(lag, obj.time_pts)
    ac = @views acov(obj.data[ch, :, :], l=(l - 1), demean=demean)

    return (ac=ac, lag=obj.epoch_time[1:l])

end
