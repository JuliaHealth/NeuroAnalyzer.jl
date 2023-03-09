export acov

"""
   acov(signal; lag, demean, norm)

Calculate autocovariance.

# Arguments

- `signal::AbstractVector`
- `lag::Int64=1`: lags range is `-lag:lag`
- `demean::Bool=false`: demean `signal` prior to calculations
- `norm::Bool=false`: normalize autocovariance

# Returns

- `acov::Vector{Float64}`
- `lags::Vector{Int64}`
"""
function acov(signal::AbstractVector; lag::Int64=1, demean::Bool=false, norm::Bool=false)

    lag < 1 && throw(ArgumentError("lag must be ≥ 1."))
    lags = collect(-lag:lag)

    if demean == true
        s_demeaned = demean(signal)
    else
        s_demeaned = signal
    end

    acov = zeros(length(lags))
    l = length(signal)

    @inbounds @simd for idx in eachindex(lags)
        # no lag
        lags[idx] == 0 && (acov[idx] = sum(s_demeaned.^2))
        # positive lag
        lags[idx] > 0 && (acov[idx] = @views sum(s_demeaned[(1 + lags[idx]):end] .* s_demeaned[1:(end - lags[idx])]))
        # negative lag
        lags[idx] < 0 && (acov[idx] = @views sum(s_demeaned[1:(end - abs(lags[idx]))] .* s_demeaned[(1 + abs(lags[idx])):end]))
    end
    norm == true && (acov ./ l)

    return acov, lags
end

"""
   acov(obj; channel, lag=1, demean=false, norm=false)

Calculate autocovariance.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `lag::Int64`: lags range is `-lag:lag`
- `demean::Bool`: demean obj prior to analysis
- `norm::Bool`: normalize autocovariance

# Returns

Named tuple containing:
- `acov::Matrix{Float64}`
- `lags::Vector{Float64}`: lags in ms
"""
function acov(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), lag::Int64=1, demean::Bool=false, norm::Bool=false)

    lag < 1 && throw(ArgumentError("lag must be ≥ 1."))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    acov = zeros(ch_n, length(-lag:lag), ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            acov[ch_idx, :, ep_idx], _ = @views acov(obj.data[channel[ch_idx], :, ep_idx], lag=lag, demean=demean, norm=norm)
        end
    end

    lags = 1/sr(obj) .* collect(-lag:lag) .* 1000

    return (acov=acov, acov_lags=lags)
end
