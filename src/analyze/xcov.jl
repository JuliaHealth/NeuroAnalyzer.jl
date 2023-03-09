export xcov

"""
   xcov(signal1, signal2; lag=1, demean=false, norm=false)

Calculate cross-covariance.

# Arguments

- `signal1::AbstractVector`
- `signal2::AbstractVector`
- `lag::Int64`: lags range is `-lag:lag`
- `demean::Bool`: demean signal prior to analysis
- `norm::Bool`: normalize cross-covariance

# Returns

- `ccov::Vector{Float64}`
- `lags::Vector{Int64}`
"""
function xcov(signal1::AbstractVector, signal2::AbstractVector; lag::Int64=1, demean::Bool=false, norm::Bool=false)

    length(signal1) == length(signal2) || throw(ArgumentError("Both signals must be of the same as length."))
    lag < 1 && throw(ArgumentError("lag must be ≥ 1."))
    lags = collect(-lag:lag)

    if demean == true
        signal1 = s_demean(signal1)
        signal2 = s_demean(signal2)
    end

    xcov = zeros(length(lags))
    l = length(signal1)

    @inbounds @fastmath @simd for idx in eachindex(lags)
        # no lag
        lags[idx] == 0 && (xcov[idx] = sum(signal1 .* signal2))
        # positive lag
        lags[idx] > 0 && (xcov[idx] = @views sum(signal1[(1 + lags[idx]):end] .* signal2[1:(end - lags[idx])]))
        # negative lag
        lags[idx] < 0 && (xcov[idx] = @views sum(signal1[1:(end - abs(lags[idx]))] .* signal2[(1 + abs(lags[idx])):end]))
    end
    norm == true && (xcov ./ l)

    return xcov, lags
end

"""
    xcov(obj; channel, lag, demean, norm)

Calculate cross-covariance.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `lag::Int64=1`: lags range is `-lag:lag`
- `demean::Bool=false`: demean signal prior to analysis
- `norm::Bool=false`: normalize cross-covariance

# Returns

Named tuple containing:
- `xcov::Matrix{Float64}`: ch1-ch1, ch1-ch2, ch1-ch3, etc.
- `lags::Vector{Float64}`: lags in ms
"""
function xcov(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), lag::Int64=1, demean::Bool=false, norm::Bool=false)

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    # create vector of lags
    lags = 1/sr(obj) .* collect(-lag:lag) .* 1000

    xcov = zeros(ch_n^2, length(lags), ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        
        # create half of the covariance matrix
        xcov_packed = Array{Vector{Float64}}(undef, ch_n, ch_n)
        Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_idx1
                xcov_packed[ch_idx1, ch_idx2], _ = @views xcov(obj.data[channel[ch_idx1], :, ep_idx], obj.data[channel[ch_idx2], :, ep_idx], lag=lag, demean=demean, norm=norm)
            end
        end
        
        # copy to the other half
        Threads.@threads for ch_idx1 in 1:(ch_n - 1)
            for ch_idx2 in (ch_idx1 + 1):ch_n
                xcov_packed[ch_idx1, ch_idx2] = @views xcov_packed[ch_idx2, ch_idx1]
            end
        end

        # unpack by channels
        Threads.@threads for ch_idx in 1:ch_n^2
            xcov[ch_idx, :, ep_idx] = @views xcov_packed[ch_idx]
        end
    end

    return (xcov=xcov, lags=lags)
end

"""
    xcov(obj1, obj2; channel1, channel2, epoch1, epoch2, lag, demean, norm)

Calculate cross-covariance between two NeuroAnalyzer NEURO objects.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type]))`: index of channels, default is all signal channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type]))`: index of channels, default is all signal channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs
- `lag::Int64=1`: lags range is `-lag:lag`
- `demean::Bool=false`: demean signal prior to analysis
- `norm::Bool=false`: normalize cross-covariance

# Returns

Named tuple containing:
- `xcov::Array{Float64, 3}`
- `lags::Vector{Float64}`: lags in ms
"""
function xcov(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)), lag::Int64=1, demean::Bool=false, norm::Bool=false)

    # check channels
    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    # check epochs
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("obj1 and obj2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)

    # create vector of lags
    lags = 1/sr(obj1) .* collect(-lag:lag) .* 1000

    xcov = zeros(length(channel1), (2 * lag + 1), length(epoch1))
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            xcov[ch_idx, :, ep_idx], _ = @views xcov(obj1.signals[channel1[ch_idx], :, epoch1[ep_idx]], obj2.signals[channel2[ch_idx], :, epoch2[ep_idx]], lag=lag, demean=demean, norm=norm)
        end
    end

    return (xcov=xcov, lags=lags)
end
