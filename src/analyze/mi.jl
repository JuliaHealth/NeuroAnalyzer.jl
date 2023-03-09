export mi

"""
    mi(signal1, signal2)

Calculate mutual information.

# Arguments

- `signal1::AbstractVector`
- `signal2::AbstractVector`

# Returns

- `mi::Float64`
"""
function mi(signal1::AbstractVector, signal2::AbstractVector)
    return get_mutual_information(signal1, signal2)
end

"""
    mi(obj; channel)

Calculate mutual information between channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `m::Array{Float64, 3}`
"""
function mi(obj::NeuroAnalyzer.NEURO; channel::Union{Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    m = zeros(ch_n, ch_n, ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        
        # create half of the matrix
        Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_idx1
                m[ch_idx1, ch_idx2, ep_idx] = @views mi(obj.data[channel[ch_idx1], :, ep_idx], obj.data[channel[ch_idx2], :, ep_idx])
            end
        end

        # copy to the other half
        Threads.@threads for ch_idx1 in 1:(ch_n - 1)
            for ch_idx2 in (ch_idx1 + 1):ch_n
                m[ch_idx1, ch_idx2, ep_idx] = @views mi[ch_idx2, ch_idx1, ep_idx]
            end
        end
    end

    return m
end

"""
    mi(obj1, obj2; channel1, channel2, epoch1, epoch2)

Calculate mutual information between two channels.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type]))`: index of channels, default is all signal channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type]))`: index of channels, default is all signal channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

- `mi::Array{Float64, 3}`
"""
function mi(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)))

    # check channels
    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    # check epochs
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("obj1 and obj2 epoch lengths must be equal."))

    ch_n = length(channel1)
    ep_n = length(epoch1)

    m = zeros(ch_n, ch_n, ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_n
                m[ch_idx1, ch_idx2, ep_idx] = @views mi(obj1.signals[channel1[ch_idx1], :, epoch1[ep_idx]], obj2.signals[channel2[ch_idx2], :, epoch2[ep_idx]])
            end
        end
    end

    return m
end