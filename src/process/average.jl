export average

"""
    average(signal)

Average all channels.

# Arguments

- `signal::AbstractArray`

# Returns

- `s_averaged::AbstractArray`
"""
function average(signal::AbstractArray)
    return mean(signal, dims=1)
end

"""
    average(signal1, signal2)

Averages `signal1` and `signal2`.

# Arguments

- `signal1::AbstractArray`
- `signal2::AbstractArray`

# Returns

- `s_averaged::Vector{Float64}`
"""
function average(signal1::AbstractArray, signal2::AbstractArray)
    return mean(hcat(signal1, signal2), dims=2)
end

"""
    average(obj; channel)

Return the average signal of channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function average(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)))

    _check_channels(obj, channel)

    obj_new = deepcopy(obj)
    keep_channel!(obj_new, channel=1)
    obj_new.data = @views s_average(obj.data[channel, :, :])
    obj_new.header.recording.labels=["averaged channel"]
    reset_components!(obj_new)
    push!(obj_new.header.history, "average(OBJ, channel=$channel)")

    return obj_new
end

"""
    average!(obj; channel)

Return the average signal of channels.  

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
"""
function average!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)))

    obj_tmp = average(obj, channel=channel)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    reset_components!(obj)

    return nothing
end

"""
    average(obj1, obj2)

Return the average signal of all `obj1` and `obj2` channels.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function average(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO)

    size(obj1.data) == size(obj2.data) || throw(ArgumentError("Both signals must have the same size."))
    ch_n = channel_n(obj1)
    ep_n = epoch_n(obj1)

    obj_new = deepcopy(obj1)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            obj_new.data[ch_idx, :, ep_idx] = @views s2_average(obj1.data[ch_idx, :, ep_idx], obj2.data[ch_idx, :, ep_idx])
        end
    end

    reset_components!(obj_new)
    push!(obj.header.history, "average(OBJ1, OBJ2)")

    return obj_new
end

