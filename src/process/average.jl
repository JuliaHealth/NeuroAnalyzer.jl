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
    average(eeg; channel)

Return the average signal of channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(eeg))`: index of channels, default is all channels

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function average(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(eeg)))

    _check_channels(eeg, channel)

    obj_new = deepcopy(eeg)
    keep_channel!(obj_new, channel=1)
    obj_new.data = @views s_average(obj.data[channel, :, :])
    obj_new.header.recording.labels=["averaged channel"]
    reset_components!(obj_new)
    push!(obj_new.header.history, "average(OBJ, channel=$channel)")

    return obj_new
end

"""
    average!(eeg; channel)

Return the average signal of channels.  

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(eeg))`: index of channels, default is all channels
"""
function average!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(eeg)))

    obj_tmp = average(eeg, channel=channel)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    reset_components!(eeg)

    return nothing
end

"""
    average(eeg1, eeg2)

Return the average signal of all `eeg1` and `eeg2` channels.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function average(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO)

    size(eeg1.data) == size(eeg2.data) || throw(ArgumentError("Both signals must have the same size."))
    ch_n = channel_n(eeg1)
    ep_n = epoch_n(eeg1)

    obj_new = deepcopy(eeg1)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            obj_new.data[ch_idx, :, ep_idx] = @views s2_average(eeg1.data[ch_idx, :, ep_idx], eeg2.data[ch_idx, :, ep_idx])
        end
    end

    reset_components!(obj_new)
    push!(obj.header[:history], "average(OBJ1, OBJ2)")

    return obj_new
end

