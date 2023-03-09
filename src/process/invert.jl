export invert_polarity
export invert_polarity!

"""
    invert_polarity(obj; channel)

Invert polarity.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function invert_polarity(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)))

    _check_channels(obj, channel)
    
    obj_new = deepcopy(obj)
    obj_new.data[channel, :, :] = .- obj_new.data[channel, :, :]
    reset_components!(obj_new)
    push!(obj_new.header.history, "invert_polarity(EEG, channel=$channel)")

    return obj_new
end

"""
    invert_polarity!(obj; channel)

Invert polarity.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: channel(s) to invert
"""
function invert_polarity!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)))

    obj_tmp = invert_polarity(obj, channel=channel)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    reset_components!(obj)

    return nothing
end
