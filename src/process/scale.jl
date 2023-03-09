export scale
export scale!

"""
    scale(eeg; channel, factor)

Multiply EEG channel(s) by `factor`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(eeg))`: index of channels, default is all channels
- `factor::Real`: channel signal is multiplied by `factor`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function scale(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(eeg)), factor::Real)

    _check_channels(eeg, channel)
    
    obj_new = deepcopy(eeg)
    obj_new.data[channel, :, :] = @views obj_new.data[channel, :, :] .* factor
    reset_components!(obj_new)
    push!(obj_new.header.history, "scale(EEG, channel=$channel)")

    return obj_new
end

"""
    scale!(eeg; channel, factor)

Multiply EEG channel(s) by `factor`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(eeg))`: index of channels, default is all channels
- `factor::Real`: channel signal is multiplied by `factor`
"""
function scale!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(eeg)), factor::Real)

    eeg_tmp = scale(eeg, channel=channel, factor=factor)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return nothing
end
