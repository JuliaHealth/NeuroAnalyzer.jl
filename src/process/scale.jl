export scale
export scale!

"""
    scale(obj; channel, factor)

Multiply channel(s) by `factor`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `factor::Real`: channel signal is multiplied by `factor`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function scale(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), factor::Real)

    _check_channels(obj, channel)
    
    obj_new = deepcopy(obj)
    obj_new.data[channel, :, :] = @views obj_new.data[channel, :, :] .* factor
    reset_components!(obj_new)
    push!(obj_new.header.history, "scale(OBJ, channel=$channel)")

    return obj_new
end

"""
    scale!(obj; channel, factor)

Multiply channel(s) by `factor`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `factor::Real`: channel signal is multiplied by `factor`
"""
function scale!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), factor::Real)

    obj_tmp = scale(obj, channel=channel, factor=factor)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return nothing
end
