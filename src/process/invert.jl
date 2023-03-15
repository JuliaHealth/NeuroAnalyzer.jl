export invert_polarity
export invert_polarity!

"""
    invert_polarity(obj; ch)

Invert polarity.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function invert_polarity(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)))

    _check_channels(obj, ch)
    
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = .- obj_new.data[ch, :, :]
    reset_components!(obj_new)
    push!(obj_new.header.history, "invert_polarity(OBJ, ch=$ch)")

    return obj_new
    
end

"""
    invert_polarity!(obj; ch)

Invert polarity.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: channel(s) to invert
"""
function invert_polarity!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)))

    obj_tmp = invert_polarity(obj, ch=ch)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return nothing

end
