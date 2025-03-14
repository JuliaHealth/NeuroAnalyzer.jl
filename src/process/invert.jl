export invert_polarity
export invert_polarity!

"""
    invert_polarity(obj; <keyword arguments>)

Invert polarity.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function invert_polarity(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::NeuroAnalyzer.NEURO

    ch = get_channel(obj, ch=ch)
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = .- obj_new.data[ch, :, :]
    reset_components!(obj_new)
    push!(obj_new.history, "invert_polarity(OBJ, ch=$ch)")

    return obj_new

end

"""
    invert_polarity!(obj; <keyword arguments>)

Invert polarity.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel(s) to invert

# Returns

Nothing
"""
function invert_polarity!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Nothing

    obj_new = invert_polarity(obj, ch=ch)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
