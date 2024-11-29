export scale
export scale!

"""
    scale(obj; <keyword arguments>)

Multiply channel(s) by `factor`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `factor::Real`: signal is multiplied by `factor`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function scale(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, factor::Real)::NeuroAnalyzer.NEURO

    ch = get_channel(obj, ch=ch)
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = @views obj_new.data[ch, :, :] .* factor
    reset_components!(obj_new)
    push!(obj_new.history, "scale(OBJ, ch=$ch, factor=$factor)")

    return obj_new

end

"""
    scale!(obj; <keyword arguments>)

Multiply channel(s) by `factor`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `factor::Real`: signal is multiplied by `factor`

# Returns

Nothing
"""
function scale!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, factor::Real)::Nothing

    obj_new = scale(obj, ch=ch, factor=factor)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end
