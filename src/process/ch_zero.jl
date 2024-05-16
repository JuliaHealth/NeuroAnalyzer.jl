export ch_zero

"""
    ch_zero(obj)

Zero channels at the beginning and at the end.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function ch_zero(obj::NeuroAnalyzer.NEURO)

    obj_new = deepcopy(obj)
    obj_new.data[:, 1, :] .= 0
    obj_new.data[:, end, :] .= 0

    reset_components!(obj_new)
    push!(obj_new.history, "zero(OBJ)")

    return obj_new

end

"""
    ch_zero!(obj)

Zero channels at the beginning and at the end.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function ch_zero!(obj::NeuroAnalyzer.NEURO)

    obj_new = ch_zero(obj)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
