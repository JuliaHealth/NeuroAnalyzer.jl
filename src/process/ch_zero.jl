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
    push!(obj_new.header.history, "zero(OBJ)")

    return obj_new
end

"""
    ch_zero!(obj)

Zero channels at the beginning and at the end.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function ch_zero!(obj::NeuroAnalyzer.NEURO)

    obj_tmp = ch_zero(obj)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return nothing
end
