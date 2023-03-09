"""
    zero(obj)

Zero channels at the beginning and at the end.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function zero(obj::NeuroAnalyzer.NEURO)

    obj_new = deepcopy(obj)
    obj_new.data[:, 1, :] .= 0
    obj_new.data[:, end, :] .= 0

    reset_components!(obj_new)
    push!(obj_new.header.history, "zero(OBJ)")

    return obj_new
end

"""
    zero!(obj)

Zero channels at the beginning and at the end.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function zero!(obj::NeuroAnalyzer.NEURO)

    obj_tmp = zero(obj)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    reset_components!(obj)

    return nothing
end
