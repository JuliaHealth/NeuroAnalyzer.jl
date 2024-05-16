export npl
export npl!

"""
    npl(obj)

Calculate non-phase-locked signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: must be ERP object

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function npl(obj::NeuroAnalyzer.NEURO)

    @assert datatype(obj) == "erp" "OBJ must be ERP."

    obj_new = deepcopy(obj)
    for ep_idx = 2:nepochs(obj_new)
        obj_new.data[:, :, ep_idx] = @views obj_new.data[:, :, ep_idx] - obj_new.data[:, :, 1]
    end
    reset_components!(obj_new)
    push!(obj_new.history, "npl(OBJ)")

    return obj_new

end

"""
    npl!(obj)

Calculate non-phase-locked signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: must be ERP object
"""
function npl!(obj::NeuroAnalyzer.NEURO)

    obj_new = npl(obj)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
