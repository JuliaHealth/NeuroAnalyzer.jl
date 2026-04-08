export npl
export npl!

"""
    npl(obj)

Calculate non-phase-locked signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be ERP/ERF object

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function npl(obj::NeuroAnalyzer.NEURO)::NeuroAnalyzer.NEURO
    # validate
    datatype(obj) in ["erp", "erf"] || throw(ArgumentError("OBJ must be ERP/ERF."))

    # create new dataset
    obj_new = deepcopy(obj)

    for ep_idx in 2:nepochs(obj_new)
        obj_new.data[:, :, ep_idx] =
            @view(obj_new.data[:, :, ep_idx]) - @view(obj_new.data[:, :, 1])
    end
    push!(obj_new.history, "npl(obj)")

    return obj_new
end

"""
    npl!(obj)

Calculate non-phase-locked signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be ERP object

# Returns

- `Nothing`
"""
function npl!(obj::NeuroAnalyzer.NEURO)::Nothing
    obj_new = npl(obj)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing
end
