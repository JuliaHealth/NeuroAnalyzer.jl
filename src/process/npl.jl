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
    obj_new = deepcopy(obj)
    npl!(obj_new)

    return obj_new
end

"""
    npl!(obj)

Calculate non-phase-locked signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be ERP/ERF object

# Returns

- `Nothing`
"""
function npl!(obj::NeuroAnalyzer.NEURO)::Nothing
    # validate
    datatype(obj) in ["erp", "erf"] || throw(ArgumentError("OBJ must be ERP/ERF."))

    for ep_idx = 2:nepochs(obj)
        obj.data[:, :, ep_idx] =
            @view(obj.data[:, :, ep_idx]) - @view(obj.data[:, :, 1])
    end
    push!(obj.history, "npl(obj)")

    return nothing
end
