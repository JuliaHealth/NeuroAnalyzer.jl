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
    obj_tmp = deepcopy(obj)

    for ep_idx = 2:nepochs(obj_tmp)
        obj_tmp.data[:, :, ep_idx] =
            @view(obj_tmp.data[:, :, ep_idx]) - @view(obj_tmp.data[:, :, 1])
    end
    push!(obj_tmp.history, "npl(obj)")

    return obj_tmp
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
    obj_tmp = npl(obj)
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj_tmp = nothing

    return nothing
end
