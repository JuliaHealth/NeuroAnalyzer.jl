export std

"""
    std(obj)

Calculate standard deviation of a NEURO object (along epochs; for ERP/ERF objects the averaged first epoch is excluded).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `s::Matrix{Float64}`
"""
function Statistics.std(obj::NeuroAnalyzer.NEURO)::Matrix{Float64}
    nepochs(obj) > 1 || throw(ArgumentError("OBJ must have > 1 epoch."))
    # ERP/ERF: the first epoch is the average
    datatype(obj) in ("erp", "erf") && nepochs(obj) < 3 &&
        throw(ArgumentError("ERP/ERF OBJ must have ≥ 2 epochs besides the average."))

    if datatype(obj) in ("erp", "erf")
        s = std(@view(obj.data[:, :, 2:end]); dims = 3)
    else
        s = std(@view(obj.data[:, :, :]); dims = 3)
    end
    s = reshape(s, size(s, 1), size(s, 2))

    return s
end
