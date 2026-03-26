export std

"""
    std(obj)

Calculate standard deviation of a NEURO object (along epochs).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `s::Matrix{Float64}`
"""
function Statistics.std(obj::NeuroAnalyzer.NEURO)::Matrix{Float64}

    nepochs(obj) > 1 || throw(ArgumentError("OBJ must have > 1 epoch."))

    if datatype(obj) == "erp"
        s = std(@view(obj.data[:, :, 2:end]), dims = 3)
    else
        s = std(@view(obj.data[:, :, :]), dims = 3)
    end
    s = reshape(s, size(s, 1), size(s, 2))

    return s

end
