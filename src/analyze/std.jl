export std

"""
    std(obj)

Calculate standard deviation of the signal data (along epochs).

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `s::Matrix{Float64}`
"""
function Statistics.std(obj::NeuroAnalyzer.NEURO)::Matrix{Float64}

    @assert nepochs(obj) > 1 "OBJ must have > 1 epoch."

    if datatype(obj) == "erp"
        s = @views std(obj.data[:, :, 2:end], dims=3)
    else
        s = @views std(obj.data[:, :, :], dims=3)
    end
    s = reshape(s, size(s, 1), size(s, 2))

    return s

end
