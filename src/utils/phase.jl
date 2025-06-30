export phases

"""
    phases(s)

Calculate phases.

# Arguments

- `s::AbstractVector`

# Returns

- `phases::Vector{Float64}`
"""
function phases(s::AbstractVector)::Vector{Float64}

    return DSP.angle.(hilbert(s))

end
