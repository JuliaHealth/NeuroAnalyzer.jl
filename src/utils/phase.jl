export phases

"""
    phases(s; pad)

Calculate phases.

# Arguments

- `s::AbstractVector`

# Returns

- `phases::Vector{Float64}`
"""
function phases(s::AbstractVector)

    return angle.(hilbert(s))

end
