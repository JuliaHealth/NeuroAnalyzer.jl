export phases

"""
    phases(s; pad)

Calculate phases.

# Arguments

- `s::AbstractVector`
- `pad::Int64=0`: number of zeros to add

# Returns

- `phases::Vector{Float64}`
"""
function phases(s::AbstractVector, pad::Int64=0)

    return angle.(fft0(s, pad))

end
