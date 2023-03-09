export phases

"""
    phases(signal; pad)

Calculate phases.

# Arguments

- `signal::AbstractArray`
- `pad::Int64=0`: number of zeros to add

# Returns

- `phases::Vector{Float64}`
"""
function phases(signal::AbstractArray, pad::Int64=0)
    return angle.(pad0(signal, pad))
end
