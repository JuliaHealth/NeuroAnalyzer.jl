export findpeaks

"""
    findpeaks(signal; d)

Find peaks.

# Arguments

- `signal::AbstractVector`
- `d::Int64=32`: distance between peeks in samples

# Returns

- `p_idx::Vector{Int64}`

"""
function findpeaks(signal::AbstractVector; d::Int64=32)
    d < 1 && throw(ArgumentError("d must be ≥ 1."))
    return findpeaks1d(signal, distance=d)[1]
end

