export findpeaks

"""
    findpeaks(signal; <keyword arguments>)

Find the indices of local peaks in a 1-D signal.

Peaks are detected using `findpeaks1d` and a minimum separation of `d` samples is enforced between consecutive peaks (non-maximum suppression).

# Arguments

- `signal::AbstractVector`: signal vector
- `d::Int64=32`: minimum distance between consecutive peaks in samples; must be ≥ 1

# Returns

- `p_idx::Vector{Int64}`: indices of detected peaks in `signal`, sorted in ascending order

# Throws

- `ArgumentError`: if `d < 1`
"""
function findpeaks(signal::AbstractVector; d::Int64 = 32)::Vector{Int64}

    !(d >= 1) && throw(ArgumentError("d must be ≥ 1."))

    return findpeaks1d(signal, distance = d)[1]

end
