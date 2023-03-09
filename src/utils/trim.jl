export trim

"""
    trim(signal; segment)

Remove segment from the signal.

# Arguments

- `signal::AbstractVector`
- `segment::Tuple{Int64, Int64}`: segment (from, to) in samples

# Returns

- `s_trimmed::Vector{Float64}`
"""
function trim(signal::AbstractVector; segment::Tuple{Int64, Int64})
    _check_segment(signal, segment[1], segment[2])
    return vcat(signal[1:segment[1] - 1], signal[(segment[2] + 1):end])
end

"""
    trim(signal; segment)

Remove segment from the signal.

# Arguments

- `signal::AbstractArray`
- `segment::Tuple{Int64, Int64}`: segment (from, to) in samples

# Returns

- `s_trimmed::Array{Float64}`
"""
function trim(signal::AbstractArray; segment::Tuple{Int64, Int64})
    _check_segment(signal[1, :, 1], segment[1], segment[2])
    return hcat(signal[:, 1:(segment[1] - 1), :], signal[:, (segment[2] + 1):end, :])
end
