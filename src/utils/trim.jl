export trim

"""
    trim(s; seg)

Remove segment from the signal.

# Arguments

- `v::AbstractVector`
- `seg::Tuple{Int64, Int64}`: segment (from, to) in samples

# Returns

- `trim::Vector{Float64}`
"""
function trim(v::AbstractVector; seg::Tuple{Int64, Int64})

    _check_segment(v, seg[1], seg[2])
    
    return vcat(v[1:seg[1] - 1], v[(seg[2] + 1):end])

end

"""
    trim(m; seg)

Remove segment from the signal.

# Arguments

- `m::AbstractMatrix`
- `seg::Tuple{Int64, Int64}`: segment (from, to) in samples

# Returns

- `trim::Array{Float64}`
"""
function trim(m::AbstractMatrix; seg::Tuple{Int64, Int64})
    
    _check_segment(m[1, :], seg[1], seg[2])

    return hcat(m[:, 1:(seg[1] - 1)], m[:, (seg[2] + 1):end])

end

"""
    trim(a; seg)

Remove segment from the signal.

# Arguments

- `a::AbstractArray`
- `seg::Tuple{Int64, Int64}`: segment (from, to) in samples

# Returns

- `trim::Array{Float64}`
"""
function trim(a::AbstractArray; seg::Tuple{Int64, Int64})

    _check_segment(a[1, :, 1], seg[1], seg[2])

    return hcat(a[:, 1:(seg[1] - 1), :], a[:, (seg[2] + 1):end, :])

end
