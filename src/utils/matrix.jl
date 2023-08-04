export m_pad0
export m_sortperm
export m_sort
export m_norm

"""
    m_pad0(m)

Pad matrix with zeros to make it square.

# Arguments

- `m::Matrix{<:Number}`

# Returns

- `m::Matrix{<:Number}`
"""
function m_pad0(m::Matrix{<:Number})

    nr, nc = size(m)

    if nr > nc
        # horizontal
        return hcat(m, repeat(zeros(eltype(m), 1), nr, nr - nc))
    elseif nr < nc
        # vertical
        return vcat(m, repeat(zeros(eltype(m), 1), nc - nr, nc))
    else
        return m
    end

end

"""
    m_sortperm(m; rev=false, dims=1)

Generates sorting index for matrix `m` by columns (`dims` = 1) or by rows (`dims` = 2).

# Arguments

- `m::AbstractMatrix`
- `rev::Bool`
- `dims::Int64`

# Returns

- `idx::Matrix{Int64}`
"""
function m_sortperm(m::AbstractMatrix; rev::Bool=false, dims::Int64=1)

    @assert dims in [1, 2] "dims must be 1 or 2."
    
    idx = zeros(Int, size(m))
    if dims == 1
        @inbounds @simd for m_idx = 1:size(m, 2)
            # sort by columns
            idx[:, m_idx] = sortperm(m[:, m_idx], rev=rev)
        end
    else
        @inbounds @simd for m_idx = 1:size(m, 1)
            # sort by rows
            idx[m_idx, :] = sortperm(m[m_idx, :], rev=rev)'
        end     
    end

    return idx

end

"""
    m_sort(m, m_idx; rev=false, dims=1)

Sorts matrix `m` using sorting index `m_idx`

# Arguments

- `m::Matrix`
- `m_idx::Vector{Int64}`
- `rev::Bool=false`
- `dims::Int64=1`: sort by columns (`dims=1`) or by rows (`dims=2`)

# Returns

- `m_sorted::Matrix`
"""
function m_sort(m::Matrix, m_idx::Vector{Int64}; rev::Bool=false, dims::Int64=1)

    @assert dims in [1, 2] "dims must be 1 or 2."

    rev == true && reverse!(m_idx)

    m_sorted = zeros(eltype(m), size(m))
    if dims == 1
        @inbounds @simd for idx = 1:size(m, 2)
            # sort by columns
            m_sorted[:, idx] = @views m[:, idx][m_idx]
        end
    else
        @inbounds @simd for idx = 1:size(m, 1)
            # sort by rows
            m_sorted[idx, :] = @views m[idx, :][m_idx]
        end
    end

    return m_sorted

end

"""
    m_norm(m)

Normalize matrix.

# Arguments

- `m::AbstractArray`

# Returns

- `m_norm::AbstractArray`
"""
function m_norm(m::AbstractArray)

    return m ./ (size(m, 2) - 1)

end
