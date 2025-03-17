export m_pad0
export m_sortperm
export m_sort
export m_norm
export vec2mat
export arr2mat

"""
    m_pad0(m)

Pad matrix with zeros to make it square.

# Arguments

- `m::AbstractMatrix`

# Returns

- `m::AbstractMatrix`
"""
function m_pad0(m::AbstractMatrix)::AbstractMatrix

    mr, mc = size(m)

    if mr > mc
        # horizontal
        return hcat(m, repeat(zeros(eltype(m), 1), mr, mr - mc))
    elseif mr < mc
        # vertical
        return vcat(m, repeat(zeros(eltype(m), 1), mc - mr, mc))
    else
        return m
    end

end

"""
    m_pad0(m, r, c)

Pad matrix with zeros to make its size r × c.

# Arguments

- `m::AbstractMatrix`
- `r::Int64`: number of rows
- `c::Int64`: number of columns

# Returns

- `m::AbstractMatrix`
"""
function m_pad0(m::AbstractMatrix, r::Int64, c::Int64)::AbstractMatrix

    @assert (r, c) > size(m) "New matrix dimensions ($r × $c) must be larger than the source matrix dimensions ($(size(m, 1)) × $(size(m, 2)))."

    mr, mc = size(m)

    # add rows
    if r > mr
        m = vcat(m, repeat(zeros(eltype(m), 1), r - mr, mc))
    end

    # add columns
    if c > mc
        m = hcat(m, repeat(zeros(eltype(m), 1), r, c - mc))
    end

    return m

end


"""
    m_sortperm(m; <keyword arguments>)

Generates matrix sorting index.

# Arguments

- `m::AbstractMatrix`
- `rev::Bool`: reverse sort
- `dims::Int64=1`: sort by columns (`dims=1`) or by rows (`dims=2`)

# Returns

- `idx::Matrix{Int64}`
"""
function m_sortperm(m::AbstractMatrix; rev::Bool=false, dims::Int64=1)::AbstractMatrix

    @assert dims in [1, 2] "dims must be 1 or 2."

    idx = zeros(Int, size(m))
    if dims == 1
        @inbounds for m_idx = axes(m, 2)
            # sort by columns
            idx[:, m_idx] = sortperm(m[:, m_idx], rev=rev)
        end
    else
        @inbounds for m_idx = axes(m, 1)
            # sort by rows
            idx[m_idx, :] = sortperm(m[m_idx, :], rev=rev)'
        end
    end

    return idx

end

"""
    m_sort(m, m_idx; <keyword arguments>)

Sorts matrix using sorting index.

# Arguments

- `m::AbstractMatrix`
- `m_idx::Vector{Int64}`: sorting index
- `rev::Bool=false`: reverse sort
- `dims::Int64=1`: sort by columns (`dims=1`) or by rows (`dims=2`)

# Returns

- `m_sorted::AbstractMatrix`
"""
function m_sort(m::AbstractMatrix, m_idx::Vector{Int64}; rev::Bool=false, dims::Int64=1)::AbstractMatrix

    @assert dims in [1, 2] "dims must be 1 or 2."

    rev && reverse!(m_idx)

    m_sorted = zeros(eltype(m), size(m))
    if dims == 1
        @inbounds for idx = axes(m, 2)
            # sort by columns
            m_sorted[:, idx] = @views m[:, idx][m_idx]
        end
    else
        @inbounds for idx = axes(m, 1)
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
function m_norm(m::AbstractArray)::AbstractArray

    return m ./ (size(m, 2) - 1)

end

"""
    vec2mat(x; <keyword arguments>)

Reshape vector into matrix using fixed segment length and overlapping.

# Arguments

- `x::AbstractVector`
- `wlen::Int64`: window length (in samples)
- `woverlap::Int64`: overlap with the previous window (in samples)

# Returns

- `m::AbstractMatrix`
"""
function vec2mat(x::AbstractVector; wlen::Int64, woverlap::Int64)::AbstractMatrix

    @assert woverlap < length(x) "woverlap must be < $(length(x))."
    @assert wlen >= 1 "wlen must be ≥ 1."
    @assert woverlap >= 0 "woverlap must be ≥ 0."

    (wlen == 1 && woverlap == 0) && return reshape(x, length(x), :)

    seg = length(x) ÷ wlen
    m = zeros(seg, wlen)
    m[1, :] = x[1:wlen]
    for idx in 2:seg
        m[idx, :] = x[(((idx - 1) * wlen + 1) - woverlap):((((idx - 1) * wlen) - woverlap) + wlen)]
    end

    return m

end

"""
    arr2mat(x)

Reshape array into matrix.

# Arguments

- `x::AbstractArray`

# Returns

- `m::AbstractMatrix`
"""
function arr2mat(x::AbstractArray)::AbstractMatrix

    @assert size(x, 1) == 1 "First dimension of x must be 1."

    m = zeros(size(x, 3), size(x, 2))
    for idx in axes(x, 3)
        m[idx, :] = x[1, :, idx]
    end

    return m

end
