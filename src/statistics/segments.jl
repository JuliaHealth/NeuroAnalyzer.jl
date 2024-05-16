export seg_mean
export seg_extract

"""
    seg_mean(seg)

Calculate mean of a segment (e.g. spectrogram).

# Arguments

- `seg::AbstractArray`

# Returns

- `seg_mean::Vector{Float64}`: averaged segment
"""
function seg_mean(seg::AbstractArray)

    return reshape(mean(mean(seg, dims=1), dims=2), size(seg, 3))

end

"""
    seg2_mean(seg1, seg2)

Calculate mean of two segments (e.g. spectrograms).

# Arguments

- `seg1::AbstractArray`
- `seg2::AbstractArray`

# Returns

Named tuple containing:
- `seg1::Vector{Float64}`: averaged segment 1
- `seg2::Vector{Float64}`: averaged segment 2
"""
function seg_mean(seg1::AbstractArray, seg2::AbstractArray)

    return (seg1=seg_mean(seg1), seg2=seg_mean(seg2))

end

"""
    seg_extract(m, r1, c1, r2, c2; v)

Extract segment from a matrix.

# Arguments

- `m::AbstractMatrix`
- `rc::NTuple{4, Int64}`: upper-left corner row and column, bottom-right corner row and column
- `c::Bool=false`: if true, use circular segment; for circular segment the segment is always returned as vector
- `v::Bool=false`: if true, return as vector (matrix m by rows over columns)

# Returns

- `seg::Union{AbstractMatrix, AbstractVector}`
"""
function seg_extract(m::AbstractMatrix, rc::NTuple{4, Int64}; v::Bool=false, c::Bool=false)

    r1 = rc[1]
    c1 = rc[2]
    r2 = rc[3]
    c2 = rc[4]

    @assert r1 > 0 "r1 must be > 0."
    @assert r2 > 0 "r2 must be > 0."
    @assert c1 > 0 "c1 must be > 0."
    @assert c2 > 0 "c2 must be > 0."
    @assert r1 <= size(m, 1) "r1 must be ≤ $(size(m, 1))."
    @assert c1 <= size(m, 2) "r2 must be ≤ $(size(m, 2))."
    @assert r2 <= size(m, 1) "c1 must be ≤ $(size(m, 1))."
    @assert c2 <= size(m, 2) "c2 must be ≤ $(size(m, 2))."

    if !c
        return !v ? m[r1:r2, c1:c2] : vec(m[r1:r2, c1:c2])
    else
        seg = zeros(Bool, size(m))
        seg_radius = distance((r1, c1), (r2, c2))
        for idx_r in 1:size(m, 1), idx_c in 1:size(m, 2)
            if distance((r1, c1), (idx_r, idx_c)) <= seg_radius
                seg[idx_r, idx_c] = true
            end
        end

        return m[seg .== true]

    end
end