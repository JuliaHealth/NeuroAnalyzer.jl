export m_pad0
export m_sortperm
export m_sort
export m_norm
export vec2mat
export arr2mat
export meshgrid

"""
    m_pad0(m)

Pad a matrix with zeros to make it square.

Adds zero columns if the matrix has more rows than columns, or zero rows if it has more columns than rows. Returns the matrix unchanged if it is already square.

# Arguments

- `m::AbstractMatrix`: input matrix

# Returns

- `AbstractMatrix`: square matrix of size `max(r, c) × max(r, c)` with the same element type as `m`.

# See also

[`m_pad0(::AbstractMatrix, ::Int64, ::Int64)`](@ref)
"""
function m_pad0(m::AbstractMatrix)::AbstractMatrix

    mr, mc = size(m)

    if mr > mc
        # append zero columns to the right to reach square size
        return hcat(m, zeros(eltype(m), mr, mr - mc))
    elseif mr < mc
        # append zero rows to the bottom to reach square size
        return vcat(m, zeros(eltype(m), mc - mr, mc))
    else
        return m
    end

end

"""
    m_pad0(m, r, c)

Pad a matrix with zeros to reach the target size `r × c`.

Rows are appended to the bottom and columns to the right as needed. At least one of `r`, `c` must be strictly larger than the corresponding dimension of `m`.

# Arguments

- `m::AbstractMatrix`: input matrix
- `r::Int64`: target number of rows; must be ≥ `size(m, 1)`
- `c::Int64`: target number of columns; must be ≥ `size(m, 2)`

# Returns

- `AbstractMatrix`: matrix of size `r × c` padded with zeros of the same element type as `m`

# Throws

- `ArgumentError`: if `r < size(m, 1)` or `c < size(m, 2)`

# See also

[`m_pad0(::AbstractMatrix)`](@ref)
"""
function m_pad0(m::AbstractMatrix, r::Int64, c::Int64)::AbstractMatrix

    # Tuple comparison in Julia is lexicographic, not element-wise; check each dimension independently to avoid a silent false-positive assertion.
    !(r >= size(m, 1)) && throw(ArgumentError("r ($r) must be ≥ size(m, 1) ($(size(m, 1)))."))
    !(c >= size(m, 2)) && throw(ArgumentError("c ($c) must be ≥ size(m, 2) ($(size(m, 2)))."))
    !(r > size(m, 1) || c > size(m, 2)) && throw(ArgumentError("At least one of r, c must exceed the current size."))

    mr, mc = size(m)
    # append zero rows to the bottom
    if r > mr
        m = vcat(m, zeros(eltype(m), r - mr, mc))
    end
    # append zero columns to the right (use updated row count r, not original mr)
    if c > mc
        m = hcat(m, zeros(eltype(m), r, c - mc))
    end

    return m

end


"""
    m_sortperm(m; <keyword arguments>)

Return the sorting permutation indices of a matrix column-wise or row-wise.

# Arguments

- `m::AbstractMatrix`: input matrix
- `rev::Bool`: if `true`, sort in descending order
- `dims::Int64=1`: sort along columns (`dims=1`) or rows (`dims=2`)

# Returns

- `Matrix{Int64}`: index matrix of the same size as `m`; each column (or row) contains the permutation that would sort that column (or row)

# Throws

- `ArgumentError`: if `dims ∉ {1, 2}`

# See also

[`m_sort`](@ref)
"""
function m_sortperm(m::AbstractMatrix; rev::Bool = false, dims::Int64 = 1)::AbstractMatrix

    !(dims in (1, 2)) && throw(ArgumentError("dims must be 1 or 2."))

    idx = zeros(Int64, size(m))
    if dims == 1
        # compute sort permutation for each column independently
        @inbounds for col in axes(m, 2)
            idx[:, col] = sortperm(m[:, col]; rev=rev)
        end
    else
        # compute sort permutation for each row independently
        @inbounds for row in axes(m, 1)
            idx[row, :] = sortperm(m[row, :]; rev=rev)'
        end
    end

    return idx

end

"""
    m_sort(m, m_idx; <keyword arguments>)

Sort a matrix using a pre-computed permutation index vector.

# Arguments

- `m::AbstractMatrix`: input matrix
- `m_idx::Vector{Int64}`: permutation index vector (e.g. from `sortperm`); this vector is **not** modified by the function
- `rev::Bool=false`: if `true`, reverse the permutation before applying it
- `dims::Int64=1`: apply permutation along columns (`dims=1`) or rows (`dims=2`)

# Returns

- `AbstractMatrix`: sorted matrix with the same size and element type as `m`

# Throws

- `ArgumentError`: if `dims ∉ {1, 2}`

# See also

[`m_sortperm`](@ref)
"""
function m_sort(
    m::AbstractMatrix,
    m_idx::Vector{Int64};
    rev::Bool = false,
    dims::Int64 = 1
)::AbstractMatrix

    !(dims in [1, 2]) && throw(ArgumentError("dims must be 1 or 2."))

    # copy to avoid mutating the caller's index vector
    perm = rev ? reverse(m_idx) : copy(m_idx)

    m_sorted = zeros(eltype(m), size(m))
    if dims == 1
        # apply permutation to each column
        @inbounds for col in axes(m, 2)
            m_sorted[:, col] = @views m[:, col][perm]
        end
    else
        # apply permutation to each row
        @inbounds for row in axes(m, 1)
            m_sorted[row, :] = @views m[row, :][perm]
        end
    end

    return m_sorted

end

"""
    m_norm(m)

Normalize an array by the number of columns minus one (`size(m, 2) - 1`).

# Arguments

- `m::AbstractArray`: input array; must have at least 2 columns

# Returns

- `AbstractArray`: array divided element-wise by `size(m, 2) - 1`

# Throws

- `ArgumentError`: if `size(m, 2) < 2` (would cause division by zero)
"""
function m_norm(m::AbstractArray)::AbstractArray

    !(size(m, 2) >= 2) && throw(ArgumentError("m must have at least 2 columns (size(m, 2) - 1 would be zero)."))
    return m ./ (size(m, 2) - 1)

end

"""
    vec2mat(x; <keyword arguments>)

Reshape a vector into a matrix using a sliding window with optional overlap.

The vector is divided into `⌊length(x) / wlen⌋` non-overlapping segments of length `wlen`. The `woverlap` parameter shifts the start of each subsequent window backward by that many samples.

# Arguments

- `x::AbstractVector`: input vector
- `wlen::Int64`: window length in samples; must be ≥ 1
- `woverlap::Int64`: number of overlapping samples between consecutive windows; must satisfy `0 ≤ woverlap < wlen`

# Returns

- `AbstractMatrix`: matrix of shape `(n_segments, wlen)`.

# Throws

- `ArgumentError`: if `wlen < 1`, `woverlap < 0`, or `woverlap ≥ wlen`
"""
function vec2mat(x::AbstractVector; wlen::Int64, woverlap::Int64)::AbstractMatrix

    !(wlen >= 1) && throw(ArgumentError("wlen must be ≥ 1."))
    !(woverlap >= 0) && throw(ArgumentError("woverlap must be ≥ 0."))
    !(woverlap < wlen) && throw(ArgumentError("woverlap must be < wlen ($wlen)."))  # was: < length(x)

    (wlen == 1 && woverlap == 0) && return reshape(x, length(x), :)

    seg = length(x) ÷ wlen
    m   = zeros(eltype(x), seg, wlen)
    m[1, :] = x[1:wlen]
    for idx in 2:seg
        start = (idx - 1) * wlen + 1 - woverlap
        m[idx, :] = x[start:(start + wlen - 1)]
    end

    return m

end

"""
    arr2mat(x)

Reshape a 3-D array of shape `(1, samples, epochs)` into a `(epochs, samples)` matrix.

# Arguments

- `x::AbstractArray`: 3-D input array; first dimension must equal 1

# Returns

- `AbstractMatrix`: matrix of shape `(size(x, 3), size(x, 2))`

# Throws

- `ArgumentError`: if `size(x, 1) ≠ 1`
"""
function arr2mat(x::AbstractArray)::AbstractMatrix

    !(size(x, 1) == 1) && throw(ArgumentError("First dimension of x must be 1; got $(size(x, 1))."))
    # equivalent to squeezing the singleton first dimension and transposing:
    # dropdims then permute, or simply index and collect row-by-row.
    return reshape(permutedims(x[1, :, :]), size(x, 3), size(x, 2))

end

"""
    meshgrid(x)

Build a 2-D mesh grid from two coordinate vectors.

Returns a pair `(mx, my)` where `mx[i]` repeats the full `x` vector for row `i`, and `my[i]` repeats `y[i]` across all columns. This matches the convention of MATLAB's `meshgrid` / NumPy's `meshgrid(..., indexing="xy")`.

# Arguments

- `x::Vector{Float64}`: x-coordinate vector of length `n`
- `y::Vector{Float64}`: y-coordinate vector of length `m`

# Returns

- `Tuple{Vector{Vector{Float64}}, Vector{Vector{Float64}}}`:
    - `mx`: `m` repetitions of `x` (one per row)
    - `my`: `m` vectors each filled with the corresponding `y[i]` value
"""
function meshgrid(
    x::Vector{Float64},
    y::Vector{Float64}
)::Tuple{
    Vector{Vector{Float64}},
    Vector{Vector{Float64}}
}

    xn = length(x)
    yn = length(y)

    # mx: each of the yn rows is the full x vector
    mx = [copy(x) for _ in 1:yn]
    # my: each of the yn rows is a constant vector equal to y[i]
    my = [fill(y[i], xn) for i in 1:yn]

    return (mx, my)

end
