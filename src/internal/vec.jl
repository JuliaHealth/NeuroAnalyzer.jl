"""
    _flipx(s)

Flip a signal along the amplitude axis (i.e. negate relative to its mean).

Equivalent to reflecting around the horizontal center line: subtract the mean, negate, then add the mean back.

# Returns

- `Vector{Float64}`: flipped signal of the same length as `s`
"""
function _flipx(s::AbstractVector)::Vector{Float64}
    m = mean(s)
    return m .- s
end

"""
    _zeros(s)

Count the number of zero crossings in signal `s`.

A zero crossing is detected wherever the sign of consecutive samples differs.
"""
_zeros(s::AbstractVector)::Int64 = count(abs.(diff(sign.(s))) .!= 0)

"""
    _window_indices(n, wlen, wstep, full)

Return start/end index pairs for a sliding window over a signal of length `n`.

# Arguments
 
- `n::Int64`: total signal length
- `wlen::Int64`: window length in samples
- `wstep::Int64`: step between window starts (the overlap between consecutive windows is `wlen - wstep`)
- `full::Bool`: if `true`, omit any trailing partial window; if `false`, include it

# Returns

- `Matrix{Int64}`: Nx2 matrix where each row is `[start, stop]` (1-based, inclusive)

# Notes

`vcat` inside a loop is O(n²) — this implementation collects start indices first then builds the matrix in one pass to avoid repeated allocation.
"""
function _window_indices(
    n::Int64,
    wlen::Int64,
    wstep::Int64,
    full::Bool,
)::Matrix{Int64}
    starts = Int64[]
    idx    = 1
    while idx - 1 + wlen <= n
        push!(starts, idx)
        idx += wstep
    end
    # optionally include a trailing partial window
    !full && idx <= n && push!(starts, idx)
    isempty(starts) && return Matrix{Int64}(undef, 0, 2)
    stops = min.(starts .+ wlen .- 1, n)
    return hcat(starts, stops)
end

"""
    _split(s; <keyword arguments>)

Split signal `s` into overlapping windows of length `wlen`.

The last window may be shorter than `wlen` if the signal length is not an exact multiple of the step.

# Arguments

- `s::AbstractVector`: input signal
- `wlen::Int64`: window length in samples
- `wstep::Int64=round(Int64, wlen * 0.9)`: step between window starts
 
# Returns
 
- `Vector{Vector{Float64}}`: vector of windows, each of length ≤ `wlen`
"""
function _split(
    s::AbstractVector;
    wlen::Int64,
    wstep::Int64 = round(Int64, wlen * 0.9),
)::Vector{Vector{Float64}}
    idx_mat = _window_indices(length(s), wlen, wstep, false)
    return [Vector{Float64}(s[idx_mat[i, 1]:idx_mat[i, 2]]) for i in axes(idx_mat, 1)]
end

"""
    _fsplit(s; <keyword arguments>)

Split signal `s` into overlapping windows of exactly `wlen` samples.

Unlike `_split`, trailing partial windows are discarded.

# Arguments

- `s::AbstractVector`: input signal
- `wlen::Int64`: window length in samples
- `wstep::Int64=round(Int64, wlen * 0.9)`: step between window starts

# Returns

- `Vector{Vector{Float64}}`: vector of complete windows, each of length `wlen`
"""
function _fsplit(
    s::AbstractVector;
    wlen::Int64,
    wstep::Int64 = round(Int64, wlen * 0.9),
)::Vector{Vector{Float64}}
    idx_mat = _window_indices(length(s), wlen, wstep, true)
    return [Vector{Float64}(s[idx_mat[i, 1]:idx_mat[i, 2]]) for i in axes(idx_mat, 1)]
end

"""
    _chunks(n; <keyword arguments>)

Return start/stop index pairs for all windows over a signal of length `n` (or `length(s)`). The last window may be shorter than `wlen`.

# Arguments

- `n::Int64`: signal length
- `wlen::Int64`: window length in samples
- `wstep::Int64=round(Int64, wlen * 0.9)`: step between window starts

# Returns

- `Matrix{Int64}`: Nx2 matrix of `[start stop]` index pairs (1-based, inclusive)
"""
function _chunks(
    n::Int64;
    wlen::Int64,
    wstep::Int64 = round(Int64, wlen * 0.9),
)::Matrix{Int64}
    return _window_indices(n, wlen, wstep, false)
end

"""
    _chunks(s; <keyword arguments>)

Return start/stop index pairs for all windows over a signal of length `n` (or `length(s)`). The last window may be shorter than `wlen`.

# Arguments

- `s::AbstractVector`: inpute signal
- `wlen::Int64`: window length in samples
- `wstep::Int64=round(Int64, wlen * 0.9)`: step between window starts

# Returns

- `Matrix{Int64}`: Nx2 matrix of `[start stop]` index pairs (1-based, inclusive)
"""
_chunks(
    s::AbstractVector;
    wlen::Int64,
    wstep::Int64 = round(Int64, wlen * 0.9),
)::Matrix{Int64} = _chunks(length(s); wlen = wlen, wstep = wstep)

"""
    _fchunks(n; <keyword arguments>)

Return start/stop index pairs for complete windows only. Trailing partial windows
are omitted.

# Arguments

- `n::Int64`: signal length
- `wlen::Int64`: window length in samples
- `wstep::Int64=round(Int64, wlen * 0.9)`: step between window starts

# Returns

- `Matrix{Int64}`: Nx2 matrix of `[start stop]` index pairs (1-based, inclusive)
"""
function _fchunks(
    n::Int64;
    wlen::Int64,
    wstep::Int64 = round(Int64, wlen * 0.9),
)::Matrix{Int64}
    return _window_indices(n, wlen, wstep, true)
end

"""
    _fchunks(s; <keyword arguments>)

Return start/stop index pairs for complete windows only. Trailing partial windows are omitted.

# Arguments

- `s::AbstractVector`: input signal
- `wlen::Int64`: window length in samples
- `wstep::Int64=round(Int64, wlen * 0.9)`: step between window starts
 
# Returns
 
- `Matrix{Int64}`: Nx2 matrix of `[start stop]` index pairs (1-based, inclusive)
"""
_fchunks(
    s::AbstractVector;
    wlen::Int64,
    wstep::Int64 = round(Int64, wlen * 0.9),
)::Matrix{Int64} = _fchunks(length(s); wlen = wlen, wstep = wstep)
