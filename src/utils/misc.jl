export linspace
export logspace
export cmax
export cmin
export cextrema
export cums
export f_nearest
export ntapers
export trtm

"""
    linspace(start, stop, length)

Return a vector of `n` evenly spaced numbers from `start` to `stop` (inclusive).

Thin wrapper around `Base.range` that always materialises the result as a `Vector{Float64}`.

# Arguments

- `start::Real`: start of the sequence
- `stop::Real`: end of the sequence (inclusive)
- `n::Int64`: number of points; must be ≥ 2

# Returns

- `Vector{Float64}`: linearly spaced sequence of length `n`

# Throws

- `ArgumentError`: if `n < 2`

# See also

[`logspace`](@ref)
"""
function linspace(start::Real, stop::Real, n::Int64)::Vector{Float64}

    @assert n >= 2 "n must be ≥ 2."
    return collect(range(start, stop, n))

end

"""
    logspace(start, stop, n)

Return a vector of `n` logarithmically spaced numbers from `start` to `stop` (inclusive).

Requires `start > 0` and `stop > 0` (logarithmic spacing is undefined for non-positive values).

# Arguments

- `start::Real`: start of the sequence; must be > 0
- `stop::Real`: end of the sequence (inclusive); must be > 0
- `n::Int64`: number of points; must be ≥ 2

# Returns
- `Vector{Float64}`: logarithmically spaced sequence of length `n`

# Throws

- `ArgumentError`: if `n < 2`, `start ≤ 0`, or `stop ≤ 0`

# See also

[`linspace`](@ref)
"""
function logspace(start::Number, stop::Number, n::Int64)::Vector{Float64}

    @assert n >= 2 "n must be ≥ 2."
    @assert start > 0 "start must be > 0."
    @assert stop > 0 "stop must be > 0."
    return Float64.(logrange(start, stop, n))

end

"""
    cmax(x)

Return the element of a complex vector with the largest magnitude.

Selects the element that maximises `|x|²` (equivalent to maximising `|x|`).

# Arguments

- `x::Vector{ComplexF64}`: input complex vector; must not be empty

# Returns

- `cmax::ComplexF64`: element of `x` with the largest absolute value

# See also

[`cmin`](@ref), [`cextrema`](@ref)
"""
function cmax(x::Vector{<:Complex})::ComplexF64

    return argmax(abs2, x)

end

"""
    cmin(x)

Return the element of a complex vector with the smallest magnitude.

Selects the element that minimises `|x|²` (equivalent to minimising `|x|`).

# Arguments

- `x::Vector{<:Complex}`: input complex vector; must not be empty

# Returns

- `ComplexF64`: element of `x` with the smallest absolute value

# See also

[`cmax`](@ref), [`cextrema`](@ref)
"""
function cmin(x::Vector{<:Complex})::ComplexF64

    return argmin(abs2, x)

end

"""
    cextrema(x)

Return the elements of a complex vector with the largest and smallest magnitudes.

# Arguments

- `x::Vector{ComplexF64}`: input complex vector; must not be empty

# Returns

Tuple containing:
- `ComplexF64`: element with the largest absolute value (`cmax`)
- `ComplexF64`: element with the smallest absolute value (`cmin`)

# See also

[`cmax`](@ref), [`cmin`](@ref)
"""
function cextrema(x::Vector{<:Complex})::Tuple{ComplexF64, ComplexF64}

    return (cmax(x), cmin(x))

end

"""
    cums(signal)

Compute the cumulative sum of a 3-dimensional array along the sample (second) axis.

# Arguments

- `s::Array{<:Real, 3}`: input array of shape `(channels, samples, epochs)`

# Returns

- `Array{Float64, 3}`: cumulative sum array of the same shape as `signal`
"""
function cums(s::Array{<:Real, 3})::Array{Float64, 3}

    ch_n, _, ep_n = size(s)

    # pre-allocate output
    # ensure Float64 output even for integer input
    csa = similar(signal, Float64)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        csa[ch_idx, :, ep_idx] = cumsum(@view(csa[ch_idx, :, ep_idx]))
    end

    return csa

end

"""
    f_nearest(m, pos)

Find the matrix position closest to a query position tuple.

Computes the Euclidean distance from every element of `m` to `p` and returns the row and column indices of the nearest element.

# Arguments

- `m::Matrix{Tuple{Float64, Float64}}`: matrix of 2-D position tuples
- `p::Tuple{Float64, Float64}`: query position

# Returns

- `Tuple{Int64, Int64}`: `(row, column)` of the nearest position in `m`
"""
function f_nearest(
    m::Matrix{Tuple{Float64, Float64}},
    p::Tuple{Float64, Float64},
)::Tuple{Int64, Int64}

    d = zeros(size(m))

    @inbounds for idx1 in axes(m, 1), idx2 in axes(m, 2)
        d[idx1, idx2] = euclidean(m[idx1, idx2], p)
    end

    # compute findmin once and reuse both the value and the CartesianIndex
    _, ci = findmin(d)

    return (ci[1], ci[2])

end

"""
    ntapers(obj; <keyword arguments>)

Return the recommended number of Slepian tapers for multi-taper spectral analysis.

The formula is `nt = floor(df × T) - 1`, where `T = epoch_len / fs` is the epoch duration in seconds. The result is clamped to a minimum of 1.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `df::Real`: desired frequency resolution (bandwidth) in Hz; must be > 0; this is the minimum separation between frequency peaks to be resolved

# Returns

- `Int64`: recommended number of Slepian tapers (≥ 1)

# Throws
- `ArgumentError`: if `df` is outside the valid range `(0, fs/2)`
"""
function ntapers(obj::NeuroAnalyzer.NEURO; df::Real)::Int64

    # validate that df lies within (0, Nyquist)
    _bin(df, (0, sr(obj) / 2))
    n = epoch_len(obj) / sr(obj) # epoch duration in seconds
    nt = round(Int64, df * n) - 1

    # guard: the formula can yield 0 or negative for very coarse resolution
    return max(nt, 1)

end

"""
    trtm(obj; <keyword arguments>)

Return a single channel's signal in trials × time format.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `ep::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj))`: epoch numbers; default use all epochs

# Returns

- `Matrix{Float64}`: matrix of shape `(n_epochs, epoch_len)`

# Throws

- `ArgumentError`: if `ch` resolves to more than one channel, or if any epoch index in `ep` is out of range
"""
function trtm(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    ep::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj)),
)::Matrix{Float64}

    _check_epochs(obj, ep)
    ch = get_channel(obj, ch = ch)
    @assert length(ch) == 1 "ch must resolve to exactly one channel."
    ch = ch[1]

    return Matrix(obj.data[ch, :, ep]')

end
