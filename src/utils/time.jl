export t2s
export s2t
export markers_s2t
export markers_s2t!
export e2t

"""
    t2s(t, fs)

Convert a time in seconds to the nearest sample number.

Sample numbering starts at 1: `t = 0` maps to sample 1, and any positive time is rounded up via `ceil`.

# Arguments

- `t::Real`: time in seconds; must be ≥ 0
- `fs::Int64`: sampling rate in Hz; must be ≥ 1

# Returns

- `Int64`: sample number (≥ 1)
"""
function t2s(t::Real, fs::Int64)::Int64
    # validate
    t >= 0 || throw(ArgumentError("t must be ≥ 0."))
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))

    return t == 0 ? 1 : ceil(Int64, t * fs)
end

"""
    s2t(s, fs)

Convert a sample number to time in seconds.

Sample numbering starts at 1: sample 1 maps to `t = 0.0`. Passing `s = 0` is invalid; a warning is emitted and `s` is silently clamped to 1.

# Arguments

- `s::Real`: sample number; must be ≥ 1 (or 0, which is clamped with a warning)
- `fs::Int64`: sampling rate in Hz; must be ≥ 1

# Returns

- `Float64`: time in seconds (rounded to 4 decimal places)
"""
function s2t(s::Real, fs::Int64)::Float64
    # validate
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))
    s >= 0 || throw(ArgumentError("s must be ≥ 0."))
    if s == 0
        _warn("Sample number 0 is invalid; clamped to 1.")
        s = 1
    end

    return round(s / fs - 1 / fs; digits = 4)
end

"""
    t2s(obj; <keyword arguments>)

Convert a time in seconds to a sample number using the object's sampling rate.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `t::Real`: time in seconds; must be ≥ 0

# Returns

- `Int64`: sample number (≥ 1)
"""
function t2s(obj::NeuroAnalyzer.NEURO; t::Real)::Int64
    return t2s(t, sr(obj))
end

"""
    s2t(obj; <keyword arguments>)

Convert a sample number to time in seconds using the object's sampling rate.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `s::Int64`: sample number; must be ≥ 1

# Returns

- `Float64`: time in seconds
"""
function s2t(obj::NeuroAnalyzer.NEURO; s::Int64)::Float64
    return s2t(s, sr(obj))
end

"""
    markers_s2t(m; <keyword arguments>)

Return a copy of a markers DataFrame with `:start` and `:length` columns converted from sample numbers to seconds.

# Arguments

- `m::DataFrame`: markers DataFrame containing integer `:start` and `:length` columns
- `fs::Int64`: sampling rate in Hz; must be ≥ 1

# Returns

- `DataFrame`: new DataFrame with `:start` and `:length` expressed in seconds
"""
function markers_s2t(m::DataFrame; fs::Int64)::DataFrame
    # validate
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))
    m_new = deepcopy(m)
    m_new[!, :start] = s2t.(m[!, :start], fs)
    m_new[!, :length] = s2t.(m[!, :length], fs)

    return m_new
end

"""
    markers_s2t!(m; <keyword arguments>)

Convert `:start` and `:length` columns of a markers DataFrame from sample numbers to seconds, in-place.

# Arguments

- `m::DataFrame`: markers DataFrame; modified in-place
- `fs::Int64`: sampling rate in Hz; must be ≥ 1

# Returns

- `Nothing`
"""
function markers_s2t!(m::DataFrame; fs::Int64)::Nothing
    # validate
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))

    m[!, :start] = s2t.(m[!, :start], fs)
    m[!, :length] = s2t.(m[!, :length], fs)

    return nothing
end

"""
    markers_s2t(obj; <keyword arguments>)

Return a copy of the object's markers DataFrame with `:start` and `:length` converted from sample numbers to seconds.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `DataFrame`: new DataFrame with `:start` and `:length` expressed in seconds
"""
function markers_s2t(obj::NeuroAnalyzer.NEURO)::DataFrame
    return markers_s2t(deepcopy(obj.markers); fs = sr(obj))
end

"""
    markers_s2t!(obj; <keyword arguments>)

Convert `:start` and `:length` columns of the object's markers DataFrame from sample numbers to seconds, in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Nothing`
"""
function markers_s2t!(obj::NeuroAnalyzer.NEURO)::Nothing
    markers_s2t!(obj.markers; fs = sr(obj))

    return nothing
end

"""
    e2t(obj; <keyword arguments>)

Return the time segment in seconds spanning a contiguous range of epoch indices. The returned segment runs from the start of `ep[1]` to the end of `ep[end]`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ep::Union{Int64, Vector{Int64}`: epoch index or vector of epoch indices

# Returns

- `Tuple{Float64, Float64}`: `(start_time, end_time)` in seconds
"""
function e2t(
    obj::NeuroAnalyzer.NEURO;
    ep::Union{Int64, AbstractUnitRange{Int64}, Vector{Int64}},
)::Tuple{Real, Real}
    ep = _n2v(ep)

    # validate
    _check_epochs(obj, ep)

    # epoch length
    el = epoch_len(obj)
    # first sample of the first epoch
    es = (ep[1] - 1) * el + 1
    # last  sample of the last  epoch
    ee = ep[end] * el

    return (obj.time_pts[es], obj.time_pts[ee])
end
