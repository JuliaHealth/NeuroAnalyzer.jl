export t2s
export s2t
export t2n
export n2t
export markers_s2t
export markers_s2t!
export e2t

# Two kinds of conversion are provided:
# - positions: `t2s` / `s2t` map a time point to a sample number and back (sample 1 is at t = 0)
# - durations: `t2n` / `n2t` map a time span to a number of samples and back
# `t2s` and `s2t` are exact inverses on the sample grid: `t2s(s2t(s, fs), fs) == s`.

"""
    t2s(t, fs)

Convert a time point in seconds to the nearest sample number.

Sample numbering starts at 1, so `t = 0` maps to sample 1 and sample `s` lies at `(s - 1) / fs` seconds. This is the inverse of [`s2t`](@ref): `t2s(s2t(s, fs), fs) == s`.

Note that `t2s(t1, fs):t2s(t2, fs)` includes both end points, i.e. it contains `t2n(t2 - t1, fs) + 1` samples. To convert a duration to a number of samples use [`t2n`](@ref).

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

    # round (not ceil) so that floating-point error, e.g. 0.3 × 10 = 3.0000000000000004, cannot shift the result by one sample
    return round(Int64, t * fs) + 1
end

"""
    s2t(s, fs)

Convert a sample number to the time point in seconds at which it lies.

Sample numbering starts at 1: sample 1 maps to `t = 0.0` and sample `s` to `(s - 1) / fs`. This is the inverse of [`t2s`](@ref). Passing `s = 0` is invalid; a warning is emitted and `s` is clamped to 1.

To convert a number of samples to a duration use [`n2t`](@ref).

# Arguments

- `s::Real`: sample number; must be ≥ 1 (or 0, which is clamped with a warning)
- `fs::Int64`: sampling rate in Hz; must be ≥ 1

# Returns

- `Float64`: time in seconds
"""
function s2t(s::Real, fs::Int64)::Float64
    # validate
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))
    s >= 0 || throw(ArgumentError("s must be ≥ 0."))
    if s == 0
        _warn("Sample number 0 is invalid; clamped to 1.")
        s = 1
    end

    # no rounding: at high sampling rates neighbouring samples are < 0.1 ms apart
    return (s - 1) / fs
end

"""
    t2s(obj; <keyword arguments>)

Convert a time point in seconds to a sample number using the object's sampling rate.

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

Convert a sample number to the time point in seconds using the object's sampling rate.

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
    t2n(t, fs)

Convert a duration in seconds to the nearest number of samples.

E.g. a 10 s window at 256 Hz contains `t2n(10, 256) == 2560` samples. To convert a time point to a sample number use [`t2s`](@ref).

# Arguments

- `t::Real`: duration in seconds; must be ≥ 0
- `fs::Int64`: sampling rate in Hz; must be ≥ 1

# Returns

- `Int64`: number of samples (≥ 0)
"""
function t2n(t::Real, fs::Int64)::Int64
    # validate
    t >= 0 || throw(ArgumentError("t must be ≥ 0."))
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))

    return round(Int64, t * fs)
end

"""
    n2t(n, fs)

Convert a number of samples to a duration in seconds.

Computes `n / fs`. This is the inverse of [`t2n`](@ref). To convert a sample number to a time point use [`s2t`](@ref).

# Arguments

- `n::Real`: number of samples; must be ≥ 0
- `fs::Int64`: sampling rate in Hz; must be ≥ 1

# Returns

- `Float64`: duration in seconds
"""
function n2t(n::Real, fs::Int64)::Float64
    # validate
    n >= 0 || throw(ArgumentError("n must be ≥ 0."))
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))

    return n / fs
end

"""
    t2n(obj; <keyword arguments>)

Convert a duration in seconds to a number of samples using the object's sampling rate.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `t::Real`: duration in seconds; must be ≥ 0

# Returns

- `Int64`: number of samples (≥ 0)
"""
function t2n(obj::NeuroAnalyzer.NEURO; t::Real)::Int64
    return t2n(t, sr(obj))
end

"""
    n2t(obj; <keyword arguments>)

Convert a number of samples to a duration in seconds using the object's sampling rate.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `n::Real`: number of samples; must be ≥ 0

# Returns

- `Float64`: duration in seconds
"""
function n2t(obj::NeuroAnalyzer.NEURO; n::Real)::Float64
    return n2t(n, sr(obj))
end

"""
    markers_s2t(m; <keyword arguments>)

Return a copy of a markers DataFrame with `:start` and `:length` columns converted from sample numbers to seconds.

`:start` (a sample number) is converted with [`s2t`](@ref), so sample 1 maps to 0 s; `:length` (a number of samples) is converted to a duration with [`n2t`](@ref).

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
    markers_s2t!(m_new; fs = fs)

    return m_new
end

"""
    markers_s2t!(m; <keyword arguments>)

Convert `:start` and `:length` columns of a markers DataFrame from sample numbers to seconds, in-place.

`:start` (a sample number) is converted with [`s2t`](@ref), so sample 1 maps to 0 s; `:length` (a number of samples) is converted to a duration with [`n2t`](@ref).

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
    # length is a duration, not a sample position: L samples last L / fs seconds
    m[!, :length] = n2t.(m[!, :length], fs)

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
    return markers_s2t(obj.markers; fs = sr(obj))
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
- `ep::Union{Int64, AbstractUnitRange{Int64}, Vector{Int64}}`: epoch index or vector of epoch indices

# Returns

- `Tuple{Float64, Float64}`: `(start_time, end_time)` in seconds
"""
function e2t(
    obj::NeuroAnalyzer.NEURO;
    ep::Union{Int64, AbstractUnitRange{Int64}, Vector{Int64}},
)::Tuple{Float64, Float64}
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
