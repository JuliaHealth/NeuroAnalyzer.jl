"""
    _time_axis(n, fs; t0 = 0.0)

Exact time axis of `n` samples at `fs` Hz starting at `t0` seconds: sample `k` lies at `t0 + (k - 1) / fs`.

Each time point is a single, correctly rounded division (no accumulated step error, no rounding to a fixed number of decimals), so `_time_axis(n, fs)[s] == s2t(s, fs)` and times that are exact multiples of `1 / fs` compare equal to their decimal literals (e.g. `1.5`, `0.1`). Values are rounded only for presentation.
"""
function _time_axis(n::Integer, fs::Real; t0::Real = 0.0)::Vector{Float64}
    n >= 0 || throw(ArgumentError("n must be ≥ 0."))
    fs > 0 || throw(ArgumentError("fs must be > 0."))
    t = collect(0:(n - 1)) ./ fs
    t0 == 0 || (t .+= t0)
    return t
end

"""
    _get_t(obj::NeuroAnalyzer.NEURO)

Compute time axes for a NEURO recording.

# Returns

- `time_pts`: global time axis spanning the full recording (all epochs concatenated), starting at 0.0 s, sampled at `sr(obj)` Hz
- `epoch_time`: time axis for a single epoch; if `obj.epoch_time` is non-empty, the axis is shifted by `obj.epoch_time[1]` (i.e. the recorded epoch-start offset)

# Notes

- Both vectors contain `n_samples × n_epochs` and `n_samples` points respectively (exact, see [`_time_axis`](@ref)).
"""
function _get_t(obj::NeuroAnalyzer.NEURO)::Tuple{Vector{Float64}, Vector{Float64}}
    fs       = sr(obj)
    n_smps   = size(obj.data, 2)   # samples per epoch
    n_epochs = size(obj.data, 3)   # number of epochs

    # global timeline: n_smps * n_epochs points, step = 1/fs seconds
    time_pts = _time_axis(n_smps * n_epochs, fs)

    # per-epoch timeline, optionally shifted by the stored epoch-start offset
    t_offset   = isempty(obj.epoch_time) ? 0.0 : obj.epoch_time[1]
    epoch_time = _time_axis(n_smps, fs; t0 = t_offset)

    return time_pts, epoch_time
end

"""
    _get_t(from::Int64, to::Int64, fs::Int64)

Build a zero-based time vector for samples `from:to` at sampling rate `fs`.

The resulting vector has `to - from + 1` elements and starts at `0.0` (exact, see [`_time_axis`](@ref)).

# Arguments

- `from`: first sample index (inclusive, 0-based)
- `to`: last sample index (inclusive, 0-based)
- `fs`: sampling rate in Hz

# Returns

- `Vector{Float64}`: time points
"""
function _get_t(from::Int64, to::Int64, fs::Int64)::Vector{Float64}
    return _time_axis(to - from + 1, fs)
end

"""
    _convert_t(t1::Float64, t2::Float64)

Format a pair of timestamps as human-readable strings

Values with absolute magnitude below 1 s are expressed in milliseconds; values ≥ 1 s are expressed in seconds.  `t1` is floored and `t2` is ceiled so the displayed interval is always at least as wide as the true interval (safe for axis-label generation).

# Returns

- `(t1, ts1, t2, ts2)` where `ts1`/`ts2` are the formatted strings.
"""
function _convert_t(t1::Float64, t2::Float64)::Tuple{Float64, String, Float64, String}
    # local helper: pick unit, apply rounding direction, build string
    function _fmt(t::Float64, round_fn::Function)::String
        if abs(t) < 1.0
            return string(round_fn(t * 1000.0; digits = 4)) * " ms"
        else
            return string(round_fn(t; digits = 4)) * " s"
        end
    end

    return t1, _fmt(t1, floor), t2, _fmt(t2, ceil)
end

"""
    _s2epoch(obj::NeuroAnalyzer.NEURO, from::Int64, to::Int64)

Map a sample range `[from, to]` to the range of epochs that contain these samples.

`from` and `to` are 1-based sample numbers; epoch `k` holds samples `(k - 1) * epoch_len + 1 : k * epoch_len`, so the result runs from the epoch containing `from` to the epoch containing `to` (partial epochs at both ends are included). A `from` value of 0 is clamped to epoch 1.

# Returns
- A single `Int64` if the sample range lies entirely within one epoch.
- An `AbstractUnitRange{Int64}` otherwise.
"""
function _s2epoch(
    obj::NeuroAnalyzer.NEURO,
    from::Int64,
    to::Int64,
)::Union{Int64, AbstractUnitRange{Int64}}
    from >= 0 || throw(ArgumentError("from must be ≥ 0."))
    to >= from || throw(ArgumentError("to must be ≥ from."))
    to <= length(obj.time_pts) ||
        throw(ArgumentError("to must be ≤ $(length(obj.time_pts))."))
    el = epoch_len(obj)

    # epoch containing sample s is ceil(s / el); clamp to 1 for from = 0
    ep_first = max(1, cld(from, el))
    ep_last  = max(1, cld(to, el))

    return ep_first == ep_last ? ep_first : ep_first:ep_last
end

"""
    _epoch2s(obj::NeuroAnalyzer.NEURO, ep::Int64) -> Tuple{Int64, Int64}

Return the first and last **1-based** sample indices of epoch `ep`.

This is the inverse of `_s2epoch`.

# Returns
- `(t1, t2)` tuple, where
    - `t1 = (ep - 1) * epoch_len(obj) + 1` — first sample (1-based)
    - `t2 =  ep      * epoch_len(obj)`     — last  sample (1-based)
"""
function _epoch2s(obj::NeuroAnalyzer.NEURO, ep::Int64)::Tuple{Int64, Int64}
    el = epoch_len(obj)
    t1 = (ep - 1) * el + 1    # first sample of this epoch (1-based)
    t2 = ep * el        # last  sample of this epoch (1-based)
    return t1, t2
end
