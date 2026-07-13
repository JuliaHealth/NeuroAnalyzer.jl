"""
    _get_t(obj::NeuroAnalyzer.NEURO)

Compute time axes for a NEURO recording.

# Returns

- `time_pts`: global time axis spanning the full recording (all epochs concatenated), starting at 0.0 s, sampled at `sr(obj)` Hz
- `epoch_time`: time axis for a single epoch; if `obj.epoch_time` is non-empty, the axis is shifted by `obj.epoch_time[1]` (i.e. the recorded epoch-start offset)

# Notes

- Both vectors contain `n_samples × n_epochs` and `n_samples` points respectively, rounded to 4 decimal places.
- Using `range(...; length=n)` instead of `collect(a:step:b)[1:end-1]` avoids the off-by-one issue and is allocation-friendlier.
"""
function _get_t(obj::NeuroAnalyzer.NEURO)::Tuple{Vector{Float64}, Vector{Float64}}
    fs       = sr(obj)
    n_smps   = size(obj.data, 2)   # samples per epoch
    n_epochs = size(obj.data, 3)   # number of epochs

    # global timeline: n_smps * n_epochs points, step = 1/fs seconds
    time_pts = round.(
        collect(range(0.0; step = 1.0 / fs, length = n_smps * n_epochs));
        digits = 4,
    )

    # per-epoch timeline, optionally shifted by the stored epoch-start offset
    t_offset   = isempty(obj.epoch_time) ? 0.0 : obj.epoch_time[1]
    epoch_time = round.(
    collect(range(t_offset; step = 1.0 / fs, length = n_smps));
    digits = 4
)

    return time_pts, epoch_time
end

"""
    _get_t(from::Int64, to::Int64, fs::Int64)

Build a zero-based time vector for samples `from:to` at sampling rate `fs`.

The resulting vector has `to - from + 1` elements, starts at `0.0`, and is rounded to 4 decimal places.

# Arguments

- `from`: first sample index (inclusive, 0-based)
- `to`: last sample index (inclusive, 0-based)
- `fs`: sampling rate in Hz

# Returns

- `Vector{Float64}`: time points
"""
function _get_t(from::Int64, to::Int64, fs::Int64)::Vector{Float64}
    # `length = to - from + 1` is exact; avoids floating-point creep in step-based ranges
    return round.(
        collect(range(0.0; step = 1.0 / fs, length = to - from + 1));
        digits = 4,
    )
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

Map a sample range `[from, to]` to the range of epoch indices that it spans.

Epoch indices are 1-based. A `from` value of 0 is clamped to epoch 1.

If `from` does not fall on an epoch boundary (i.e. it is in the interior of an epoch), that leading partial epoch is excluded and the range starts with the next full epoch.

# Returns
- A single `Int64` if the sample range lies entirely within one epoch.
- An `AbstractUnitRange{Int64}` otherwise.

# Notes
- `from` and `to` are 0-based sample indices.
- The exclusion of a leading partial epoch may be intentional (e.g. when the caller wants only epochs that start at or after `from`). Verify this matches your use-case before relying on it.
"""
function _s2epoch(
    obj::NeuroAnalyzer.NEURO,
    from::Int64,
    to::Int64,
)::Union{Int64, AbstractUnitRange{Int64}}
    from >= 0 || throw(ArgumentError("from must be ≥ 0."))
    to <= length(obj.time_pts) || throw(ArgumentError("to must be ≤ $(length(obj.time_pts))."))
    el = epoch_len(obj)

    ep_first = floor(Int64, from / el)
    ep_last  = ceil(Int64, to / el)

    # if `from` is strictly inside an epoch (not on a boundary), the bounding
    # epoch starts one step later — the partial epoch is excluded
    if from > 0 && !iszero(mod(from, el))
        ep_first += 1
    end

    # clamp to 1-based indexing (handles `from = 0`)
    ep_first = max(1, ep_first)

    return ep_first == ep_last ? ep_first : ep_first:ep_last
end

"""
    _epoch2s(obj::NeuroAnalyzer.NEURO, ep::Int64) -> Tuple{Float64, Float64}

Return the first and last **1-based** sample indices of epoch `ep`.

This is the inverse of `_s2epoch`.

# Returns
- `(t1, t2)` Tupe, where
    - `t1 = (ep - 1) * epoch_len(obj) + 1` — first sample (1-based)
    - `t2 =  ep      * epoch_len(obj)`     — last  sample (1-based)
"""
function _epoch2s(obj::NeuroAnalyzer.NEURO, ep::Int64)::Tuple{Int64, Int64}
    el = epoch_len(obj)
    t1 = (ep - 1) * el + 1    # first sample of this epoch (1-based)
    t2 = ep * el        # last  sample of this epoch (1-based)
    return t1, t2
end
