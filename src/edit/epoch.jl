export epoch
export epoch!
export epoch_ts
export epoch_ts!
export subepoch
export subepoch!

"""
    epoch(obj; <keyword arguments>)

Return a copy of `obj` split into epochs.

Epochs are created either by splitting at marker positions (when `marker` is specified) or by dividing the signal into fixed-length windows (`ep_len`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be continuous (1 epoch)
- `marker::String=""`: marker value to split at; if empty, split by `ep_len`
- `offset::Real=0`: time offset in seconds for marker-based epoching (each epoch begins at `marker_time - offset`)
- `ep_len::Union{Real, Nothing}=nothing`: epoch length in seconds

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function epoch(
    obj::NeuroAnalyzer.NEURO;
    marker::String = "",
    offset::Real = 0,
    ep_len::Union{Real, Nothing} = nothing,
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    epoch!(obj_new; marker = marker, offset = offset, ep_len = ep_len)

    return obj_new
end

"""
    epoch!(obj; <keyword arguments>)

Split `obj` into epochs in-place.

Epochs are created either by splitting at marker positions (when `marker` is specified) or by dividing the signal into fixed-length windows (`ep_len`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be continuous (1 epoch)
- `marker::String=""`: marker value to split at; if empty, split by `ep_len`
- `offset::Real=0`: time offset in seconds for marker-based epoching (each epoch begins at `marker_time - offset`)
- `ep_len::Union{Real, Nothing}=nothing`: epoch length in seconds

# Returns

- `Nothing`
"""
function epoch!(
    obj::NeuroAnalyzer.NEURO;
    marker::String = "",
    offset::Real = 0,
    ep_len::Union{Real, Nothing} = nothing,
)::Nothing
    # validate
    nepochs(obj) == 1 ||
        throw(ArgumentError("epoch() must be applied to continuous object."))

    # create ID for epochs
    epoch_id = if marker != ""
        marker
    elseif !isnothing(ep_len)
        "length_$(ep_len)s"
    else
        "full"
    end

    if marker != ""
        # marker-based epoching
        _has_markers(obj) || throw(ArgumentError("OBJ does not contain markers."))
        _check_markers(obj, marker)
        isnothing(ep_len) &&
            throw(ArgumentError("ep_len must be specified for marker-based epoching."))

        mrk_idx   = findall(obj.markers.value .== marker)
        mrk_start = obj.markers[mrk_idx, :start]
        mrk_len   = obj.markers[mrk_idx, :length]

        # remove markers that would begin before the signal start
        for idx = length(mrk_start):-1:1
            if mrk_start[idx] - offset < obj.time_pts[1]
                deleteat!(mrk_start, idx)
                deleteat!(mrk_len, idx)
            end
        end

        isempty(mrk_start) && throw(
            ArgumentError(
                "No markers remain after applying offset; all markers fall before signal start.",
            ),
        )

        offset + ep_len >= maximum(mrk_len) || throw(
            ArgumentError(
                "offset + ep_len must be ≥ $(maximum(mrk_len)) (maximum marker length).",
            ),
        )

        epochs, obj.markers = _make_epochs_bymarkers(
            obj.data;
            marker  = marker,
            markers = deepcopy(obj.markers),
            # 1-based sample number of each marker (a marker at t = 0 s is sample 1)
            marker_start = t2s.(mrk_start, sr(obj)),
            offset       = t2n(offset, sr(obj)),
            ep_len       = t2n(ep_len, sr(obj)),
            fs           = sr(obj),
        )

    else
        # fixed-length epoching
        if !isnothing(ep_len)
            ep_len <= signal_len(obj) / sr(obj) || throw(
                ArgumentError(
                    "ep_len must be ≤ signal length ($(signal_len(obj) / sr(obj)) s).",
                ),
            )
            ep_len = round(Int64, ep_len * sr(obj))
        end

        epochs = _make_epochs(obj.data; ep_len = ep_len)

        # remove markers that fall outside the new epoch grid (trailing samples are dropped)
        n_total = size(epochs, 2) * size(epochs, 3)
        obj.markers = _remap_markers(
            obj.markers, sr(obj), s -> 1 <= s <= n_total ? s : nothing,
        )
    end

    obj.data = epochs
    obj.header.recording[:epoch_id] = epoch_id
    obj.header.recording[:bad_channel] = obj.header.recording[:bad_channel]
    obj.time_pts, obj.epoch_time = _get_t(obj)
    obj.epoch_time .-= offset

    push!(obj.history, "epoch(obj; marker=$marker, offset=$offset, ep_len=$ep_len)")

    return nothing
end

"""
    epoch_ts(obj; <keyword arguments>)

Return a copy of `obj` with the epoch time axis shifted by `ts` seconds.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ts::Real`: time shift in seconds (positive → later start, negative → earlier start)

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function epoch_ts(obj::NeuroAnalyzer.NEURO; ts::Real)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    epoch_ts!(obj_new; ts = ts)

    return obj_new
end

"""
    epoch_ts!(obj; <keyword arguments>)

Shift the epoch time axis by `ts` seconds in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ts::Real`: time shift in seconds (positive → later start, negative → earlier start)

# Returns

- `Nothing`
"""
function epoch_ts!(obj::NeuroAnalyzer.NEURO; ts::Real)::Nothing

    obj.epoch_time .+= ts
    push!(obj.history, "epoch_ts(obj, ts=$ts)")

    return nothing
end

"""
    subepoch(obj; <keyword arguments>)

Return a copy of `obj` trimmed to a sub-range within each epoch.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ep_start::Real`: sub-epoch start in seconds (relative to epoch start)
- `ep_end::Real`: sub-epoch end in seconds (relative to epoch start)

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function subepoch(
    obj::NeuroAnalyzer.NEURO;
    ep_start::Real,
    ep_end::Real,
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    subepoch!(obj_new; ep_start = ep_start, ep_end = ep_end)

    return obj_new
end

"""
    subepoch!(obj; <keyword arguments>)

Trim each epoch to a sub-range in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ep_start::Real`: sub-epoch start in seconds (relative to epoch start)
- `ep_end::Real`: sub-epoch end in seconds (relative to epoch start)

# Returns

- `Nothing`
"""
function subepoch!(
    obj::NeuroAnalyzer.NEURO;
    ep_start::Real,
    ep_end::Real,
)::Nothing
    # validate
    ep_time = obj.epoch_time
    ep_start >= ep_time[1] || throw(ArgumentError("ep_start must be ≥ $(ep_time[1])."))
    ep_end <= ep_time[end] || throw(ArgumentError("ep_end must be ≤ $(ep_time[end])."))
    ep_start < ep_end || throw(ArgumentError("ep_start must be < ep_end."))

    ep_start_idx = vsearch(ep_start, ep_time)
    ep_end_idx   = vsearch(ep_end, ep_time)

    obj.data = obj.data[:, ep_start_idx:ep_end_idx, :]
    obj.epoch_time = ep_time[ep_start_idx:ep_end_idx]
    obj.time_pts, _ = _get_t(obj)

    # markers: keep those within the retained window of their epoch and move
    # them to the corresponding sample of the trimmed epoch
    l_old = epoch_len(obj)
    l_new = ep_end_idx - ep_start_idx + 1
    obj.markers = _remap_markers(
        obj.markers,
        sr(obj),
        s -> begin
            ep = (s - 1) ÷ l_old
            w = (s - 1) % l_old + 1
            ep_start_idx <= w <= ep_end_idx ? ep * l_new + w - ep_start_idx + 1 : nothing
        end,
    )

    push!(obj.history, "subepoch(obj; ep_start=$ep_start, ep_end=$ep_end)")

    return nothing
end
