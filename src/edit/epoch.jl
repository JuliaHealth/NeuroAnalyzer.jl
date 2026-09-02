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
    # validate
    nepochs(obj) == 1 ||
        throw(ArgumentError("epoch() must be applied to continuous object."))

    # create new dataset
    obj_tmp = deepcopy(obj)

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

        mrk_idx   = findall(obj_tmp.markers.value .== marker)
        mrk_start = obj_tmp.markers[mrk_idx, :start]
        mrk_len   = obj_tmp.markers[mrk_idx, :length]

        # remove markers that would begin before the signal start
        for idx in length(mrk_start):-1:1
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

        epochs, obj_tmp.markers = _make_epochs_bymarkers(
            obj_tmp.data;
            marker       = marker,
            markers      = deepcopy(obj_tmp.markers),
            marker_start = round.(Int64, mrk_start .* sr(obj)),
            offset       = round(Int64, offset * sr(obj)),
            ep_len       = round(Int64, ep_len * sr(obj)),
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

        # remove markers that fall outside the new epoch grid
        for marker_idx in DataFrames.nrow(obj_tmp.markers):-1:1
            round(Int64, sr(obj) * obj_tmp.markers[marker_idx, :start]) in
            0:(size(epochs, 2) * size(epochs, 3)) ||
                deleteat!(obj_tmp.markers, marker_idx)
        end
    end

    obj_tmp.data = epochs
    obj_tmp.header.recording[:epoch_id] = epoch_id
    obj_tmp.header.recording[:bad_channel] = zeros(Bool, size(obj_tmp.data, 1))
    obj_tmp.time_pts, obj_tmp.epoch_time = _get_t(obj_tmp)
    obj_tmp.epoch_time .-= offset

    push!(obj_tmp.history, "epoch(obj; marker=$marker, offset=$offset, ep_len=$ep_len)")

    return obj_tmp
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
    obj_tmp = epoch(obj; marker = marker, offset = offset, ep_len = ep_len)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj.time_pts = obj_tmp.time_pts
    obj.epoch_time = obj_tmp.epoch_time
    obj.markers = obj_tmp.markers
    obj_tmp = nothing

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

    # create new dataset
    obj_tmp = deepcopy(obj)

    obj_tmp.epoch_time .+= ts
    push!(obj_tmp.history, "epoch_ts(obj, ts=$ts)")

    return obj_tmp
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
    obj_tmp = epoch_ts(obj; ts = ts)
    obj.history = obj_tmp.history
    obj.time_pts = obj_tmp.time_pts
    obj.epoch_time = obj_tmp.epoch_time
    obj_tmp = nothing

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
    # validate
    ep_time = obj.epoch_time
    ep_start >= ep_time[1] || throw(ArgumentError("ep_start must be ≥ $(ep_time[1])."))
    ep_end <= ep_time[end] || throw(ArgumentError("ep_end must be ≤ $(ep_time[end])."))
    ep_start < ep_end || throw(ArgumentError("ep_start must be < ep_end."))

    # create new dataset
    obj_tmp = deepcopy(obj)

    ep_start_idx = vsearch(ep_start, ep_time)
    ep_end_idx   = vsearch(ep_end, ep_time)

    obj_tmp.data = obj.data[:, ep_start_idx:ep_end_idx, :]
    obj_tmp.epoch_time = ep_time[ep_start_idx:ep_end_idx]
    obj_tmp.time_pts, _ = _get_t(obj_tmp)

    # compute per-epoch time windows in absolute signal time
    ep_tps = _epochs_tps(obj)
    ep_tps[2, :] = ep_tps[1, :] .+ ep_end
    ep_tps[1, :] .+= ep_start

    mrk_start = obj.markers.start
    mrk_epoch = _markers_epochs(obj)

    # remove markers outside the retained window
    for mrk_idx in length(mrk_start):-1:1
        ep = mrk_epoch[mrk_idx]
        if mrk_start[mrk_idx] < ep_tps[1, ep] || mrk_start[mrk_idx] > ep_tps[2, ep]
            deleteat!(obj_tmp.markers, mrk_idx)
            deleteat!(mrk_epoch, mrk_idx)
        end
    end

    # shift remaining marker timestamps to align with the trimmed epochs
    mrk_start_new = deepcopy(obj_tmp.markers.start)
    for mrk_idx in eachindex(mrk_start_new)
        ep = mrk_epoch[mrk_idx]
        mrk_start_new[mrk_idx] -= (
            ep_start + (ep - 1) * (ep_start + (obj.time_pts[epoch_len(obj)] - ep_end))
        )
    end
    obj_tmp.markers.start = round.(mrk_start_new; digits = 3)

    push!(obj_tmp.history, "subepoch(obj; ep_start=$ep_start, ep_end=$ep_end)")

    return obj_tmp
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
function subepoch!(obj::NeuroAnalyzer.NEURO; ep_start::Real, ep_end::Real)::Nothing
    obj_tmp = subepoch(obj; ep_start = ep_start, ep_end = ep_end)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj.time_pts = obj_tmp.time_pts
    obj.epoch_time = obj_tmp.epoch_time
    obj.markers = obj_tmp.markers
    obj_tmp = nothing

    return nothing
end
