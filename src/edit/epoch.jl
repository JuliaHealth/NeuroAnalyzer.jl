export epoch
export epoch!
export epoch_ts
export epoch_ts!
export subepoch
export subepoch!

"""
    epoch(obj; <keyword arguments>)

Return a copy of `obj` split into epochs.

Epochs are created either by splitting at marker positions (when `marker_value` is specified) or by dividing the signal into fixed-length windows (`ep_len`). Marker-based epochs get a condition label (epoch metadata, `obj.epochs.condition`; see `get_epochs()`, `epoch_conditions()`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be continuous (1 epoch)
- `marker_value::Union{String, Vector{String}, Regex}=""`: value(s) of the markers to split at: one value, a list of values or a regular expression (e.g. `r"^s7"`); if empty, split by `ep_len`
- `marker_type::String=""`: if not empty, only markers of this type are used (e.g. `"stim"`, so that an annotation with the same value is not used)
- `offset::Real=0`: time offset in seconds for marker-based epoching (each epoch begins at `marker_time - offset`)
- `ep_len::Union{Real, Nothing}=nothing`: epoch length in seconds
- `exclude_bad_segments::Bool=NeuroAnalyzer.exclude_bad_segments`: if `true`, epochs that overlap a bad segment (marker of type `"BAD"`, also one that begins before the epoch) of all channels are skipped, epochs that overlap a bad segment of some channels are marked as rejected (`obj.epochs.reject`) (default: the `exclude_bad_segments` preference, see `set_exclude_bad_segments()`)
- `condition::Union{String, AbstractDict}="value"`: condition label of each marker-based epoch (`obj.epochs.condition`): the name of a marker column (`"value"`, `"type"` or an extra column, e.g. `"trial_type"` of BIDS events) or a dictionary of marker values and condition labels (e.g. `Dict("s70" => "face", "s71" => "face", "s130" => "house")`; values not in the dictionary are their own labels)
- `overlap::Real=0`: overlap of consecutive fixed-length epochs in seconds (0 ≤ `overlap` < `ep_len`); markers inside the overlapping parts appear in each epoch; not used for marker-based epoching

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function epoch(
    obj::NeuroAnalyzer.NEURO;
    marker_value::Union{String, Vector{String}, Regex} = "",
    marker_type::String = "",
    offset::Real = 0,
    ep_len::Union{Real, Nothing} = nothing,
    overlap::Real = 0,
    condition::Union{String, AbstractDict} = "value",
    exclude_bad_segments::Bool = NeuroAnalyzer.exclude_bad_segments,
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    epoch!(
        obj_new;
        marker_value = marker_value,
        marker_type = marker_type,
        offset = offset,
        ep_len = ep_len,
        overlap = overlap,
        condition = condition,
        exclude_bad_segments = exclude_bad_segments,
    )

    return obj_new
end

"""
    epoch!(obj; <keyword arguments>)

Split `obj` into epochs in-place.

Epochs are created either by splitting at marker positions (when `marker_value` is specified) or by dividing the signal into fixed-length windows (`ep_len`). Marker-based epochs get a condition label (epoch metadata, `obj.epochs.condition`; see `get_epochs()`, `epoch_conditions()`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be continuous (1 epoch)
- `marker_value::Union{String, Vector{String}, Regex}=""`: value(s) of the markers to split at: one value, a list of values or a regular expression (e.g. `r"^s7"`); if empty, split by `ep_len`
- `marker_type::String=""`: if not empty, only markers of this type are used (e.g. `"stim"`, so that an annotation with the same value is not used)
- `offset::Real=0`: time offset in seconds for marker-based epoching (each epoch begins at `marker_time - offset`)
- `ep_len::Union{Real, Nothing}=nothing`: epoch length in seconds
- `exclude_bad_segments::Bool=NeuroAnalyzer.exclude_bad_segments`: if `true`, epochs that overlap a bad segment (marker of type `"BAD"`, also one that begins before the epoch) of all channels are skipped, epochs that overlap a bad segment of some channels are marked as rejected (`obj.epochs.reject`) (default: the `exclude_bad_segments` preference, see `set_exclude_bad_segments()`)
- `condition::Union{String, AbstractDict}="value"`: condition label of each marker-based epoch (`obj.epochs.condition`): the name of a marker column (`"value"`, `"type"` or an extra column, e.g. `"trial_type"` of BIDS events) or a dictionary of marker values and condition labels (e.g. `Dict("s70" => "face", "s71" => "face", "s130" => "house")`; values not in the dictionary are their own labels)
- `overlap::Real=0`: overlap of consecutive fixed-length epochs in seconds (0 ≤ `overlap` < `ep_len`); markers inside the overlapping parts appear in each epoch; not used for marker-based epoching

# Returns

- `Nothing`
"""
function epoch!(
    obj::NeuroAnalyzer.NEURO;
    marker_value::Union{String, Vector{String}, Regex} = "",
    marker_type::String = "",
    offset::Real = 0,
    ep_len::Union{Real, Nothing} = nothing,
    overlap::Real = 0,
    condition::Union{String, AbstractDict} = "value",
    exclude_bad_segments::Bool = NeuroAnalyzer.exclude_bad_segments,
)::Nothing
    # validate
    nepochs(obj) == 1 ||
        throw(ArgumentError("epoch() must be applied to continuous object."))
    overlap >= 0 || throw(ArgumentError("overlap must be ≥ 0."))
    # marker-based epoching
    by_marker = !(marker_value isa String && marker_value == "")
    marker_value isa Vector && isempty(marker_value) &&
        throw(ArgumentError("marker_value must not be an empty list."))
    overlap == 0 || (!by_marker && !isnothing(ep_len)) || throw(
        ArgumentError("overlap can only be used for fixed-length epoching (ep_len)."),
    )
    overlap == 0 || overlap < ep_len ||
        throw(ArgumentError("overlap must be < ep_len."))

    # epoch metadata of the continuous signal: onset of its first sample in the original
    # recording and its source
    ep_tbl   = _epochs_table(obj)
    # markers of the continuous signal (bad segments are checked against them)
    mrk0     = copy(obj.markers)
    onset0   = ep_tbl[1, :onset]
    source   = ep_tbl[1, :source]
    fs       = sr(obj)

    if by_marker
        # marker-based epoching
        _has_markers(obj) || throw(ArgumentError("OBJ does not contain markers."))
        isnothing(ep_len) &&
            throw(ArgumentError("ep_len must be specified for marker-based epoching."))

        mrk_sel = _marker_select(obj.markers, marker_value, marker_type)
        any(mrk_sel) || throw(
            ArgumentError(
                "Marker $(repr(marker_value))" *
                (marker_type == "" ? "" : " of type '$marker_type'") * " not found.",
            ),
        )
        # condition label of each selected marker
        mrk_cond = _marker_conditions(obj.markers, condition)
        mrk_idx  = findall(mrk_sel)
        # remove markers that would begin before the signal start
        mrk_idx   = [idx for idx in mrk_idx if obj.markers[idx, :start] - offset >= obj.time_pts[1]]
        mrk_start = obj.markers[mrk_idx, :start]
        mrk_len   = obj.markers[mrk_idx, :length]

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

        # 1-based sample number of each marker (a marker at t = 0 s is sample 1)
        mrk_smp = t2s.(mrk_start, fs)
        epochs, obj.markers = _make_epochs_bymarkers(
            obj.data;
            marker  = "",
            markers = deepcopy(obj.markers),
            marker_start = mrk_smp,
            offset       = t2n(offset, fs),
            ep_len       = t2n(ep_len, fs),
            fs           = fs,
        )
        # first samples of the epochs kept by _make_epochs_bymarkers (epochs that fall,
        # even partially, outside the signal are dropped)
        ep_starts = mrk_smp .- t2n(offset, fs)
        ep_keep = (ep_starts .>= 1) .& (ep_starts .+ t2n(ep_len, fs) .- 1 .<= signal_len(obj))
        ep_starts = ep_starts[ep_keep]
        ep_cond = mrk_cond[mrk_idx][ep_keep]

    else
        # fixed-length epoching
        ep_n_smp = signal_len(obj)
        if !isnothing(ep_len)
            ep_len <= signal_len(obj) / sr(obj) || throw(
                ArgumentError(
                    "ep_len must be ≤ signal length ($(signal_len(obj) / sr(obj)) s).",
                ),
            )
            ep_n_smp = round(Int64, ep_len * sr(obj))
        end

        ovl_n = round(Int64, overlap * sr(obj))
        if ovl_n == 0
            epochs = _make_epochs(obj.data; ep_len = ep_n_smp)
            ep_starts = collect(1:ep_n_smp:(ep_n_smp * (size(epochs, 3) - 1) + 1))

            # remove markers that fall outside the new epoch grid (trailing samples are dropped)
            n_total = size(epochs, 2) * size(epochs, 3)
            obj.markers = _remap_markers(
                obj.markers, sr(obj), s -> 1 <= s <= n_total ? s : nothing,
            )
        else
            ovl_n < ep_n_smp || throw(
                ArgumentError("overlap must be < ep_len ($(ep_n_smp) samples after rounding)."),
            )
            # windows start every (ep_len - overlap) samples; trailing samples that do not
            # fill a complete window are dropped
            step = ep_n_smp - ovl_n
            ep_starts = collect(1:step:(signal_len(obj) - ep_n_smp + 1))
            epochs, obj.markers = _make_epochs_bymarkers(
                obj.data;
                marker       = "",
                markers      = deepcopy(obj.markers),
                marker_start = ep_starts,
                offset       = 0,
                ep_len       = ep_n_smp,
                fs           = sr(obj),
            )
        end
        ep_cond = fill("", length(ep_starts))
    end

    size(epochs, 3) == length(ep_starts) ||
        throw(ArgumentError("Internal error: epochs and their onsets do not match."))
    obj.data = epochs
    obj.epochs = DataFrame(
        :condition => ep_cond,
        :source    => fill(source, length(ep_starts)),
        :reject    => falses(length(ep_starts)),
        :onset     => onset0 .+ (ep_starts .- 1) ./ fs,
    )
    obj.time_pts, obj.epoch_time = _get_t(obj)
    obj.epoch_time .-= offset

    push!(
        obj.history,
        "epoch(obj; marker_value=$(repr(marker_value)), marker_type=$marker_type, offset=$offset, ep_len=$ep_len, overlap=$overlap, condition=$(repr(condition)), exclude_bad_segments=$exclude_bad_segments)",
    )

    # bad segments of the continuous signal (markers of type "BAD"; also those that begin
    # before an epoch, whose markers are not kept in the epochs): epochs overlapping a bad
    # segment of all channels are skipped, epochs overlapping a bad segment of some
    # channels are marked as rejected
    if exclude_bad_segments
        n_ep = size(obj.data, 2)
        drop = Int64[]
        rej = Int64[]
        for m in eachrow(mrk0)
            m[:type] == "BAD" || continue
            m1 = t2s(m[:start], fs)
            m2 = max(m1, round(Int64, (m[:start] + m[:length]) * fs))
            for (k, s0) in enumerate(ep_starts)
                (m1 <= s0 + n_ep - 1 && m2 >= s0) || continue
                push!(isempty(m[:channel]) ? drop : rej, k)
            end
        end
        drop = sort(unique(drop))
        rej = setdiff(unique(rej), drop)
        obj.epochs[rej, :reject] .= true
        !isempty(rej) && _info("$(length(rej)) epoch(s) with bad segments of some channels marked as rejected")
        if !isempty(drop)
            length(drop) < nepochs(obj) ||
                throw(ArgumentError("All epochs overlap bad segments."))
            # (delete_epoch!() keeps the epoch time axis, i.e. the offset)
            delete_epoch!(obj; ep = drop)
            _info("$(length(drop)) epoch(s) overlapping bad segments skipped")
        end
    end

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
    # epoch onsets: first retained sample of each epoch
    obj.epochs = copy(_epochs_table(obj))
    obj.epochs[!, :onset] .+= (ep_start_idx - 1) / sr(obj)
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
