export delete_channel
export delete_channel!
export keep_channel
export keep_channel!

"""
    delete_channel(obj; <keyword arguments>)

Return a copy of `obj` with the specified channel(s) removed.

All per-channel header fields are updated, including MEG sensor index lists (`:magnetometers`, `:gradiometers`), SSP channel mask and projection vectors, sensor geometry (`:coil_loc`), and MEP stimulation intensities and coil types.

For NIRS objects, a deleted intensity or OD channel is deleted together with the channels of the same type and optode pair at the other wavelengths (the conversion to concentration requires all wavelengths of a pair); use `delete_pair()` to remove all channels of a source-detector pair. The optode labels and locations are kept (use `delete_optode()` to remove optodes).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s) to remove

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function delete_channel(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
)::NeuroAnalyzer.NEURO
    # resolve channel names to integer indices
    ch_idx = get_channel(obj; ch = ch)
    ch_idx = ch_idx isa Integer ? [Int64(ch_idx)] : Vector{Int64}(ch_idx)
    isempty(ch_idx) && throw(ArgumentError("No channels selected."))

    return _delete_channel(obj, ch_idx)
end

# NIRS: add to `ch_idx` the intensity / OD channels of the same type and optode pair at the
# other wavelengths (channels described by :wavelength_index / :optode_pairs come first)
function _nirs_pair_wl(
    obj::NeuroAnalyzer.NEURO,
    ch_idx::Vector{Int64};
    msg::String,
)::Vector{Int64}
    rec = obj.header.recording
    n_opt = length(rec[:wavelength_index])
    ct = rec[:channel_type]
    chp = rec[:optode_pairs]
    add = Int64[]
    for i in ch_idx
        (i <= n_opt && ct[i] in ("nirs_int", "nirs_od")) || continue
        for j = 1:n_opt
            (j in ch_idx || j in add) && continue
            ct[j] == ct[i] && chp[j, :] == chp[i, :] && push!(add, j)
        end
    end
    isempty(add) && return ch_idx
    _info("$msg the other wavelengths of the same optode pairs: $(labels(obj)[sort(add)])")

    return vcat(ch_idx, add)
end

# index-based core of delete_channel(); used directly by keep_channel() so
# that channels are never re-resolved by (possibly duplicated) labels
function _delete_channel!(
    obj::NeuroAnalyzer.NEURO,
    ch_idx::Vector{Int64},
)::Nothing
    ch_n = nchannels(obj)
    ch_idx = unique(ch_idx)
    isempty(ch_idx) && throw(ArgumentError("No channels selected."))
    all(i -> 1 <= i <= ch_n, ch_idx) ||
        throw(ArgumentError("Channel indices must be in [1, $ch_n]."))
    # NIRS: an intensity or OD channel is deleted together with the channels of the same
    # type and optode pair at the other wavelengths (OD -> concentration needs all of them)
    datatype(obj) == "nirs" &&
        (ch_idx = _nirs_pair_wl(obj, ch_idx; msg = "Also deleting"))
    length(ch_idx) < ch_n || throw(
        ArgumentError(
            "Number of channels to delete ($(length(ch_idx))) must be smaller than number of all channels ($ch_n).",
        ),
    )

    # channels that remain, in original order
    keep = setdiff(_c(ch_n), ch_idx)
    # history entry: labels before the header is modified
    del_lbl = labels(obj)[sort(ch_idx)]

    # locations: remove entries of deleted channels (matched by label);
    # labels also used by a kept channel are left untouched
    del_labels = setdiff(labels(obj)[ch_idx], labels(obj)[keep])
    obj.locs = obj.locs[.!in(del_labels).(obj.locs[!, :label]), :]

    _select_channels!(obj, keep)

    # markers refer to channels by label (empty = all channels): deleted channels are
    # removed from them, markers of deleted channels only are removed
    n_del = _markers_delete_channels!(obj.markers, del_labels)
    n_del > 0 && _info("$n_del marker(s) of deleted channels removed")

    _history!(obj, "delete_channel", (; ch = del_lbl))

    return nothing
end

# keep the channels `idx` (indices, in the new order: a subset and/or a permutation of the
# channels): data and all per-channel header fields; locations and markers refer to the
# channels by label and are not changed
function _select_channels!(obj::NeuroAnalyzer.NEURO, idx::Vector{Int64})::Nothing
    ch_n = nchannels(obj)
    allunique(idx) && all(i -> 1 <= i <= ch_n, idx) ||
        throw(ArgumentError("Channel indices must be unique and in [1, $ch_n]."))
    # new position of each old channel (0: removed)
    new_pos = zeros(Int64, ch_n)
    new_pos[idx] = 1:length(idx)

    rec = obj.header.recording
    dt = rec[:data_type]

    # universal per-channel header vectors
    _rec_filters!(rec)
    for field in (:label, :channel_type, :bad_channel, :unit, :lowpass, :highpass, :notch)
        rec[field] = rec[field][idx]
    end

    # type-specific per-channel header vectors
    if dt in _EEG_LIKE_DATA
        for field in (:transducers, :gain)
            rec[field] = rec[field][idx]
        end

    elseif dt == "mep"
        for field in (:stimulation_intensity, :coil_type)
            rec[field] = rec[field][idx]
        end

    elseif dt in _MEG_LIKE_DATA
        rec[:coil_type] = rec[:coil_type][idx]

        # sensor index lists: removed channels dropped, the others renumbered
        for field in (:gradiometers, :magnetometers)
            rec[field] = sort(Int64[new_pos[i] for i in rec[field] if new_pos[i] > 0])
        end

        # sensor geometry (absent in objects created before it was introduced)
        if haskey(rec, :coil_loc) && !isempty(rec[:coil_loc])
            rec[:coil_loc] = rec[:coil_loc][idx, :]
        end
        if haskey(rec, :coil_id) && !isempty(rec[:coil_id])
            rec[:coil_id] = rec[:coil_id][idx]
        end

        # SSP: the columns of ssp_data are the channels marked in ssp_channels (in the
        # channel order)
        if haskey(rec, :ssp_channels) && any(rec[:ssp_channels])
            ssp_mask = rec[:ssp_channels]
            ssp_col = cumsum(ssp_mask)
            cols = Int64[ssp_col[i] for i in idx if ssp_mask[i]]
            if length(cols) < count(ssp_mask)
                _warn(
                    "$(count(ssp_mask) - length(cols)) deleted channel(s) are part of SSP projectors; their weights are removed from the projection vectors.",
                )
            end
            rec[:ssp_channels] = ssp_mask[idx]
            !isempty(rec[:ssp_data]) && (rec[:ssp_data] = rec[:ssp_data][:, cols])
        end

    elseif dt == "nirs"
        # the first entries of wavelength_index / optode_pairs describe the first channels
        # (intensity or OD); they must stay the first channels
        n_opt = length(rec[:wavelength_index])
        opt = Base.filter(<=(n_opt), idx)
        idx[1:length(opt)] == opt || throw(
            ArgumentError(
                "The NIRS channels described by optode pairs (channels 1:$n_opt) must stay the first channels.",
            ),
        )
        rec[:wavelength_index] = rec[:wavelength_index][opt]
        rec[:optode_pairs] = rec[:optode_pairs][opt, :]
    end

    rec[:channel_order] = _sort_channels(rec[:channel_type])
    obj.data = obj.data[idx, :, :]

    return nothing
end

# out-of-place version of _delete_channel!()
function _delete_channel(
    obj::NeuroAnalyzer.NEURO,
    ch_idx::Vector{Int64},
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    _delete_channel!(obj_new, ch_idx)

    return obj_new
end

"""
    delete_channel!(obj; <keyword arguments>)

Delete channel(s) in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s) to remove

# Returns

- `Nothing`
"""
function delete_channel!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
)::Nothing
    # resolve channel names to integer indices
    ch_idx = get_channel(obj; ch = ch)
    ch_idx = ch_idx isa Integer ? [Int64(ch_idx)] : Vector{Int64}(ch_idx)
    isempty(ch_idx) && throw(ArgumentError("No channels selected."))
    _delete_channel!(obj, ch_idx)

    return nothing
end

"""
    keep_channel(obj; <keyword arguments>)

Return a copy of `obj` retaining only the specified channel(s).

For NIRS objects, keeping an intensity or OD channel keeps the channels of the same type and optode pair at the other wavelengths.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s) to keep

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function keep_channel(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    keep_channel!(obj_new; ch = ch)

    return obj_new
end

"""
    keep_channel!(obj; <keyword arguments>)

Keep only the specified channel(s) in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s) to keep

# Returns

- `Nothing`
"""
function keep_channel!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
)::Nothing
    ch_n = nchannels(obj)
    # resolve channel names to integer indices
    ch_idx = get_channel(obj; ch = ch)
    ch_idx = unique(ch_idx isa Integer ? [Int64(ch_idx)] : Vector{Int64}(ch_idx))
    isempty(ch_idx) && throw(ArgumentError("No channels selected."))
    # NIRS: keeping an intensity or OD channel keeps the channels of the same type and
    # optode pair at the other wavelengths
    datatype(obj) == "nirs" && (ch_idx = _nirs_pair_wl(obj, ch_idx; msg = "Also keeping"))
    # nothing to remove
    length(ch_idx) == ch_n && return nothing
    # labels of the kept channels (history), before the header is modified
    keep_lbl = labels(obj)[sort(ch_idx)]
    # delete by index: re-resolving by label would select wrong channels
    # when labels are duplicated
    n0 = length(obj.history)
    _delete_channel!(obj, setdiff(_c(ch_n), ch_idx))
    # record the user's call rather than the internal deletion
    _history!(obj, "keep_channel", (; ch = keep_lbl); from = n0)

    return nothing
end
