export delete_channel
export delete_channel!
export keep_channel
export keep_channel!

"""
    delete_channel(obj; <keyword arguments>)

Return a copy of `obj` with the specified channel(s) removed.

All per-channel header fields are updated, including MEG sensor index lists (`:magnetometers`, `:gradiometers`), SSP channel mask and projection vectors, sensor geometry (`:coil_loc`), and MEP stimulation intensities and coil types.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s) to remove
- `del_opt::Bool=false`: set `true` only when called from `delete_optode()`; bypasses the NIRS guard that prevents direct signal-channel deletion

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function delete_channel(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    del_opt::Bool = false,
)::NeuroAnalyzer.NEURO
    # resolve channel names to integer indices
    ch_idx = get_channel(obj; ch = ch)
    ch_idx = ch_idx isa Integer ? [Int64(ch_idx)] : Vector{Int64}(ch_idx)
    isempty(ch_idx) && throw(ArgumentError("No channels selected."))

    return _delete_channel(obj, ch_idx; del_opt = del_opt)
end

# index-based core of delete_channel(); used directly by keep_channel() so
# that channels are never re-resolved by (possibly duplicated) labels
function _delete_channel(
    obj::NeuroAnalyzer.NEURO,
    ch_idx::Vector{Int64};
    del_opt::Bool = false,
)::NeuroAnalyzer.NEURO
    ch_n = nchannels(obj)
    ch_idx = unique(ch_idx)
    isempty(ch_idx) && throw(ArgumentError("No channels selected."))
    all(i -> 1 <= i <= ch_n, ch_idx) ||
        throw(ArgumentError("Channel indices must be in [1, $ch_n]."))
    length(ch_idx) < ch_n || throw(
        ArgumentError(
            "Number of channels to delete ($(length(ch_idx))) must be smaller than number of all channels ($ch_n).",
        ),
    )

    # channels that remain, in original order
    keep = setdiff(_c(ch_n), ch_idx)
    keep_mask = falses(ch_n)
    keep_mask[keep] .= true

    rec = obj.header.recording
    dt = rec[:data_type]

    # NIRS signal channels must be removed with delete_optode() (auxiliary channels may be deleted)
    if dt == "nirs" && !del_opt
        any(idx -> startswith(rec[:channel_type][idx], "nirs") && rec[:channel_type][idx] != "nirs_aux", ch_idx) &&
            throw(ArgumentError("NIRS signal channels must be deleted using delete_optode()."))
    end

    # create new dataset
    obj_tmp = deepcopy(obj)
    rec_tmp = obj_tmp.header.recording

    # locations: remove entries of deleted channels (matched by label);
    # labels also used by a kept channel are left untouched
    del_labels = setdiff(labels(obj)[ch_idx], labels(obj)[keep])
    obj_tmp.locs = obj_tmp.locs[.!in(del_labels).(obj_tmp.locs[!, :label]), :]

    # universal per-channel header vectors
    for field in (:label, :channel_type, :bad_channel, :unit)
        rec_tmp[field] = rec[field][keep]
    end

    # type-specific per-channel header vectors
    if dt in _EEG_LIKE_DATA
        for field in (:prefiltering, :transducers, :gain)
            rec_tmp[field] = rec[field][keep]
        end

    elseif dt in ("sensors", "eda", "tpt")
        rec_tmp[:prefiltering] = rec[:prefiltering][keep]

    elseif dt == "mep"
        for field in (:stimulation_intensity, :coil_type)
            rec_tmp[field] = rec[field][keep]
        end

    elseif dt in _MEG_LIKE_DATA
        for field in (:prefiltering, :coil_type)
            rec_tmp[field] = rec[field][keep]
        end

        # sensor index lists: drop deleted channels, renumber the rest
        new_idx = cumsum(keep_mask)
        for field in (:gradiometers, :magnetometers)
            rec_tmp[field] = Int64[new_idx[i] for i in rec[field] if keep_mask[i]]
        end

        # sensor geometry (absent in objects created before it was introduced)
        if haskey(rec, :coil_loc) && !isempty(rec[:coil_loc])
            rec_tmp[:coil_loc] = rec[:coil_loc][keep, :]
        end
        if haskey(rec, :coil_id) && !isempty(rec[:coil_id])
            rec_tmp[:coil_id] = rec[:coil_id][keep]
        end

        # SSP: columns of ssp_data correspond to channels marked in ssp_channels
        if haskey(rec, :ssp_channels) && any(rec[:ssp_channels])
            ssp_mask = rec[:ssp_channels]
            ssp_cols_keep = keep_mask[ssp_mask]
            if !all(ssp_cols_keep)
                _warn(
                    "$(count(.!ssp_cols_keep)) deleted channel(s) are part of SSP projectors; their weights are removed from the projection vectors.",
                )
            end
            rec_tmp[:ssp_channels] = ssp_mask[keep]
            !isempty(rec[:ssp_data]) && (rec_tmp[:ssp_data] = rec[:ssp_data][:, ssp_cols_keep])
        end

    elseif dt == "nirs"
        # only the first entries of wavelength_index / optode_pairs correspond
        # to channels; delete from the highest index down so lower indices stay valid
        for idx in sort(ch_idx; rev = true)
            if idx in eachindex(rec_tmp[:wavelength_index])
                deleteat!(rec_tmp[:wavelength_index], idx)
            end
            if idx in axes(rec_tmp[:optode_pairs], 1)
                rec_tmp[:optode_pairs] = rec_tmp[:optode_pairs][setdiff(axes(rec_tmp[:optode_pairs], 1), idx), :]
            end
        end
    end

    rec_tmp[:channel_order] = _sort_channels(rec_tmp[:channel_type])
    obj_tmp.data = obj_tmp.data[keep, :, :]

    # markers refer to channels by number (0 = all channels): renumber them,
    # markers of deleted channels are removed
    if DataFrames.nrow(obj_tmp.markers) > 0
        new_idx = cumsum(keep_mask)
        mrk_ch = obj_tmp.markers[!, :channel]
        on_deleted = [c != 0 && !keep_mask[c] for c in mrk_ch]
        any(on_deleted) && _info("$(count(on_deleted)) marker(s) of deleted channels removed.")
        obj_tmp.markers = obj_tmp.markers[.!on_deleted, :]
        obj_tmp.markers[!, :channel] = [c == 0 ? 0 : new_idx[c] for c in obj_tmp.markers[!, :channel]]
    end

    push!(obj_tmp.history, "delete_channel(obj; ch=$(labels(obj)[sort(ch_idx)]))")

    return obj_tmp
end

"""
    delete_channel!(obj; <keyword arguments>)

Delete channel(s) in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s) to remove
- `del_opt::Bool=false`: set `true` only when called from `delete_optode()`; bypasses the NIRS guard that prevents direct signal-channel deletion

# Returns

- `Nothing`
"""
function delete_channel!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    del_opt::Bool = false,
)::Nothing
    obj_tmp = delete_channel(obj; ch = ch, del_opt = del_opt)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj.locs = obj_tmp.locs
    obj.markers = obj_tmp.markers

    return nothing
end

"""
    keep_channel(obj; <keyword arguments>)

Return a copy of `obj` retaining only the specified channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s) to keep
- `del_opt::Bool=false`: passed to `delete_channel()`; set `true` only when called from optode-handling functions

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function keep_channel(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    del_opt::Bool = false,
)::NeuroAnalyzer.NEURO
    ch_n = nchannels(obj)

    # resolve channel names to integer indices
    ch_idx = get_channel(obj; ch = ch)
    ch_idx = unique(ch_idx isa Integer ? [Int64(ch_idx)] : Vector{Int64}(ch_idx))
    isempty(ch_idx) && throw(ArgumentError("No channels selected."))

    # nothing to remove: still return a copy, never the input object
    length(ch_idx) == ch_n && return deepcopy(obj)

    # delete by index: re-resolving by label would select wrong channels
    # when labels are duplicated
    obj_tmp = _delete_channel(obj, setdiff(_c(ch_n), ch_idx); del_opt = del_opt)
    # record the user's call rather than the internal deletion
    obj_tmp.history[end] = "keep_channel(obj; ch=$(labels(obj)[sort(ch_idx)]))"

    return obj_tmp
end

"""
    keep_channel!(obj; <keyword arguments>)

Keep only the specified channel(s) in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s) to keep
- `del_opt::Bool=false`: passed to `delete_channel()`; set `true` only when called from optode-handling functions

# Returns

- `Nothing`
"""
function keep_channel!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    del_opt::Bool = false,
)::Nothing
    obj_tmp = keep_channel(obj; ch = ch, del_opt = del_opt)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj.locs = obj_tmp.locs
    obj.markers = obj_tmp.markers

    return nothing
end
