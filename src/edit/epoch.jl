export epoch
export epoch!
export epoch_ts
export epoch_ts!
export subepoch
export subepoch!

"""
    epoch(obj; <keyword arguments>)

Split into epochs. Return signal that is split either by markers (if specified), by epoch length or by number of epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `marker::String=""`: marker description to split at
- `offset::Real=0`: time offset (in seconds) for marker-based epoching (each epoch time will start at `marker time - offset`)
- `ep_n::Union{Int64, Nothing}=nothing`: number of epochs
- `ep_len::Union{Real, Nothing}=nothing`: epoch length in seconds

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function epoch(obj::NeuroAnalyzer.NEURO; marker::String="", offset::Real=0, ep_n::Union{Int64, Nothing}=nothing, ep_len::Union{Real, Nothing}=nothing)

    obj_new = deepcopy(obj)

    # create ID for epochs
    if marker != ""
        epoch_id = marker
    elseif marker == ""
        if ep_len !== nothing
            epoch_id = "length_$(ep_len)s"
        else
            epoch_id = "n_$(ep_n)"
        end
    end

    if marker != ""
        # split by markers
        @assert _has_markers(obj) "OBJ does not contain markers."
        @assert ep_len !== nothing "ep_len must be specified."
        _check_markers(obj, marker)

        # get marker positions
        @assert any(obj_new.markers[!, :description] .== marker) "OBJ does not contain marker $marker."
        mrk_idx = getindex.(findall(obj_new.markers[!, :description] .== marker))
        mrk_start = obj_new.markers[mrk_idx, :start]
        mrk_len = obj_new.markers[mrk_idx, :length]

        # remove markers that would be before signal start
        for idx in length(mrk_start):-1:1
            if mrk_start[idx] - offset < obj.time_pts[1]
                deleteat!(mrk_start, idx)
                deleteat!(mrk_len, idx)
            end
        end

        @assert offset + ep_len >= maximum(mrk_len) "offset + ep_len must be ≥ $(maximum(mrk_len)) (maximum marker length)."

        # split into epochs
        epochs, obj_new.markers = _make_epochs_bymarkers(obj_new.data, marker=marker, markers=deepcopy(obj_new.markers), marker_start=round.(Int64, mrk_start * sr(obj)), offset=round(Int64, offset * sr(obj)), ep_len=round(Int64, ep_len * sr(obj)), fs=sr(obj))

    else
        if ep_len !== nothing
            @assert ep_len <= signal_len(obj) / sr(obj) "ep_len must be ≤ signal length ($(signal_len(obj) / sr(obj)))."
            ep_len = round(Int64, ep_len * sr(obj))
        end
        # split by ep_len or ep_n
        epochs = _make_epochs(obj.data, ep_n=ep_n, ep_len=ep_len)

        # delete markers outside epochs
        for marker_idx in nrow(obj_new.markers):-1:1
            round(Int64, sr(obj) * obj_new.markers[marker_idx, :start]) in 0:size(epochs, 2) * size(epochs, 3) || deleteat!(obj_new.markers, marker_idx)
        end
    end

    # update signal
    obj_new.data = epochs

    # create ID for epochs
    obj_new.header.recording[:epoch_id] = epoch_id

    # bad channels
    obj_new.header.recording[:bad_channel] = zeros(Bool, size(epochs, 1), size(epochs, 3))

    # update time
    obj_new.time_pts, obj_new.epoch_time = _get_t(obj_new)
    obj_new.epoch_time .-= offset

    reset_components!(obj_new)
    push!(obj_new.history, "epoch(OBJ, marker=$marker, offset=$offset, ep_n=$ep_n, ep_len=$ep_len)")

    return obj_new

end

"""
    epoch!(obj; <keyword arguments>)

Split into epochs. Return signal that is split either by markers (if specified), by epoch length or by number of epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `marker::String=""`: marker description to split at
- `offset::Real=0`: time offset (in seconds) for marker-based epoching (each epoch time will start at `marker time - offset`)
- `ep_n::Union{Int64, Nothing}=nothing`: number of epochs
- `ep_len::Union{Real, Nothing}=nothing`: epoch length in seconds
"""
function epoch!(obj::NeuroAnalyzer.NEURO; marker::String="", offset::Real=0, ep_n::Union{Int64, Nothing}=nothing, ep_len::Union{Real, Nothing}=nothing)

    obj_new = epoch(obj, marker=marker, offset=offset, ep_n=ep_n, ep_len=ep_len)
    obj.header = obj_new.header
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time

    return nothing

end

"""
    epoch_ts(obj; <keyword arguments>)

Edit epochs time start.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ts::Real`: time start in seconds

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function epoch_ts(obj::NeuroAnalyzer.NEURO; ts::Real)

    obj_new = deepcopy(obj)
    obj_new.epoch_time .+= ts

    push!(obj_new.history, "epoch_ts(OBJ, ts=$ts)")

    return obj_new

end

"""
    epoch_ts!(obj; <keyword arguments>)

Edit OBJ epochs time start.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ts::Real`: time start in seconds

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function epoch_ts!(obj::NeuroAnalyzer.NEURO; ts::Real)

    obj_new = epoch_ts(obj, ts=ts)
    obj.history = obj_new.history
    obj.epoch_time = obj_new.epoch_time

    return nothing

end

"""
    subepoch(obj; <keyword arguments>)

Extract sub-epochs with a reduced time range.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ep_start::Real`: sub-epoch start
- `ep_end::Real`: sub-epoch end

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function subepoch(obj::NeuroAnalyzer.NEURO; ep_start::Real, ep_end::Real)

    obj_new = deepcopy(obj)
    ep_time = obj.epoch_time

    @assert ep_start >= ep_time[1] "ep_start must be ≥ $(ep_time[1])."
    @assert ep_end <= ep_time[end] "ep_end must be ≤ $(ep_time[end])."

    ep_start_idx = vsearch(ep_start, ep_time)
    ep_end_idx = vsearch(ep_end, ep_time)

    # update signal
    obj_new.data = obj.data[:, ep_start_idx:ep_end_idx, :]

    # update time points
    obj_new.epoch_time = ep_time[ep_start_idx:ep_end_idx]
    obj_new.time_pts, _ = NeuroAnalyzer._get_t(obj_new)

    # delete markers outside epochs
    ep_tps = NeuroAnalyzer._epochs_tps(obj)
    ep_tps[2, :] = ep_tps[1, :] .+ ep_end
    ep_tps[1, :] .+= ep_start
    mrk_start = obj.markers[!, :start]
    mrk_epoch = NeuroAnalyzer._markers_epochs(obj)
    for mrk_idx in length(mrk_start):-1:1
        if mrk_start[mrk_idx] < ep_tps[1, mrk_epoch[mrk_idx]] || mrk_start[mrk_idx] > ep_tps[2, mrk_epoch[mrk_idx]]
            deleteat!(obj_new.markers, mrk_idx)
            deleteat!(mrk_epoch, mrk_idx)
        end
    end

    # shift markers time points
    mrk_start = deepcopy(obj_new.markers[!, :start])
    for mrk_idx in eachindex(mrk_start)
        mrk_start[mrk_idx] -= (ep_start + ((mrk_epoch[mrk_idx] - 1) * (ep_start + (obj.time_pts[(epoch_len(obj))] - ep_end))))
    end
    mrk_start = round.(mrk_start, digits=3)

    obj_new.markers[!, :start] = mrk_start

    reset_components!(obj_new)
    push!(obj_new.history, "subepoch(OBJ, ep_start=$ep_start, ep_start=$ep_end)")

    return obj_new

end

"""
    subepoch!(obj; <keyword arguments>)

Extract sub-epochs with a reduced time range.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ep_start::Real`: sub-epoch start
- `ep_end::Real`: sub-epoch end
"""
function subepoch!(obj::NeuroAnalyzer.NEURO; ep_start::Real, ep_end::Real)

    obj_new = subepoch(obj, ep_start=ep_start, ep_end=ep_end)
    obj.header = obj_new.header
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time

    return nothing

end
