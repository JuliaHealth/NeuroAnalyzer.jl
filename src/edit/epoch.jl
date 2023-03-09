export epoch
export epoch!
export epoch_time
export epoch_time!
export extract_epoch

"""
    epoch(obj; marker, ep_offset, ep_n, ep_len)

Split OBJ into epochs. Return signal that is split either by markers (if specified), by epoch length or by number of epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `marker::String="": marker name to split at
- `ep_offset::Int64=0": time offset (in samples) for marker-based epoching (each epoch time will start at marker time - ep_offset)
- `ep_n::Union{Int64, Nothing}=nothing`: number of epochs
- `ep_len::Union{Int64, Nothing}`=nothing: epoch length in samples

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function epoch(obj::NeuroAnalyzer.NEURO; marker::String="", ep_offset::Real=0, ep_n::Union{Int64, Nothing}=nothing, ep_len::Union{Int64, Nothing}=nothing)

    obj_new = deepcopy(obj)

    if marker != ""
        # split by markers
        if obj.header.has_markers == true
            ep_len === nothing && throw(ArgumentError("ep_len must be specified."))
            ep_offset == 0 && throw(ArgumentError("ep_offset must be specified."))
            _check_markers(obj, marker)

            # get marker positions
            marker_idx = []
            for idx in 1:length(obj.markers[!, :description])
                obj_new.markers[idx, :description] == marker && push!(marker_idx, idx)
            end
            marker_start = obj_new.markers[!, :start][marker_idx]

            # split into epochs
            epochs, obj_new.markers = _make_epochs_bymarkers(obj.data, markers=obj_new.markers, marker_start=marker_start, ep_offset=ep_offset, ep_len=ep_len)
        else
            throw(ArgumentError("OBJ does not contain markers."))
        end
    else
        # split by ep_len or ep_n
        epochs = _make_epochs(obj.data, ep_n=ep_n, ep_len=ep_len)

        # delete markers outside epochs
        for marker_idx in nrow(obj_new.markers):-1:1
            obj_new.markers[marker_idx, :start] in 1:size(epochs, 2) * size(epochs, 3) || deleteat!(obj_new.markers, marker_idx)
        end
    end

    # create new dataset
    ep_n = size(epochs, 3)
    epoch_duration_samples = size(epochs, 2)
    epoch_duration_seconds = size(epochs, 2) / obj.header[:sampling_rate]
    duration_samples = size(epochs, 2) * size(epochs, 3)
    duration_seconds = duration_samples / obj.header[:sampling_rate]
    time_pts = collect(0:(1 / obj.header[:sampling_rate]):duration_seconds)
    time_pts = time_pts[1:(end - 1)]

    # update signal
    obj_new.data = epochs

    # update time
    obj_new.time_pts = time_pts

    # update epochs time
    fs = sr(obj_new)
    new_epochs_time = linspace(-s2t(ep_offset, fs), epoch_duration_seconds - s2t(ep_offset, fs), epoch_duration_samples)
    obj_new.epoch_time = new_epochs_time

    # update header
    obj_new.header.recording[:duration_samples] = duration_samples
    obj_new.header.recording[:duration_seconds] = duration_seconds
    obj_new.header.recording[:epoch_n] = ep_n
    obj_new.header.recording[:epoch_duration_samples] = epoch_duration_samples
    obj_new.header.recording[:epoch_duration_seconds] = epoch_duration_seconds

    reset_components!(obj_new)
    push!(obj_new.header.recording.history, "epoch(OBJ, ep_n=$ep_n, ep_len=$ep_len)")

    return obj_new
end

"""
    epoch!(obj; marker, ep_offset, ep_n, ep_len)

Split OBJ into epochs. Return signal that is split either by markers (if specified), by epoch length or by number of epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `marker::String="": marker name to split at
- `ep_offset::Int64=0": time offset (in samples) for marker-based epoching (each epoch time will start at marker time - ep_offset)
- `ep_n::Union{Int64, Nothing}=nothing`: number of epochs
- `ep_len::Union{Int64, Nothing}`=nothing: epoch length in samples
"""
function epoch!(obj::NeuroAnalyzer.NEURO; marker::String="", ep_offset::Real=0, ep_n::Union{Int64, Nothing}=nothing, ep_len::Union{Int64, Nothing}=nothing)

    obj_tmp = epoch(obj, marker=marker, ep_offset=ep_offset, ep_n=ep_n, ep_len=ep_len)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    obj.time_pts = obj_tmp.time_pts
    obj.epoch_time = obj_tmp.epoch_time
    reset_components!(obj)

    return nothing
end

"""
    epoch_time(obj; ts)

Edit epochs time start.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ts::Real`: time start in seconds

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function epoch_time(obj::NeuroAnalyzer.NEURO; ts::Real)

    ep_len = epoch_len(obj)
    fs = sr(obj)
    new_epochs_time = linspace(ts, ts + (ep_len / fs), ep_len)
    obj_new = deepcopy(obj)
    obj_new.epoch_time = new_epochs_time

    push!(obj_new.header.recording.history, "epoch_time(OBJ, ts=$ts)")

    return obj_new
end

"""
    epoch_time!(obj; ts)

Edit OBJ epochs time start.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ts::Real`: time start in seconds

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function epoch_time!(obj::NeuroAnalyzer.NEURO; ts::Real)

    obj_tmp = epoch_time(obj, ts=ts)
    obj.header = obj_tmp.header
    obj.epoch_time = obj_tmp.epoch_time

    return nothing
end

"""
    extract_epoch(obj; epoch)

Extract OBJ epoch.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `epoch::Int64`: epoch index

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function extract_epoch(obj::NeuroAnalyzer.NEURO; epoch::Int64)

    _check_epochs(obj, epoch)

    s_new = reshape(obj.data[:, :, epoch], channel_n(obj), signal_len(obj), 1)
    obj_new = deepcopy(obj)
    obj_new.data = s_new
    obj_new.epoch_time = obj.epoch_time
    obj_new.header.recording[:epoch_n] = 1
    obj_new.header.recording[:duration_samples] = obj_new.header.recording[:epoch_duration_samples]
    obj_new.header.recording[:duration_seconds] = obj_new.header.recording[:epoch_duration_seconds]

    reset_components!(obj_new)
    push!(obj_new.header.recording.history, "extract_epoch(OBJ, epoch=$epoch)")

    return obj_new
end

