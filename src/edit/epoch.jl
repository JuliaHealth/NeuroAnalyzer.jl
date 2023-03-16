export epoch
export epoch!
export epoch_time
export epoch_time!

"""
    epoch(obj; marker, ep_offset, ep_n, ep_len)

Split OBJ into epochs. Return signal that is split either by markers (if specified), by epoch length or by number of epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `marker::String="": marker name to split at
- `ep_offset::Int64=0": time offset (in samples) for marker-based epoching (each epoch time will start at `marker time - ep_offset`)
- `ep_n::Union{Int64, Nothing}=nothing`: number of epochs
- `ep_len::Union{Int64, Nothing}`=nothing: epoch length in samples

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function epoch(obj::NeuroAnalyzer.NEURO; marker::String="", ep_offset::Real=0, ep_n::Union{Int64, Nothing}=nothing, ep_len::Union{Int64, Nothing}=nothing)

    obj_new = deepcopy(obj)

    if marker != ""
        # split by markers
        if _has_markers(obj) == true
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

    # update signal
    obj_new.data = epochs

    # update time
    obj_new.time_pts = round.(collect(obj_new.time_pts[1]:(1 / sr(obj_new)):((size(obj_new.data, 2) * size(obj_new.data, 3)) / sr(obj_new))) .- (1 / sr(obj_new)), digits=3)
    obj_new.epoch_time = round.(collect(obj_new.time_pts[1]:(1 / sr(obj_new)):(size(obj_new.data, 2) / sr(obj_new))) .- (1 / sr(obj_new)), digits=3)

    reset_components!(obj_new)
    push!(obj_new.history, "epoch(OBJ, marker=$marker, ep_offset=$ep_offset, ep_n=$ep_n, ep_len=$ep_len)")

    return obj_new

end

"""
    epoch!(obj; marker, ep_offset, ep_n, ep_len)

Split OBJ into epochs. Return signal that is split either by markers (if specified), by epoch length or by number of epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `marker::String="": marker name to split at
- `ep_offset::Int64=0": time offset (in samples) for marker-based epoching (each epoch time will start at `marker time - ep_offset`)
- `ep_n::Union{Int64, Nothing}=nothing`: number of epochs
- `ep_len::Union{Int64, Nothing}`=nothing: epoch length in samples
"""
function epoch!(obj::NeuroAnalyzer.NEURO; marker::String="", ep_offset::Real=0, ep_n::Union{Int64, Nothing}=nothing, ep_len::Union{Int64, Nothing}=nothing)

    obj_new = epoch(obj, marker=marker, ep_offset=ep_offset, ep_n=ep_n, ep_len=ep_len)
    obj.header = obj_new.header
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time

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

    push!(obj_new.history, "epoch_time(OBJ, ts=$ts)")

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

    obj_new = epoch_time(obj, ts=ts)
    obj.history = obj_new.history
    obj.epoch_time = obj_new.epoch_time

    return nothing
    
end
