export trim
export trim!

"""
    trim(obj; segment, remove_epochs)

Trim OBJ signal by removing parts of the signal.

# Arguments

- OBJ
- `segment::Tuple{Int64, Int64}`: segment to be removed (from, to) in samples
- `remove_epochs::Bool=true`: if true, remove epochs containing signal to trim or remove signal and re-epoch trimmed signal

# Returns

- OBJ
"""
function trim(obj::NeuroAnalyzer.NEURO; segment::Tuple{Int64, Int64}, remove_epochs::Bool=true)

    if remove_epochs == false
        obj_new = deepcopy(obj)
        epoch_n(obj) > 1 && (epoch!(obj_new, ep_n=1))
        _check_segment(obj_new, segment[1], segment[2])
        obj_new.data = trim(obj_new.data, segment=segment)
        t_trimmed = collect(0:(1 / sr(obj)):(size(obj_new.data, 2) / sr(obj)))[1:(end - 1)]
        obj_new.time_pts = t_trimmed
        obj_new.epoch_time = t_trimmed .+ obj.epoch_time[1]
        if epoch_n(obj) > 1
            if epoch_len(obj) <= signal_len(obj_new)
                epoch!(obj_new, epoch_len=epoch_len(obj))
            else
                _info("Cannot apply original epoch length, returning single-epoch OBJ.")
            end
        end
        obj_new.markers = _delete_markers(obj_new.markers, segment)
        obj_new.markers = _shift_markers(obj_new.markers, segment[1], length(segment[1]:segment[2]))
    else
        epoch_n(obj) == 1 && throw(ArgumentError("OBJ has only one epoch, cannot use remove_epochs=true."))
        epochs = _s2epoch(obj, segment[1], segment[2])
        _info("Removing epochs: $epochs.")
        obj_new = delete_epoch(obj, epoch=epochs)
    end

    reset_components!(obj_new)
    push!(obj_new.history, "trim(OBJ, segment=$segment, remove_epochs=$remove_epochs)")

    return obj_new
end

"""
    trim!(obj; segment, remove_epochs)

Trim OBJ signal by removing parts of the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `segment::Tuple{Int64, Int64}`: segment to be removed (from, to) in samples
- `remove_epochs::Bool=true`: remove epochs containing signal to trim (remove_epochs=true) or remove signal and remove epoching
"""
function trim!(obj::NeuroAnalyzer.NEURO; segment::Tuple{Int64, Int64}, remove_epochs::Bool=true)

    obj_new = trim(obj, segment=segment, remove_epochs=remove_epochs)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time

    return nothing

end
