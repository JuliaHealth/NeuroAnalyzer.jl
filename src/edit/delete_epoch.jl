export delete_epoch
export delete_epoch!
export keep_epoch
export keep_epoch!

"""
    delete_epoch(obj; epoch)

Remove epoch(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) to be removed

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function delete_epoch(obj::NeuroAnalyzer.NEURO; epoch::Union{Int64, Vector{Int64}, AbstractRange})

    epoch_n(obj) == 1 && throw(ArgumentError("You cannot delete the last epoch."))
    typeof(epoch) <: AbstractRange && (epoch = collect(epoch))
    length(epoch) == epoch_n(obj) && throw(ArgumentError("You cannot delete all epochs."))
    length(epoch) > 1 && (epoch = sort!(epoch, rev=true))
    _check_epochs(obj, epoch)

    obj_new = deepcopy(obj)

    # remove epoch
    obj_new.data = obj_new.data[:, :, setdiff(1:end, (epoch))]

    # remove markers within deleted epochs and shift markers after the deleted epoch
    for ep_idx in epoch
        t1, t2 = _epoch2s(obj, ep_idx)
        obj_new.markers = _delete_markers(obj_new.markers, (t1, t2))
        obj_new.markers = _shift_markers(obj_new.markers, t1, length(t1:t2))
    end

    # update headers
    obj_new.header.recording[:epoch_n] -= length(epoch)
    ep_n = obj_new.header.recording[:epoch_n]
    obj_new.header.recording[:duration_samples] = ep_n * size(obj.data, 2)
    obj_new.header.recording[:duration_seconds] = round((ep_n * size(obj.data, 2)) / sr(obj), digits=2)
    obj_new.header.recording[:epoch_duration_samples] = size(obj.data, 2)
    obj_new.header.recording[:epoch_duration_seconds] = round(size(obj.data, 2) / sr(obj), digits=2)

    reset_components!(obj_new)
    push!(obj_new.header.history, "delete_epoch(OBJ, $epoch)")
    
    return obj_new
end

"""
    delete_epoch!(obj; epoch)

Remove epoch(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) to be removed
"""
function delete_epoch!(obj::NeuroAnalyzer.NEURO; epoch::Union{Int64, Vector{Int64}, AbstractRange})

    obj_tmp = delete_epoch(obj, epoch=epoch)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    obj.time_pts = obj_tmp.time_pts
    obj.markers = obj_tmp.markers
    reset_components!(obj)

    return nothing
end

"""
    keep_epoch(obj; epoch)

Keep epoch(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) to keep

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function keep_epoch(obj::NeuroAnalyzer.NEURO; epoch::Union{Int64, Vector{Int64}, AbstractRange})

    epoch_n(obj) == 1 && throw(ArgumentError("contains only one epoch."))
    typeof(epoch) <: AbstractRange && (epoch = collect(epoch))
    length(epoch) > 1 && (epoch = sort!(epoch, rev=true))
    _check_epochs(obj, epoch)

    epoch_list = collect(1:epoch_n(obj))
    epoch_to_remove = setdiff(epoch_list, epoch)

    length(epoch_to_remove) > 1 && (epoch_to_remove = sort!(epoch_to_remove, rev=true))

    obj_new = delete_epoch(obj, epoch=epoch_to_remove)
    reset_components!(obj_new)
    push!(obj_new.header.history, "keep_epoch(OBJ, $epoch)")    

    return obj_new
end

"""
    keep_epoch!(obj; epoch)

Keep OBJ epoch(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) to keep
"""
function keep_epoch!(obj::NeuroAnalyzer.NEURO; epoch::Union{Int64, Vector{Int64}, AbstractRange})

    obj_tmp = keep_epoch(obj, epoch=epoch)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    obj.time_pts = obj_tmp.time_pts
    obj.markers = obj_tmp.markers
    reset_components!(obj)

    return nothing
end
