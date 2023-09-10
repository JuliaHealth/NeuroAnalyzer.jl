export delete_epoch
export delete_epoch!
export keep_epoch
export keep_epoch!

"""
    delete_epoch(obj; ep)

Remove epoch(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ep::Union{Int64, Vector{Int64}, <:AbstractRange}`: epoch number(s) to be removed

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function delete_epoch(obj::NeuroAnalyzer.NEURO; ep::Union{Int64, Vector{Int64}, <:AbstractRange})

    @assert nepochs(obj) > 1 "You cannot delete the last epoch."
    typeof(ep) <: AbstractRange && (ep = collect(ep))
    @assert length(ep) < nepochs(obj) "Number of epochs to delete ($(length(ep))) must be smaller than number of all epochs."
    length(ep) > 1 && (ep = sort!(ep, rev=true))
    _check_epochs(obj, ep)

    obj_new = deepcopy(obj)

    # remove epoch
    obj_new.data = obj_new.data[:, :, setdiff(1:end, (ep))]

    # remove markers within deleted epochs and shift markers after the deleted epoch
    for ep_idx in ep
        t1, t2 = _epoch2s(obj, ep_idx)
        obj_new.markers = _delete_markers(obj_new.markers, (t1, t2), sr(obj))
        obj_new.markers = _shift_markers(obj_new.markers, t1, length(t1:t2), sr(obj))
    end

    # update time
    obj_new.time_pts, obj_new.epoch_time = _get_t(obj_new)

    reset_components!(obj_new)
    push!(obj_new.history, "delete_epoch(OBJ, $ep)")
    
    return obj_new

end

"""
    delete_epoch!(obj; ep)

Remove epoch(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ep::Union{Int64, Vector{Int64}, <:AbstractRange}`: epoch number(s) to be removed
"""
function delete_epoch!(obj::NeuroAnalyzer.NEURO; ep::Union{Int64, Vector{Int64}, <:AbstractRange})

    obj_new = delete_epoch(obj, ep=ep)
    obj.header = obj_new.header
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.time_pts = obj_new.time_pts
    obj.markers = obj_new.markers

    return nothing
end

"""
    keep_epoch(obj; ep)

Keep epoch(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ep::Union{Int64, Vector{Int64}, <:AbstractRange}`: epoch number(s) to keep

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function keep_epoch(obj::NeuroAnalyzer.NEURO; ep::Union{Int64, Vector{Int64}, <:AbstractRange})

    @assert nepochs(obj) > 1 "OBJ contains only one epoch."
    typeof(ep) <: AbstractRange && (ep = collect(ep))
    length(ep) > 1 && (ep = sort!(ep, rev=true))
    _check_epochs(obj, ep)

    ep_list = collect(1:nepochs(obj))
    ep_to_remove = setdiff(ep_list, ep)

    length(ep_to_remove) > 1 && (ep_to_remove = sort!(ep_to_remove, rev=true))

    obj_new = delete_epoch(obj, ep=ep_to_remove)
    reset_components!(obj_new)
    push!(obj_new.history, "keep_epoch(OBJ, $ep)")    

    return obj_new

end

"""
    keep_epoch!(obj; ep)

Keep OBJ epoch(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ep::Union{Int64, Vector{Int64}, <:AbstractRange}`: epoch number(s) to keep
"""
function keep_epoch!(obj::NeuroAnalyzer.NEURO; ep::Union{Int64, Vector{Int64}, <:AbstractRange})

    obj_new = keep_epoch(obj, ep=ep)
    obj.header = obj_new.header
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.time_pts = obj_new.time_pts
    obj.markers = obj_new.markers

    return nothing

end
