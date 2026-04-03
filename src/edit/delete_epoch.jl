export delete_epoch
export delete_epoch!
export keep_epoch
export keep_epoch!

"""
    delete_epoch(obj; <keyword arguments>)

Return a copy of `obj` with the specified epoch(s) removed.

Markers within deleted epochs are dropped; markers after deleted epochs are shifted to remain aligned with the new time axis.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ep::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}`: epoch numbers to remove

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function delete_epoch(
    obj::NeuroAnalyzer.NEURO;
    ep::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}},
)::NeuroAnalyzer.NEURO
    # validate
    nepochs(obj) > 1 || throw(ArgumentError("You cannot delete the last epoch."))
    ep = _n2v(ep)
    ep_sorted = sort(collect(ep); rev = true)
    length(ep_sorted) < nepochs(obj) || throw(
        ArgumentError(
            "Number of epochs to delete ($(length(ep_sorted))) must be less than " *
            "the total number of epochs ($(nepochs(obj))).",
        ),
    )
    _check_epochs(obj, ep_sorted)

    # create new dataset
    obj_new = deepcopy(obj)

    # remove epoch
    obj_new = deepcopy(obj)
    obj_new.data = obj_new.data[:, :, setdiff(1:nepochs(obj), ep_sorted)]

    epoch_ranges = [_epoch2s(obj, e) for e in sort(collect(ep_sorted))]
    for (t1, t2) in reverse(epoch_ranges) # process latest epochs first to preserve offsets
        obj_new.markers = _delete_markers(obj_new.markers, (t1, t2))
        obj_new.markers = _shift_markers(obj_new.markers, (t1, t2))
    end

    obj_new.time_pts, obj_new.epoch_time = _get_t(obj_new)
    push!(obj_new.history, "delete_epoch(obj; ep=$ep)")

    return obj_new
end

"""
    delete_epoch!(obj; <keyword arguments>)

Delete epoch(s) in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ep::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}`: epoch numbers to remove

# Returns

- `Nothing`
"""
function delete_epoch!(
    obj::NeuroAnalyzer.NEURO; ep::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}},
)::Nothing
    obj_new = delete_epoch(obj; ep = ep)
    obj.header = obj_new.header
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time
    obj.markers = obj_new.markers

    return nothing
end

"""
    keep_epoch(obj; <keyword arguments>)

Return a copy of `obj` retaining only the specified epoch(s).

Implemented by computing the complement set and delegating to `delete_epoch`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ep::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}`: epoch numbers to keep

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function keep_epoch(
    obj::NeuroAnalyzer.NEURO;
    ep::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}},
)::NeuroAnalyzer.NEURO
    # validate
    nepochs(obj) > 1 || throw(ArgumentError("OBJ contains only one epoch."))
    _check_epochs(obj, ep)

    ep_to_remove = setdiff(1:nepochs(obj), ep)
    isempty(ep_to_remove) && return deepcopy(obj) # nothing to remove

    obj_new = delete_epoch(obj; ep = ep_to_remove)
    push!(obj_new.history, "keep_epoch(obj; ep=$ep)")

    return obj_new
end

"""
    keep_epoch!(obj; <keyword arguments>)

Keep only the specified epoch(s) in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ep::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}`: epoch numbers to keep

# Returns

- `Nothing`
"""
function keep_epoch!(
    obj::NeuroAnalyzer.NEURO; ep::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}},
)::Nothing
    obj_new = keep_epoch(obj; ep = ep)
    obj.header = obj_new.header
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time
    obj.markers = obj_new.markers

    return nothing
end
