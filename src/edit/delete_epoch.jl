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
    obj_tmp = deepcopy(obj)

    # remove epoch
    obj_tmp = deepcopy(obj)
    obj_tmp.data = obj_tmp.data[:, :, setdiff(1:nepochs(obj), ep_sorted)]

    epoch_ranges = [_epoch2s(obj, e) for e in sort(collect(ep_sorted))]
    for (t1, t2) in reverse(epoch_ranges) # process latest epochs first to preserve offsets
        obj_tmp.markers = _delete_markers(obj_tmp.markers, (t1, t2))
        obj_tmp.markers = _shift_markers(obj_tmp.markers, (t1, t2))
    end

    obj_tmp.time_pts, obj_tmp.epoch_time = _get_t(obj_tmp)
    push!(obj_tmp.history, "delete_epoch(obj; ep=$ep)")

    return obj_tmp
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
    obj_tmp = delete_epoch(obj; ep = ep)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj.time_pts = obj_tmp.time_pts
    obj.epoch_time = obj_tmp.epoch_time
    obj.markers = obj_tmp.markers
    obj_tmp = nothing

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

    obj_tmp = delete_epoch(obj; ep = ep_to_remove)
    push!(obj_tmp.history, "keep_epoch(obj; ep=$ep)")

    return obj_tmp
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
    obj_tmp = keep_epoch(obj; ep = ep)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj.time_pts = obj_tmp.time_pts
    obj.epoch_time = obj_tmp.epoch_time
    obj.markers = obj_tmp.markers
    obj_tmp = nothing

    return nothing
end
