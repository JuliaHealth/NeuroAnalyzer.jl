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
    obj_new = deepcopy(obj)
    delete_epoch!(obj_new; ep = ep)

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
    obj::NeuroAnalyzer.NEURO;
    ep::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}},
)::Nothing
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

    # remove epoch
    obj.data = obj.data[:, :, setdiff(1:nepochs(obj), ep_sorted)]

    # marker times are in seconds; epoch e spans [(e - 1) * ep_dur, e * ep_dur)
    ep_dur = epoch_duration(obj)
    # process latest epochs first to preserve offsets
    for e in ep_sorted
        t1 = (e - 1) * ep_dur
        t2 = e * ep_dur
        deleteat!(obj.markers, t1 .<= obj.markers.start .< t2)
        obj.markers[obj.markers.start .>= t2, :start] .-= ep_dur
    end

    obj.time_pts, obj.epoch_time = _get_t(obj)
    push!(obj.history, "delete_epoch(obj; ep=$ep)")

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
    obj_new = deepcopy(obj)
    keep_epoch!(obj_new; ep = ep)

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
    # validate
    nepochs(obj) > 1 || throw(ArgumentError("OBJ contains only one epoch."))
    _check_epochs(obj, ep)

    ep_to_remove = setdiff(1:nepochs(obj), ep)
    isempty(ep_to_remove) && return nothing # nothing to remove

    delete_epoch!(obj; ep = ep_to_remove)
    push!(obj.history, "keep_epoch(obj; ep=$ep)")

    return nothing
end
