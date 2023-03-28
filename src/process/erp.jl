export erp
export erp!

"""
    erp(obj)

Average epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function erp(obj::NeuroAnalyzer.NEURO)

    obj_new = deepcopy(obj)
    obj_new.data = mean(obj_new.data, dims=3)[:, :, :]
    obj_new.time_pts, obj_new.epoch_time = _get_t(obj_new)

    # remove markers of deleted epochs
    for marker_idx in nrow(obj_new.markers):-1:1
        obj_new.markers[marker_idx, :start] > size(obj_new.data, 2) && deleteat!(obj_new.markers, marker_idx)
    end

    reset_components!(obj_new)
    push!(obj_new.history, "erp(OBJ)")

    return obj_new
end

"""
    erp!(obj)

Average epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function erp!(obj::NeuroAnalyzer.NEURO)

    obj_new = erp(obj)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history
    obj.header = obj_new.header
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time
    obj.markers = obj_new.markers

    return nothing

end
