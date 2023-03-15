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
    time_pts = collect(0:(1 / sr(obj)):(size(obj_new.data, 2) / sr(obj)))
    time_pts = round.(time_pts[1:(end - 1)], digits=3)
    obj_new.time_pts = time_pts
    obj_new.epoch_time = obj_new.epoch_time[1] .+ time_pts

    # remove markers of deleted epochs
    for marker_idx in nrow(obj_new.markers):-1:1
        obj_new.markers[marker_idx, :start] > size(obj_new.data, 2) && deleteat!(obj_new.markers, marker_idx)
    end

    reset_components!(obj_new)
    push!(obj_new.header.history, "erp(OBJ)")

    return obj_new
end

"""
    erp!(obj)

Average epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function erp!(obj::NeuroAnalyzer.NEURO)

    obj_tmp = erp(obj)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    obj.time_pts = obj_tmp.time_pts
    obj.markers = obj_tmp.markers
    obj.components = obj_tmp.components

    return nothing

end
