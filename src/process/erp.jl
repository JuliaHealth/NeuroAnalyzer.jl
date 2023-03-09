"""
    erp(obj)

Average epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function erp(obj::NeuroAnalyzer.NEURO)

    obj_new = deepcopy(obj)
    obj_new.data = mean(obj_new.data, dims=3)[:, :, :]
    duration_samples = size(obj_new.data, 2)
    duration_seconds = duration_samples / sr(obj)
    time_pts = collect(0:(1 / sr(obj)):duration_seconds)
    obj_new.time_pts = time_pts[1:(end - 1)]
    obj_new.header.recording[:duration_samples] = duration_samples
    obj_new.header.recording[:duration_seconds] = duration_seconds
    obj_new.header.recording[:epoch_n] = 1

    # remove markers of deleted epochs
    for marker_idx in nrow(obj_new.markers):-1:1
        obj_new.markers[marker_idx, :start] > duration_samples && deleteat!(obj_new.markers, marker_idx)
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
- `epoch::Union{Vector{Int64}, AbstractRange}=1:epoch_n(obj)`: epochs to average; default is all epochs
"""
function erp!(obj::NeuroAnalyzer.NEURO)

    obj_tmp = erp(obj)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    obj.time_pts = obj_tmp.time_pts
    obj.markers = obj_tmp.markers
    reset_components!(obj)

    return nothing
end
