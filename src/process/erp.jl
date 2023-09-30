export erp
export erp!

"""
    erp(obj; bl)

Average epochs. Non-signal channels are removed. `OBJ.header.recording[:data_type]` becomes `erp`. First epoch is the ERP.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `bl::Real=0`: baseline is the first `bl` seconds; if `bl` is greater than 0, DC value is calculated as mean of the first `bl` seconds and subtracted from the signal

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function erp(obj::NeuroAnalyzer.NEURO; bl::Real=0)

    nchannels(obj) > length(signal_channels(obj)) && _warn("Non-signal channels will be removed.")

    obj_new = keep_channel(obj, ch=signal_channels(obj))
    obj_new.data = cat(mean(obj_new.data, dims=3)[:, :, :], obj_new.data, dims=3)
    obj_new.header.recording[:data_type] = "erp"
    obj_new.time_pts, obj_new.epoch_time = _get_t(obj_new)

    # remove DC
    if bl != 0
        obj_new.data[:, :, 1] = remove_dc(obj_new.data[:, :, 1], round(Int64, bl * sr(obj)))
    end

    # remove markers of deleted epochs
    for marker_idx in nrow(obj_new.markers):-1:1
        obj_new.markers[marker_idx, :start] > size(obj_new.data, 2) && deleteat!(obj_new.markers, marker_idx)
    end
    obj_new.markers[!, :start] .+= (obj_new.epoch_time[1] * sr(obj_new))

    reset_components!(obj_new)
    push!(obj_new.history, "erp(OBJ, bl=$bl)")

    return obj_new
    
end

"""
    erp!(obj; bl)

Average epochs. Non-signal channels are removed. `OBJ.header.recording[:data_type]` becomes `erp`. First epoch is the ERP.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `bl::Real=0`: baseline is the first `bl` seconds; if `bl` is greater than 0, DC value is calculated as mean of the first `n` samples and subtracted from the signal
"""
function erp!(obj::NeuroAnalyzer.NEURO; bl::Real=0)

    obj_new = erp(obj, bl=bl)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history
    obj.header = obj_new.header
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time
    obj.markers = obj_new.markers

    return nothing

end
