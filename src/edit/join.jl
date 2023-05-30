export join
export join!

"""
    join(obj1, obj2)

Join two NeuroAnalyzer objects. Both objects must have the same data type, number of channels, epochs and sampling rate.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function join(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO)

    obj1.header.recording[:data_type] == obj1.header.recording[:data_type] || throw(ArgumentError("Both objects must have the same data type."))
    sr(obj1) == sr(obj2) || throw(ArgumentError("Both objects must have the same sampling rate."))
    channel_n(obj1) == channel_n(obj2) || throw(ArgumentError("Both objects must have the same number of channels."))
    epoch_n(obj1) == epoch_n(obj2) || throw(ArgumentError("Both objects must have the same number of epochs."))

    obj_new = deepcopy(obj1)

    # merge data
    obj_new.data = hcat(obj1.data, obj2.data)
    
    # regenerate time points
    obj_new.time_pts, obj_new.epoch_time = NeuroAnalyzer._get_t(obj_new)

    # merge markers
    if nrow(obj2.markers) > 0
        obj_new.markers = vcat(obj1.markers, obj2.markers)
    end
    if nrow(obj1.markers) > 0
        obj_new.markers[(nrow(obj1.markers) + 1):end, :start] .+= (signal_len(obj1) / sr(obj1))
    end
    
    reset_components!(obj_new)

    push!(obj_new.history, "join(OBJ1, OBJ2)")

    return obj_new

end

"""
    join!(obj1, obj2)

Join two NeuroAnalyzer objects. Both objects must have the same data type, number of channels, epochs and sampling rate.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
"""
function join!(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO)

    obj_new = NeuroAnalyzer.join(obj1, obj2)
    obj1.data = obj_new.data
    obj1.history = obj_new.history
    obj1.components = obj_new.components
    obj1.time_pts = obj_new.time_pts
    obj1.epoch_time = obj_new.epoch_time
    obj1.markers = obj_new.markers

    return obj_new

end