export join
export join!

"""
    join(obj1, obj2)

Join two NeuroAnalyzer objects. Each `obj2` epoch are horizontally concatenated (along time) with respective `obj1` epoch. Both objects must have the same data type, number of channels, epochs and sampling rate, but may differ in epoch lengths.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function join(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO)::NeuroAnalyzer.NEURO

    @assert datatype(obj1) == obj1.header.recording[:data_type] "OBJ1 and OBJ2 must have the same data type."
    @assert sr(obj1) == sr(obj2) "OBJ1 and OBJ2 must have the same sampling rate."
    @assert nchannels(obj1) == nchannels(obj2) "OBJ1 and OBJ2 must have the same number of channels."
    @assert nepochs(obj1) == nepochs(obj2) "OBJ1 and OBJ2 must have the same number of epochs."

    obj_new = deepcopy(obj1)

    # merge data
    obj_new.data = hcat(obj1.data, obj2.data)

    # regenerate time points
    obj_new.time_pts, obj_new.epoch_time = _get_t(obj_new)

    # merge markers
    DataFrame.nrow(obj2.markers) > 0 && (obj_new.markers = vcat(obj1.markers, obj2.markers))
    DataFrame.nrow(obj1.markers) > 0 && (obj_new.markers[(DataFrame.nrow(obj1.markers) + 1):end, :start] .+= (signal_len(obj1) / sr(obj1)))

    reset_components!(obj_new)

    push!(obj_new.history, "join(OBJ1, OBJ2)")

    return obj_new

end

"""
    join!(obj1, obj2)

Join two NeuroAnalyzer objects into the first object. Each `obj2` epoch are horizontally concatenated (along time) with respective `obj1` epoch. Both objects must have the same data type, number of channels, epochs and sampling rate, but may differ in epoch lengths.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`

# Returns

Nothing
"""
function join!(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO)::Nothing

    obj_new = NeuroAnalyzer.join(obj1, obj2)
    obj1.data = obj_new.data
    obj1.history = obj_new.history
    obj1.components = obj_new.components
    obj1.time_pts = obj_new.time_pts
    obj1.epoch_time = obj_new.epoch_time
    obj1.markers = obj_new.markers

    return nothing

end