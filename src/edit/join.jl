export join
export join!

"""
    join(obj1, obj2)

Join two NeuroAnalyzer objects. Each `obj2` epoch are horizontally concatenated (along time) with respective `obj1` epoch. Both objects must have the same data type, number of channels, epochs and sampling rate, but may differ in epoch lengths.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: input NEURO object
- `obj2::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function join(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO)::NeuroAnalyzer.NEURO
    # validate
    datatype(obj1) == obj1.header.recording[:data_type] ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same data type."))
    sr(obj1) == sr(obj2) ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same sampling rate."))
    nchannels(obj1) == nchannels(obj2) ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same number of channels."))
    nepochs(obj1) == nepochs(obj2) ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same number of epochs."))

    obj_tmp = deepcopy(obj1)

    # merge data
    obj_tmp.data = hcat(obj1.data, obj2.data)

    # regenerate time points
    obj_tmp.time_pts, obj_tmp.epoch_time = _get_t(obj_tmp)

    # merge markers
    DataFrames.nrow(obj2.markers) > 0 &&
        (obj_tmp.markers = vcat(obj1.markers, obj2.markers))
    DataFrames.nrow(obj1.markers) > 0 && (
        obj_tmp.markers[(DataFrames.nrow(obj1.markers) + 1):end, :start] .+= (
            signal_len(obj1) / sr(obj1)
        )
    )

    push!(obj_tmp.history, "join(obj1, obj2)")

    return obj_tmp
end

"""
    join!(obj1, obj2)

Join two NeuroAnalyzer objects into the first object. Each `obj2` epoch are horizontally concatenated (along time) with respective `obj1` epoch. Both objects must have the same data type, number of channels, epochs and sampling rate, but may differ in epoch lengths.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: input NEURO object
- `obj2::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Nothing`
"""
function join!(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO)::Nothing
    obj_tmp = NeuroAnalyzer.join(obj1, obj2)
    obj1.data = obj_tmp.data
    obj1.history = obj_tmp.history
    obj1.time_pts = obj_tmp.time_pts
    obj1.epoch_time = obj_tmp.epoch_time
    obj1.markers = obj_tmp.markers
    obj_tmp = nothing

    return nothing
end
