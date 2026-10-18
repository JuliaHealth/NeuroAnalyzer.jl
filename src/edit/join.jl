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
    datatype(obj1) == datatype(obj2) ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same data type."))
    sr(obj1) == sr(obj2) ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same sampling rate."))
    nchannels(obj1) == nchannels(obj2) ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same number of channels."))
    labels(obj1) == labels(obj2) ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same channel labels (in the same order)."))
    nepochs(obj1) == nepochs(obj2) ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same number of epochs."))

    obj_tmp = deepcopy(obj1)

    # merge data: each obj2 epoch is appended to the respective obj1 epoch
    obj_tmp.data = hcat(obj1.data, obj2.data)

    # regenerate time points
    obj_tmp.time_pts, obj_tmp.epoch_time = _get_t(obj_tmp)

    # merge markers: positions within each joined epoch
    # (obj1 samples first, then obj2 samples)
    fs = sr(obj1)
    l1, l2 = epoch_len(obj1), epoch_len(obj2)
    m1 = _remap_markers(obj1.markers, fs, s -> ((s - 1) ÷ l1) * (l1 + l2) + (s - 1) % l1 + 1)
    m2 = _remap_markers(obj2.markers, fs, s -> ((s - 1) ÷ l2) * (l1 + l2) + l1 + (s - 1) % l2 + 1)
    obj_tmp.markers = sort(vcat(m1, m2), :start)

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
