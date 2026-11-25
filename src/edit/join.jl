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
    obj_new = deepcopy(obj1)
    join!(obj_new, obj2)

    return obj_new
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
    # validate
    datatype(obj1) == datatype(obj2) ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same data type."))
    sr(obj1) == sr(obj2) ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same sampling rate."))
    nchannels(obj1) == nchannels(obj2) ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same number of channels."))
    labels(obj1) == labels(obj2) ||
        throw(
            ArgumentError(
                "OBJ1 and OBJ2 must have the same channel labels (in the same order).",
            ),
        )
    nepochs(obj1) == nepochs(obj2) ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same number of epochs."))

    # merge markers: positions within each joined epoch (obj1 samples first, then obj2
    # samples); calculated before obj1 is modified (obj2 may be obj1)
    fs = sr(obj1)
    l1, l2 = epoch_len(obj1), epoch_len(obj2)
    m1 =
        _remap_markers(obj1.markers, fs, s -> ((s - 1) ÷ l1) * (l1 + l2) + (s - 1) % l1 + 1)
    m2 = _remap_markers(
        obj2.markers,
        fs,
        s -> ((s - 1) ÷ l2) * (l1 + l2) + l1 + (s - 1) % l2 + 1,
    )

    # merge data: each obj2 epoch is appended to the respective obj1 epoch
    # epoch metadata of obj1 is kept (each obj2 epoch continues an obj1 epoch)
    obj1.epochs = _epochs_table(obj1)
    obj1.data = hcat(obj1.data, obj2.data)
    # regenerate time points
    obj1.time_pts, obj1.epoch_time = _get_t(obj1)
    # extra marker columns of one object only get missing in the markers of the other
    obj1.markers = sort(vcat(m1, m2; cols = :union), :start)
    _history!(obj1, "join", NamedTuple())

    return nothing
end
