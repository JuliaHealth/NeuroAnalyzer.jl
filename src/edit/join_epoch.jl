export join_epoch
export join_epoch!
export append_epoch
export append_epoch!

"""
    join_epoch(obj1, obj2)

Join two NeuroAnalyzer objects in time: each `obj2` epoch is appended (along time) to the respective `obj1` epoch, e.g. two parts of a continuous recording. Both objects must have the same data type, channel labels, number of epochs and sampling rate, but may differ in epoch lengths.

A marker of type `"boundary"` (value `"join"`) is placed at each joint (the first sample of the `obj2` signal), so that filtering does not run across it (see `filter()`).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: input NEURO object
- `obj2::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object

# See also

`append_epoch()` to append the epochs of other objects as new epochs.
"""
function join_epoch(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj1)
    join_epoch!(obj_new, obj2)

    return obj_new
end

"""
    join_epoch!(obj1, obj2)

Join two NeuroAnalyzer objects in time into the first object: each `obj2` epoch is appended (along time) to the respective `obj1` epoch, e.g. two parts of a continuous recording. Both objects must have the same data type, channel labels, number of epochs and sampling rate, but may differ in epoch lengths.

A marker of type `"boundary"` (value `"join"`) is placed at each joint (the first sample of the `obj2` signal), so that filtering does not run across it (see `filter()`).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `obj2::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Nothing`
"""
function join_epoch!(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO)::Nothing
    # validate
    _check_compatible(obj1, obj2)
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
    obj1.markers = vcat(m1, m2; cols = :union)
    # boundary markers at the joints (first sample of the obj2 signal in each epoch)
    for ep = 1:nepochs(obj1)
        _markers_push!(
            obj1.markers,
            (
                type = "boundary",
                start = s2t((ep - 1) * (l1 + l2) + l1 + 1, fs),
                length = 0.0,
                value = "join",
                channel = String[],
            ),
        )
    end
    sort!(obj1.markers, :start)
    _history!(obj1, "join_epoch", NamedTuple(); args = (obj2,))

    return nothing
end

"""
    append_epoch(obj, objs...)

Append the epochs of other objects as new epochs, e.g. epochs of several recordings or sessions of one subject. All objects must have the same data type, channel labels, sampling rate and epoch length (continuous objects of the same length become one epoch each).

The epoch metadata (condition, source, rejection, onset) and the markers of all objects are kept; the source of epochs without one is the file name of their object. A channel marked as bad in any object is marked as bad.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `objs::NeuroAnalyzer.NEURO...`: objects whose epochs are appended

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object

# Notes

The epochs are separate segments of the signal: filters do not run across epochs, so no boundary markers are needed. The header (except bad channels) is that of `obj`; a warning is shown if the low-pass, high-pass or notch frequencies of the objects differ.

# See also

`join_epoch()` to join objects in time.
"""
function append_epoch(obj::NeuroAnalyzer.NEURO, objs::NeuroAnalyzer.NEURO...)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    append_epoch!(obj_new, objs...)

    return obj_new
end

"""
    append_epoch!(obj, objs...)

Append the epochs of other objects as new epochs in-place, e.g. epochs of several recordings or sessions of one subject. All objects must have the same data type, channel labels, sampling rate and epoch length (continuous objects of the same length become one epoch each).

The epoch metadata (condition, source, rejection, onset) and the markers of all objects are kept; the source of epochs without one is the file name of their object. A channel marked as bad in any object is marked as bad.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `objs::NeuroAnalyzer.NEURO...`: objects whose epochs are appended

# Returns

- `Nothing`

# Notes

The epochs are separate segments of the signal: filters do not run across epochs, so no boundary markers are needed. The header (except bad channels) is that of `obj`; a warning is shown if the low-pass, high-pass or notch frequencies of the objects differ.
"""
function append_epoch!(obj::NeuroAnalyzer.NEURO, objs::NeuroAnalyzer.NEURO...)::Nothing
    # validate
    isempty(objs) && throw(ArgumentError("At least one object must be appended."))
    for o in objs
        _check_compatible(obj, o)
        epoch_len(o) == epoch_len(obj) || throw(
            ArgumentError(
                "Epoch lengths differ ($(epoch_len(obj)) and $(epoch_len(o)) samples); use join_epoch() to join objects in time.",
            ),
        )
    end

    fs = sr(obj)
    l = epoch_len(obj)
    r = obj.header.recording
    _rec_filters!(r)

    # epoch metadata and markers (shifted by the epochs before them); computed before obj
    # is modified (an object may be appended to itself)
    ep_tab = [_epochs_meta(o) for o in (obj, objs...)]
    n_ep = [nepochs(o) for o in (obj, objs...)]
    mrk = DataFrame[]
    for (k, o) in enumerate((obj, objs...))
        shift = sum(n_ep[1:(k - 1)]; init = 0) * l
        push!(mrk, _remap_markers(o.markers, fs, s -> s + shift))
    end
    filt_diff = false
    for o in objs
        ro = o.header.recording
        _rec_filters!(ro)
        filt_diff |= ro[:lowpass] != r[:lowpass] || ro[:highpass] != r[:highpass] ||
                     ro[:notch] != r[:notch]
        r[:bad_channel] = r[:bad_channel] .| ro[:bad_channel]
    end
    any(o -> !(o.epoch_time ≈ obj.epoch_time), objs) &&
        _warn("The epoch time axes of the objects differ; that of the first object is kept.")
    filt_diff && _warn("The low-pass, high-pass or notch frequencies of the objects differ; those of the first object are kept.")
    data = cat(obj.data, (o.data for o in objs)...; dims = 3)

    obj.data = data
    obj.epochs = vcat(ep_tab...)
    obj.markers = sort(vcat(mrk...; cols = :union), :start)
    # the epoch time axis (e.g. the pre-stimulus offset) of obj is kept
    et = copy(obj.epoch_time)
    obj.time_pts, _ = _get_t(obj)
    obj.epoch_time = et

    _history!(obj, "append_epoch", NamedTuple(); args = objs)

    return nothing
end

# data type, sampling rate and channel labels of two objects match
function _check_compatible(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO)::Nothing
    datatype(obj1) == datatype(obj2) ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same data type."))
    sr(obj1) == sr(obj2) ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same sampling rate."))
    nchannels(obj1) == nchannels(obj2) ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same number of channels."))
    labels(obj1) == labels(obj2) || throw(
        ArgumentError("OBJ1 and OBJ2 must have the same channel labels (in the same order)."),
    )
    return nothing
end

# epoch metadata of an object; epochs without a source get the file name of the object
function _epochs_meta(obj::NeuroAnalyzer.NEURO)::DataFrame
    t = copy(_epochs_table(obj))
    src = _epochs_source(obj.header)
    t[!, :source] = String[isempty(s) ? src : s for s in t[!, :source]]
    return t
end
