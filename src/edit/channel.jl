export channel_type
export set_channel_type
export set_channel_type!
export rename_channel
export rename_channel!
export edit_channel
export edit_channel!
export replace_channel
export replace_channel!
export add_label
export add_label!
export append_channel
export append_channel!
export insert_channel
export insert_channel!

# ---------------------------------------------------------------------------
# internal helpers
# ---------------------------------------------------------------------------

# resolve a channel name that must match exactly one channel
function _single_channel(obj::NeuroAnalyzer.NEURO, ch::String)::Int64
    idx = get_channel(obj; ch = ch)
    idx isa Integer && return Int64(idx)
    isempty(idx) && throw(ArgumentError("No channels selected."))
    length(idx) == 1 || throw(ArgumentError("ch must resolve to exactly one channel."))
    return Int64(idx[1])
end

# true if the object holds (not yet applied) SSP projection vectors
_has_ssp(obj::NeuroAnalyzer.NEURO)::Bool =
    haskey(obj.header.recording, :ssp_data) &&
    !isempty(obj.header.recording[:ssp_data])

# rename a label in OBJ.locs (all matching rows)
function _rename_locs!(locs::DataFrame, old::String, new::String)::Nothing
    nrow(locs) == 0 && return nothing
    locs[locs[!, :label] .== old, :label] .= new
    return nothing
end

# keep MEG sensor index lists consistent with a channel's type
function _update_meg_sensor_lists!(rec::Dict, ch::Int64, type::String)::Nothing
    for field in (:magnetometers, :gradiometers)
        haskey(rec, field) || continue
        rec[field] = Base.filter(!isequal(ch), rec[field])
    end
    if type == "mag"
        rec[:magnetometers] = sort(push!(get(rec, :magnetometers, Int64[]), ch))
    elseif type == "grad"
        rec[:gradiometers] = sort(push!(get(rec, :gradiometers, Int64[]), ch))
    end
    return nothing
end

"""
    channel_type(obj; <keyword arguments>)

Return the type string of a single channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel

# Returns

- `String`: channel type (e.g. `"eeg"`, `"eog"`)
"""
function channel_type(obj::NeuroAnalyzer.NEURO; ch::String)::String
    return obj.header.recording[:channel_type][_single_channel(obj, ch)]
end

"""
    set_channel_type(obj; <keyword arguments>)

Return a copy of `obj` with the type of one channel changed. Channel order is recomputed; for MEG data the magnetometer/gradiometer index lists are updated.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `type::String`: new channel type (must be a recognized type)

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function set_channel_type(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    type::String,
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    set_channel_type!(obj_new; ch = ch, type = type)

    return obj_new
end

"""
    set_channel_type!(obj; <keyword arguments>)

Set the type of one channel in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `type::String`: new channel type (must be a recognized type)

# Returns

- `Nothing`
"""
function set_channel_type!(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    type::String,
)::Nothing
    # validate
    type = lowercase(type)
    _check_var(type, string.(channel_types), "type")
    ch_idx = _single_channel(obj, ch)

    rec = obj.header.recording

    rec[:channel_type][ch_idx] = type
    rec[:channel_order] = _sort_channels(rec[:channel_type])
    (type in ("mag", "grad") || haskey(rec, :magnetometers)) &&
        _update_meg_sensor_lists!(rec, ch_idx, type)
    # per-channel fields needed by the new type (e.g. :gain of an EEG channel)
    for f in _rec_fields_for_type(type)
        _rec_ensure_ch_field!(rec, f, nchannels(obj))
    end

    _history!(obj, "set_channel_type", (; ch = labels(obj)[ch_idx], type))

    return nothing
end

"""
    rename_channel(obj; <keyword arguments>)

Return a copy of `obj` with one channel renamed. Also updates the matching entries in `obj.locs`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `name::String`: new name (must not be empty or already exist)

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function rename_channel(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    name::String,
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    rename_channel!(obj_new; ch = ch, name = name)

    return obj_new
end

"""
    rename_channel!(obj; <keyword arguments>)

Rename one channel in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `name::String`: new name (must not be empty or already exist)

# Returns

- `Nothing`
"""
function rename_channel!(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    name::String,
)::Nothing
    isempty(name) && throw(ArgumentError("name cannot be empty."))
    name in labels(obj) && throw(ArgumentError("Channel \"$name\" already exists."))
    ch_idx = _single_channel(obj, ch)
    old = labels(obj)[ch_idx]

    obj.header.recording[:label][ch_idx] = name
    _rename_locs!(obj.locs, old, name)
    _markers_relabel!(obj.markers, [old], [name])

    _history!(obj, "rename_channel", (; ch = old, name))

    return nothing
end

"""
    edit_channel(obj; <keyword arguments>)

Return a copy of `obj` with one field of one channel's header entry changed. Editing `:label` behaves as `rename_channel()`, editing `:channel_type` as `set_channel_type()`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `field::Symbol`: field to edit (`:channel_type` or `:label`)
- `value::String`: new value

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function edit_channel(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    field::Symbol,
    value::String,
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    edit_channel!(obj_new; ch = ch, field = field, value = value)

    return obj_new
end

"""
    edit_channel!(obj; <keyword arguments>)

Edit one field of one channel's header entry in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `field::Symbol`: field to edit (`:channel_type` or `:label`)
- `value::String`: new value

# Returns

- `Nothing`
"""
function edit_channel!(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    field::Symbol,
    value::String,
)::Nothing
    # validate
    isempty(value) && throw(ArgumentError("value cannot be empty."))
    _check_var(field, [:channel_type, :label], "field")
    ch_idx = _single_channel(obj, ch)

    # label before renaming (history)
    lbl = labels(obj)[ch_idx]
    n0 = length(obj.history)
    if field == :label
        rename_channel!(obj; ch = ch, name = value)
    else
        set_channel_type!(obj; ch = ch, type = value)
    end

    # record the user's call instead of the delegated one
    _history!(obj, "edit_channel", (; ch = lbl, field, value); from = n0)

    return nothing
end

"""
    replace_channel(obj; <keyword arguments>)

Return a copy of `obj` with one channel's data replaced.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `s::AbstractArray`: replacement data, shape (1, epoch_len, nepochs)

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function replace_channel(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    s::AbstractArray,
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    replace_channel!(obj_new; ch = ch, s = s)

    return obj_new
end

"""
    replace_channel!(obj; <keyword arguments>)

Replace one channel's data in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `s::AbstractArray`: replacement data, shape (1, epoch_len, nepochs)

# Returns

- `Nothing`
"""
function replace_channel!(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    s::AbstractArray,
)::Nothing
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)
    size(s) == (1, epoch_len(obj), nepochs(obj)) || throw(
        ArgumentError(
            "data size ($(size(s, 1)) × $(size(s, 2)) × $(size(s, 3))) must be 1 × $(epoch_len(obj)) × $(nepochs(obj)).",
        ),
    )
    ch_idx = _single_channel(obj, ch)

    _has_ssp(obj) && obj.header.recording[:ssp_channels][ch_idx] &&
        _warn("Channel is part of SSP projections; apply them before modifying data.")

    obj.data[ch_idx, :, :] = s

    _history!(obj, "replace_channel", (; ch = labels(obj)[ch_idx], s))

    return nothing
end

"""
    add_label(obj; <keyword arguments>)

Return a copy of `obj` with channel labels replaced. Labels in `obj.locs` are renamed accordingly.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `clabels::Vector{String}`: new labels; must be unique, non-empty and have length equal to `nchannels(obj)`

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function add_label(obj::NeuroAnalyzer.NEURO; clabels::Vector{String})::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    add_label!(obj_new; clabels = clabels)

    return obj_new
end

"""
    add_label!(obj; <keyword arguments>)

Replace channel labels in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `clabels::Vector{String}`: new labels; must be unique, non-empty and have length equal to `nchannels(obj)`

# Returns

- `Nothing`
"""
function add_label!(obj::NeuroAnalyzer.NEURO; clabels::Vector{String})::Nothing
    # validate
    length(clabels) == nchannels(obj) ||
        throw(ArgumentError("clabels length must be $(nchannels(obj))."))
    any(isempty, clabels) && throw(ArgumentError("clabels cannot contain empty labels."))
    allunique(clabels) || throw(ArgumentError("clabels must be unique."))

    old = labels(obj)
    obj.header.recording[:label] = copy(clabels)

    # rename locs in one pass (avoids clashes when labels are permuted)
    if nrow(obj.locs) > 0
        mapping = Dict(zip(old, clabels))
        obj.locs[!, :label] = [get(mapping, l, l) for l in obj.locs[!, :label]]
    end
    # markers refer to channels by label (renamed in one pass, like the locations)
    _markers_relabel!(obj.markers, old, clabels)

    _history!(obj, "add_label", (; clabels))

    return nothing
end

"""
    append_channel(obj; <keyword arguments>)

Return a copy of `obj` with new channel(s) appended.

If `obj.data` is empty the object is initialized with the provided data (time axes are created from the sampling rate); otherwise the new channels are concatenated along the first (channel) dimension.

All per-channel header fields are extended with defaults: no filter for the low-pass, high-pass and notch frequencies (Nyquist frequency, 0 Hz, none), `""` for transducers and coil type, `1.0` for gain, a zero row for MEG sensor geometry (`:coil_loc`) and `0` for its coil code (`:coil_id`), `false` for SSP channel mask, `0` for MEP stimulation intensity. New MEG channels of type `"mag"`/`"grad"` are added to the magnetometer/gradiometer lists. NIRS channels cannot be appended (they are described by optode pairs); other channels (e.g. EEG, auxiliary) can be appended to NIRS objects.

Channels of another object are appended with `append_channel(obj1, obj2)`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `data::Array{<:Number, 3}`: channel data, shape (ch_n, epoch_len, epoch_n)
- `label::Union{String, Vector{String}}`: channel label(s)
- `type::Union{String, Vector{String}}`: channel type(s)
- `unit::Union{String, Vector{String}}`: channel unit(s)

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function append_channel(
    obj::NeuroAnalyzer.NEURO;
    data::Array{<:Number, 3},
    label::Union{String, Vector{String}},
    type::Union{String, Vector{String}},
    unit::Union{String, Vector{String}},
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    append_channel!(obj_new; data = data, label = label, type = type, unit = unit)

    return obj_new
end

"""
    append_channel!(obj; <keyword arguments>)

Add new channel(s) in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `data::Array{<:Number, 3}`: channel data, shape (ch_n, epoch_len, epoch_n)
- `label::Union{String, Vector{String}}`: channel label(s)
- `type::Union{String, Vector{String}}`: channel type(s)
- `unit::Union{String, Vector{String}}`: channel unit(s)

# Returns

- `Nothing`
"""
function append_channel!(
    obj::NeuroAnalyzer.NEURO;
    data::Array{<:Number, 3},
    label::Union{String, Vector{String}},
    type::Union{String, Vector{String}},
    unit::Union{String, Vector{String}},
)::Nothing
    _append_channel!(obj; data = data, label = label, type = type, unit = unit)
    _history!(
        obj, "append_channel",
        (; label = label isa String ? [label] : label,
         type = lowercase.(type isa String ? [type] : type),
         unit = unit isa String ? [unit] : unit, data),
    )

    return nothing
end

# core of append_channel!(): data and header fields of the new channels (no history);
# `nirs_opt = true`: NIRS optical channels may be added (their optode description is set by
# the caller)
function _append_channel!(
    obj::NeuroAnalyzer.NEURO;
    data::Array{<:Number, 3},
    label::Union{String, Vector{String}},
    type::Union{String, Vector{String}},
    unit::Union{String, Vector{String}},
    nirs_opt::Bool = false,
)::Nothing
    ch_n = size(data, 1)
    ch_n > 0 || throw(ArgumentError("data must contain at least one channel."))

    # normalize to vectors for uniform handling below
    label_v = label isa String ? [label] : copy(label)
    type_v = lowercase.(type isa String ? [type] : type)
    unit_v = unit isa String ? [unit] : copy(unit)

    rec = obj.header.recording
    dt = rec[:data_type]
    # objects with NIRS channels: NIRS channels are described by optode pairs, other channels
    # (e.g. EEG, auxiliary) can be appended
    !nirs_opt && (dt == "nirs" || _nirs_n_opt(rec) > 0) &&
        any(t -> startswith(t, "nirs") && t != "nirs_aux", type_v) &&
        throw(ArgumentError("NIRS channels must be added using optode functions."))

    # validate
    has_data = length(obj.data) > 0
    if has_data
        epoch_len(obj) == size(data, 2) || throw(
            ArgumentError(
                "Epoch length mismatch: new data has $(size(data, 2)), object has $(epoch_len(obj)).",
            ),
        )
        nepochs(obj) == size(data, 3) || throw(
            ArgumentError(
                "Epoch count mismatch: new data has $(size(data, 3)), object has $(nepochs(obj)).",
            ),
        )
    end
    for (name, v) in (("labels", label_v), ("types", type_v), ("units", unit_v))
        length(v) == ch_n || throw(
            ArgumentError(
                "Number of $name ($(length(v))) must equal number of new channels ($ch_n).",
            ),
        )
    end
    any(isempty, label_v) && throw(ArgumentError("Labels cannot be empty."))
    allunique(label_v) || throw(ArgumentError("New labels must be unique."))
    old_labels = has_data ? labels(obj) : String[]
    dup = intersect(label_v, old_labels)
    isempty(dup) || throw(ArgumentError("Channel(s) $dup already exist."))
    for t in type_v
        t in string.(channel_types) || throw(ArgumentError("Unknown channel type \"$t\"."))
    end

    r = obj.header.recording
    n_old = has_data ? nchannels(obj) : 0

    obj.data = has_data ? cat(obj.data, data; dims = 1) : Array{Float64}(data)
    r[:label] = vcat(old_labels, label_v)
    r[:channel_type] = vcat(has_data ? r[:channel_type] : String[], type_v)
    r[:unit] = vcat(has_data ? r[:unit] : String[], unit_v)
    _vext(field, default) = vcat(
        (has_data && haskey(r, field) && length(r[field]) == n_old) ? r[field] :
        typeof(default)[],
        fill(default, ch_n),
    )
    r[:bad_channel] = _vext(:bad_channel, false)
    # low-pass and high-pass frequency: no filter
    r[:lowpass] = _vext(:lowpass, Float64(r[:sampling_rate]) / 2)
    r[:highpass] = _vext(:highpass, 0.0)
    # notch frequencies: none (one vector per channel, not shared)
    r[:notch] = vcat(
        (has_data && haskey(r, :notch)) ? r[:notch] : Vector{Float64}[],
        [Float64[] for _ = 1:ch_n],
    )

    # per-channel fields present in the header (of any modality): default values for the
    # new channels
    for (f, d) in _REC_CH_FIELDS
        _rec_has_ch_field(r, f, n_old) && (r[f] = vcat(r[f], fill(d, ch_n)))
    end
    _rec_has_ch_field(r, :coil_loc, n_old) && (r[:coil_loc] = vcat(r[:coil_loc], zeros(ch_n, 12)))
    # fields needed by the new channel types (e.g. :gain of EEG channels appended to a NIRS
    # object) and, for an empty object, by its data type
    need = unique(vcat(_rec_fields_for_type.(type_v)...))
    if !has_data
        dt in _EEG_LIKE_DATA && append!(need, [:transducers, :gain])
        dt in _MEG_LIKE_DATA && append!(need, [:coil_type, :ssp_channels])
        dt == "mep" && append!(need, [:stimulation_intensity, :coil_type])
    end
    for f in unique(need)
        _rec_ensure_ch_field!(r, f, n_old + ch_n)
    end
    for t in unique(type_v), (k, v) in _rec_scalar_defaults(t)
        haskey(r, k) || (r[k] = deepcopy(v))
    end
    # MEG sensor index lists
    if any(in(("mag", "grad")), type_v) || haskey(r, :magnetometers)
        for field in (:magnetometers, :gradiometers)
            has_data && haskey(r, field) || (r[field] = Int64[])
        end
        for (i, t) in enumerate(type_v)
            t == "mag" && push!(r[:magnetometers], n_old + i)
            t == "grad" && push!(r[:gradiometers], n_old + i)
        end
    end

    r[:channel_order] = _sort_channels(r[:channel_type])

    # time axes for a previously empty object
    if !has_data
        fs = r[:sampling_rate]
        fs > 0 ||
            throw(ArgumentError("Sampling rate must be set to initialize an empty object."))
        n_s = size(data, 2)
        n_e = size(data, 3)
        obj.time_pts = _time_axis(n_s * n_e, fs)
        obj.epoch_time = _time_axis(n_s, fs)
    end

    return nothing
end

"""
    append_channel(obj1, obj2; <keyword arguments>)

Return a copy of `obj1` with channels of `obj2` appended, e.g. EEG with auxiliary, EOG or ECG channels recorded separately, or EEG recorded simultaneously with NIRS.

The objects must have the same sampling rate, epoch length and number of epochs (the signals must be synchronized). The objects may contain channels of different types (e.g. EEG and NIRS): the header fields of all channel types are merged:

- per-channel fields: types, units, low-pass, high-pass and notch frequencies, bad channel flags, gain and transducers (EEG-like channels), coil type, geometry and code (MEG), stimulation intensity (MEP); fields missing in `obj1` are added (with default values for its channels)
- recording fields missing in `obj1` (e.g. line frequency, reference, NIRS wavelengths and optode labels, MEG device-to-head transform); if both objects have different values, those of `obj1` are kept (with a warning)
- NIRS optical channels (described by optode pairs) of `obj2`, with the other wavelengths of their optode pairs: placed before the channels of `obj1` (optical channels are the first channels), with their optode description; `obj1` must not contain optical channels
- channel locations

Markers and SSP projectors of `obj2` are not copied.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: object the channels are appended to
- `obj2::NeuroAnalyzer.NEURO`: object the channels are taken from
- `ch::Union{String, Vector{String}, Regex}="all"`: channel name(s) of `obj2`

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object

# Notes

The channel labels of `obj2` must not exist in `obj1` (rename them with `rename_channel()`). The data type of `obj1` is kept; functions for a channel type (e.g. referencing of EEG channels, NIRS processing) use the channels of that type.
"""
function append_channel(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex} = "all",
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj1)
    append_channel!(obj_new, obj2; ch = ch)

    return obj_new
end

"""
    append_channel!(obj1, obj2; <keyword arguments>)

Append channels of `obj2` to `obj1` in-place, e.g. EEG with auxiliary, EOG or ECG channels recorded separately, or EEG recorded simultaneously with NIRS.

The objects must have the same sampling rate, epoch length and number of epochs (the signals must be synchronized). The objects may contain channels of different types (e.g. EEG and NIRS): the header fields of all channel types are merged:

- per-channel fields: types, units, low-pass, high-pass and notch frequencies, bad channel flags, gain and transducers (EEG-like channels), coil type, geometry and code (MEG), stimulation intensity (MEP); fields missing in `obj1` are added (with default values for its channels)
- recording fields missing in `obj1` (e.g. line frequency, reference, NIRS wavelengths and optode labels, MEG device-to-head transform); if both objects have different values, those of `obj1` are kept (with a warning)
- NIRS optical channels (described by optode pairs) of `obj2`, with the other wavelengths of their optode pairs: placed before the channels of `obj1` (optical channels are the first channels), with their optode description; `obj1` must not contain optical channels
- channel locations

Markers and SSP projectors of `obj2` are not copied.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: object the channels are appended to; modified in-place
- `obj2::NeuroAnalyzer.NEURO`: object the channels are taken from
- `ch::Union{String, Vector{String}, Regex}="all"`: channel name(s) of `obj2`

# Returns

- `Nothing`

# Notes

The channel labels of `obj2` must not exist in `obj1` (rename them with `rename_channel()`). The data type of `obj1` is kept; functions for a channel type (e.g. referencing of EEG channels, NIRS processing) use the channels of that type.
"""
function append_channel!(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex} = "all",
)::Nothing
    # validate
    sr(obj1) == sr(obj2) || throw(
        ArgumentError(
            "Sampling rates differ ($(sr(obj1)) and $(sr(obj2)) Hz); resample one of the objects first.",
        ),
    )
    epoch_len(obj1) == epoch_len(obj2) && nepochs(obj1) == nepochs(obj2) || throw(
        ArgumentError(
            "Signal lengths differ: $(epoch_len(obj1)) × $(nepochs(obj1)) and $(epoch_len(obj2)) × $(nepochs(obj2)) samples × epochs.",
        ),
    )
    idx = get_channel(obj2; ch = ch)
    idx = idx isa Integer ? [Int64(idx)] : Vector{Int64}(idx)
    isempty(idx) && throw(ArgumentError("No channels selected."))

    r2 = obj2.header.recording
    _rec_filters!(r2)
    lbl = labels(obj2)[idx]
    n_old = nchannels(obj1)
    n0 = length(obj1.history)
    r1 = obj1.header.recording

    # NIRS optical channels of obj2 (described by optode pairs): added before the channels of
    # obj1 (optical channels are the first channels), with their optode description
    n_opt2 = _nirs_n_opt(r2)
    opt = Base.filter(<=(n_opt2), idx)
    if !isempty(opt)
        _nirs_n_opt(r1) == 0 || throw(
            ArgumentError(
                "OBJ1 already contains NIRS optical channels; NIRS channels are merged with optode functions.",
            ),
        )
        # (the other wavelengths of the same optode pairs are needed by od2conc())
        opt2 = _nirs_pair_wl(obj2, opt; msg = "Also appending")
        if length(opt2) > length(opt)
            idx = sort(unique(vcat(idx, opt2)))
            lbl = labels(obj2)[idx]
            opt = Base.filter(<=(n_opt2), idx)
        end
        idx = vcat(opt, setdiff(idx, opt))
        lbl = labels(obj2)[idx]
    end

    # recording-level fields of obj2 missing (or empty) in obj1, before the channels are
    # appended (which adds default values of some of them)
    miss1 = Set(k for k in keys(r2) if !haskey(r1, k) || _rec_empty(r1[k]))

    _append_channel!(
        obj1;
        data = obj2.data[idx, :, :],
        label = lbl,
        type = r2[:channel_type][idx],
        unit = r2[:unit][idx],
        nirs_opt = !isempty(opt),
    )

    # per-channel header fields of obj2: copied (fields missing in obj1 are added, with
    # default values for its channels)
    new = (n_old + 1):(n_old + length(idx))
    ch_n = nchannels(obj1)
    for f in (:lowpass, :highpass, :bad_channel)
        r1[f][new] = r2[f][idx]
    end
    r1[:notch][new] = [copy(n) for n in r2[:notch][idx]]
    for f in keys(_REC_CH_FIELDS)
        _rec_has_ch_field(r2, f, nchannels(obj2)) || continue
        _rec_ensure_ch_field!(r1, f, ch_n)
        r1[f][new] = r2[f][idx]
    end
    if _rec_has_ch_field(r2, :coil_loc, nchannels(obj2))
        _rec_has_ch_field(r1, :coil_loc, ch_n) || (r1[:coil_loc] = zeros(ch_n, 12))
        r1[:coil_loc][new, :] = r2[:coil_loc][idx, :]
    end
    # SSP projectors of obj2 are not copied (their vectors refer to the channels of obj2)
    if haskey(r1, :ssp_channels) && any(r1[:ssp_channels][new])
        r1[:ssp_channels][new] .= false
        _warn("SSP projectors of OBJ2 are not copied.")
    end

    # recording-level fields of obj2 missing in obj1 (e.g. line frequency, reference, NIRS
    # wavelengths and optode labels, MEG device-to-head transform); differences are reported
    skip = Set([
        :label, :channel_type, :unit, :bad_channel, :lowpass, :highpass, :notch,
        :channel_order, :coil_loc, :magnetometers, :gradiometers, :ssp_data, :ssp_labels,
        :wavelength_index, :optode_pairs, :data_type, :file_name, :file_size_mb,
        :file_type, :sampling_rate, keys(_REC_CH_FIELDS)...,
    ])
    differ = Symbol[]
    for (k, v) in r2
        k in skip && continue
        if k in miss1
            r1[k] = deepcopy(v)
        elseif !_rec_empty(v) && !isequal(r1[k], v) &&
               k in (:line_frequency, :reference, :wavelengths, :coord_system, :length_unit)
            push!(differ, k)
        end
    end
    isempty(differ) || _warn(
        "Header fields differ, the values of OBJ1 are kept: $(join(string.(differ), ", ")).",
    )

    # NIRS optical channels: moved to the front, optode description of obj2
    if !isempty(opt)
        n_o = length(opt)
        _select_channels!(obj1, vcat(new[1:n_o], 1:n_old, new[(n_o + 1):end]))
        r1[:wavelengths] = copy(r2[:wavelengths])
        r1[:wavelength_index] = r2[:wavelength_index][opt]
        r1[:optode_pairs] = r2[:optode_pairs][opt, :]
        for f in (:src_labels, :det_labels, :optode_labels)
            haskey(r2, f) && (r1[f] = copy(r2[f]))
        end
        _info("NIRS optical channels of OBJ2 are placed before the channels of OBJ1.")
    end

    # locations of the appended channels (NIRS optical channels: of their optodes)
    if DataFrames.nrow(obj2.locs) > 0
        ll = isempty(opt) ? lbl :
             vcat(lbl, get(r2, :optode_labels, String[]), get(r2, :src_labels, String[]),
                  get(r2, :det_labels, String[]))
        l2 = obj2.locs[in(ll).(obj2.locs[!, :label]), :]
        l2 = l2[.!in(obj1.locs[!, :label]).(l2[!, :label]), :]
        if DataFrames.nrow(l2) > 0
            obj1.locs = DataFrames.nrow(obj1.locs) == 0 ? copy(l2) :
                vcat(obj1.locs, l2; cols = :union)
        end
    end

    _history!(obj1, "append_channel", (; ch = lbl); args = (obj2,), from = n0)

    return nothing
end

"""
    insert_channel(obj; <keyword arguments>)

Return a copy of `obj` with new channel(s) inserted at a given position.

The per-channel header fields of the new channels are set as by `append_channel()`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `data::Array{<:Number, 3}`: channel data, shape (ch_n, epoch_len, epoch_n)
- `label::Union{String, Vector{String}}`: channel label(s)
- `type::Union{String, Vector{String}}`: channel type(s)
- `unit::Union{String, Vector{String}}`: channel unit(s)
- `pos::Int64`: position of the (first) new channel: `1` inserts before the first channel, `nchannels(obj) + 1` appends (as `append_channel()`)

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function insert_channel(
    obj::NeuroAnalyzer.NEURO;
    data::Array{<:Number, 3},
    label::Union{String, Vector{String}},
    type::Union{String, Vector{String}},
    unit::Union{String, Vector{String}},
    pos::Int64,
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    insert_channel!(obj_new; data = data, label = label, type = type, unit = unit, pos = pos)

    return obj_new
end

"""
    insert_channel!(obj; <keyword arguments>)

Insert new channel(s) at a given position in-place.

The per-channel header fields of the new channels are set as by `append_channel()`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `data::Array{<:Number, 3}`: channel data, shape (ch_n, epoch_len, epoch_n)
- `label::Union{String, Vector{String}}`: channel label(s)
- `type::Union{String, Vector{String}}`: channel type(s)
- `unit::Union{String, Vector{String}}`: channel unit(s)
- `pos::Int64`: position of the (first) new channel: `1` inserts before the first channel, `nchannels(obj) + 1` appends (as `append_channel()`)

# Returns

- `Nothing`
"""
function insert_channel!(
    obj::NeuroAnalyzer.NEURO;
    data::Array{<:Number, 3},
    label::Union{String, Vector{String}},
    type::Union{String, Vector{String}},
    unit::Union{String, Vector{String}},
    pos::Int64,
)::Nothing
    # validate
    ch_n = length(obj.data) > 0 ? nchannels(obj) : 0
    _in(pos, (1, ch_n + 1), "pos")

    # appended, then moved to pos (all per-channel fields follow the channels)
    n0 = length(obj.history)
    append_channel!(obj; data = data, label = label, type = type, unit = unit)
    n_new = nchannels(obj) - ch_n
    if pos <= ch_n
        _select_channels!(obj, vcat(1:(pos - 1), (ch_n + 1):(ch_n + n_new), pos:ch_n))
    end

    _history!(
        obj,
        "insert_channel",
        (; label = labels(obj)[pos:(pos + n_new - 1)], type, unit, pos, data);
        from = n0,
    )

    return nothing
end
