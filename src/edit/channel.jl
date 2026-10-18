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
export add_channel
export add_channel!

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
    datatype(obj) in _MEG_LIKE_DATA &&
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
    # validate
    type = lowercase(type)
    _check_var(type, string.(channel_types), "type")
    ch_idx = _single_channel(obj, ch)

    # create new dataset
    obj_tmp = deepcopy(obj)
    rec = obj_tmp.header.recording

    rec[:channel_type][ch_idx] = type
    rec[:channel_order] = _sort_channels(rec[:channel_type])
    datatype(obj) in _MEG_LIKE_DATA && _update_meg_sensor_lists!(rec, ch_idx, type)

    push!(obj_tmp.history, "set_channel_type(obj; ch=$(labels(obj)[ch_idx]), type=$type)")

    return obj_tmp
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
function set_channel_type!(obj::NeuroAnalyzer.NEURO; ch::String, type::String)::Nothing
    obj_tmp = set_channel_type(obj; ch = ch, type = type)
    obj.header = obj_tmp.header
    obj.history = obj_tmp.history

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
    isempty(name) && throw(ArgumentError("name cannot be empty."))
    name in labels(obj) && throw(ArgumentError("Channel \"$name\" already exists."))
    ch_idx = _single_channel(obj, ch)
    old = labels(obj)[ch_idx]

    obj_tmp = deepcopy(obj)
    obj_tmp.header.recording[:label][ch_idx] = name
    _rename_locs!(obj_tmp.locs, old, name)

    push!(obj_tmp.history, "rename_channel(obj; ch=$old, name=$name)")

    return obj_tmp
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
function rename_channel!(obj::NeuroAnalyzer.NEURO; ch::String, name::String)::Nothing
    obj_tmp = rename_channel(obj; ch = ch, name = name)
    obj.header = obj_tmp.header
    obj.history = obj_tmp.history
    obj.locs = obj_tmp.locs

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
    # validate
    isempty(value) && throw(ArgumentError("value cannot be empty."))
    _check_var(field, [:channel_type, :label], "field")
    ch_idx = _single_channel(obj, ch)

    obj_tmp = if field == :label
        rename_channel(obj; ch = ch, name = value)
    else
        set_channel_type(obj; ch = ch, type = value)
    end

    # record the user's call instead of the delegated one
    obj_tmp.history[end] =
        "edit_channel(obj; ch=$(labels(obj)[ch_idx]), field=$field, value=$value)"

    return obj_tmp
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
    obj_tmp = edit_channel(obj; ch = ch, field = field, value = value)
    obj.header = obj_tmp.header
    obj.history = obj_tmp.history
    obj.locs = obj_tmp.locs

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

    # create new dataset
    obj_tmp = deepcopy(obj)
    obj_tmp.data[ch_idx, :, :] = s

    push!(obj_tmp.history, "replace_channel(obj; ch=$(labels(obj)[ch_idx]), s)")

    return obj_tmp
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
function replace_channel!(obj::NeuroAnalyzer.NEURO; ch::String, s::AbstractArray)::Nothing
    obj_tmp = replace_channel(obj; ch = ch, s = s)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history

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
    # validate
    length(clabels) == nchannels(obj) ||
        throw(ArgumentError("clabels length must be $(nchannels(obj))."))
    any(isempty, clabels) && throw(ArgumentError("clabels cannot contain empty labels."))
    allunique(clabels) || throw(ArgumentError("clabels must be unique."))

    # create new dataset
    obj_tmp = deepcopy(obj)
    old = labels(obj)
    obj_tmp.header.recording[:label] = copy(clabels)

    # rename locs in one pass (avoids clashes when labels are permuted)
    if nrow(obj_tmp.locs) > 0
        mapping = Dict(zip(old, clabels))
        obj_tmp.locs[!, :label] = [get(mapping, l, l) for l in obj_tmp.locs[!, :label]]
    end

    push!(obj_tmp.history, "add_label(obj, clabels=$clabels)")

    return obj_tmp
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
    obj_tmp = add_label(obj; clabels = clabels)
    obj.header = obj_tmp.header
    obj.history = obj_tmp.history
    obj.locs = obj_tmp.locs

    return nothing
end

"""
    add_channel(obj; <keyword arguments>)

Return a copy of `obj` with new channel(s) appended.

If `obj.data` is empty the object is initialized with the provided data (time axes are created from the sampling rate); otherwise the new channels are concatenated along the first (channel) dimension.

All per-channel header fields are extended with defaults: `""` for prefiltering, transducers and coil type, `1.0` for gain, a zero row for MEG sensor geometry (`:coil_loc`) and `0` for its coil code (`:coil_id`), `false` for SSP channel mask, `0` for MEP stimulation intensity. New MEG channels of type `"mag"`/`"grad"` are added to the magnetometer/gradiometer lists. NIRS channels must be added with optode functions.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `data::Array{<:Number, 3}`: channel data, shape (ch_n, epoch_len, epoch_n)
- `label::Union{String, Vector{String}}`: channel label(s)
- `type::Union{String, Vector{String}}`: channel type(s)
- `unit::Union{String, Vector{String}}`: channel unit(s)

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function add_channel(
    obj::NeuroAnalyzer.NEURO;
    data::Array{<:Number, 3},
    label::Union{String, Vector{String}},
    type::Union{String, Vector{String}},
    unit::Union{String, Vector{String}},
)::NeuroAnalyzer.NEURO
    ch_n = size(data, 1)
    ch_n > 0 || throw(ArgumentError("data must contain at least one channel."))

    # normalize to vectors for uniform handling below
    label_v = label isa String ? [label] : copy(label)
    type_v = lowercase.(type isa String ? [type] : type)
    unit_v = unit isa String ? [unit] : copy(unit)

    rec = obj.header.recording
    dt = rec[:data_type]
    dt == "nirs" &&
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

    # create new dataset
    obj_tmp = deepcopy(obj)
    r = obj_tmp.header.recording
    n_old = has_data ? nchannels(obj) : 0

    # extend a per-channel field (or create it for an empty object)
    function _extend!(field::Symbol, default)
        old = (has_data && haskey(r, field)) ? r[field] : typeof(default)[]
        r[field] = vcat(old, fill(default, ch_n))
        return nothing
    end

    obj_tmp.data = has_data ? cat(obj.data, data; dims = 1) : Array{Float64}(data)
    r[:label] = vcat(old_labels, label_v)
    r[:channel_type] = vcat(has_data ? r[:channel_type] : String[], type_v)
    r[:unit] = vcat(has_data ? r[:unit] : String[], unit_v)
    _extend!(:bad_channel, false)

    if dt in _EEG_LIKE_DATA
        _extend!(:prefiltering, "")
        _extend!(:transducers, "")
        _extend!(:gain, 1.0)

    elseif dt in ("sensors", "eda", "tpt")
        _extend!(:prefiltering, "")

    elseif dt == "mep"
        _extend!(:stimulation_intensity, 0)
        _extend!(:coil_type, "")

    elseif dt in _MEG_LIKE_DATA
        _extend!(:prefiltering, "")
        _extend!(:coil_type, "")
        _extend!(:ssp_channels, false)
        if haskey(r, :coil_loc) && !isempty(r[:coil_loc])
            r[:coil_loc] = vcat(r[:coil_loc], zeros(ch_n, 12))
        end
        if haskey(r, :coil_id) && !isempty(r[:coil_id])
            r[:coil_id] = vcat(r[:coil_id], zeros(Int64, ch_n))
        end
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
        fs > 0 || throw(ArgumentError("Sampling rate must be set to initialize an empty object."))
        n_s = size(data, 2)
        n_e = size(data, 3)
        obj_tmp.time_pts = _time_axis(n_s * n_e, fs)
        obj_tmp.epoch_time = _time_axis(n_s, fs)
    end

    push!(
        obj_tmp.history,
        "add_channel(obj; data, label=$label_v, type=$type_v, unit=$unit_v)",
    )

    return obj_tmp
end

"""
    add_channel!(obj; <keyword arguments>)

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
function add_channel!(
    obj::NeuroAnalyzer.NEURO;
    data::Array{<:Number, 3},
    label::Union{String, Vector{String}},
    type::Union{String, Vector{String}},
    unit::Union{String, Vector{String}},
)::Nothing
    obj_tmp = add_channel(obj; data = data, label = label, type = type, unit = unit)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.time_pts = obj_tmp.time_pts
    obj.epoch_time = obj_tmp.epoch_time
    obj.history = obj_tmp.history

    return nothing
end
