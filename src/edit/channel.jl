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
    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    length(ch) == 1 || throw(ArgumentError("ch must resolve to exactly one channel."))

    return obj.header.recording[:channel_type][ch[1]]
end

"""
    set_channel_type(obj; <keyword arguments>)

Return a copy of `obj` with the type of one channel changed.

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

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    length(ch) == 1 || throw(ArgumentError("ch must resolve to exactly one channel."))
    ch = ch[1]

    # create new dataset
    obj_new = deepcopy(obj)

    obj_new.header.recording[:channel_type][ch] = type
    push!(obj_new.history, "set_channel_type(obj; ch=$ch, type=$type)")

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
function set_channel_type!(obj::NeuroAnalyzer.NEURO; ch::String, type::String)::Nothing
    obj_new = set_channel_type(obj; ch = ch, type = type)
    obj.header = obj_new.header
    obj.history = obj_new.history

    return nothing
end

"""
    rename_channel(obj; <keyword arguments>)

Return a copy of `obj` with one channel renamed. Also updates the matching entry in `obj.locs` if one exists.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `name::String`: new name (must not already exist)

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function rename_channel(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    name::String,
)::NeuroAnalyzer.NEURO
    clabels = obj.header.recording[:label]
    name in clabels && throw(ArgumentError("Channel \"$name\" already exists."))

    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    length(ch) == 1 || throw(ArgumentError("ch must resolve to exactly one channel."))
    ch = ch[1]

    obj_new = deepcopy(obj)
    obj_new.header.recording[:label][ch] = name

    # update matching locs entry if present.
    l_result = _find_bylabel(obj_new.locs, labels(obj)[ch])
    if !isempty(l_result)
        l_idx = l_result isa Int64 ? l_result : l_result[1]
        obj_new.locs[l_idx, :label] = name
    end

    push!(obj_new.history, "rename_channel(obj; ch=$ch, name=$name)")

    return obj_new
end

"""
    rename_channel!(obj; <keyword arguments>)

Rename one channel in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `name::String`: new name (must not already exist)

# Returns

- `Nothing`
"""
function rename_channel!(obj::NeuroAnalyzer.NEURO; ch::String, name::String)::Nothing
    obj_new = rename_channel(obj; ch = ch, name = name)
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.locs = obj_new.locs

    return nothing
end

"""
    edit_channel(obj; <keyword arguments>)

Return a copy of `obj` with one field of one channel's header entry changed.

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

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    length(ch) == 1 || throw(ArgumentError("ch must resolve to exactly one channel."))
    ch = ch[1]

    # create new dataset
    obj_new = deepcopy(obj)
    obj_new.header.recording[field][ch] = value

    push!(obj_new.history, "edit_channel(obj; ch=$ch, field=$field, value=$value)")

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
    obj_new = edit_channel(obj; ch = ch, field = field, value = value)
    obj.header = obj_new.header
    obj.history = obj_new.history

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

    size(s) == (1, epoch_len(obj), nepochs(obj)) ||
        throw(
            ArgumentError(
                "data size ($(size(s, 2)) × $(size(s, 3))) must be the same as channel size ($(size(obj, 2)) × $(size(obj, 3))).",
            ),
        )
    datatype(obj) == "meg" && size(obj.header.recording[:ssp_data]) != (0,) &&
        _warn("OBJ contains SSP projections data; apply them before modifying data.")

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    length(ch) == 1 || throw(ArgumentError("ch must resolve to exactly one channel."))
    ch = ch[1]

    # create new dataset
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = s

    push!(obj_new.history, "replace_channel(obj; ch=$ch, s)")

    return obj_new
end

"""
    replace_channel!(obj; <keyword arguments>)

Replace one channel's data in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `s::Array{Float64, 3}`: replacement data, shape (1, epoch_len, nepochs)

# Returns

- `Nothing`
"""
function replace_channel!(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    s::Array{Float64, 3},
)::Nothing
    obj_new = replace_channel(obj; ch = ch, s = s)
    obj.header = obj_new.header
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing
end

"""
    add_label(obj; <keyword arguments>)

Return a copy of `obj` with channel labels replaced.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `clabels::Vector{String}`: new labels; must have length equal to `nchannels(obj)`

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function add_label(obj::NeuroAnalyzer.NEURO; clabels::Vector{String})::NeuroAnalyzer.NEURO
    # validate
    length(clabels) == nchannels(obj) ||
        throw(ArgumentError("clabels length must be $(nchannels(obj))."))

    # create new dataset
    obj_new = deepcopy(obj)

    obj_new.header.recording[:label] = clabels
    push!(obj_new.history, "add_label(obj, clabels=$clabels)")

    return obj_new
end

"""
    add_label!(obj; <keyword arguments>)

Replace channel labels in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `clabels::Vector{String}`: new labels; must have length equal to `nchannels(obj)`

# Returns

- `Nothing`
"""
function add_label!(obj::NeuroAnalyzer.NEURO; clabels::Vector{String})::Nothing
    obj_new = add_label(obj; clabels = clabels)
    obj.header = obj_new.header
    obj.history = obj_new.history

    return nothing
end

"""
    add_channel(obj; <keyword arguments>)

Return a copy of `obj` with new channel(s) appended.

If `obj.data` is empty the object is initialized with the provided data; otherwise the new channels are concatenated along the first (channel) dimension.

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

    # normalize to vectors for uniform handling below
    label_v = label isa String ? [label] : label
    type_v  = type isa String ? [type] : type
    unit_v  = unit isa String ? [unit] : unit

    # validate
    if length(obj.data) > 0
        signal_len(obj) == size(data, 2) || throw(
            ArgumentError(
                "Epoch length mismatch: new data has $(size(data, 2)), object has $(signal_len(obj)).",
            ),
        )
        nepochs(obj) == size(data, 3) || throw(
            ArgumentError(
                "Epoch count mismatch: new data has $(size(data, 3)), object has $(nepochs(obj)).",
            ),
        )
    end
    length(label_v) == ch_n || throw(
        ArgumentError(
            "Number of labels ($(length(label_v))) must equal number of new channels ($ch_n).",
        ),
    )
    length(type_v) == ch_n || throw(
        ArgumentError(
            "Number of types ($(length(type_v))) must equal number of new channels ($ch_n).",
        ),
    )
    length(unit_v) == ch_n || throw(
        ArgumentError(
            "Number of units ($(length(unit_v))) must equal number of new channels ($ch_n).",
        ),
    )

    for t in type_v
        t in channel_types || throw(ArgumentError("Unknown channel type \"$t\"."))
    end

    datatype(obj) == "meg" && size(obj.header.recording[:ssp_data]) != (0,) &&
        _warn("OBJ contains SSP projections data; apply them before modifying data.")

    # create new dataset
    obj_new = deepcopy(obj)

    if length(obj.data) > 0
        obj_new.data                            = cat(obj.data, data; dims = 1)
        obj_new.header.recording[:label]        = vcat(obj.header.recording[:label], label_v)
        obj_new.header.recording[:channel_type] = vcat(obj.header.recording[:channel_type], string.(type_v))
        obj_new.header.recording[:unit]         = vcat(obj.header.recording[:unit], unit_v)
        obj_new.header.recording[:bad_channel]  = vcat(obj.header.recording[:bad_channel], zeros(Bool, ch_n))

        max_ord = maximum(obj_new.header.recording[:channel_order])
        obj_new.header.recording[:channel_order] = vcat(
            obj.header.recording[:channel_order],
            collect((max_ord + 1):(max_ord + ch_n)),
        )
    else
        obj_new.data                             = data
        obj_new.header.recording[:label]         = label_v
        obj_new.header.recording[:channel_type]  = string.(type_v)
        obj_new.header.recording[:unit]          = unit_v
        obj_new.header.recording[:channel_order] = collect(1:ch_n)
        obj_new.header.recording[:bad_channel]   = zeros(Bool, ch_n)
    end

    push!(
        obj_new.history,
        "add_channel(obj; data, label=$label_v, type=$type_v, unit=$unit_v)",
    )

    return obj_new
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
    obj_new = add_channel(obj; data = data, label = label, type = type, unit = unit)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time
    obj.history = obj_new.history

    return nothing
end
