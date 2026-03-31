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

Get channel type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel

# Returns

- `String`
"""
function channel_type(obj::NeuroAnalyzer.NEURO; ch::String)::String

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    length(ch) == 1 || throw(ArgumentError("ch must resolve to exactly one channel."))
    ch = ch[1]

    return obj.header.recording[:channel_type][ch]
end

"""
    set_channel_type(obj; <keyword arguments>)

Set channel type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `type::String`: new type

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
    ch = get_channel(ch; ch = ch)
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

Set channel type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `type::String`

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

Rename channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `name::String`: new name

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function rename_channel(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    name::String,
)::NeuroAnalyzer.NEURO

    # create new dataset
    obj_new = deepcopy(obj)

    clabels = obj_new.header.recording[:label]
    name in clabels && throw(ArgumentError("Channel $name already exist."))

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)[1]
    obj_new.header.recording[:label][ch] = name

    # rename label in locs
    l_idx = _find_bylabel(obj_new.locs, labels(obj)[ch])[1]
    !isnothing(l_idx) && (obj_new.locs[l_idx, :label] = name)

    push!(obj_new.history, "rename_channel(obj; ch=$ch, name=$name)")

    return obj_new
end

"""
    rename_channel!(obj; <keyword arguments>)

Rename channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `name::String`: new name

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

Edit channel properties (`:channel_type` or `:label`) in `OBJ.header.recording`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `field::Symbol`
- `value::String`

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
    isnothing(value) && throw(ArgumentError("value cannot be empty."))
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

Edit channel properties (`:channel_type` or `:label`) in `OBJ.header.recording`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `field::Symbol`
- `value::String`

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

Replace channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `s::AbstractArray`: signal array, shape (channels, samples, epochs)

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
                "signal size ($(size(s))) must be the same as channel size ($(size(obj.data[ch, :, :])).",
            ),
        )
    datatype(obj) == "meg" && size(obj.header.recording[:ssp_data]) != (0,) ||
        _warn(
            "OBJ contains SSP projections data, you should apply them before modifying OBJ data.",
        )

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    length(ch) == 1 || throw(ArgumentError("ch must resolve to exactly one channel."))
    ch = ch[1]

    # create new dataset
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = s

    push!(obj_new.history, "replace_channel(obj; ch=$ch, s")

    return obj_new
end

"""
    replace_channel!(obj; <keyword arguments>)

Replace channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `s::Array{Float64, 3}`: signal to replace with

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

Add channel labels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `clabels::Vector{String}`

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
    push!(obj_new.history, "add_label(OBJ, clabels=$clabels")

    return obj_new
end

"""
    add_label!(obj; <keyword arguments>)

Add channel labels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `clabels::Vector{String}`

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

Add channels data to an empty `NeuroAnalyzer.NEURO` object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `data::Array{<:Number, 3}`: channels data
- `label::Union{String, Vector{String}}`: channels labels
- `type::Union{String, Vector{String}}`: channels types

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

    # validate
    if length(obj.data) > 0
        signal_len(obj) == size(data, 2) ||
            throw(
                ArgumentError(
                    "Epoch length of the new data ($(size(data, 2))) and the object data ($(signal_len(obj)))must be equal.",
                ),
            )
        nepochs(obj) == size(data, 3) ||
            throw(
                ArgumentError(
                    "Number of epochs of the new data ($(size(data, 3))) and the object data ($(nepochs(obj))) must be equal.",
                ),
            )
    end
    length(label) == size(data, 1) ||
        throw(
            ArgumentError(
                "Number of labels ($(length(label))) and number of data channels ($(size(data, 1))) must be equal.",
            ),
        )
    length(type) == size(data, 1) ||
        throw(
            ArgumentError(
                "Number of channel types ($(length(type))) and number of data channels ($(size(data, 1))) must be equal.",
            ),
        )
    length(unit) == size(data, 1) ||
        throw(
            ArgumentError(
                "Number of channel units ($(length(unit))) and number of data channels ($(size(data, 1))) must be equal.",
            ),
        )

    for idx in eachindex(type)
        type[idx] in channel_types ||
            throw(ArgumentError("Unknown channel type $(type[idx])."))
    end

    datatype(obj) == "meg" && size(obj.header.recording[:ssp_data]) != (0,) ||
        _warn(
            "OBJ contains SSP projections data, you should apply them before modifying OBJ data.",
        )

    # create new dataset
    obj_new = deepcopy(obj)

    if length(obj.data) > 0
        obj_new.data = [obj.data; data]
        obj_new.header.recording[:label] = [obj.header.recording[:label], label]
        obj_new.header.recording[:channel_type] = [
            obj.header.recording[:channel_type],
            string.(type),
        ]
        obj_new.header.recording[:unit] = [obj.header.recording[:unit], unit]
        obj_new.header.recording[:channel_order] = [
            obj_new.header.recording[:channel_order],
            collect(
                maximum(
                    obj_new.header.recording[:channel_order],
                ):(
                    maximum(obj_new.header.recording[:channel_order]) + size(
                    data, 1,
                )
                ),
            ),
        ]
        obj_new.header.recording[:bad_channel] = [
            obj_new.header.recording[:bad_channel],
            zeros(Bool, size(data, 1)),
        ]
    else
        obj_new.data = data
        obj_new.header.recording[:label] = label
        obj_new.header.recording[:channel_type] = string.(type)
        obj_new.header.recording[:unit] = unit
        obj_new.header.recording[:channel_order] = collect(1:size(data, 1))
        obj_new.header.recording[:bad_channel] = zeros(Bool, size(data, 1))
    end

    push!(obj_new.history, "add_channel(OBJ, data, label=$label, type=$type, unit=$unit)")

    return obj_new
end

"""
    add_channel!(obj; <keyword arguments>)

Add channels data to an empty `NeuroAnalyzer.NEURO` object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `data::Array{<:Number, 3}`: channels data
- `label::Union{String, Vector{String}}`: channels labels
- `type::Union{String, Vector{String}}`: channels types

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
