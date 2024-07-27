export get_channel
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
    get_channel(obj; <keyword arguments>)

Return list of channel names of specified type or their numbers if names are specified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}=""`: channel name or list of channel names
- `type::Union{String, Vector{String}}="all"`: channels types
- `wl::Real`: return NIRS channels for wavelength (in nm)

# Returns

- `ch::Vector{String}`
"""
function get_channel(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}="", type::Union{String, Vector{String}}="all", wl::Real=0)

    # return physical channel numbers
    ch != "" && return _ch_idx(obj, ch)

    # return channel names
    ch = String[]
    isa(type, String) && (type = [type])
    for idx in type
        _check_var(idx, channel_types, "type")
    end

    l = labels(obj)
    if wl == 0        
        if type == ["all"]
            ch = l
        else
            for type_idx in eachindex(type)
                for ch_idx in eachindex(obj.header.recording[:channel_type])
                    obj.header.recording[:channel_type][ch_idx] == type[type_idx] && push!(ch, l[ch_idx])
                end
            end
        end
    else
        _check_datatype(obj, ["nirs"])
        @assert wl in obj.header.recording[:wavelengths] "OBJ does not contain data for $wl wavelength. Available wavelengths: $(obj.header.recording[:wavelengths])."
        wl_idx = findfirst(isequal(wl), obj.header.recording[:wavelengths])
        for ch_idx in eachindex(obj.header.recording[:wavelength_index])
            obj.header.recording[:wavelength_index][ch_idx] == wl_idx && push!(ch, l[ch_idx])
        end
    end

    return ch

end

"""
    channel_type(obj; <keyword arguments>)

Get channel type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name

# Returns

- `cht::String`
"""
function channel_type(obj::NeuroAnalyzer.NEURO; ch::String)

    ch = get_channel(obj, ch=ch)
    cht = obj.header.recording[:channel_type][ch]

    return cht[1]

end

"""
    set_channel_type(obj; <keyword arguments>)

Set channel type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `type::String`: new type

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function set_channel_type(obj::NeuroAnalyzer.NEURO; ch::String, type::String)

    type = lowercase(type)
    _check_var(type, string.(channel_types), "type")

    # create new dataset
    obj_new = deepcopy(obj)
    ch = get_channel(obj_new, ch=ch)[1]
    obj_new.header.recording[:channel_type][ch] = type

    # add entry to :history field
    push!(obj_new.history, "set_channel_type(OBJ, ch=$ch, type=$type)")

    return obj_new

end

"""
    set_channel_type!(obj; <keyword arguments>)

Set channel type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `type::String`
"""
function set_channel_type!(obj::NeuroAnalyzer.NEURO; ch::String, type::String)

    obj_new = set_channel_type(obj, ch=ch, type=type)
    obj.header = obj_new.header
    obj.history = obj_new.history

    return nothing

end

"""
    rename_channel(obj; <keyword arguments>)

Rename channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `name::String`: new name

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function rename_channel(obj::NeuroAnalyzer.NEURO; ch::String, name::String)

    # create new dataset
    obj_new = deepcopy(obj)
    clabels = obj_new.header.recording[:label]
    @assert !(name in clabels) "Channel $name already exist."

    ch = get_channel(obj, ch=ch)[1]
    obj_new.header.recording[:label][ch] = name

    # rename label in locs
    l_idx = _find_bylabel(obj_new.locs, labels(obj)[ch])[1]
    !isnothing(l_idx) && (obj_new.locs[l_idx, :label] = name)

    push!(obj_new.history, "rename_channel(OBJ, ch=$ch, name=$name)")

    return obj_new

end

"""
    rename_channel!(obj; <keyword arguments>)

Rename channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `name::String`: new name
"""
function rename_channel!(obj::NeuroAnalyzer.NEURO; ch::String, name::String)

    obj_new = rename_channel(obj, ch=ch, name=name)
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.locs = obj_new.locs

    return nothing

end

"""
    edit_channel(obj; <keyword arguments>)

Edit channel properties (`:channel_type` or `:label`) in `OBJ.header.recording`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `field::Symbol`
- `value::String`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function edit_channel(obj::NeuroAnalyzer.NEURO; ch::String, field::Symbol, value::String)

    @assert value !== nothing "value cannot be empty."
    ch = get_channel(obj, ch=ch)[1]
    _check_var(field, [:channel_type, :label], "field")

    obj_new = deepcopy(obj)
    obj_new.header.recording[field][ch] = value

    push!(obj_new.history, "edit_channel(OBJ, ch=$ch, field=$field, value=$value)")

    return obj_new

end

"""
    edit_channel!(obj; <keyword arguments>)

Edit channel properties (`:channel_type` or `:label`) in `OBJ.header.recording`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `field::Symbol`
- `value::String`
"""
function edit_channel!(obj::NeuroAnalyzer.NEURO; ch::String, field::Symbol, value::String)

    obj_new = edit_channel(obj, ch=ch, field=field, value=value)
    obj.header = obj_new.header
    obj.history = obj_new.history

    return nothing

end

"""
    replace_channel(obj; <keyword arguments>)

Replace channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `s::AbstractArray`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function replace_channel(obj::NeuroAnalyzer.NEURO; ch::String, s::AbstractArray)

    @assert size(s) == (1, epoch_len(obj), nepochs(obj)) "signal size ($(size(s))) must be the same as channel size ($(size(obj.data[ch, :, :]))."

    ch = get_channel(obj, ch=ch)[1]
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = s

    reset_components!(obj_new)
    push!(obj_new.history, "replace_channel(OBJ, ch=$ch, s")

    return obj_new

end

"""
    replace_channel!(obj; <keyword arguments>)

Replace channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `s::Array{Float64, 3}`: signal to replace with
"""
function replace_channel!(obj::NeuroAnalyzer.NEURO; ch::String, s::Array{Float64, 3})

    obj_new = replace_channel(obj, ch=ch, s=s)
    obj.header = obj_new.header
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end

"""
    add_label(obj; <keyword arguments>)

Add channel labels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `clabels::Vector{String}`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function add_label(obj::NeuroAnalyzer.NEURO; clabels::Vector{String})

    @assert length(clabels) == nchannels(obj) "clabels length must be $(nchannels(obj))."

    obj_new = deepcopy(obj)
    obj_new.header.recording[:label] = clabels

    push!(obj_new.history, "add_label(OBJ, clabels=$clabels")

    return obj_new

end

"""
    add_label!(obj; <keyword arguments>)

Add channel labels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `clabels::Vector{String}`
"""
function add_label!(obj::NeuroAnalyzer.NEURO; clabels::Vector{String})

    obj_new = add_label(obj, clabels=clabels)
    obj.header = obj_new.header
    obj.history = obj_new.history

    return nothing

end

"""
    add_channel(obj; <keyword arguments>)

Add channels data to an empty `NeuroAnalyzer.NEURO` object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `data::Array{<:Number, 3}`: channels data
- `label::Union{String, Vector{String}}`: channels labels
- `type::Union{String, Vector{String}}`: channels types

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function add_channel(obj::NeuroAnalyzer.NEURO; data::Array{<:Number, 3}, label::Union{String, Vector{String}}, type::Union{String, Vector{String}}, unit::Union{String, Vector{String}})

    if length(obj.data) > 0
        @assert signal_len(obj) == size(data, 2) "Epoch length of the new data and the object data must be equal."
        @assert nepochs(obj) == size(data, 3) "Number of epochs of the new data and the object data must be equal."
    end
    @assert length(label) == size(data, 1) "Number of labels and number of data channels must be equal."
    @assert length(type) == size(data, 1) "Number of channel types and number of data channels must be equal."
    @assert length(unit) == size(data, 1) "Number of channel units and number of data channels must be equal."

    for idx in eachindex(type)
        @assert type[idx] in channel_types "Unknown channel type $(type[idx])."
    end

    obj_new = deepcopy(obj)
    if length(obj.data) > 0
        obj_new.data = [obj.data; data]
        obj_new.header.recording[:label] = [obj.header.recording[:label]; label]
        obj_new.header.recording[:channel_type] = [obj.header.recording[:channel_type]; string.(type)]
        obj_new.header.recording[:unit] = [obj.header.recording[:unit]; unit]
    else
        obj_new.data = data
        obj_new.header.recording[:label] = label
        obj_new.header.recording[:channel_type] = string.(type)
        obj_new.header.recording[:unit] = unit
    end

    push!(obj_new.history, "add_channel(OBJ, data, label=$label, type=$type, unit=$unit)")

    return obj_new

end

"""
    add_channel!(obj; <keyword arguments>)

Add channels data to an empty `NeuroAnalyzer.NEURO` object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `data::Array{<:Number, 3}`: channels data
- `label::Union{String, Vector{String}}`: channels labels
- `type::Union{String, Vector{String}}`: channels types
"""
function add_channel!(obj::NeuroAnalyzer.NEURO; data::Array{<:Number, 3}, label::Union{String, Vector{String}}, type::Union{String, Vector{String}}, unit::Union{String, Vector{String}})

    obj_new = add_channel(obj, data=data, label=label, type=type, unit=unit)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time
    obj.history = obj_new.history

    return nothing

end
