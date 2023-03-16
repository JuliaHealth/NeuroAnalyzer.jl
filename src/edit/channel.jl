export channel_type
export channel_type!
export get_channel
export rename_channel
export rename_channel!
export edit_channel
export edit_channel!
export replace_channel
export replace_channel!
export add_labels
export add_labels!

"""
    channel_type(obj; ch, type)

Change channel type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, String}`
- `type::String`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function channel_type(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, String}, type::String)

    type = lowercase(type)
    clabels = labels(obj)

    # create new dataset
    obj_new = deepcopy(obj)
    types = obj_new.header.recording[:channel_type]
    
    if typeof(ch) == String
        ch_found = nothing
        for idx in eachindex(clabels)
            if ch == clabels[idx]
                types[idx] = type
                ch_found = idx
            end
        end
        if ch_found === nothing
            throw(ArgumentError("Channel name ($ch) does not match signal labels."))
        end
    else
        _check_channels(obj, ch)
        types[ch] = type
    end
    obj_new.header.recording[:channel_type] = types
    
    # add entry to :history field
    push!(obj_new.header.history, "channel_type(OBJ, ch=$ch, type=$type)")

    return obj_new

end

"""
    channel_type!(obj; ch, new_name)

Change channel type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, String}`
- `type::String`
"""
function channel_type!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, String}, type::String)

    obj_new = channel_type(obj, ch=ch, type=type)
    obj.header = obj_new.header
    obj.history = obj_new.history

    return nothing

end

"""
    get_channel(obj; ch)

Return channel number (if provided by name) or name (if provided by number).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, String}`: channel number or name

# Returns

- `ch_idx::Union{Int64, String}`: channel number or name
"""
function get_channel(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, String})

    clabels = labels(obj)

    if typeof(ch) == String
        # get channel by name
        ch_idx = nothing
        for idx in eachindex(clabels)
            if lowercase(ch) == lowercase(clabels[idx])
                ch_idx = idx
            end
        end
        if ch_idx === nothing
            throw(ArgumentError("Channel name ($ch) does not match signal labels."))
        end
        return ch_idx
    else
        # get channel by number
        _check_channels(obj, ch)
        return clabels[ch]
    end

end

"""
    rename_channel(obj; ch, name)

Rename channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, String}`: channel number or name
- `name::String`: new name

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function rename_channel(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, String}, name::String)

    # create new dataset
    obj_new = deepcopy(obj)
    clabels = labels(obj_new)
    name in clabels && throw(ArgumentError("Channel $name already exist."))

    if typeof(ch) == String
        # get channel by name
        ch_found = nothing
        for idx in eachindex(clabels)
            if ch == clabels[idx]
                clabels[idx] = name
                ch_found = idx
            end
        end
        if ch_found === nothing
            throw(ArgumentError("Channel name ($ch )does not match channel labels."))
        end
    else
        # get channel by number
        _check_channels(obj, ch)
        clabels[ch] = name
    end
    obj_new.header.recording[:labels] = clabels
    
    push!(obj_new.header.history, "rename_channel(OBJ, ch=$ch, name=$name)")

    return obj_new

end

"""
    rename_channel!(obj; ch, name)

Rename channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, String}`: channel number or name
- `name::String`: new name
"""
function rename_channel!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, String}, name::String)

    obj_new = rename_channel(obj, ch=ch, name=name)
    obj.header = obj_new.header
    obj.history = obj_new.history

    return nothing

end

"""
    edit_channel(obj; ch, field, value)

Edit channel properties (`:channel_type` or `:labels`) in `OBJ.header.recording`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Int64`
- `field::Symbol`
- `value::Any`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function edit_channel(obj::NeuroAnalyzer.NEURO; ch::Int64, field::Symbol, value::Any)
    
    value === nothing && throw(ArgumentError("value cannot be empty."))
    _check_channels(obj, ch)
    _check_var(field, [:channel_type, :labels], "field")    

    obj_new = deepcopy(obj)
    typeof(obj_new.header.recording[field][ch]) == typeof(value) || throw(ArgumentError("field type ($(eltype(obj_new.header.recording[field]))) does not mach value type ($(typeof(value)))."))
    obj_new.header.recording[field][ch] = value

    push!(obj_new.header.history, "edit_channel(OBJ, ch=$ch, field=$field, value=$value)")   

    return obj_new

end

"""
    edit_channel!(obj; ch, field, value)

Edit channel properties (`:channel_type` or `:labels`) in `OBJ.header.recording`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Int64`
- `field::Symbol`
- `value::Any`
"""
function edit_channel!(obj::NeuroAnalyzer.NEURO; ch::Int64, field::Symbol, value::Any)
    
    obj_new = edit_channel(obj, ch=ch, field=field, value=value)
    obj.header = obj_new.header
    obj.history = obj_new.history

    return nothing

end

"""
    replace_channel(obj; ch, s)

Replace channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, String}`: channel number or name
- `s::Array{Float64, 3}`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function replace_channel(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, String}, s::Array{Float64, 3})

    ch_idx = nothing
    clabels = labels(obj)

    if typeof(ch) == String
        for idx in eachindex(clabels)
            if ch == clabels[idx]
                ch_idx = idx
            end
        end
        ch_idx === nothing && throw(ArgumentError("Channel name ($ch) does not match OBJ labels."))
    else
        _check_channels(obj, ch)
        ch_idx = ch
    end

    obj_new = deepcopy(obj)
    size(s) == (1, epoch_len(obj_new), epoch_n(obj_new)) || throw(ArgumentError("signal size ($(size(s))) must be the same as channel size ($(size(obj_new.data[ch_idx, :, :]))."))

    obj_new.data[ch_idx, :, :] = s

    reset_components!(obj_new)
    push!(obj_new.header.history, "replace_channel(OBJ, ch=$ch, s=$s")

    return obj_new

end

"""
    replace_channel!(obj; ch, s)

Replace channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, String}`: channel number or name
- `s::Array{Float64, 3}`: signal to replace with
"""
function replace_channel!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, String}, s::Array{Float64, 3})

    obj_new = replace_channel(obj, ch=ch, s=s)
    obj.header = obj_new.header
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end

"""
    add_labels(obj; labels)

Add channel labels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `clabels::Vector{String}`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function add_labels(obj::NeuroAnalyzer.NEURO; clabels::Vector{String})

    length(clabels) == channel_n(obj) || throw(ArgumentError("clabels length must be $(channel_n(obj))."))
    
    obj_new = deepcopy(obj)
    obj_new.header.recording[:labels] = clabels

    push!(obj_new.header.history, "add_labels(OBJ, clabels=$clabels")
 
    return obj_new
end

"""
    add_labels!(obj::NeuroAnalyzer.NEURO; clabels::Vector{String})

Add OBJ channel labels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `clabels::Vector{String}`
"""
function add_labels!(obj::NeuroAnalyzer.NEURO; clabels::Vector{String})

    obj_new = add_labels(obj, clabels=clabels)
    obj.header = obj_new.header
    obj.history = obj_new.history

    return nothing

end
