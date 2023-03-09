export channel_type
export channel_type!
export get_channel
export rename_channel
export rename_channel!
export extract_channel
export edit_channel
export edit_channel!
export replace_channel
export replace_channel!
export add_labels
export add_labels!

"""
    channel_type(obj; channel, type)

Change channel type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, String}`
- `type::String`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function channel_type(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, String}, type::String)

    type = lowercase(type)
    labels = labels(obj)

    # create new dataset
    obj_new = deepcopy(obj)
    types = obj_new.header.recording[:channel_type]
    
    if typeof(channel) == String
        channel_found = nothing
        for idx in eachindex(labels)
            if channel == labels[idx]
                types[idx] = type
                channel_found = idx
            end
        end
        if channel_found === nothing
            throw(ArgumentError("Channel name ($channel) does not match signal labels."))
        end
    else
        _check_channels(obj, channel)
        types[channel] = type
    end
    obj_new.header.recording[:channel_type] = types
    
    # add entry to :history field
    push!(obj_new.header.recording.history, "channel_type(OBJ, channel=$channel, type=$type)")

    return obj_new
end

"""
    channel_type!(obj; channel, new_name)

Change channel type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, String}`
- `type::String`
"""
function channel_type!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, String}, type::String)

    obj_tmp = channel_type(obj, channel=channel, type=type)
    obj.header = obj_tmp.header

    return nothing
end

"""
    get_channel(obj; channel)

Return channel number (if provided by name) or name (if provided by number).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, String}`: channel number or name

# Returns

- `channel_idx::Union{Int64, String}`: channel number or name
"""
function get_channel(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, String})

    labels = labels(obj)
    if typeof(channel) == String
        # get channel by name
        channel_idx = nothing
        for idx in eachindex(labels)
            if lowercase(channel) == lowercase(labels[idx])
                channel_idx = idx
            end
        end
        if channel_idx === nothing
            throw(ArgumentError("Channel name ($channel) does not match signal labels."))
        end
        return channel_idx
    else
        # get channel by number
        _check_channels(obj, channel)
        return labels[channel]
    end
end

"""
    rename_channel(obj; channel, name)

Rename channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, String}`: channel number or name
- `name::String`: new name

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function rename_channel(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, String}, name::String)

    # create new dataset
    obj_new = deepcopy(obj)
    labels = labels(obj_new)
    name in labels && throw(ArgumentError("Channel $name already exist."))

    if typeof(channel) == String
        # get channel by name
        channel_found = nothing
        for idx in eachindex(labels)
            if channel == labels[idx]
                labels[idx] = name
                channel_found = idx
            end
        end
        if channel_found === nothing
            throw(ArgumentError("Channel name ($channel )does not match channel labels."))
        end
    else
        # get channel by number
        _check_channels(obj, channel)
        labels[channel] = name
    end
    obj_new.header.recording[:labels] = labels
    
    # add entry to :history field
    push!(obj_new.header.recording.history, "rename_channel(OBJ, channel=$channel, name=$name)")

    return obj_new
end

"""
    rename_channel!(obj; channel, name)

Rename channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, String}`: channel number or name
- `name::String`: new name
"""
function rename_channel!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, String}, name::String)

    obj.header[:labels] = rename_channel(obj, channel=channel, name=name).header.recording[:labels]
    push!(obj.header[:history], "rename_channel!(OBJ, channel=$channel, name=$name)")

    return nothing
end

"""
    extract_channel(obj; channel)

Extract channel data.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, String}`: channel number or name

# Returns

- `channel::Vector{Float64}`
"""
function extract_channel(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, String})

    labels = labels(obj)
    if typeof(channel) == String
        # get channel by name
        channel_idx = nothing
        for idx in eachindex(labels)
            if channel == labels[idx]
                channel_idx = idx
            end
        end
        if channel_idx === nothing
            throw(ArgumentError("Channel name ($channel )does not match channel labels."))
        end
    else
        # get channel by number
        _check_channels(obj, channel)
        channel_idx = channel
    end    
    channel = reshape(obj.data[channel_idx, :, :], 1, epoch_len(obj), epoch_n(obj))

    return channel
end

"""
    edit_channel(obj; channel, field, value)

Edit channel properties.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Int64`
- `field::Symbol`
- `value::Any`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function edit_channel(obj::NeuroAnalyzer.NEURO; channel::Int64, field::Symbol, value::Any)
    
    value === nothing && throw(ArgumentError("value cannot be empty."))
    _check_channels(obj, channel)
    _check_var(field, [:channel_type, :labels], "field")    

    obj_new = deepcopy(obj)
    typeof(obj_new.header.recording[field][channel]) == typeof(value) || throw(ArgumentError("field type ($(eltype(obj_new.header.recording[field]))) does not mach value type ($(typeof(value)))."))
    obj_new.header.recording[field][channel] = value

    # add entry to :history field
    push!(obj_new.header.recording.history, "edit_channel(OBJ, channel=$channel, field=$field, value=$value)")   

    return obj_new
end

"""
    edit_channel!(obj; channel, field, value)

Edit channel properties.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Int64`
- `field::Symbol`
- `value::Any`
"""
function edit_channel!(obj::NeuroAnalyzer.NEURO; channel::Int64, field::Symbol, value::Any)
    
    obj_tmp = edit_channel(obj, channel=channel, field=field, value=value)
    obj.header = obj_tmp.header

    return nothing
end

"""
    replace_channel(obj; channel, signal)

Replace channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, String}`: channel number or name
- `signal::Array{Float64, 3}`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function replace_channel(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, String}, signal::Array{Float64, 3})

    channel_idx = nothing
    labels = labels(obj)
    if typeof(channel) == String
        for idx in eachindex(labels)
            if channel == labels[idx]
                channel_idx = idx
            end
        end
        channel_idx === nothing && throw(ArgumentError("Channel name ($channel) does not match signal labels."))
    else
        _check_channels(obj, channel)
        channel_idx = channel
    end

    obj_new = deepcopy(obj)
    size(signal) == (1, epoch_len(obj_new), epoch_n(obj_new)) || throw(ArgumentError("signal size ($(size(signal))) must be the same as channel size ($(size(obj_new.data[channel_idx, :, :]))."))
    obj_new.data[channel_idx, :, :] = signal
    reset_components!(obj_new)

    # add entry to :history field
    push!(obj_new.header.recording.history, "replace_channel(OBJ, channel=$channel, signal=$signal")

    return obj_new
end

"""
    replace_channel!(obj; channel, signal)

Replace channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, String}`: channel number or name
- `signal::Array{Float64, 3}`
"""
function replace_channel!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, String}, signal::Array{Float64, 3})

    obj_tmp = replace_channel(obj, channel=chanel, signal=signal)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    reset_components!(obj)

    return nothing
end

"""
    add_labels(obj; labels)

Add channel labels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `labels::Vector{String}`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function add_labels(obj::NeuroAnalyzer.NEURO; labels::Vector{String})

    length(labels) == channel_n(obj) || throw(ArgumentError("labels length must be $(channel_n(obj))."))
    
    obj_new = deepcopy(obj)
    obj_new.header.recording[:labels] = labels

    push!(obj_new.header.recording.history, "add_labels(OBJ, labels=$labels")
 
    return obj_new
end

"""
    add_labels!(obj::NeuroAnalyzer.NEURO; labels::Vector{String})

Add OBJ channel labels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `labels::Vector{String}`
"""
function add_labels!(obj::NeuroAnalyzer.NEURO; labels::Vector{String})

    length(labels) == channel_n(obj) || throw(ArgumentError("labels length must be $(channel_n(obj))."))
    obj_tmp = add_labels(obj, labels=labels)
    obj.header = obj_tmp.header

    return nothing
end
