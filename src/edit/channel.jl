export signal_channels
export get_channel_bytype
export get_channel_bywl
export get_channel_type
export set_channel_type
export set_channel_type!
export get_channel
export rename_channel
export rename_channel!
export edit_channel
export edit_channel!
export replace_channel
export replace_channel!
export add_labels
export add_labels!
export add_channel
export add_channel!

"""
    signal_channels(obj)

Return all signal (e.g. EEG or MEG) channels; signal is determined by `:data_type` variable in `obj.header.recording`). For MEG data type, 'meg', `grad` and `mag` channels are returned.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns
 
- `chs::Vector{Int64}`
"""
function signal_channels(obj::NeuroAnalyzer.NEURO)

    dt = Symbol(obj.header.recording[:data_type])

    if dt === :eeg
        chs = union(get_channel_bytype(obj, type=[:eeg, :ecog, :seeg, :ref, :eog]))
    elseif dt === :seeg
        chs = union(get_channel_bytype(obj, type=[:eeg, :ecog, :seeg, :ref, :eog]))
    elseif dt === :ecog
        chs = union(get_channel_bytype(obj, type=[:eeg, :ecog, :seeg, :ref, :eog]))
    elseif dt === :meg
        chs = union(get_channel_bytype(obj, type=[:meg, :mag, :grad]))
    elseif dt === :nirs
        chs = union(get_channel_bytype(obj, type=[:nirs_int, :nirs_od, :nirs_dmean, :nirs_dvar, :nirs_dskew, :nirs_mua, :nirs_musp, :nirs_hbo, :nirs_hbr, :nirs_hbt, :nirs_h2o, :nirs_lipid, :nirs_bfi, :nirs_hrf_dod, :nirs_hrf_dmean, :nirs_hrf_dvar, :nirs_hrf_dskew, :nirs_hrf_hbo, :nirs_hrf_hbr, :nirs_hrf_hbt, :nirs_hrf_bfi, :nirs_aux]))
    elseif dt === :erp
        chs = union(get_channel_bytype(obj, type=[:eeg, :ecog, :seeg, :meg, :mag, :grad, :nirs_int, :nirs_od, :nirs_dmean, :nirs_dvar, :nirs_dskew, :nirs_mua, :nirs_musp, :nirs_hbo, :nirs_hbr, :nirs_hbt, :nirs_h2o, :nirs_lipid, :nirs_bfi, :nirs_hrf_dod, :nirs_hrf_dmean, :nirs_hrf_dvar, :nirs_hrf_dskew, :nirs_hrf_hbo, :nirs_hrf_hbr, :nirs_hrf_hbt, :nirs_hrf_bfi, :nirs_aux]))
    else
        chs = get_channel_bytype(obj, type=dt)
    end

    return chs

end

"""
    signal_channels(dt)

Return all signal (e.g. EEG or MEG) channels for a given data type. For MEG data type, 'meg', `grad` and `mag` channels are returned.

# Arguments

- `dt::String`: data type of an recording object
- `ct::Vector{String}`: list of channel types

# Returns
 
- `chs::Vector{Int64}`
"""
function signal_channels(dt::Symbol, ct::Vector{String})

    if dt === :meg
        chs = union(get_channel_bytype(ct, type=[:meg, :mag, :grad]))
    elseif dt === :nirs
        chs = union(get_channel_bytype(ct, type=[:nirs_int, :nirs_od, :nirs_dmean, :nirs_dvar, :nirs_dskew, :nirs_mua, :nirs_musp, :nirs_hbo, :nirs_hbr, :nirs_hbt, :nirs_h2o, :nirs_lipid, :nirs_bfi, :nirs_hrf_dod, :nirs_hrf_dmean, :nirs_hrf_dvar, :nirs_hrf_dskew, :nirs_hrf_hbo, :nirs_hrf_hbr, :nirs_hrf_hbt, :nirs_hrf_bfi, :nirs_aux]))
    elseif dt === :erp
        chs = union(get_channel_bytype(ct, type=[:eeg, :meg, :mag, :grad, :nirs_int, :nirs_od, :nirs_dmean, :nirs_dvar, :nirs_dskew, :nirs_mua, :nirs_musp, :nirs_hbo, :nirs_hbr, :nirs_hbt, :nirs_h2o, :nirs_lipid, :nirs_bfi, :nirs_hrf_dod, :nirs_hrf_dmean, :nirs_hrf_dvar, :nirs_hrf_dskew, :nirs_hrf_hbo, :nirs_hrf_hbr, :nirs_hrf_hbt, :nirs_hrf_bfi, :nirs_aux]))
    else
        chs = get_channel_bytype(ct, type=dt)
    end

    return chs

end

"""
    get_channel_bytype(obj; type)

Return channel number(s) for channel of `type` type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `type::Union{Symbol, Vector{Symbol}}=:all`: channel type

# Returns

- `ch_idx::Vector{Int64}`
"""
function get_channel_bytype(obj::NeuroAnalyzer.NEURO; type::Union{Symbol, Vector{Symbol}}=:all)

    if type isa Symbol
        _check_var(type, NeuroAnalyzer.channel_types, "type")
    else
        for idx in type
            _check_var(idx, NeuroAnalyzer.channel_types, "type")
        end
    end
        
    if type === :all
        ch_idx = _c(nchannels(obj))
    elseif type isa Symbol
        ch_idx = Vector{Int64}()
        for idx in 1:nchannels(obj)
            lowercase(obj.header.recording[:channel_type][idx]) == string(type) && push!(ch_idx, idx)
        end
    else
        ch_idx = Vector{Int64}()
        for idx1 in 1:nchannels(obj)
            for idx2 in type
                lowercase(obj.header.recording[:channel_type][idx1]) == string(idx2) && push!(ch_idx, idx1)
            end
        end        
    end

    return ch_idx

end

"""
    get_channel_bytype(ct; type)

Return channel number(s) for channel of `type` type.

# Arguments

- `ct::Vector{String}`: list of channel types
- `type::Union{Symbol, Vector{Symbol}}=:all`: channel type

# Returns

- `ch_idx::Vector{Int64}`
"""
function get_channel_bytype(ct::Vector{String}; type::Union{Symbol, Vector{Symbol}}=:all)

    if type isa Symbol
        _check_var(type, channel_types, "type")
    else
        for idx in eachindex(type)
            _check_var(type[idx], channel_types, "type")
        end
    end
        
    if type === :all
        ch_idx = _c(length(ct))
    elseif type isa Symbol
        ch_idx = Vector{Int64}()
        for idx in 1:length(ct)
            lowercase(ct[idx]) == string(type) && (push!(ch_idx, idx))
        end
    else
        ch_idx = Vector{Int64}()
        for idx1 in 1:length(ct)
            for idx2 in eachindex(type)
                lowercase(ct[idx]) == string(type[idx2]) && (push!(ch_idx, idx1))
            end
        end        
    end

    length(ch_idx) == 1 && (ch_idx = ch_idx[1])

    return ch_idx

end

"""
    get_channel_bywl(obj; wl)

Return NIRS channel number(s) for wavelength `wl`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `wl::Real`: wavelength (in nm)

# Returns

- `ch_idx::Vector{Int64}`
"""
function get_channel_bywl(obj::NeuroAnalyzer.NEURO; wl::Real)

    _check_datatype(obj, [:nirs])
    @assert wl in obj.header.recording[:wavelengths] "OBJ does not contain data for $wl wavelength. Available wavelengths: $(obj.header.recording[:wavelengths])."

    wl_idx = findfirst(isequal(wl), obj.header.recording[:wavelengths])
    ch_idx = Int64[]
    for idx in eachindex(obj.header.recording[:wavelength_index])
        obj.header.recording[:wavelength_index][idx] == wl_idx && push!(ch_idx, idx)
    end

    return ch_idx

end

"""
    get_channel_type(obj; ch, type)

Get channel type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, String}`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function get_channel_type(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, String})

    # create new dataset
    types = obj.header.recording[:channel_type]
    
    if ch isa String
        ch_found = nothing
        for idx in eachindex(clabels)
            if ch == clabels[idx]
                types[idx] = type
                ch_found = idx
            end
        end
        @assert ch_found !== nothing "Channel name ($ch) does not match signal labels."
    else
        _check_channels(obj, ch)
    end

    return obj.header.recording[:channel_type][ch]

end

"""
    set_channel_type(obj; ch, type)

Set channel type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, String}`
- `type::String`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function set_channel_type(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, String}, type::String)

    type = lowercase(type)
    NeuroAnalyzer._check_var(type, string.(NeuroAnalyzer.channel_types), "type")

    clabels = labels(obj)

    # create new dataset
    obj_new = deepcopy(obj)
    types = obj_new.header.recording[:channel_type]
    
    if ch isa String
        ch_found = nothing
        for idx in eachindex(clabels)
            if ch == clabels[idx]
                types[idx] = type
                ch_found = idx
            end
        end
        @assert ch_found !== nothing "Channel name ($ch) does not match signal labels."
    else
        _check_channels(obj, ch)
        types[ch] = type
    end
    obj_new.header.recording[:channel_type] = types
    
    # add entry to :history field
    push!(obj_new.history, "set_channel_type(OBJ, ch=$ch, type=$type)")

    return obj_new

end

"""
    set_channel_type!(obj; ch, new_name)

Set channel type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, String}`
- `type::String`
"""
function set_channel_type!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, String}, type::String)

    obj_new = set_channel_type(obj, ch=ch, type=type)
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

    if ch isa String
        # get channel by name
        ch_idx = nothing
        for idx in eachindex(clabels)
            if lowercase(ch) == lowercase(clabels[idx])
                ch_idx = idx
            end
        end
        @assert ch_idx !== nothing "Channel name ($ch) does not match signal labels."
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
    @assert !(name in clabels) "Channel $name already exist."

    if ch isa String
        # get channel by name
        ch_found = nothing
        for idx in eachindex(clabels)
            if ch == clabels[idx]
                clabels[idx] = name
                ch_found = idx
            end
        end
        @assert ch_found !== nothing "Channel name ($ch )does not match channel labels."
    else
        # get channel by number
        _check_channels(obj, ch)
        clabels[ch] = name
    end
    
    # rename channel locs label
    if ch isa String
        l_idx = NeuroAnalyzer._find_bylabel(obj_new.locs, ch)
    else
        l_idx = NeuroAnalyzer._find_bylabel(obj_new.locs, obj.header.recording[:labels][ch])
    end
    if length(l_idx) == 1
        l_idx = l_idx[1]
        obj_new.locs[l_idx, :labels] = name
    end

    push!(obj_new.history, "rename_channel(OBJ, ch=$ch, name=$name)")

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
    obj.locs = obj_new.locs

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
    
    @assert value !== nothing "value cannot be empty."
    _check_channels(obj, ch)
    _check_var(field, [:channel_type, :labels], "field")    

    obj_new = deepcopy(obj)
    @assert obj_new.header.recording[field][ch] isa typeof(value) "field type ($(eltype(obj_new.header.recording[field]))) does not mach value type ($(typeof(value)))."
    obj_new.header.recording[field][ch] = value

    push!(obj_new.history, "edit_channel(OBJ, ch=$ch, field=$field, value=$value)")   

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

    if ch isa String
        for idx in eachindex(clabels)
            if ch == clabels[idx]
                ch_idx = idx
            end
        end
        @assert ch_idx !== nothing "Channel name ($ch) does not match OBJ labels."
    else
        _check_channels(obj, ch)
        ch_idx = ch
    end

    obj_new = deepcopy(obj)
    @assert size(s) == (1, epoch_len(obj_new), nepochs(obj_new)) "signal size ($(size(s))) must be the same as channel size ($(size(obj_new.data[ch_idx, :, :]))."

    obj_new.data[ch_idx, :, :] = s

    reset_components!(obj_new)
    push!(obj_new.history, "replace_channel(OBJ, ch=$ch, s=$s")

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

    @assert length(clabels) == nchannels(obj) "clabels length must be $(nchannels(obj))."
    
    obj_new = deepcopy(obj)
    obj_new.header.recording[:labels] = clabels

    push!(obj_new.history, "add_labels(OBJ, clabels=$clabels")
 
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

"""
    add_channel(obj; data, label, type)

Add channel(s) data to empty `NeuroAnalyzer.NEURO` object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `data::Array{<:Number, 3}`: channel(s) data
- `label::Union{String, Vector{String}}=string.(_c(size(data, 1)))`: channel(s) label(s)
- `type::Union{Symbol, Vector{Symbol}}`: channel(s) type(s)

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function add_channel(obj::NeuroAnalyzer.NEURO; data::Array{<:Number, 3}, label::Union{String, Vector{String}}=string.(_c(size(data, 1))), type::Union{Symbol, Vector{Symbol}}, unit::Union{String, Vector{String}}=repeat([""], size(data, 1)))

    @assert length(obj.data) == 0 "OBJ already contains data."
    @assert length(label) == size(data, 1) "Number of labels and number of data channels must be equal."
    @assert length(type) == size(data, 1) "Number of channel types and number of data channels must be equal."
    @assert length(unit) == size(data, 1) "Number of channel units and number of data channels must be equal."

    for idx in eachindex(type)
        @assert type[idx] in channel_types "Unknown channel type $(type[idx])."
    end

    obj_new = deepcopy(obj)
    obj_new.data = data
    obj_new.header.recording[:labels] = label
    obj_new.header.recording[:channel_type] = string.(type)
    obj_new.header.recording[:units] = unit

    push!(obj_new.history, "add_channel(OBJ, data, label=$label, type=$type, unit=$unit)")

    return obj_new

end

"""
    add_channel!(obj; data, label, type)

Add channel(s) data to empty `NeuroAnalyzer.NEURO` object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `data::Array{<:Number, 3}`: channel(s) data
- `label::Union{String, Vector{String}}=string.(_c(size(data, 1)))`: channel(s) label(s)
- `type::Union{Symbol, Vector{Symbol}}`: channel(s) type(s)
"""
function add_channel!(obj::NeuroAnalyzer.NEURO; data::Array{<:Number, 3}, label::Union{String, Vector{String}}=string.(_c(size(data, 1))), type::Union{Symbol, Vector{Symbol}}, unit::Union{String, Vector{String}}=repeat([""], size(data, 1)))

    obj_new = add_channel(obj, data=data, label=label, type=type, unit=unit)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time
    obj.history = obj_new.history

    return nothing

end