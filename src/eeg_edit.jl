"""
    eeg_add_component(eeg; c, v)

Add component.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `c::Symbol`: component name
- `v::Any`: component value

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_add_component(eeg::NeuroAnalyzer.EEG; c::Symbol, v::Any)

    eeg_new = deepcopy(eeg)
    c in eeg_new.eeg_header[:components] && throw(ArgumentError("Component $c already exists. Use eeg_delete_component() to remove it prior the operation."))
    push!(eeg_new.eeg_header[:components], c)
    push!(eeg_new.eeg_components, v)
    push!(eeg_new.eeg_header[:history], "eeg_add_component(EEG, c=$c, v=$v)")

    return eeg_new
end

"""
    eeg_add_component!(eeg; c, v)

Add component.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `c::Symbol`: component name
- `v::Any`: component value
"""
function eeg_add_component!(eeg::NeuroAnalyzer.EEG; c::Symbol, v::Any)

    c in eeg.eeg_header[:components] && throw(ArgumentError("Component $c already exists. Use eeg_delete_component!() to remove it prior the operation."))
    push!(eeg.eeg_header[:components], c)
    push!(eeg.eeg_components, v)
    push!(eeg.eeg_header[:history], "eeg_add_component!(EEG, c=$c, v=$v)")

    return nothing
end

"""
    eeg_list_components(eeg)

List component names.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `components::Vector{Symbol}`
"""
function eeg_list_components(eeg::NeuroAnalyzer.EEG)
    return eeg.eeg_header[:components]
end

"""
    eeg_extract_component(eeg, c)

Extract component values.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `c::Symbol`: component name

# Returns

- `component::Any`
"""
function eeg_extract_component(eeg::NeuroAnalyzer.EEG; c::Symbol)

    c in eeg.eeg_header[:components] || throw(ArgumentError("Component $c does not exist. Use eeg_list_component() to view existing components."))
    
    for idx in eachindex(eeg.eeg_header[:components])
        if c == eeg.eeg_header[:components][idx]
            return eeg.eeg_components[idx]
        end
    end
end

"""
    eeg_delete_component(eeg; c)

Delete component. 

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `c::Symbol`: component name

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_delete_component(eeg::NeuroAnalyzer.EEG; c::Symbol)

    c in eeg.eeg_header[:components] || throw(ArgumentError("Component $c does not exist. Use eeg_list_component() to view existing components."))
    
    eeg_new = deepcopy(eeg)
    for idx in eachindex(eeg.eeg_header[:components])
        if c == eeg_new.eeg_header[:components][idx]
            deleteat!(eeg_new.eeg_components, idx)
            deleteat!(eeg_new.eeg_header[:components], idx)
            push!(eeg_new.eeg_header[:history], "eeg_delete_component(EEG, c=$c)")
            return eeg_new
        end
    end
end

"""
    eeg_delete_component!(eeg; c)

Delete component.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `c::Symbol`: component name
"""
function eeg_delete_component!(eeg::NeuroAnalyzer.EEG; c::Symbol)

    c in eeg.eeg_header[:components] || throw(ArgumentError("Component $c does not exist. Use eeg_list_component() to view existing components."))
    
    for idx in length(eeg.eeg_header[:components]):-1:1
        if c == eeg.eeg_header[:components][idx]
            deleteat!(eeg.eeg_components, idx)
            deleteat!(eeg.eeg_header[:components], idx)
            push!(eeg.eeg_header[:history], "eeg_delete_component(EEG, c=$c)")
        end
    end

    return nothing
end

"""
    eeg_reset_components(eeg)

Remove all components.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_reset_components(eeg::NeuroAnalyzer.EEG)
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_header[:components] = []
    eeg_new.eeg_components = []
    return eeg_new
end

"""
    eeg_reset_components!(eeg)

Remove all components.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_reset_components!(eeg::NeuroAnalyzer.EEG)
    eeg.eeg_header[:components] = []
    eeg.eeg_components = []
    return nothing
end

"""
    eeg_component_idx(eeg, c)

Return component index.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `c::Symbol`: component name

# Return

- `c_idx::Int64`
"""
function eeg_component_idx(eeg::NeuroAnalyzer.EEG; c::Symbol)

    c in eeg.eeg_header[:components] || throw(ArgumentError("Component $c does not exist. Use eeg_list_component() to view existing components."))
    c_idx = findfirst(isequal(c), eeg.eeg_header[:components])

    return c_idx
end

"""
    eeg_component_type(eeg, c)

Return component data type.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `c::Symbol`: component name

# Return

- `c_type::DataType`
"""
function eeg_component_type(eeg::NeuroAnalyzer.EEG; c::Symbol)

    c in eeg.eeg_header[:components] || throw(ArgumentError("Component $c does not exist. Use eeg_list_component() to view existing components."))
    c_idx = eeg_component_idx(eeg; c=c)
    c_type = typeof(eeg.eeg_components[c_idx])

    return c_type
end

"""
    eeg_rename_component(eeg, c_old, c_new)

Rename component.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `c_old::Symbol`: old component name
- `c_new::Symbol`: new component name

# Return

- `eeg_new:EEG`
"""
function eeg_rename_component(eeg::NeuroAnalyzer.EEG; c_old::Symbol, c_new::Symbol)

    c_old in eeg.eeg_header[:components] || throw(ArgumentError("Component $c_old does not exist. Use eeg_list_component() to view existing components."))
    c_new in eeg.eeg_header[:components] && throw(ArgumentError("Component $c_new already exists. Use eeg_list_component() to view existing components."))

    eeg_new = deepcopy(eeg)
    c_idx = eeg_component_idx(eeg, c=c_old)
    eeg_new.eeg_header[:components][c_idx] = c_new

    push!(eeg_new.eeg_header[:history], "eeg_rename_component(EEG, c_old=$c_old, c_new=$c_new)")

    return eeg_new
end

"""
    eeg_rename_component!(eeg, c_old, c_new)

Rename component.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `c_old::Symbol`: old component name
- `c_new::Symbol`: new component name
"""
function eeg_rename_component!(eeg::NeuroAnalyzer.EEG; c_old::Symbol, c_new::Symbol)

    c_old in eeg.eeg_header[:components] || throw(ArgumentError("Component $c_old does not exist. Use eeg_list_component() to view existing components."))
    c_new in eeg.eeg_header[:components] && throw(ArgumentError("Component $c_new already exists. Use eeg_list_component() to view existing components."))

    c_idx = eeg_component_idx(eeg, c=c_old)
    eeg.eeg_header[:components][c_idx] = c_new

    push!(eeg.eeg_header[:history], "eeg_rename_component!(EEG, c_old=$c_old, c_new=$c_new)")

    return nothing
end

"""
    eeg_delete_channel(eeg; channel)

Delete EEG channel(s).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to be removed

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_delete_channel(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})

    typeof(channel) <: AbstractRange && (channel = collect(channel))
    channel_n = eeg_channel_n(eeg)
    length(channel) > 1 && (channel = sort!(channel, rev=true))
    length(channel) == channel_n && throw(ArgumentError("You cannot delete all channels."))

    _check_channels(eeg, channel)

    eeg_new = deepcopy(eeg)

    # update headers
    for idx in channel
        loc = findfirst(isequal(lowercase(eeg_new.eeg_header[:labels][idx])), lowercase.(string.(eeg_new.eeg_locs[!, :labels])))
        loc !== nothing && deleteat!(eeg_new.eeg_locs, loc)
        deleteat!(eeg_new.eeg_header[:labels], idx)
        deleteat!(eeg_new.eeg_header[:channel_type], idx)
        deleteat!(eeg_new.eeg_header[:transducers], idx)
        deleteat!(eeg_new.eeg_header[:physical_dimension], idx)
        deleteat!(eeg_new.eeg_header[:prefiltering], idx)
        deleteat!(eeg_new.eeg_header[:gain], idx)
    end
    eeg_new.eeg_header[:channel_n] -= length(channel)

    # remove channel
    eeg_new.eeg_signals = eeg_new.eeg_signals[setdiff(1:end, (channel)), :, :]

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_delete_channel(EEG, $channel)")

    return eeg_new
end

"""
    eeg_delete_channel!(eeg; channel)

Delete EEG channel(s).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to be removed
"""
function eeg_delete_channel!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})

    typeof(channel) <: AbstractRange && (channel = collect(channel))
    channel_n = eeg_channel_n(eeg)
    length(channel) > 1 && (channel = sort!(channel, rev=true))
    length(channel) == channel_n && throw(ArgumentError("You cannot delete all channels."))

    _check_channels(eeg, channel)

    # update headers
    for idx in channel
        loc = findfirst(isequal(lowercase(eeg.eeg_header[:labels][idx])), lowercase.(string.(eeg.eeg_locs[!, :labels])))
        loc !== nothing && deleteat!(eeg.eeg_locs, loc)
        deleteat!(eeg.eeg_header[:labels], idx)
        deleteat!(eeg.eeg_header[:channel_type], idx)
        deleteat!(eeg.eeg_header[:transducers], idx)
        deleteat!(eeg.eeg_header[:physical_dimension], idx)
        deleteat!(eeg.eeg_header[:prefiltering], idx)
        deleteat!(eeg.eeg_header[:gain], idx)
    end
    eeg.eeg_header[:channel_n] -= length(channel)

    # remove channel
    eeg.eeg_signals = eeg.eeg_signals[setdiff(1:end, (channel)), :, :]

    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_delete_channel!(EEG, $channel)")

    return nothing
end

"""
    eeg_keep_channel(eeg; channel)

Keep EEG channel(s).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to keep

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_keep_channel(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})

    typeof(channel) <: AbstractRange && (channel = collect(channel))
    _check_channels(eeg, channel)

    channel_n = eeg_channel_n(eeg)
    channels_to_remove = setdiff(collect(1:channel_n), channel)
    length(channels_to_remove) == channel_n && throw(ArgumentError("You cannot delete all channels."))

    return eeg_delete_channel(eeg, channel=channels_to_remove)
end

"""
    eeg_keep_channel!(eeg; channel)

Keep EEG channel(s).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to keep
"""
function eeg_keep_channel!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})

    typeof(channel) <: AbstractRange && (channel = collect(channel))
    _check_channels(eeg, channel)

    channel_n = eeg_channel_n(eeg)
    channels_to_remove = setdiff(collect(1:channel_n), channel)
    # length(channels_to_remove) > 1 && sort!(channels_to_remove, rev=true)
    length(channels_to_remove) == channel_n && throw(ArgumentError("You cannot delete all channels."))

    eeg_delete_channel!(eeg, channel=channels_to_remove)
end

"""
    eeg_get_channel(eeg; channel)

Return EEG channel number (if provided by name) or name (if provided by number).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, String}`: channel number or name

# Returns

- `channel_idx::Union{Int64, String}`: channel number or name
"""
function eeg_get_channel(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, String})

    labels = eeg_labels(eeg)
    if typeof(channel) == String
        channel_idx = nothing
        for idx in eachindex(labels)
            if channel == labels[idx]
                channel_idx = idx
            end
        end
        if channel_idx === nothing
            throw(ArgumentError("Channel name ($channel) does not match signal labels."))
        end
        return channel_idx
    else
        _check_channels(eeg, channel)
        return labels[channel]
    end
end

"""
    eeg_rename_channel(eeg; channel, name)

Rename EEG channel.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, String}`: channel number or name
- `name::String`: new name

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_rename_channel(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, String}, name::String)

    # create new dataset
    eeg_new = deepcopy(eeg)
    labels = eeg_labels(eeg_new)
    name in labels && throw(ArgumentError("Channel $name already exist."))

    if typeof(channel) == String
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
        _check_channels(eeg, channel)
        labels[channel] = name
    end
    eeg_new.eeg_header[:labels] = labels
    
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_rename_channel(EEG, channel=$channel, name=$name)")

    return eeg_new
end

"""
    eeg_rename_channel!(eeg; channel, name)

Rename EEG channel.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, String}`: channel number or name
- `name::String`: new name
"""
function eeg_rename_channel!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, String}, name::String)

    eeg.eeg_header[:labels] = eeg_rename_channel(eeg, channel=channel, name=name).eeg_header[:labels]
    push!(eeg.eeg_header[:history], "eeg_rename_channel!(EEG, channel=$channel, name=$name)")

    return nothing
end

"""
    eeg_extract_channel(eeg; channel)

Extract EEG channel data.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, String}`: channel number or name

# Returns

- `channel::Vector{Float64}`
"""
function eeg_extract_channel(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, String})

    labels = eeg_labels(eeg)
    if typeof(channel) == String
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
        _check_channels(eeg, channel)
        channel_idx = channel
    end    
    eeg_channel = reshape(eeg.eeg_signals[channel_idx, :, :], 1, eeg_epoch_len(eeg), eeg_epoch_n(eeg))

    return eeg_channel
end

"""
    eeg_history(eeg)

Show EEG processing history.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_history(eeg::NeuroAnalyzer.EEG)
    return eeg.eeg_header[:history]
end

"""
    eeg_labels(eeg)

Return EEG channel labels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `labels::Vector{String}`
"""
function eeg_labels(eeg::NeuroAnalyzer.EEG)
    length(eeg.eeg_header[:labels]) == 0 && throw(ArgumentError("EEG has no labels."))
    return eeg.eeg_header[:labels]
end

"""
    eeg_sr(eeg)

Return EEG sampling rate.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `sr::Int64`
"""
function eeg_sr(eeg::NeuroAnalyzer.EEG)
    return eeg.eeg_header[:sampling_rate]
end

"""
    eeg_channel_n(eeg; type=:eeg)

Return number of EEG channels of `type`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `type::Vector{Symbol}=:all`: channel type: `:all`, `:eeg`, `:meg`, `:ecg`, `:eog`, `:emg`, `:ref`, `:mrk`

# Returns

- `channel_n::Int64`
"""
function eeg_channel_n(eeg::NeuroAnalyzer.EEG; type::Symbol=:all)

    _check_var(type, [:all, :eeg, :meg, :ecg, :eog, :emg, :ref, :mrk], "type")
    length(eeg.eeg_header[:channel_type]) == 0 && throw(ArgumentError("EEG has no channel types."))
    channel_n = 0
    for idx in 1:eeg.eeg_header[:channel_n]
        eeg.eeg_header[:channel_type][idx] == string(type) && (channel_n += 1)
    end
    type === :all && (channel_n = size(eeg.eeg_signals, 1))

    return channel_n
end

"""
    eeg_epoch_n(eeg)

Return number of EEG epochs.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `epoch_n::Int64`
"""
function eeg_epoch_n(eeg::NeuroAnalyzer.EEG)
    epoch_n = eeg.eeg_header[:epoch_n]
    return epoch_n
end

"""
    eeg_signal_len(eeg)

Return length of EEG signal.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `signal_len::Int64`
"""
function eeg_signal_len(eeg::NeuroAnalyzer.EEG)
    return eeg.eeg_header[:eeg_duration_samples]
end

"""
    eeg_epoch_len(eeg)

Return length of EEG epoch.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `epoch_len::Int64`
"""
function eeg_epoch_len(eeg::NeuroAnalyzer.EEG)
    epoch_len = eeg.eeg_header[:epoch_duration_samples]
    return epoch_len
end

"""
    eeg_info(eeg)

Show EEG info.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_info(eeg::NeuroAnalyzer.EEG)

    println("            Signal type: $(uppercase(eeg.eeg_header[:signal_type]))")
    println("              File name: $(eeg.eeg_header[:eeg_filename])")
    println("            File format: $(eeg.eeg_header[:eeg_filetype])")
    println("         File size [MB]: $(eeg.eeg_header[:eeg_filesize_mb])")
    println("       Memory size [MB]: $(round(Base.summarysize(eeg) / 1024^2, digits=2))")
    println("     Sampling rate (Hz): $(eeg_sr(eeg))")
    if eeg.eeg_header[:markers] == false
        println("            markers: no")
    else
        println("            markers: yes")
    end
    println("Signal length (samples): $(eeg_signal_len(eeg))")
    println("Signal length (seconds): $(round(eeg.eeg_header[:eeg_duration_seconds], digits=2))")
    println("     Number of channels: $(eeg_channel_n(eeg))")
    println("       Number of epochs: $(eeg_epoch_n(eeg))")
    println(" Epoch length (samples): $(eeg_epoch_len(eeg))")
    println(" Epoch length (seconds): $(round(eeg.eeg_header[:epoch_duration_seconds], digits=2))")
    if eeg.eeg_header[:reference] == ""
        println("         Reference type: unknown")
    else
        println("         Reference type: $(eeg.eeg_header[:reference])")
    end
    if length(eeg_labels(eeg)) == 0
        println("                 Labels: no")
    else
        println("                 Labels: yes")
    end
    if eeg.eeg_header[:channel_locations] == false
        println("      Channel locations: no")
    else
        println("      Channel locations: yes")
    end
    if eeg.eeg_header[:components] != []
        print("             Components: ")
        c = eeg_list_components(eeg)
        if length(c) == 1
            println(c[1])
        else
            for idx in 1:(length(c) - 1)
                print(c[idx], ", ")
            end
            println(c[end])
        end
    else
        println("             Components: no")
    end
    println("Channels:")
    for idx in eachindex(eeg.eeg_header[:labels])
        println("\tchannel: $idx\tlabel: $(rpad(eeg.eeg_header[:labels][idx], 16, " "))\ttype: $(uppercase(eeg.eeg_header[:channel_type][idx]))")
    end
end

"""
    eeg_epochs(eeg; epoch_n=nothing, epoch_len=nothing, average=false)

Split EEG into epochs. Return signal that is split either by epoch length or by number of epochs. 

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `epoch_n::Union{Int64, Nothing}=nothing`: number of epochs
- `epoch_len::Union{Int64, Nothing}`=nothing: epoch length in samples
- `average::Bool=false`: average all epochs, return one averaged epoch

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_epochs(eeg::NeuroAnalyzer.EEG; epoch_n::Union{Int64, Nothing}=nothing, epoch_len::Union{Int64, Nothing}=nothing, average::Bool=false)

    # split into epochs
    epochs = _make_epochs(eeg.eeg_signals, epoch_n=epoch_n, epoch_len=epoch_len, average=average)

    # create new dataset
    epoch_n = size(epochs, 3)
    epoch_duration_samples = size(epochs, 2)
    epoch_duration_seconds = size(epochs, 2) / eeg.eeg_header[:sampling_rate]
    eeg_duration_samples = size(epochs, 2) * size(epochs, 3)
    eeg_duration_seconds = eeg_duration_samples / eeg.eeg_header[:sampling_rate]
    eeg_time = collect(0:(1 / eeg.eeg_header[:sampling_rate]):epoch_duration_seconds)
    eeg_time = eeg_time[1:(end - 1)]
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals = epochs
    eeg_new.eeg_time = eeg_time

    # update epochs time
    fs = eeg_sr(eeg_new)
    ts = eeg.eeg_epochs_time[1]
    new_epochs_time = linspace(ts, ts + (epoch_duration_samples / fs), epoch_duration_samples)
    eeg_new.eeg_epochs_time = new_epochs_time

    eeg_new.eeg_header[:eeg_duration_samples] = eeg_duration_samples
    eeg_new.eeg_header[:eeg_duration_seconds] = eeg_duration_seconds
    eeg_new.eeg_header[:epoch_n] = epoch_n
    eeg_new.eeg_header[:epoch_duration_samples] = epoch_duration_samples
    eeg_new.eeg_header[:epoch_duration_seconds] = epoch_duration_seconds

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_epochs(EEG, epoch_n=$epoch_n, epoch_len=$epoch_len, average=$average)")

    return eeg_new
end

"""
    eeg_epochs!(eeg; epoch_n=nothing, epoch_len=nothing, average=false)

Split EEG into epochs. Return signal that is split either by epoch length or by number of epochs. 

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `epoch_n::Union{Int64, Nothing}=nothing`: number of epochs
- `epoch_len::Union{Int64, Nothing}`=nothing: epoch length in samples
- `average::Bool=false`: average all epochs, return one averaged epoch
"""
function eeg_epochs!(eeg::NeuroAnalyzer.EEG; epoch_n::Union{Int64, Nothing}=nothing, epoch_len::Union{Int64, Nothing}=nothing, average::Bool=false)

    eeg_tmp = eeg_epochs(eeg, epoch_n=epoch_n, epoch_len=epoch_len, average=average)
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg.eeg_signals = eeg_tmp.eeg_signals

    eeg_reset_components!(eeg)
    return nothing
end

"""
    eeg_extract_epoch(eeg; epoch)

Extract EEG epoch.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `epoch::Int64`: epoch index

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_extract_epoch(eeg::NeuroAnalyzer.EEG; epoch::Int64)

    _check_epochs(eeg, epoch)

    s_new = reshape(eeg.eeg_signals[:, :, epoch], eeg_channel_n(eeg), eeg_signal_len(eeg), 1)
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals = s_new
    eeg_new.eeg_epochs_time = eeg.eeg_epochs_time
    eeg_new.eeg_header[:epoch_n] = 1
    eeg_new.eeg_header[:eeg_duration_samples] = eeg_new.eeg_header[:epoch_duration_samples]
    eeg_new.eeg_header[:eeg_duration_seconds] = eeg_new.eeg_header[:epoch_duration_seconds]

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_extract_epoch(EEG, epoch=$epoch)")

    return eeg_new
end

"""
    eeg_trim(eeg; segment, remove_epochs)

Trim EEG signal by removing parts of the signal.

# Arguments

- EEG
- `segment::Tuple{Int64, Int64}`: segment to be removed (from, to) in samples
- `remove_epochs::Bool=true`: if true, remove epochs containing signal to trim or remove signal and re-epoch trimmed signal

# Returns

- EEG
"""
function eeg_trim(eeg::NeuroAnalyzer.EEG; segment::Tuple{Int64, Int64}, remove_epochs::Bool=true)

    if remove_epochs == false
        eeg_new = deepcopy(eeg)
        eeg_epoch_n(eeg) > 1 && (eeg_epochs!(eeg_new, epoch_n=1))
        _check_segment(eeg_new, segment[1], segment[2])
        eeg_new.eeg_signals = s_trim(eeg_new.eeg_signals, segment=segment)
        t_trimmed = collect(0:(1 / eeg_sr(eeg)):(size(eeg_new.eeg_signals, 2) / eeg_sr(eeg)))[1:(end - 1)]
        eeg_new.eeg_time = t_trimmed
        eeg_new.eeg_epochs_time = t_trimmed
        eeg_new.eeg_header[:eeg_duration_samples] -= length(segment[1]:segment[2])
        eeg_new.eeg_header[:eeg_duration_seconds] -= length(segment[1]:segment[2]) * (1 / eeg_sr(eeg))
        eeg_new.eeg_header[:epoch_duration_samples] -= length(segment[1]:segment[2])
        eeg_new.eeg_header[:epoch_duration_seconds] -= length(segment[1]:segment[2]) * (1 / eeg_sr(eeg))
        if eeg_epoch_n(eeg) > 1
            if eeg_epoch_len(eeg) <= eeg_signal_len(eeg_new)
                eeg_epochs!(eeg_new, epoch_len=eeg_epoch_len(eeg))
            else
                verbose == true && @info "Cannot apply original epoch length, returning single-epoch EEG."
            end
        end
    else
        eeg_epoch_n(eeg) == 1 && throw(ArgumentError("EEG has only one epoch, cannot use remove_epochs=true."))
        epochs = _t2epoch(eeg, segment[1], segment[2])
        verbose == true && @info "Removing $epochs epochs."
        eeg_new = eeg_delete_epoch(eeg, epoch=epochs)
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_trim(EEG, segment=$segment, remove_epochs=$remove_epochs)")

    return eeg_new
end

"""
    eeg_trim!(eeg; segment, remove_epochs)

Trim EEG signal by removing parts of the signal.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `segment::Tuple{Int64, Int64}`: segment to be removed (from, to) in samples
- `remove_epochs::Bool=true`: remove epochs containing signal to trim (remove_epochs=true) or remove signal and remove epoching
"""
function eeg_trim!(eeg::NeuroAnalyzer.EEG; segment::Tuple{Int64, Int64}, remove_epochs::Bool=true)

    eeg_tmp = eeg_trim(eeg, segment=segment, remove_epochs=remove_epochs)
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg.eeg_signals = eeg_tmp.eeg_signals

    eeg_reset_components!(eeg)
    return nothing
end

"""
    eeg_edit_header(eeg; field, value)

Change value of EEG header.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `field::Symbol`
- `value::Any`

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_edit_header(eeg::NeuroAnalyzer.EEG; field::Symbol, value::Any)

    value === nothing && throw(ArgumentError("value cannot be empty."))

    eeg_new = deepcopy(eeg)
    field in keys(eeg_new.eeg_header) || throw(ArgumentError("$field does not exist."))
    typeof(eeg_new.eeg_header[field]) == typeof(value) || throw(ArgumentError("field type ($(typeof(eeg_new.eeg_header[field]))) does not mach value type ($(typeof(value)))."))
    eeg_new.eeg_header[field] = value
    push!(eeg_new.eeg_header[:history], "eeg_edit(EEG, field=$field, value=$value)")    

    return eeg_new
end

"""
    eeg_edit_header!(eeg; field, value)

Change value of EEG header.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `field::Symbol`
- `value::Any`
"""
function eeg_edit_header!(eeg::NeuroAnalyzer.EEG; field::Symbol, value::Any)
    eeg.eeg_header = eeg_edit_header(eeg, field=field, value=value).eeg_header
    return nothing
end

"""
    eeg_show_header(eeg)

Show keys and values of EEG header.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_show_header(eeg::NeuroAnalyzer.EEG)
    for (key, value) in eeg.eeg_header
        println("$key: $value")
    end
end

"""
    eeg_delete_epoch(eeg; epoch)

Remove EEG epoch(s).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) to be removed

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_delete_epoch(eeg::NeuroAnalyzer.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_epoch_n(eeg) == 1 && throw(ArgumentError("You cannot delete the last epoch."))
    typeof(epoch) <: AbstractRange && (epoch = collect(epoch))
    length(epoch) == eeg_epoch_n(eeg) && throw(ArgumentError("You cannot delete all epochs."))
    length(epoch) > 1 && (epoch = sort!(epoch, rev=true))
    _check_epochs(eeg, epoch)

    eeg_new = deepcopy(eeg)

    # remove epoch
    eeg_new.eeg_signals = eeg_new.eeg_signals[:, :, setdiff(1:end, (epoch))]

    # update headers
    eeg_new.eeg_header[:epoch_n] -= length(epoch)
    epoch_n = eeg_new.eeg_header[:epoch_n]
    eeg_new.eeg_header[:eeg_duration_samples] = epoch_n * size(eeg.eeg_signals, 2)
    eeg_new.eeg_header[:eeg_duration_seconds] = round((epoch_n * size(eeg.eeg_signals, 2)) / eeg_sr(eeg), digits=2)
    eeg_new.eeg_header[:epoch_duration_samples] = size(eeg.eeg_signals, 2)
    eeg_new.eeg_header[:epoch_duration_seconds] = round(size(eeg.eeg_signals, 2) / eeg_sr(eeg), digits=2)

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_delete_epoch(EEG, $epoch)")
    
    return eeg_new
end

"""
    eeg_delete_epoch!(eeg; epoch)

Remove EEG epoch(s).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) to be removed
"""
function eeg_delete_epoch!(eeg::NeuroAnalyzer.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_tmp = eeg_delete_epoch(eeg, epoch=epoch)
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg.eeg_signals = eeg_tmp.eeg_signals

    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_delete_epoch!(EEG, $epoch)")

    return nothing
end

"""
    eeg_keep_epoch(eeg; epoch)

Keep EEG epoch(s).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) to keep

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_keep_epoch(eeg::NeuroAnalyzer.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_epoch_n(eeg) == 1 && throw(ArgumentError("EEG contains only one epoch."))
    typeof(epoch) <: AbstractRange && (epoch = collect(epoch))
    length(epoch) > 1 && (epoch = sort!(epoch, rev=true))
    _check_epochs(eeg, epoch)

    epoch_list = collect(1:eeg_epoch_n(eeg))
    epoch_to_remove = setdiff(epoch_list, epoch)

    length(epoch_to_remove) > 1 && (epoch_to_remove = sort!(epoch_to_remove, rev=true))

    eeg_new = eeg_delete_epoch(eeg, epoch=epoch_to_remove)
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_keep_epoch(EEG, $epoch)")    

    return eeg_new
end

"""
    eeg_keep_epoch!(eeg; epoch)

Keep EEG epoch(s).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) to keep
"""
function eeg_keep_epoch!(eeg::NeuroAnalyzer.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_tmp = eeg_keep_epoch(eeg, epoch=epoch)
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg.eeg_signals = eeg_tmp.eeg_signals

    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_keep_epoch(EEG, $epoch)")

    return nothing
end

"""
    eeg_detect_bad(eeg; method, ch_t)

Detect bad EEG channels and epochs.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))`: index of channels, default is all EEG/MEG channels
- `method::Vector{Symbol}=[:flat, :rmse, :rmsd, :euclid, :p2p, :var]`: detection method:
    - `:flat`: flat channel(s)
    - `:p2p`: peak-to-peak amplitude; good for detecting transient artifacts
    - `:var`: mean signal variance outside of 95%CI and variance inter-quartile outliers
    - `:rmse`: RMSE vs average channel outside of 95%CI
    - `:rmsd`: RMSD
    - `:euclid`: Euclidean distance
- `w::Int64=10`: window width in samples (signal is averaged within `w`-width window)
- `ftol::Float64=0.1`: tolerance (signal is flat within `-tol` to `+tol`), `eps()` gives very low tolerance
- `fr::Float64=0.3`: acceptable ratio (0.0 to 1.0) of flat segments within a channel before marking it as flat
- `p::Float64=0.95`: probability threshold (0.0 to 1.0) for marking channel as bad; also threshold for `:p2p` detection: above `mean + p * std` and below `mean - p * std`, here p (as percentile) will be converted to z-score (0.9 (90th percentile): 1.282, 0.95 (95th percentile): 1.645, 0.975 (97.5th percentile): 1.960, 0.99 (99th percentile): 2.326) 
- `tc::Float64=0.3`: threshold (0.0 to 1.0) of bad channels ratio to mark the epoch as bad

# Returns

Named tuple containing:
- `bad_m::Matrix{Bool}`: matrix of bad channels × epochs
- `bad_epochs::Vector{Int64}`: list of bad epochs
"""
function eeg_detect_bad(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type])), method::Vector{Symbol}=[:flat, :rmse, :rmsd, :euclid, :p2p], w::Int64=10, ftol::Float64=0.1, fr::Float64=0.3, p::Float64=0.95, tc::Float64=0.2)

    for idx in method
        _check_var(idx, [:flat, :rmse, :rmsd, :euclid, :var, :p2p], "method")
    end

    p < 0 || p > 1 && throw(ArgumentError("p must be ≥ 0.0 and ≤ 1.0"))
    tc < 0 || tc > 1 && throw(ArgumentError("t must be ≥ 0.0 and ≤ 1.0"))

    _check_channels(eeg, channel)
    channel_n = length(channel)
    epoch_n = eeg_epoch_n(eeg)    

    #signal = eeg_demean(eeg, channel=channel).eeg_signals
    signal = eeg.eeg_signals

    bad_m = zeros(Bool, channel_n, epoch_n)
    bad_epochs = Vector{Int64}()

    if :flat in method
        @inbounds @simd for epoch_idx in 1:epoch_n
            bad_channels_score = 0
            bad_channels = zeros(Bool, channel_n)
            Threads.@threads for channel_idx in 1:channel_n
                r = @views s_detect_channel_flat(signal[channel[channel_idx], :, epoch_idx], w=w, tol=ftol)
                if r > fr
                    bad_channels_score += 1
                    bad_channels[channel_idx] = true
                end
            end
            [bad_channels[channel_idx] == true && (bad_m[channel_idx, epoch_idx] = true) for channel_idx in 1:channel_n]
            (bad_channels_score / channel_n) > tc && push!(bad_epochs, epoch_idx)
        end
    end
    
    if :rmse in method
        @inbounds @simd for epoch_idx in 1:epoch_n
            ch_m = @views vec(median(signal[channel, :, epoch_idx], dims=1))
            bad_channels_score = 0
            bad_channels = zeros(Bool, channel_n)

            rmse_ch = zeros(channel_n)
            Threads.@threads for channel_idx in 1:channel_n
                rmse_ch[channel_idx] = @views s2_rmse(signal[channel[channel_idx], :, epoch_idx], ch_m)
            end
            Threads.@threads for channel_idx in 1:channel_n
                if rmse_ch[channel_idx] < HypothesisTests.confint(OneSampleTTest(rmse_ch))[1] || rmse_ch[channel_idx] > HypothesisTests.confint(OneSampleTTest(rmse_ch))[2]
                    bad_channels_score += 1
                    bad_channels[channel_idx] = true
                end
            end
            [bad_channels[channel_idx] == true && (bad_m[channel_idx, epoch_idx] = true) for channel_idx in 1:channel_n]
            (bad_channels_score / channel_n) > tc && push!(bad_epochs, epoch_idx)
        end
    end

    if :rmsd in method
        @inbounds @simd for epoch_idx in 1:epoch_n
            ch_m = @views vec(median(signal[channel, :, epoch_idx], dims=1))
            bad_channels_score = 0
            bad_channels = zeros(Bool, channel_n)

            rmsd_ch = zeros(channel_n)
            Threads.@threads for channel_idx in 1:channel_n
                rmsd_ch[channel_idx] = @views Distances.rmsd(signal[channel[channel_idx], :, epoch_idx], ch_m)
            end
            Threads.@threads for channel_idx in 1:channel_n
                if rmsd_ch[channel_idx] < HypothesisTests.confint(OneSampleTTest(rmsd_ch))[1] || rmsd_ch[channel_idx] > HypothesisTests.confint(OneSampleTTest(rmsd_ch))[2]
                    bad_channels_score += 1
                    bad_channels[channel_idx] = true
                end
            end
            [bad_channels[channel_idx] == true && (bad_m[channel_idx, epoch_idx] = true) for channel_idx in 1:channel_n]
            (bad_channels_score / channel_n) > tc && push!(bad_epochs, epoch_idx)
        end
    end

    if :euclid in method
        @inbounds @simd for epoch_idx in 1:epoch_n
            ch_m = @views vec(median(signal[channel, :, epoch_idx], dims=1))
            bad_channels_score = 0
            bad_channels = zeros(Bool, channel_n)

            ed_ch = zeros(channel_n)
            Threads.@threads for channel_idx in 1:channel_n
                ed_ch[channel_idx] = @views Distances.euclidean(signal[channel[channel_idx], :, epoch_idx], ch_m)
            end
            Threads.@threads for channel_idx in 1:channel_n
                if ed_ch[channel_idx] < HypothesisTests.confint(OneSampleTTest(ed_ch))[1] || ed_ch[channel_idx] > HypothesisTests.confint(OneSampleTTest(ed_ch))[2]
                    bad_channels_score += 1
                    bad_channels[channel_idx] = true
                end
            end
            [bad_channels[channel_idx] == true && (bad_m[channel_idx, epoch_idx] = true) for channel_idx in 1:channel_n]
            (bad_channels_score / channel_n) > tc && push!(bad_epochs, epoch_idx)
        end
    end

    if :var in method
        s_v = @views var(signal[channel, :, :], dims=2)
        # mean variance
        s_mv = @views vec(mean(s_v, dims=3))
        # variance outliers
        o = reshape(outlier_detect(vec(s_v), method=:iqr), channel_n, epoch_n)

        @inbounds @simd for epoch_idx in 1:epoch_n
            bad_channels_score = 0
            bad_channels = zeros(Bool, channel_n)
            ch_v = @views vec(var(signal[channel, :, epoch_idx], dims=2))
            s_mv = vcat(s_mv, ch_v)
            Threads.@threads for channel_idx in 1:channel_n
                #if ch_v[channel_idx] > HypothesisTests.confint(OneSampleTTest(s_mv))[2] || o[channel_idx, epoch_idx]
                if o[channel_idx, epoch_idx]
                    bad_channels_score += 1
                    bad_channels[channel_idx] = true
                end
            end
            [bad_channels[channel_idx] == true && (bad_m[channel_idx, epoch_idx] = true) for channel_idx in 1:channel_n]
            (bad_channels_score / channel_n) > tc && push!(bad_epochs, epoch_idx)
        end
    end

    if :p2p in method
        @inbounds @simd for epoch_idx in 1:epoch_n
            bad_channels_score = 0
            bad_channels = zeros(Bool, channel_n)
            Threads.@threads for channel_idx in 1:channel_n
                s = Vector{Float64}()
                for length_idx in 1:w:(length(signal[channel[channel_idx], :, epoch_idx]) - w)
                    @views push!(s, median(signal[channel[channel_idx], length_idx:(length_idx + w), epoch_idx]))
                end
                p2p = @views round.(diff(s), digits=-2)
                s_m = @views mean(signal[channel[channel_idx], :, epoch_idx])
                s_s = @views std(signal[channel[channel_idx], :, epoch_idx])
                s_u = s_m + quantile.(Normal(), p) * s_s
                s_l = s_m - quantile.(Normal(), p) * s_s
                p2p_p = zeros(Bool, length(p2p))
                p2p_p[p2p .> s_u] .= true
                p2p_p[p2p .< s_l] .= true
                if sum(p2p_p) > 0
                    bad_channels_score += 1
                    bad_channels[channel_idx] = true
                end
            end
            for channel_idx in 1:channel_n
                bad_channels[channel_idx] == true && (bad_m[channel_idx, epoch_idx] = true)
            end 
            (bad_channels_score / channel_n) > tc && push!(bad_epochs, epoch_idx)
        end
    end

    return (bad_m=bad_m, bad_epochs=sort(unique(bad_epochs)))
end

"""
    eeg_add_labels(eeg, labels)

Add EEG channel labels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `labels::Vector{String}`

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_add_labels(eeg::NeuroAnalyzer.EEG, labels::Vector{String})

    length(labels) == eeg_channel_n(eeg) || throw(ArgumentError("labels length must be $(eeg_channel_n(eeg))."))
    
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_header[:labels] = labels

    push!(eeg_new.eeg_header[:history], "eeg_add_labels(EEG, labels=$labels")
 
    return eeg_new
end

"""
    eeg_add_labels!(eeg::NeuroAnalyzer.EEG, labels::Vector{String})

Add EEG channel labels.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `labels::Vector{String}`
"""
function eeg_add_labels!(eeg::NeuroAnalyzer.EEG, labels::Vector{String})

    length(labels) == eeg_channel_n(eeg) || throw(ArgumentError("labels length must be $(eeg_channel_n(eeg))."))
    eeg.eeg_header[:labels] = labels
    push!(eeg.eeg_header[:history], "eeg_add_labels(EEG, labels=$labels")

    return nothing
end

"""
    eeg_edit_channel(eeg; channel, field, value)

Edit EEG channel properties.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64`
- `field::Symbol`
- `value::Any`

# Returns

- `eeg_new::NeuroAnalyzer.EEG`
"""
function eeg_edit_channel(eeg::NeuroAnalyzer.EEG; channel::Int64, field::Symbol, value::Any)
    
    value === nothing && throw(ArgumentError("value cannot be empty."))
    _check_channels(eeg, channel)
    _check_var(field, [:channel_type, :labels], "field")    

    eeg_new = deepcopy(eeg)
    typeof(eeg_new.eeg_header[field][channel]) == typeof(value) || throw(ArgumentError("field type ($(eltype(eeg_new.eeg_header[field]))) does not mach value type ($(typeof(value)))."))
    eeg_new.eeg_header[field][channel] = value

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_edit_channel(EEG, channel=$channel, field=$field, value=$value)")   

    return eeg_new
end

"""
    eeg_edit_channel!(eeg; channel, field, value)

Edit EEG channel properties.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Int64`
- `field::Symbol`
- `value::Any`
"""
function eeg_edit_channel!(eeg::NeuroAnalyzer.EEG; channel::Int64, field::Symbol, value::Any)
    
    eeg_tmp = eeg_edit_channel(eeg, channel=channel, field=field, value=value)
    eeg.eeg_header = eeg_tmp.eeg_header

    push!(eeg.eeg_header[:history], "eeg_edit_channel(EEG, channel=$channel, field=$field, value=$value)")

    return nothing
end

"""
    eeg_keep_channel_type(eeg; type)

Keep EEG channels of `type` type.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `type::Symbol=:eeg`: type of channels to keep

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_keep_channel_type(eeg::NeuroAnalyzer.EEG; type::Symbol=:eeg)

    _check_var(type, [:all, :eeg, :meg, :ecg, :eog, :emg, :ref, :mrk], "type")

    channels_idx = Vector{Int64}()
    for idx in 1:eeg_channel_n(eeg, type=:all)
        eeg.eeg_header[:channel_type][idx] == string(type) && push!(channels_idx, idx)
    end
    eeg_new = eeg_keep_channel(eeg, channel=channels_idx)
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_keep_channel_type(EEG, type=$type")

    return eeg_new
end

"""
    eeg_keep_channel_type!(eeg; type)

Keep EEG channels of `type` type.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `type::Symbol=:eeg`: type of channels to keep
"""
function eeg_keep_channel_type!(eeg::NeuroAnalyzer.EEG; type::Symbol=:eeg)

    eeg_tmp = eeg_keep_channel_type(eeg, type=type)
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg.eeg_signals = eeg_tmp.eeg_signals

    push!(eeg.eeg_header[:history], "eeg_keep_channel_type!(EEG, type=$type")

    return nothing
end

"""
    eeg_view_note(eeg)

Return EEG note.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_view_note(eeg::NeuroAnalyzer.EEG)
    return eeg.eeg_header[:note]
end

"""
    eeg_copy(eeg)

Make copy of EEG.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `eeg_copy::NeuroAnalyzer.EEG`
"""
function eeg_copy(eeg::NeuroAnalyzer.EEG)
    return deepcopy(eeg)
end

"""
    eeg_epochs_time(eeg; ts)

Edit EEG epochs time start.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `ts::Real`: time start in seconds

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_epochs_time(eeg::NeuroAnalyzer.EEG; ts::Real)

    epoch_len = eeg_epoch_len(eeg)
    fs = eeg_sr(eeg)
    new_epochs_time = linspace(ts, ts + (epoch_len / fs), epoch_len)
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_epochs_time = new_epochs_time

    push!(eeg_new.eeg_header[:history], "eeg_epochs_time(EEG, ts=$ts)")

    return eeg_new
end

"""
    eeg_epochs_time!(eeg; ts)

Edit EEG epochs time start.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `ts::Real`: time start in seconds

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_epochs_time!(eeg::NeuroAnalyzer.EEG; ts::Real)

    eeg_tmp = eeg_epochs_time(eeg, ts=ts)
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg.eeg_epochs_time = eeg_tmp.eeg_epochs_time

    return nothing
end

"""
    eeg_add_note(eeg; note)

Add EEG note.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `note::String`

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_add_note(eeg::NeuroAnalyzer.EEG; note::String)
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_header[:note] = note
    return eeg_new
end

"""
    eeg_add_note!(eeg; note)

Add EEG note.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `note::String`
"""
function eeg_add_note!(eeg::NeuroAnalyzer.EEG; note::String)
    eeg.eeg_header[:note] = note
    return nothing
end

"""
    eeg_delete_note(eeg)

Delete EEG note.

# Arguments

- `eeg::NeuroAnalyzer.EEG`

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_delete_note(eeg::NeuroAnalyzer.EEG)
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_header[:note] = ""
    return eeg_new
end

"""
    eeg_delete_note!(eeg)

Delete EEG note.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_delete_note!(eeg::NeuroAnalyzer.EEG)
    eeg.eeg_header[:note] = ""
    return nothing
end

"""
    eeg_replace_channel(eeg; channel, signal)

Replace EEG channel.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, String}`: channel number or name
- `signal::Array{Float64, 3}

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_replace_channel(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, String}, signal::Array{Float64, 3})

    channel_idx = nothing
    labels = eeg_labels(eeg)
    if typeof(channel) == String
        for idx in eachindex(labels)
            if channel == labels[idx]
                channel_idx = idx
            end
        end
        channel_idx === nothing && throw(ArgumentError("Channel name ($channel) does not match signal labels."))
    else
        _check_channels(eeg, channel)
        channel_idx = channel
    end

    eeg_new = deepcopy(eeg)
    size(signal) == (1, eeg_epoch_len(eeg_new), eeg_epoch_n(eeg_new)) || throw(ArgumentError("signal size ($(size(signal))) must be the same as EEG channel size ($(size(eeg_new.eeg_signals[channel_idx, :, :]))."))
    eeg_new.eeg_signals[channel_idx, :, :] = signal
    eeg_reset_components!(eeg_new)

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_replace_channel(EEG, channel=$channel, signal=$signal")

    return eeg_new
end

"""
    eeg_replace_channel!(eeg; channel, signal)

Replace EEG channel.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, String}`: channel number or name
- `signal::Array{Float64, 3}
"""
function eeg_replace_channel!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, String}, signal::Array{Float64, 3})

    eeg_tmp = eeg_replace_channel(eeg, channel=chanel, signal=signal)
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg.eeg_signals = eeg_tmp.eeg_signals

    eeg_reset_components!(eeg)
    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_replace_channel(EEG, channel=$channel, signal=$signal")

    return nothing
end

"""
    eeg_interpolate_channel(eeg; channel, epoch, m, q)

Interpolate EEG channel(s) using planar interpolation.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to interpolate
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) within to interpolate
- `imethod::Symbol=:sh`: interpolation method Shepard (`:sh`), Multiquadratic (`:mq`), InverseMultiquadratic (`:imq`), ThinPlate (`:tp`), NearestNeighbour (`:nn`), Gaussian (`:ga`)
- `interpolation_factor::Int64=100`: interpolation quality

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_interpolate_channel(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange}, epoch::Union{Int64, Vector{Int64}, AbstractRange}, imethod::Symbol=:sh, interpolation_factor::Int64=100)

    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")
    eeg.eeg_header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available, use eeg_load_electrodes() or eeg_add_electrodes() first."))

    typeof(channel) == Vector{Int64} && sort!(channel, rev=true)

    eeg_new = deepcopy(eeg)
    _check_channels(eeg, channel)
    _check_epochs(eeg, epoch)

    loc_x1 = eeg.eeg_locs[!, :loc_x]
    loc_y1 = eeg.eeg_locs[!, :loc_y]
    
    eeg_delete_channel!(eeg, channel=channel)
    loc_x2 = eeg.eeg_locs[!, :loc_x]
    loc_y2 = eeg.eeg_locs[!, :loc_y]
    channels = eeg_channel_idx(eeg, type=Symbol(eeg.eeg_header[:signal_type]))

    epoch_n = length(epoch)
    epoch_len = eeg_epoch_len(eeg)

    s_interpolated = zeros(Float64, length(channel), epoch_len, epoch_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(epoch_n * epoch_len, 1))

    Threads.@threads for epoch_idx in eachindex(epoch)
        @inbounds @simd for length_idx in 1:epoch_len
            s_tmp, x, y = @views _interpolate(eeg.eeg_signals[channels, length_idx, epoch[epoch_idx]], loc_x2, loc_y2, interpolation_factor, imethod, :none)
            for channel_idx in eachindex(channel)
                x_idx = vsearch(loc_x1[channel[channel_idx]], x)
                y_idx = vsearch(loc_y1[channel[channel_idx]], y)
                s_interpolated[channel_idx, length_idx, epoch_idx] = s_tmp[x_idx, y_idx]
            end

            # update progress bar
            progress_bar == true && next!(p)
        end
    end

    eeg_new.eeg_signals[channel, :, epoch] = s_interpolated

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_interpolate_channel(EEG, channel=$channel, epoch=$epoch, imethod=$imethod, interpolation_factor=$interpolation_factor)")

    return eeg_new
end

"""
    eeg_interpolate_channel!(eeg; channel, epoch, imethod, interpolation_factor)

Interpolate EEG channel(s) using planar interpolation.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, Vector{Int64}, Abstractrange}`: channel number(s) to interpolate
- `epoch::Union{Int64, Vector{Int64}, Abstractrange}`: epoch number(s) within to interpolate
- `imethod::Symbol=:sh`: interpolation method Shepard (`:sh`), Multiquadratic (`:mq`), InverseMultiquadratic (`:imq`), ThinPlate (`:tp`), NearestNeighbour (`:nn`), Gaussian (`:ga`)
- `interpolation_factor::Int64=100`: interpolation quality
"""
function eeg_interpolate_channel!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, Vector{Int64}}, m::Symbol=:shepard, interpolation_factor::Int64=100)

    eeg.eeg_signals = eeg_interpolate_channel(eeg, channel=channel, epoch=epoch, imethod=imethod, interpolation_factor=interpolation_factor).eeg_signals
    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_interpolate_channel!(EEG, channel=$channel, epoch=$epoch, imethod=$imethod, interpolation_factor=$interpolation_factor)")

    return nothing
end

"""
    loc_flipy(locs; planar, spherical)

Flip channel locations along y axis.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function loc_flipy(locs::DataFrame; planar::Bool=true, spherical::Bool=true)

    locs_new = deepcopy(locs)

    for idx in eachindex(locs[!, :labels])
        if planar == true
            t = locs_new[!, :loc_theta][idx]
            q = _angle_quadrant(t)
            q == 1 && (t = 360 - t)
            q == 2 && (t = 180 + (180 - t))
            q == 3 && (t = 180 - (t - 180))
            q == 4 && (t = 360 - t)
            t = mod(t, 360)
            locs_new[!, :loc_theta][idx] = t
        end
        spherical == true && (locs_new[!, :loc_y][idx] = -locs_new[!, :loc_y][idx])
    end
    loc_cart2sph!(locs_new)

    return locs_new
end

"""
    loc_flipy!(locs; planar, spherical)

Flip channel locations along y axis.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates
"""
function loc_flipy!(locs::DataFrame; planar::Bool=true, spherical::Bool=true)

    for idx in eachindex(locs[!, :labels])
        if planar == true
            t = locs[!, :loc_theta][idx]
            q = _angle_quadrant(t)
            q == 1 && (t = 360 - t)
            q == 2 && (t = 180 + (180 - t))
            q == 3 && (t = 180 - (t - 180))
            q == 4 && (t = 360 - t)
            t = mod(t, 360)
            locs[!, :loc_theta][idx] = t
        end
        spherical == true && (locs[!, :loc_y][idx] = -locs[!, :loc_y][idx])
    end
    loc_cart2sph!(locs)

    return nothing
end

"""
    loc_flipx(locs; planar, spherical)

Flip channel locations along x axis.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function loc_flipx(locs::DataFrame; planar::Bool=true, spherical::Bool=true)

    locs_new = deepcopy(locs)

    for idx in eachindex(locs[!, :labels])
        if planar == true
            t = locs_new[!, :loc_theta][idx]
            q = _angle_quadrant(t)
            q == 1 && (t = 90 + (90 - t))
            q == 2 && (t = 90 - (t - 90))
            q == 3 && (t = 270 + (270 - t))
            q == 4 && (t = 270 - (t - 270))
            t = mod(t, 360)
            locs_new[!, :loc_theta][idx] = t
        end
        spherical == true && (locs_new[!, :loc_x][idx] = -locs_new[!, :loc_x][idx])
    end
    loc_cart2sph!(locs_new)

    return locs_new
end

"""
    loc_flipx!(locs; planar, spherical)

Flip channel locations along x axis.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates
"""
function loc_flipx!(locs::DataFrame; planar::Bool=true, spherical::Bool=true)

    for idx in eachindex(locs[!, :labels])
        if planar == true
            t = locs[!, :loc_theta][idx]
            q = _angle_quadrant(t)
            q == 1 && (t = 90 + (90 - t))
            q == 2 && (t = 90 - (t - 90))
            q == 3 && (t = 270 + (270 - t))
            q == 4 && (t = 270 - (t - 270))
            t = mod(t, 360)
            locs[!, :loc_theta][idx] = t
        end
        spherical == true && (locs[!, :loc_x][idx] = -locs[!, :loc_x][idx])
    end
    loc_cart2sph!(locs)

    return nothing
end

"""
    loc_flipz(locs)

Flip channel locations along z axis.

# Arguments

- `locs::DataFrame`

# Returns

- `locs_new::DataFrame`
"""
function loc_flipz(locs::DataFrame)

    locs_new = deepcopy(locs)

    for idx in eachindex(locs[!, :labels])
        locs_new[!, :loc_z][idx] = -locs_new[!, :loc_z][idx]
    end
    loc_cart2sph!(locs_new)

    return locs_new
end

"""
    loc_flipz!(eeg)

Flip channel locations along z axis.

# Arguments

- `locs::DataFrame`
"""
function loc_flipz!(locs::DataFrame)

    for idx in eachindex(locs[!, :labels])
        locs[!, :loc_z][idx] = -locs[!, :loc_z][idx]
    end
    loc_cart2sph!(locs)

    return nothing
end

"""
    eeg_channel_type(eeg; channel, type)

Change EEG channel type.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, String}`
- `type::String`

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_channel_type(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, String}, type::String)

    type = lowercase(type)
    labels = eeg_labels(eeg)

    # create new dataset
    eeg_new = deepcopy(eeg)
    types = eeg_new.eeg_header[:channel_type]
    
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
        _check_channels(eeg, channel)
        types[channel] = type
    end
    eeg_new.eeg_header[:channel_type] = types
    
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_channel_type(EEG, channel=$channel, type=$type)")

    return eeg_new
end

"""
    eeg_channel_type!(eeg; channel, new_name)

Change EEG channel type.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, String}`
- `type::String`
"""
function eeg_channel_type!(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, String}, type::String)

    eeg.eeg_header[:channel_type] = eeg_channel_idx(eeg, channel=channel, type=type)[:channel_type]
    push!(eeg.eeg_header[:history], "eeg_channel_type(EEG, channel=$channel, type=$type)")

    return nothing
end

"""
    eeg_edit_electrode(eeg; <keyword arguments>)

Edit EEG electrode.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{String, Int64}`: channel number or name
- `x::Union{Real, Nothing}=nothing`: Cartesian X spherical coordinate
- `y::Union{Real, Nothing}=nothing`: Cartesian Y spherical coordinate
- `z::Union{Real, Nothing}=nothing`: Cartesian Z spherical coordinate
- `theta::Union{Real, Nothing}=nothing`: polar planar theta coordinate
- `radius::Union{Real, Nothing}=nothing`: polar planar radius coordinate
- `theta_sph::Union{Real, Nothing}=nothing`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
- `radius_sph::Union{Real, Nothing}=nothing`: spherical radius, the distance from the origin to the point
- `phi_sph::Union{Real, Nothing}=nothing`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees
- `name::String=""`: channel name
- `type::String=""`: channel type

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_edit_electrode(eeg::NeuroAnalyzer.EEG; channel::Union{String, Int64}, x::Union{Real, Nothing}=nothing, y::Union{Real, Nothing}=nothing, z::Union{Real, Nothing}=nothing, theta::Union{Real, Nothing}=nothing, radius::Union{Real, Nothing}=nothing, theta_sph::Union{Real, Nothing}=nothing, radius_sph::Union{Real, Nothing}=nothing, phi_sph::Union{Real, Nothing}=nothing, name::String="", type::String="")

    eeg_new = deepcopy(eeg)
    channel = _get_channel_idx(eeg_labels(eeg_new), channel)

    name != "" && eeg_rename_channel!(eeg_new, channel=channel, name=name)
    type != "" && eeg_channel_type!(eeg_new, channel=channel, type=type)

    x !== nothing && (eeg_new.eeg_locs[!, :loc_x][channel] = x)
    y !== nothing && (eeg_new.eeg_locs[!, :loc_y][channel] = y)
    z !== nothing && (eeg_new.eeg_locs[!, :loc_z][channel] = z)
    theta !== nothing && (eeg_new.eeg_locs[!, :loc_theta][channel] = theta)
    radius !== nothing && (eeg_new.eeg_locs[!, :loc_radius][channel] = radius)
    theta_sph !== nothing && (eeg_new.eeg_locs[!, :loc_theta_sph][channel] = theta_sph)
    radius_sph !== nothing && (eeg_new.eeg_locs[!, :loc_radius_sph][channel] = radius_sph)
    phi_sph !== nothing && (eeg_new.eeg_locs[!, :loc_phi_sph][channel] = phi_sph)

    (x !== nothing || y !== nothing || z !== nothing || theta !== nothing || radius !== nothing || theta_sph !== nothing  || radius_sph !== nothing || phi_sph !== nothing) && (eeg_new.eeg_header[:channel_locations] == true)

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_edit_electrode(EEG; channel=$channel, x=$x, y=$y, z=$z, theta=$theta, radius=$radius, theta_sph=$theta_sph, radius_sph=$radius_sph, phi_sph=$phi_sph, name=$name, type=$type)")

    return eeg_new
end

"""
    eeg_edit_electrode!(eeg; <keyword arguments>)

Edit EEG electrode.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{String, Int64}`: channel number or name
- `x::Union{Real, Nothing}=nothing`: Cartesian X spherical coordinate
- `y::Union{Real, Nothing}=nothing`: Cartesian Y spherical coordinate
- `z::Union{Real, Nothing}=nothing`: Cartesian Z spherical coordinate
- `theta::Union{Real, Nothing}=nothing`: polar planar theta coordinate
- `radius::Union{Real, Nothing}=nothing`: polar planar radius coordinate
- `theta_sph::Union{Real, Nothing}=nothing`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
- `radius_sph::Union{Real, Nothing}=nothing`: spherical radius, the distance from the origin to the point
- `phi_sph::Union{Real, Nothing}=nothing`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees
- `name::String=""`: channel name
- `type::String=""`: channel type
"""
function eeg_edit_electrode!(eeg::NeuroAnalyzer.EEG; channel::Union{String, Int64}, x::Union{Real, Nothing}=nothing, y::Union{Real, Nothing}=nothing, z::Union{Real, Nothing}=nothing, theta::Union{Real, Nothing}=nothing, radius::Union{Real, Nothing}=nothing, theta_sph::Union{Real, Nothing}=nothing, radius_sph::Union{Real, Nothing}=nothing, phi_sph::Union{Real, Nothing}=nothing, name::String="", type::String="")

    eeg_tmp = eeg_edit_electrode(eeg, channel=channel, x=x, y=y, z=z, theta=theta, radius=radius, theta_sph=theta_sph, radius_sph=radius_sph, phi_sph=phi_sph, name=name, type=type)
    eeg.eeg_header = eeg_tmp.eeg_header
    eeg.eeg_locs = eeg_tmp.eeg_locs

    eeg_reset_components!(eeg)
    return nothing
end

"""
    eeg_electrode_loc(eeg; channel, output)

Return locations of EEG channel electrode.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `channel::Union{Int64, String}`: channel number or name
- `output::Bool=true`: print output if true

# Returns

Named tuple containing:
- `theta::Float64`: polar planar theta coordinate
- `radius::Float64`: polar planar radius coordinate
- `x::Float64`: Cartesian X spherical coordinate
- `y::Float64`: Cartesian Y spherical coordinate
- `z::Float64`: Cartesian Z spherical coordinate
- `theta_sph::Float64`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
- `radius_sph::Float64`: spherical radius, the distance from the origin to the point
- `phi_sph::Float64`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees
"""
function eeg_electrode_loc(eeg::NeuroAnalyzer.EEG; channel::Union{Int64, String}, output::Bool=true)

    eeg.eeg_header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available, use eeg_load_electrodes() or eeg_add_electrodes() first."))

    channel = _get_channel_idx(eeg_labels(eeg), channel)

    x = eeg.eeg_locs[!, :loc_x][channel]
    y = eeg.eeg_locs[!, :loc_y][channel]
    z = eeg.eeg_locs[!, :loc_z][channel]
    theta = eeg.eeg_locs[!, :loc_theta][channel]
    radius = eeg.eeg_locs[!, :loc_radius][channel]
    theta_sph = eeg.eeg_locs[!, :loc_theta_sph][channel]
    radius_sph = eeg.eeg_locs[!, :loc_radius_sph][channel]
    phi_sph = eeg.eeg_locs[!, :loc_phi_sph][channel]

    if output
        println("Channel: $channel")
        println("  Label: $(eeg_labels(eeg)[channel])")
        println("  theta: $theta (planar)")
        println(" radius: $radius (planar)")
        println("      X: $x (spherical)")
        println("      Y: $y (spherical)")
        println("      Z: $z (spherical)")
        println(" radius: $radius_sph (spherical)")
        println("  theta: $theta_sph (spherical)")
        println("    phi: $phi_sph (spherical)")
    end
    
    return (theta=theta, radius=radius, x=x, y=y, z=z, theta_sph=theta_sph, radius_sph=radius_sph, phi_sph=phi_sph)
end

"""
    loc_swapxy(locs; planar, spherical)

Swap channel locations x and y axes.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function loc_swapxy(locs::DataFrame; planar::Bool=true, spherical::Bool=true)

    locs_new = deepcopy(locs)

    for idx in eachindex(locs[!, :labels])
        if planar == true
            t = deg2rad(locs_new[!, :loc_theta][idx])
            t += pi / 2
            locs_new[!, :loc_theta][idx] = rad2deg(t)
        end
        if spherical == true
            locs_new[!, :loc_x][idx], locs_new[!, :loc_y][idx] = locs_new[!, :loc_y][idx], locs_new[!, :loc_x][idx]
        end
    end
    loc_cart2sph!(locs_new)

    return locs_new
end

"""
    loc_swapxy!(locs; planar, spherical)

Swap channel locations x and y axes.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates
"""
function loc_swapxy!(locs::DataFrame; planar::Bool=true, spherical::Bool=true)

    for idx in eachindex(locs[!, :labels])
        if planar == true
            t = deg2rad(locs[!, :loc_theta][idx])
            t += pi / 2
            locs[!, :loc_theta][idx] = rad2deg(t)
        end
        if spherical == true
            locs[!, :loc_x][idx], locs[!, :loc_y][idx] = locs[!, :loc_y][idx], locs[!, :loc_x][idx]
        end
    end
    loc_cart2sph!(locs)

    return nothing
end

"""
    loc_sph2cart(locs)

Convert spherical locations to Cartesian.

# Arguments

- `locs::DataFrame`

# Returns

- `locs_new::DataFrame`
"""
function loc_sph2cart(locs::DataFrame)

    locs_new = deepcopy(locs)

    for idx in eachindex(locs[!, :labels])
        r = locs_new[!, :loc_radius_sph][idx]
        t = locs_new[!, :loc_theta_sph][idx]
        p = locs_new[!, :loc_phi_sph][idx]
        x, y, z = sph2cart(r, t, p)
        locs_new[!, :loc_x][idx] = x
        locs_new[!, :loc_y][idx] = y
        locs_new[!, :loc_z][idx] = z
    end

    return locs_new
end

"""
    loc_sph2cart!(locs)

Convert spherical locations to Cartesian.

# Arguments

- `locs::DataFrame`
"""
function loc_sph2cart!(locs::DataFrame)

    for idx in eachindex(locs[!, :labels])
        r = locs[!, :loc_radius_sph][idx]
        t = locs[!, :loc_theta_sph][idx]
        p = locs[!, :loc_phi_sph][idx]
        x, y, z = sph2cart(r, t, p)
        locs[!, :loc_x][idx] = x
        locs[!, :loc_y][idx] = y
        locs[!, :loc_z][idx] = z
    end

    return nothing
end

"""
    loc_cart2sph(locs)

Convert Cartesian locations to spherical.

# Arguments

- `locs::DataFrame`

# Returns

- `locs_new::DataFrame`
"""
function loc_cart2sph(locs::DataFrame)

    locs_new = deepcopy(locs)

    for idx in eachindex(locs[!, :labels])
        x = locs_new[!, :loc_x][idx]
        y = locs_new[!, :loc_y][idx]
        z = locs_new[!, :loc_z][idx]
        r, t, p = round.(cart2sph(x, y, z), digits=2)
        locs_new[!, :loc_radius_sph][idx] = r
        locs_new[!, :loc_theta_sph][idx] = t
        locs_new[!, :loc_phi_sph][idx] = p
    end

    return locs_new
end

"""
    loc_cart2sph!(locs)

Convert Cartesian locations to spherical.

# Arguments

- `locs::DataFrame`
"""
function loc_cart2sph!(locs::DataFrame)

    for idx in eachindex(locs[!, :labels])
        x = locs[!, :loc_x][idx]
        y = locs[!, :loc_y][idx]
        z = locs[!, :loc_z][idx]
        r, t, p = cart2sph(x, y, z)
        locs[!, :loc_radius_sph][idx] = r
        locs[!, :loc_theta_sph][idx] = t
        locs[!, :loc_phi_sph][idx] = p
    end

    return nothing
end

"""
    eeg_view_markers(eeg)

Show markers.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_view_markers(eeg::NeuroAnalyzer.EEG)
    eeg.eeg_header[:markers] == true || throw(ArgumentError("EEG has no markers."))
    for marker_idx in 1:size(eeg.eeg_markers, 1)
        println("ID: $(rpad(("'" * eeg.eeg_markers[!, :id][marker_idx] * "'"), 24, " ")) start [sample]: $(rpad(eeg.eeg_markers[!, :start][marker_idx], 8, " ")) length [samples]: $(rpad(eeg.eeg_markers[!, :length][marker_idx], 8, " ")) description: $(rpad(("'" * eeg.eeg_markers[!, :description][marker_idx] * "'"), 24, " ")) channel: $(eeg.eeg_markers[!, :channel][marker_idx])")
    end
end

"""
    eeg_delete_marker(eeg; n)

Delete marker.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `n::Int64`: marker number

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_delete_marker(eeg::NeuroAnalyzer.EEG; n::Int64)
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_header[:markers] == true || throw(ArgumentError("EEG has no markers."))
    nn = size(eeg_new.eeg_markers, 1)
    (n < 1 || n > nn) && throw(ArgumentError("n has to be ≥ 1 and ≤ $nn."))
    deleteat!(eeg_new.eeg_markers, n)
    size(eeg_new.eeg_markers, 1) == 0 && (eeg_new.eeg_header[:markers] = false)
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_delete_marker(EEG; n=$n)")
    
    return eeg_new
end

"""
    eeg_delete_marker!(eeg; n)

Delete marker.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `n::Int64`: marker number
"""
function eeg_delete_marker!(eeg::NeuroAnalyzer.EEG; n::Int64)

    eeg_tmp = eeg_delete_marker(eeg, n=n)
    eeg.eeg_header[:markers] = eeg_tmp.eeg_header[:markers]
    eeg.eeg_markers = eeg_tmp.eeg_markers

    eeg_reset_components!(eeg)
    return nothing
end

"""
    eeg_add_marker(eeg; id, start, len, desc)

Add marker.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `id::String`: marker ID
- `start::Int64`: marker time in samples
- `len::Int64`: marker length in samples
- `desc::String`: marker description
- `channel::Int64`: channel number, if 0 then marker is related to all channels

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_add_marker(eeg::NeuroAnalyzer.EEG; id::String, start::Int64, len::Int64, desc::String, channel::Int64)

    start < 1 && throw(ArgumentError("start must be > 0."))
    len < 1 && throw(ArgumentError("len must be > 0."))
    start > eeg_epoch_len(eeg) && throw(ArgumentError("start must be ≤ $(eeg_epoch_len(eeg))."))
    len > eeg_epoch_len(eeg) && throw(ArgumentError("len must be ≤ $(eeg_epoch_len(eeg))."))
    start + len > eeg_epoch_len(eeg) + 1 && throw(ArgumentError("start + len must be ≤ $(eeg_epoch_len(eeg) + 1)."))

    eeg_new = deepcopy(eeg)
    eeg_new.eeg_header[:markers] = true
    append!(eeg_new.eeg_markers, DataFrame(:id => id, :start => start, :length => len, :description => desc, :channel => channel))
    sort!(eeg_new.eeg_markers)
    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_add_marker(EEG; id=$id, start=$start, len=$len, desc=$desc, channel=$channel)")

    return eeg_new
end

"""
    eeg_add_marker!(eeg; id, start, len, desc)

Add marker.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `id::String`: marker ID
- `start::Int64`: marker time in samples
- `len::Int64`: marker length in samples
- `desc::String`: marker description
- `channel::Int64`: channel number, if 0 then marker is related to all channels
"""
function eeg_add_marker!(eeg::NeuroAnalyzer.EEG; id::String, start::Int64, len::Int64, desc::String, channel::Int64)

    eeg_tmp = eeg_add_marker(eeg, id=id, start=start, len=len, desc=desc, channel=channel)
    eeg.eeg_header[:markers] = eeg_tmp.eeg_header[:markers]
    eeg.eeg_markers = eeg_tmp.eeg_markers

    eeg_reset_components!(eeg)
    return nothing
end

"""
    eeg_channel_idx(eeg; type=:eeg)

Return index of EEG channels of `type` type.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `type::Vector{Symbol}=:all`: channel type

# Returns

- `channel_n::Int64`
"""
function eeg_channel_idx(eeg::NeuroAnalyzer.EEG; type::Symbol=:all)

    _check_var(type, [:all, :eeg, :meg, :ecg, :eog, :emg, :ref, :mrk], "type")
    channel_idx = Vector{Int64}()
    for idx in 1:eeg_channel_n(eeg)
        eeg.eeg_header[:channel_type][idx] == string(type) && (push!(channel_idx, idx))
    end

    return channel_idx
end

"""
    eeg_vch(eeg; f)

Calculate virtual channel using formula `f`.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `f::String`: channel calculation formula, e.g. `"cz / mean(fp1 + fp2)"`; case of labels in the formula is ignored, all standard Julia math operators are available, channel labels must be the same as of the EEG object

# Returns
 
- `vc::Array{Float64, 3}`: single channel × time × epochs
"""
function eeg_vch(eeg::NeuroAnalyzer.EEG; f::String)

    epoch_n = eeg_epoch_n(eeg)
    f = lowercase(f)
    labels = lowercase.(eeg_labels(eeg))
    vc = zeros(1, eeg_epoch_len(eeg), epoch_n)
    Threads.@threads for epoch_idx in 1:epoch_n
        f_tmp = f
        for channel_idx in eachindex(labels)
            occursin(labels[channel_idx], f) == true && (f_tmp = replace(f_tmp, labels[channel_idx] => "$(eeg.eeg_signals[channel_idx, :, epoch_idx])"))
        end
        try
            @inbounds vc[1, :, epoch_idx] = eval(Meta.parse("@. " * f_tmp))
        catch
            @error "Formula is incorrect, check channel labels and operators."
        end
    end

    return vc
end

"""
    eeg_edit_marker(eeg; n, id, start, len, desc)

Edit EEG marker.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `n::Int64`: marker number
- `id::String`: marker ID
- `start::Int64`: marker time in samples
- `len::Int64`: marker length in samples
- `desc::String`: marker description
- `channel::Int64`: channel number, if 0 then marker is related to all channels

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_edit_marker(eeg::NeuroAnalyzer.EEG; n::Int64, id::String, start::Int64, len::Int64, desc::String, channel::Int64)

    eeg.eeg_header[:markers] == true || throw(ArgumentError("EEG has no markers."))
    start < 1 && throw(ArgumentError("start must be > 0."))
    len < 1 && throw(ArgumentError("len must be > 0."))
    start > eeg_epoch_len(eeg) && throw(ArgumentError("start must be ≤ $(eeg_epoch_len(eeg))."))
    len > eeg_epoch_len(eeg) && throw(ArgumentError("len must be ≤ $(eeg_epoch_len(eeg))."))
    start + len + 1 > eeg_epoch_len(eeg) && throw(ArgumentError("start + len must be ≤ $(eeg_epoch_len(eeg) + 1)."))

    nn = size(eeg.eeg_markers, 1)
    n < 1 || n > nn && throw(ArgumentError("n has to be ≥ 1 and ≤ $nn."))
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_markers[n, :] = Dict(:id => id, :start => start, :length => len, :description => desc, :channel => channel)
     eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_edit_marker(EEG; id=$id, start=$start, len=$len, desc=$desc, channel=$channel)")

    return eeg_new
end

"""
    eeg_edit_marker!(eeg; n, id, start, len, desc)

Edit EEG marker.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `n::Int64`: marker number
- `id::String`: marker ID
- `start::Int64`: marker time in samples
- `len::Int64`: marker length in samples
- `desc::String`: marker description
- `channel::Int64`: channel number, if 0 then marker is related to all channels
"""
function eeg_edit_marker!(eeg::NeuroAnalyzer.EEG; n::Int64, id::String, start::Int64, len::Int64, desc::String, channel::Int64)

    eeg_tmp = eeg_edit_marker(eeg, n=n, id=id, start=start, len=len, desc=desc, channel=channel)
    eeg.eeg_header[:markers] = eeg_tmp.eeg_header[:markers]
    eeg.eeg_markers = eeg_tmp.eeg_markers

    eeg_reset_components!(eeg)
    return nothing
end