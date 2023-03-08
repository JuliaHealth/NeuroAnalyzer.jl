"""
    eeg_delete_channel(eeg; channel)

Delete EEG channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to be removed

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_delete_channel(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange})

    typeof(channel) <: AbstractRange && (channel = collect(channel))
    ch_n = channel_n(eeg)
    length(channel) > 1 && (channel = sort!(channel, rev=true))
    length(channel) == ch_n && throw(ArgumentError("You cannot delete all channels."))

    _check_channels(eeg, channel)

    eeg_new = deepcopy(eeg)

    # update headers
    for idx in channel
        loc = findfirst(isequal(lowercase(eeg_new.header.recording[:labels][idx])), lowercase.(string.(eeg_new.eeg_locs[!, :labels])))
        loc !== nothing && deleteat!(eeg_new.eeg_locs, loc)
        deleteat!(eeg_new.header.recording[:labels], idx)
        deleteat!(eeg_new.header.recording[:channel_type], idx)
        deleteat!(eeg_new.header.recording[:transducers], idx)
        deleteat!(eeg_new.header.recording[:units], idx)
        deleteat!(eeg_new.header.recording[:prefiltering], idx)
        deleteat!(eeg_new.header.recording[:gain], idx)
    end
    eeg_new.header.recording[:ch_n] -= length(channel)

    # remove channel
    eeg_new.eeg_signals = eeg_new.eeg_signals[setdiff(1:end, (channel)), :, :]

    eeg_reset_components!(eeg_new)
    push!(eeg_new.header.recording[:history], "eeg_delete_channel(EEG, $channel)")

    return eeg_new
end

"""
    eeg_delete_channel!(eeg; channel)

Delete EEG channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to be removed
"""
function eeg_delete_channel!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange})

    typeof(channel) <: AbstractRange && (channel = collect(channel))
    ch_n = channel_n(eeg)
    length(channel) > 1 && (channel = sort!(channel, rev=true))
    length(channel) == ch_n && throw(ArgumentError("You cannot delete all channels."))

    _check_channels(eeg, channel)

    # update headers
    for idx in channel
        loc = findfirst(isequal(lowercase(obj.header[:labels][idx])), lowercase.(string.(obj.locs[!, :labels])))
        loc !== nothing && deleteat!(obj.locs, loc)
        deleteat!(obj.header[:labels], idx)
        deleteat!(obj.header[:channel_type], idx)
        deleteat!(obj.header[:transducers], idx)
        deleteat!(obj.header[:units], idx)
        deleteat!(obj.header[:prefiltering], idx)
        deleteat!(obj.header[:gain], idx)
    end
    obj.header[:ch_n] -= length(channel)

    # remove channel
    obj.data = obj.data[setdiff(1:end, (channel)), :, :]

    eeg_reset_components!(eeg)
    push!(obj.header[:history], "eeg_delete_channel!(EEG, $channel)")

    return nothing
end

"""
    eeg_keep_channel(eeg; channel)

Keep EEG channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to keep

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_keep_channel(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange})

    typeof(channel) <: AbstractRange && (channel = collect(channel))
    _check_channels(eeg, channel)

    ch_n = channel_n(eeg)
    channels_to_remove = setdiff(collect(1:ch_n), channel)
    length(channels_to_remove) == ch_n && throw(ArgumentError("You cannot delete all channels."))

    return eeg_delete_channel(eeg, channel=channels_to_remove)
end

"""
    eeg_keep_channel!(eeg; channel)

Keep EEG channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to keep
"""
function eeg_keep_channel!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange})

    typeof(channel) <: AbstractRange && (channel = collect(channel))
    _check_channels(eeg, channel)

    ch_n = channel_n(eeg)
    channels_to_remove = setdiff(collect(1:ch_n), channel)
    length(channels_to_remove) == ch_n && throw(ArgumentError("You cannot delete all channels."))

    eeg_delete_channel!(eeg, channel=channels_to_remove)
end

"""
    eeg_get_channel(eeg; channel)

Return EEG channel number (if provided by name) or name (if provided by number).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, String}`: channel number or name

# Returns

- `channel_idx::Union{Int64, String}`: channel number or name
"""
function eeg_get_channel(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, String})

    labels = eeg_labels(eeg)
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
        _check_channels(eeg, channel)
        return labels[channel]
    end
end

"""
    eeg_rename_channel(eeg; channel, name)

Rename EEG channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, String}`: channel number or name
- `name::String`: new name

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_rename_channel(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, String}, name::String)

    # create new dataset
    eeg_new = deepcopy(eeg)
    labels = eeg_labels(eeg_new)
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
        _check_channels(eeg, channel)
        labels[channel] = name
    end
    eeg_new.header.recording[:labels] = labels
    
    # add entry to :history field
    push!(eeg_new.header.recording[:history], "eeg_rename_channel(EEG, channel=$channel, name=$name)")

    return eeg_new
end

"""
    eeg_rename_channel!(eeg; channel, name)

Rename EEG channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, String}`: channel number or name
- `name::String`: new name
"""
function eeg_rename_channel!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, String}, name::String)

    obj.header[:labels] = eeg_rename_channel(eeg, channel=channel, name=name).header.recording[:labels]
    push!(obj.header[:history], "eeg_rename_channel!(EEG, channel=$channel, name=$name)")

    return nothing
end

"""
    eeg_extract_channel(eeg; channel)

Extract EEG channel data.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, String}`: channel number or name

# Returns

- `channel::Vector{Float64}`
"""
function eeg_extract_channel(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, String})

    labels = eeg_labels(eeg)
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
        _check_channels(eeg, channel)
        channel_idx = channel
    end    
    eeg_channel = reshape(obj.data[channel_idx, :, :], 1, epoch_len(eeg), eeg_epoch_n(eeg))

    return eeg_channel
end

"""
    eeg_epoch(eeg; marker, ep_offset, ep_n, epoch_len)

Split EEG into epochs. Return signal that is split either by markers (if specified), by epoch length or by number of epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `marker::String="": marker name to split at
- `ep_offset::Int64=0": time offset (in samples) for marker-based epoching (each epoch time will start at marker time - ep_offset)
- `ep_n::Union{Int64, Nothing}=nothing`: number of epochs
- `epoch_len::Union{Int64, Nothing}`=nothing: epoch length in samples

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_epoch(obj::NeuroAnalyzer.NEURO; marker::String="", ep_offset::Real=0, ep_n::Union{Int64, Nothing}=nothing, epoch_len::Union{Int64, Nothing}=nothing)

    eeg_new = deepcopy(eeg)

    if marker != ""
        # split by markers
        if obj.header.markers == true
            epoch_len === nothing && throw(ArgumentError("epoch_len must be specified."))
            ep_offset == 0 && throw(ArgumentError("ep_offset must be specified."))
            _check_markers(eeg, marker)

            # get marker positions
            marker_idx = []
            for idx in 1:length(obj.markers[!, :description])
                eeg_new.eeg_markers[idx, :description] == marker && push!(marker_idx, idx)
            end
            marker_start = eeg_new.eeg_markers[!, :start][marker_idx]

            # split into epochs
            epochs, eeg_new.eeg_markers = _make_epochs_bymarkers(obj.data, markers=eeg_new.eeg_markers, marker_start=marker_start, ep_offset=ep_offset, epoch_len=epoch_len)
        else
            throw(ArgumentError("EEG does not contain markers."))
        end
    else
        # split by epoch_len or ep_n
        epochs = _make_epochs(obj.data, ep_n=ep_n, epoch_len=epoch_len)

        # delete markers outside epochs
        for marker_idx in nrow(eeg_new.eeg_markers):-1:1
            eeg_new.eeg_markers[marker_idx, :start] in 1:size(epochs, 2) * size(epochs, 3) || deleteat!(eeg_new.eeg_markers, marker_idx)
        end
    end

    # create new dataset
    ep_n = size(epochs, 3)
    epoch_duration_samples = size(epochs, 2)
    epoch_duration_seconds = size(epochs, 2) / obj.header[:sampling_rate]
    eeg_duration_samples = size(epochs, 2) * size(epochs, 3)
    eeg_duration_seconds = eeg_duration_samples / obj.header[:sampling_rate]
    eeg_time = collect(0:(1 / obj.header[:sampling_rate]):eeg_duration_seconds)
    eeg_time = eeg_time[1:(end - 1)]

    # update signal
    eeg_new.eeg_signals = epochs

    # update time
    eeg_new.eeg_time = eeg_time

    # update epochs time
    fs = eeg_sr(eeg_new)
    new_epochs_time = linspace(-s2t(ep_offset, fs), epoch_duration_seconds - s2t(ep_offset, fs), epoch_duration_samples)
    eeg_new.eeg_epoch_time = new_epochs_time

    # update header
    eeg_new.header.recording[:eeg_duration_samples] = eeg_duration_samples
    eeg_new.header.recording[:eeg_duration_seconds] = eeg_duration_seconds
    eeg_new.header.recording[:epoch_n] = ep_n
    eeg_new.header.recording[:epoch_duration_samples] = epoch_duration_samples
    eeg_new.header.recording[:epoch_duration_seconds] = epoch_duration_seconds

    eeg_reset_components!(eeg_new)
    push!(eeg_new.header.recording[:history], "eeg_epoch(EEG, ep_n=$ep_n, epoch_len=$epoch_len)")

    return eeg_new
end

"""
    eeg_epoch!(eeg; marker, ep_offset, ep_n, epoch_len)

Split EEG into epochs. Return signal that is split either by markers (if specified), by epoch length or by number of epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `marker::String="": marker name to split at
- `ep_offset::Int64=0": time offset (in samples) for marker-based epoching (each epoch time will start at marker time - ep_offset)
- `ep_n::Union{Int64, Nothing}=nothing`: number of epochs
- `epoch_len::Union{Int64, Nothing}`=nothing: epoch length in samples
"""
function eeg_epoch!(obj::NeuroAnalyzer.NEURO; marker::String="", ep_offset::Real=0, ep_n::Union{Int64, Nothing}=nothing, epoch_len::Union{Int64, Nothing}=nothing)

    eeg_tmp = eeg_epoch(eeg, marker=marker, ep_offset=ep_offset, ep_n=ep_n, epoch_len=epoch_len)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    obj.time_pts = eeg_tmp.eeg_time
    obj.epoch_time = eeg_tmp.eeg_epoch_time
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_erp(eeg)

Average EEG epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_erp(obj::NeuroAnalyzer.NEURO)

    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals = mean(eeg_new.eeg_signals, dims=3)[:, :, :]
    eeg_duration_samples = size(eeg_new.eeg_signals, 2)
    eeg_duration_seconds = eeg_duration_samples / eeg_sr(eeg)
    eeg_time = collect(0:(1 / eeg_sr(eeg)):eeg_duration_seconds)
    eeg_new.eeg_time = eeg_time[1:(end - 1)]
    eeg_new.header.recording[:eeg_duration_samples] = eeg_duration_samples
    eeg_new.header.recording[:eeg_duration_seconds] = eeg_duration_seconds
    eeg_new.header.recording[:epoch_n] = 1

    # remove markers of deleted epochs
    for marker_idx in nrow(eeg_new.eeg_markers):-1:1
        eeg_new.eeg_markers[marker_idx, :start] > eeg_duration_samples && deleteat!(eeg_new.eeg_markers, marker_idx)
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.header.recording[:history], "eeg_erp(EEG)")

    return eeg_new
end

"""
    eeg_erp!(eeg)

Average EEG epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `epoch::Union{Vector{Int64}, AbstractRange}=1:eeg_epoch_n(eeg)`: epochs to average; default is all epochs
"""
function eeg_erp!(obj::NeuroAnalyzer.NEURO)

    eeg_tmp = eeg_erp(eeg)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    obj.time_pts = eeg_tmp.eeg_time
    obj.markers = eeg_tmp.eeg_markers
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_extract_epoch(eeg; epoch)

Extract EEG epoch.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `epoch::Int64`: epoch index

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_extract_epoch(obj::NeuroAnalyzer.NEURO; epoch::Int64)

    _check_epochs(eeg, epoch)

    s_new = reshape(obj.data[:, :, epoch], channel_n(eeg), eeg_signal_len(eeg), 1)
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals = s_new
    eeg_new.eeg_epoch_time = obj.epoch_time
    eeg_new.header.recording[:epoch_n] = 1
    eeg_new.header.recording[:eeg_duration_samples] = eeg_new.header.recording[:epoch_duration_samples]
    eeg_new.header.recording[:eeg_duration_seconds] = eeg_new.header.recording[:epoch_duration_seconds]

    eeg_reset_components!(eeg_new)
    push!(eeg_new.header.recording[:history], "eeg_extract_epoch(EEG, epoch=$epoch)")

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
function eeg_trim(obj::NeuroAnalyzer.NEURO; segment::Tuple{Int64, Int64}, remove_epochs::Bool=true)

    if remove_epochs == false
        eeg_new = deepcopy(eeg)
        eeg_epoch_n(eeg) > 1 && (eeg_epoch!(eeg_new, ep_n=1))
        _check_segment(eeg_new, segment[1], segment[2])
        eeg_new.eeg_signals = s_trim(eeg_new.eeg_signals, segment=segment)
        t_trimmed = collect(0:(1 / eeg_sr(eeg)):(size(eeg_new.eeg_signals, 2) / eeg_sr(eeg)))[1:(end - 1)]
        eeg_new.eeg_time = t_trimmed
        eeg_new.eeg_epoch_time = t_trimmed .+ obj.epoch_time[1]
        eeg_new.header.recording[:eeg_duration_samples] -= length(segment[1]:segment[2])
        eeg_new.header.recording[:eeg_duration_seconds] -= length(segment[1]:segment[2]) * (1 / eeg_sr(eeg))
        eeg_new.header.recording[:epoch_duration_samples] -= length(segment[1]:segment[2])
        eeg_new.header.recording[:epoch_duration_seconds] -= length(segment[1]:segment[2]) * (1 / eeg_sr(eeg))
        if eeg_epoch_n(eeg) > 1
            if epoch_len(eeg) <= eeg_signal_len(eeg_new)
                eeg_epoch!(eeg_new, epoch_len=epoch_len(eeg))
            else
                _info("Cannot apply original epoch length, returning single-epoch EEG.")
            end
        end
        eeg_new.eeg_markers = _delete_markers(eeg_new.eeg_markers, segment)
        eeg_new.eeg_markers = _shift_markers(eeg_new.eeg_markers, segment[1], length(segment[1]:segment[2]))
    else
        eeg_epoch_n(eeg) == 1 && throw(ArgumentError("EEG has only one epoch, cannot use remove_epochs=true."))
        epochs = _s2epoch(eeg, segment[1], segment[2])
        _info("Removing epochs: $epochs.")
        eeg_new = eeg_delete_epoch(eeg, epoch=epochs)
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.header.recording[:history], "eeg_trim(EEG, segment=$segment, remove_epochs=$remove_epochs)")

    return eeg_new
end

"""
    eeg_trim!(eeg; segment, remove_epochs)

Trim EEG signal by removing parts of the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `segment::Tuple{Int64, Int64}`: segment to be removed (from, to) in samples
- `remove_epochs::Bool=true`: remove epochs containing signal to trim (remove_epochs=true) or remove signal and remove epoching
"""
function eeg_trim!(obj::NeuroAnalyzer.NEURO; segment::Tuple{Int64, Int64}, remove_epochs::Bool=true)

    eeg_tmp = eeg_trim(eeg, segment=segment, remove_epochs=remove_epochs)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data

    eeg_reset_components!(eeg)
    return nothing
end

"""
    eeg_edit_header(eeg; field, value)

Change value of EEG header.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `field::Symbol`
- `value::Any`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_edit_header(obj::NeuroAnalyzer.NEURO; field::Symbol, value::Any)

    value === nothing && throw(ArgumentError("value cannot be empty."))

    eeg_new = deepcopy(eeg)
    field in keys(eeg_new.header.recording) || throw(ArgumentError("$field does not exist."))
    typeof(eeg_new.header.recording[field]) == typeof(value) || throw(ArgumentError("field type ($(typeof(eeg_new.header.recording[field]))) does not mach value type ($(typeof(value)))."))
    eeg_new.header.recording[field] = value
    push!(eeg_new.header.recording[:history], "eeg_edit(EEG, field=$field, value=$value)")    

    return eeg_new
end

"""
    eeg_edit_header!(eeg; field, value)

Change value of EEG header.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `field::Symbol`
- `value::Any`
"""
function eeg_edit_header!(obj::NeuroAnalyzer.NEURO; field::Symbol, value::Any)
    obj.header = eeg_edit_header(eeg, field=field, value=value).header.recording
    return nothing
end

"""
    eeg_show_header(eeg)

Show keys and values of EEG header.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_show_header(obj::NeuroAnalyzer.NEURO)
    for (key, value) in obj.header
        println("$key: $value")
    end
end

"""
    eeg_delete_epoch(eeg; epoch)

Remove EEG epoch(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) to be removed

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_delete_epoch(obj::NeuroAnalyzer.NEURO; epoch::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_epoch_n(eeg) == 1 && throw(ArgumentError("You cannot delete the last epoch."))
    typeof(epoch) <: AbstractRange && (epoch = collect(epoch))
    length(epoch) == eeg_epoch_n(eeg) && throw(ArgumentError("You cannot delete all epochs."))
    length(epoch) > 1 && (epoch = sort!(epoch, rev=true))
    _check_epochs(eeg, epoch)

    eeg_new = deepcopy(eeg)

    # remove epoch
    eeg_new.eeg_signals = eeg_new.eeg_signals[:, :, setdiff(1:end, (epoch))]

    # remove markers within deleted epochs and shift markers after the deleted epoch
    for epoch_idx in epoch
        t1, t2 = _epoch2s(eeg, epoch_idx)
        eeg_new.eeg_markers = _delete_markers(eeg_new.eeg_markers, (t1, t2))
        eeg_new.eeg_markers = _shift_markers(eeg_new.eeg_markers, t1, length(t1:t2))
    end

    # update headers
    eeg_new.header.recording[:epoch_n] -= length(epoch)
    ep_n = eeg_new.header.recording[:epoch_n]
    eeg_new.header.recording[:eeg_duration_samples] = ep_n * size(obj.data, 2)
    eeg_new.header.recording[:eeg_duration_seconds] = round((ep_n * size(obj.data, 2)) / eeg_sr(eeg), digits=2)
    eeg_new.header.recording[:epoch_duration_samples] = size(obj.data, 2)
    eeg_new.header.recording[:epoch_duration_seconds] = round(size(obj.data, 2) / eeg_sr(eeg), digits=2)

    eeg_reset_components!(eeg_new)
    push!(eeg_new.header.recording[:history], "eeg_delete_epoch(EEG, $epoch)")
    
    return eeg_new
end

"""
    eeg_delete_epoch!(eeg; epoch)

Remove EEG epoch(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) to be removed
"""
function eeg_delete_epoch!(obj::NeuroAnalyzer.NEURO; epoch::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_tmp = eeg_delete_epoch(eeg, epoch=epoch)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    obj.time_pts = eeg_tmp.eeg_time
    obj.markers = eeg_tmp.eeg_markers
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_keep_epoch(eeg; epoch)

Keep EEG epoch(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) to keep

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_keep_epoch(obj::NeuroAnalyzer.NEURO; epoch::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_epoch_n(eeg) == 1 && throw(ArgumentError("EEG contains only one epoch."))
    typeof(epoch) <: AbstractRange && (epoch = collect(epoch))
    length(epoch) > 1 && (epoch = sort!(epoch, rev=true))
    _check_epochs(eeg, epoch)

    epoch_list = collect(1:eeg_epoch_n(eeg))
    epoch_to_remove = setdiff(epoch_list, epoch)

    length(epoch_to_remove) > 1 && (epoch_to_remove = sort!(epoch_to_remove, rev=true))

    eeg_new = eeg_delete_epoch(eeg, epoch=epoch_to_remove)
    eeg_reset_components!(eeg_new)
    push!(eeg_new.header.recording[:history], "eeg_keep_epoch(EEG, $epoch)")    

    return eeg_new
end

"""
    eeg_keep_epoch!(eeg; epoch)

Keep EEG epoch(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) to keep
"""
function eeg_keep_epoch!(obj::NeuroAnalyzer.NEURO; epoch::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_tmp = eeg_keep_epoch(eeg, epoch=epoch)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    obj.time_pts = eeg_tmp.eeg_time
    obj.markers = eeg_tmp.eeg_markers
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_detect_bad(eeg; method, ch_t)

Detect bad EEG channels and epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg, type=Symbol(obj.header.recording[:data_type]))`: index of channels, default is all EEG channels
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
function eeg_detect_bad(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_get_channel_bytype(eeg, type=Symbol(obj.header.recording[:data_type])), method::Vector{Symbol}=[:flat, :rmse, :rmsd, :euclid, :p2p], w::Int64=10, ftol::Float64=0.1, fr::Float64=0.3, p::Float64=0.95, tc::Float64=0.2)

    for idx in method
        _check_var(idx, [:flat, :rmse, :rmsd, :euclid, :var, :p2p], "method")
    end

    p < 0 || p > 1 && throw(ArgumentError("p must be ≥ 0.0 and ≤ 1.0"))
    tc < 0 || tc > 1 && throw(ArgumentError("t must be ≥ 0.0 and ≤ 1.0"))

    _check_channels(eeg, channel)
    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)    

    signal = obj.data

    bad_m = zeros(Bool, ch_n, ep_n)
    bad_epochs = Vector{Int64}()

    if :flat in method
        @inbounds @simd for epoch_idx in 1:ep_n
            bad_channels_score = 0
            bad_channels = zeros(Bool, ch_n)
            Threads.@threads for channel_idx in 1:ch_n
                r = @views s_detect_channel_flat(signal[channel[channel_idx], :, epoch_idx], w=w, tol=ftol)
                if r > fr
                    bad_channels_score += 1
                    bad_channels[channel_idx] = true
                end
            end
            [bad_channels[channel_idx] == true && (bad_m[channel_idx, epoch_idx] = true) for channel_idx in 1:ch_n]
            (bad_channels_score / ch_n) > tc && push!(bad_epochs, epoch_idx)
        end
    end
    
    if :rmse in method
        @inbounds @simd for epoch_idx in 1:ep_n
            ch_m = @views vec(median(signal[channel, :, epoch_idx], dims=1))
            bad_channels_score = 0
            bad_channels = zeros(Bool, ch_n)

            rmse_ch = zeros(ch_n)
            Threads.@threads for channel_idx in 1:ch_n
                rmse_ch[channel_idx] = @views s2_rmse(signal[channel[channel_idx], :, epoch_idx], ch_m)
            end
            Threads.@threads for channel_idx in 1:ch_n
                if rmse_ch[channel_idx] < HypothesisTests.confint(OneSampleTTest(rmse_ch))[1] || rmse_ch[channel_idx] > HypothesisTests.confint(OneSampleTTest(rmse_ch))[2]
                    bad_channels_score += 1
                    bad_channels[channel_idx] = true
                end
            end
            [bad_channels[channel_idx] == true && (bad_m[channel_idx, epoch_idx] = true) for channel_idx in 1:ch_n]
            (bad_channels_score / ch_n) > tc && push!(bad_epochs, epoch_idx)
        end
    end

    if :rmsd in method
        @inbounds @simd for epoch_idx in 1:ep_n
            ch_m = @views vec(median(signal[channel, :, epoch_idx], dims=1))
            bad_channels_score = 0
            bad_channels = zeros(Bool, ch_n)

            rmsd_ch = zeros(ch_n)
            Threads.@threads for channel_idx in 1:ch_n
                rmsd_ch[channel_idx] = @views Distances.rmsd(signal[channel[channel_idx], :, epoch_idx], ch_m)
            end
            Threads.@threads for channel_idx in 1:ch_n
                if rmsd_ch[channel_idx] < HypothesisTests.confint(OneSampleTTest(rmsd_ch))[1] || rmsd_ch[channel_idx] > HypothesisTests.confint(OneSampleTTest(rmsd_ch))[2]
                    bad_channels_score += 1
                    bad_channels[channel_idx] = true
                end
            end
            [bad_channels[channel_idx] == true && (bad_m[channel_idx, epoch_idx] = true) for channel_idx in 1:ch_n]
            (bad_channels_score / ch_n) > tc && push!(bad_epochs, epoch_idx)
        end
    end

    if :euclid in method
        @inbounds @simd for epoch_idx in 1:ep_n
            ch_m = @views vec(median(signal[channel, :, epoch_idx], dims=1))
            bad_channels_score = 0
            bad_channels = zeros(Bool, ch_n)

            ed_ch = zeros(ch_n)
            Threads.@threads for channel_idx in 1:ch_n
                ed_ch[channel_idx] = @views Distances.euclidean(signal[channel[channel_idx], :, epoch_idx], ch_m)
            end
            Threads.@threads for channel_idx in 1:ch_n
                if ed_ch[channel_idx] < HypothesisTests.confint(OneSampleTTest(ed_ch))[1] || ed_ch[channel_idx] > HypothesisTests.confint(OneSampleTTest(ed_ch))[2]
                    bad_channels_score += 1
                    bad_channels[channel_idx] = true
                end
            end
            [bad_channels[channel_idx] == true && (bad_m[channel_idx, epoch_idx] = true) for channel_idx in 1:ch_n]
            (bad_channels_score / ch_n) > tc && push!(bad_epochs, epoch_idx)
        end
    end

    if :var in method
        s_v = @views var(signal[channel, :, :], dims=2)
        # mean variance
        s_mv = @views vec(mean(s_v, dims=3))
        # variance outliers
        o = reshape(outlier_detect(vec(s_v), method=:iqr), ch_n, ep_n)

        @inbounds @simd for epoch_idx in 1:ep_n
            bad_channels_score = 0
            bad_channels = zeros(Bool, ch_n)
            ch_v = @views vec(var(signal[channel, :, epoch_idx], dims=2))
            s_mv = vcat(s_mv, ch_v)
            Threads.@threads for channel_idx in 1:ch_n
                #if ch_v[channel_idx] > HypothesisTests.confint(OneSampleTTest(s_mv))[2] || o[channel_idx, epoch_idx]
                if o[channel_idx, epoch_idx]
                    bad_channels_score += 1
                    bad_channels[channel_idx] = true
                end
            end
            [bad_channels[channel_idx] == true && (bad_m[channel_idx, epoch_idx] = true) for channel_idx in 1:ch_n]
            (bad_channels_score / ch_n) > tc && push!(bad_epochs, epoch_idx)
        end
    end

    if :p2p in method
        @inbounds @simd for epoch_idx in 1:ep_n
            bad_channels_score = 0
            bad_channels = zeros(Bool, ch_n)
            Threads.@threads for channel_idx in 1:ch_n
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
            for channel_idx in 1:ch_n
                bad_channels[channel_idx] == true && (bad_m[channel_idx, epoch_idx] = true)
            end 
            (bad_channels_score / ch_n) > tc && push!(bad_epochs, epoch_idx)
        end
    end

    return (bad_m=bad_m, bad_epochs=sort(unique(bad_epochs)))
end

"""
    eeg_add_labels(eeg; labels)

Add EEG channel labels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `labels::Vector{String}`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_add_labels(obj::NeuroAnalyzer.NEURO; labels::Vector{String})

    length(labels) == channel_n(eeg) || throw(ArgumentError("labels length must be $(channel_n(eeg))."))
    
    eeg_new = deepcopy(eeg)
    eeg_new.header.recording[:labels] = labels

    push!(eeg_new.header.recording[:history], "eeg_add_labels(EEG, labels=$labels")
 
    return eeg_new
end

"""
    eeg_add_labels!(obj::NeuroAnalyzer.NEURO; labels::Vector{String})

Add EEG channel labels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `labels::Vector{String}`
"""
function eeg_add_labels!(obj::NeuroAnalyzer.NEURO; labels::Vector{String})

    length(labels) == channel_n(eeg) || throw(ArgumentError("labels length must be $(channel_n(eeg))."))
    eeg_tmp = eeg_add_labels(eeg, labels=labels)
    obj.header = obj_tmp.header

    return nothing
end

"""
    eeg_edit_channel(eeg; channel, field, value)

Edit EEG channel properties.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Int64`
- `field::Symbol`
- `value::Any`

# Returns

- `eeg_new::NeuroAnalyzer.NEURO`
"""
function eeg_edit_channel(obj::NeuroAnalyzer.NEURO; channel::Int64, field::Symbol, value::Any)
    
    value === nothing && throw(ArgumentError("value cannot be empty."))
    _check_channels(eeg, channel)
    _check_var(field, [:channel_type, :labels], "field")    

    eeg_new = deepcopy(eeg)
    typeof(eeg_new.header.recording[field][channel]) == typeof(value) || throw(ArgumentError("field type ($(eltype(eeg_new.header.recording[field]))) does not mach value type ($(typeof(value)))."))
    eeg_new.header.recording[field][channel] = value

    # add entry to :history field
    push!(eeg_new.header.recording[:history], "eeg_edit_channel(EEG, channel=$channel, field=$field, value=$value)")   

    return eeg_new
end

"""
    eeg_edit_channel!(eeg; channel, field, value)

Edit EEG channel properties.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Int64`
- `field::Symbol`
- `value::Any`
"""
function eeg_edit_channel!(obj::NeuroAnalyzer.NEURO; channel::Int64, field::Symbol, value::Any)
    
    eeg_tmp = eeg_edit_channel(eeg, channel=channel, field=field, value=value)
    obj.header = obj_tmp.header

    return nothing
end

"""
    eeg_keep_channel_type(eeg; type)

Keep EEG channels of `type` type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `type::Symbol=:eeg`: type of channels to keep

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_keep_channel_type(obj::NeuroAnalyzer.NEURO; type::Symbol=:eeg)

    _check_var(type, [:all, :eeg, :meg, :ecg, :eog, :emg, :ref, :mrk], "type")

    channels_idx = Vector{Int64}()
    for idx in 1:channel_n(eeg, type=:all)
        obj.header[:channel_type][idx] == string(type) && push!(channels_idx, idx)
    end
    eeg_new = eeg_keep_channel(eeg, channel=channels_idx)
    eeg_reset_components!(eeg_new)
    push!(eeg_new.header.recording[:history], "eeg_keep_channel_type(EEG, type=$type")

    return eeg_new
end

"""
    eeg_keep_channel_type!(eeg; type)

Keep EEG channels of `type` type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `type::Symbol=:eeg`: type of channels to keep
"""
function eeg_keep_channel_type!(obj::NeuroAnalyzer.NEURO; type::Symbol=:eeg)

    eeg_tmp = eeg_keep_channel_type(eeg, type=type)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_view_note(eeg)

Return EEG note.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_view_note(obj::NeuroAnalyzer.NEURO)
    return obj.header[:note]
end

"""
    eeg_copy(eeg)

Make copy of EEG.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `eeg_copy::NeuroAnalyzer.NEURO`
"""
function eeg_copy(obj::NeuroAnalyzer.NEURO)
    return deepcopy(eeg)
end

"""
    eeg_epoch_time(eeg; ts)

Edit EEG epochs time start.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ts::Real`: time start in seconds

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_epoch_time(obj::NeuroAnalyzer.NEURO; ts::Real)

    epoch_len = epoch_len(eeg)
    fs = eeg_sr(eeg)
    new_epochs_time = linspace(ts, ts + (epoch_len / fs), epoch_len)
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_epoch_time = new_epochs_time

    push!(eeg_new.header.recording[:history], "eeg_epoch_time(EEG, ts=$ts)")

    return eeg_new
end

"""
    eeg_epoch_time!(eeg; ts)

Edit EEG epochs time start.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ts::Real`: time start in seconds

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_epoch_time!(obj::NeuroAnalyzer.NEURO; ts::Real)

    eeg_tmp = eeg_epoch_time(eeg, ts=ts)
    obj.header = obj_tmp.header
    obj.epoch_time = eeg_tmp.eeg_epoch_time

    return nothing
end

"""
    eeg_add_note(eeg; note)

Add EEG note.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `note::String`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_add_note(obj::NeuroAnalyzer.NEURO; note::String)
    eeg_new = deepcopy(eeg)
    eeg_new.header.recording[:note] = note
    return eeg_new
end

"""
    eeg_add_note!(eeg; note)

Add EEG note.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `note::String`
"""
function eeg_add_note!(obj::NeuroAnalyzer.NEURO; note::String)
    obj.header[:note] = note
    return nothing
end

"""
    eeg_delete_note(eeg)

Delete EEG note.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_delete_note(obj::NeuroAnalyzer.NEURO)
    eeg_new = deepcopy(eeg)
    eeg_new.header.recording[:note] = ""
    return eeg_new
end

"""
    eeg_delete_note!(eeg)

Delete EEG note.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_delete_note!(obj::NeuroAnalyzer.NEURO)
    obj.header[:note] = ""
    return nothing
end

"""
    eeg_replace_channel(eeg; channel, signal)

Replace EEG channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, String}`: channel number or name
- `signal::Array{Float64, 3}`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_replace_channel(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, String}, signal::Array{Float64, 3})

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
    size(signal) == (1, epoch_len(eeg_new), eeg_epoch_n(eeg_new)) || throw(ArgumentError("signal size ($(size(signal))) must be the same as EEG channel size ($(size(eeg_new.eeg_signals[channel_idx, :, :]))."))
    eeg_new.eeg_signals[channel_idx, :, :] = signal
    eeg_reset_components!(eeg_new)

    # add entry to :history field
    push!(eeg_new.header.recording[:history], "eeg_replace_channel(EEG, channel=$channel, signal=$signal")

    return eeg_new
end

"""
    eeg_replace_channel!(eeg; channel, signal)

Replace EEG channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, String}`: channel number or name
- `signal::Array{Float64, 3}`
"""
function eeg_replace_channel!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, String}, signal::Array{Float64, 3})

    eeg_tmp = eeg_replace_channel(eeg, channel=chanel, signal=signal)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_plinterpolate_channel(eeg; channel, epoch, m, q)

Interpolate EEG channel(s) using planar interpolation.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to interpolate
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) within to interpolate
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `interpolation_factor::Int64=100`: interpolation quality

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_plinterpolate_channel(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}, epoch::Union{Int64, Vector{Int64}, AbstractRange}, imethod::Symbol=:sh, interpolation_factor::Int64=100)

    for idx in channel
        idx in eeg_get_channel_bytype(eeg, type=Symbol(obj.header.recording[:data_type])) || throw(ArgumentError("channel must be EEG/MEG signal channel(s); cannot interpolate non-EEG/MEG channels."))
    end

    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")
    obj.header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available, use eeg_load_electrodes() or eeg_add_electrodes() first."))

    typeof(channel) == Vector{Int64} && sort!(channel, rev=true)

    eeg_new = deepcopy(eeg)
    eeg_tmp = deepcopy(eeg)
    _check_channels(eeg, channel)
    _check_epochs(eeg, epoch)

    locs_x1 = obj.locs[!, :loc_x]
    locs_y1 = obj.locs[!, :loc_y]
    
    eeg_tmp = eeg_delete_channel(eeg, channel=channel)
    locs_x2 = eeg_tmp.eeg_locs[!, :loc_x]
    locs_y2 = eeg_tmp.eeg_locs[!, :loc_y]
    channels = eeg_get_channel_bytype(eeg_tmp, type=Symbol(obj.header.recording[:data_type]))

    ep_n = length(epoch)
    epoch_len = epoch_len(eeg_tmp)

    s_interpolated = zeros(Float64, length(channel), epoch_len, ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * epoch_len, 1))

    @inbounds @simd for epoch_idx in eachindex(epoch)
        Threads.@threads for length_idx in 1:epoch_len
            s_tmp, x, y = @views _interpolate(obj_tmp.data[channels, length_idx, epoch[epoch_idx]], locs_x2, locs_y2, interpolation_factor, imethod, :none)
            for channel_idx in eachindex(channel)
                x_idx = vsearch(locs_x1[channel[channel_idx]], x)
                y_idx = vsearch(locs_y1[channel[channel_idx]], y)
                s_interpolated[channel_idx, length_idx, epoch_idx] = s_tmp[x_idx, y_idx]
            end

            # update progress bar
            progress_bar == true && next!(p)
        end
    end

    eeg_new.eeg_signals[channel, :, epoch] = s_interpolated

    eeg_reset_components!(eeg_new)
    push!(eeg_new.header.recording[:history], "eeg_plinterpolate_channel(EEG, channel=$channel, epoch=$epoch, imethod=$imethod, interpolation_factor=$interpolation_factor)")

    return eeg_new
end

"""
    eeg_plinterpolate_channel!(eeg; channel, epoch, imethod, interpolation_factor)

Interpolate EEG channel(s) using planar interpolation.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to interpolate
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) within to interpolate
- `imethod::Symbol=:sh`: interpolation method Shepard (`:sh`), Multiquadratic (`:mq`), InverseMultiquadratic (`:imq`), ThinPlate (`:tp`), NearestNeighbour (`:nn`), Gaussian (`:ga`)
- `interpolation_factor::Int64=100`: interpolation quality
"""
function eeg_plinterpolate_channel!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}}, epoch::Union{Int64, Vector{Int64}, AbstractRange}, imethod::Symbol=:shepard, interpolation_factor::Int64=100)

    eeg_tmp = eeg_plinterpolate_channel(eeg, channel=channel, epoch=epoch, imethod=imethod, interpolation_factor=interpolation_factor)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_channel_type(eeg; channel, type)

Change EEG channel type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, String}`
- `type::String`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_channel_type(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, String}, type::String)

    type = lowercase(type)
    labels = eeg_labels(eeg)

    # create new dataset
    eeg_new = deepcopy(eeg)
    types = eeg_new.header.recording[:channel_type]
    
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
    eeg_new.header.recording[:channel_type] = types
    
    # add entry to :history field
    push!(eeg_new.header.recording[:history], "eeg_channel_type(EEG, channel=$channel, type=$type)")

    return eeg_new
end

"""
    eeg_channel_type!(eeg; channel, new_name)

Change EEG channel type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, String}`
- `type::String`
"""
function eeg_channel_type!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, String}, type::String)

    eeg_tmp = eeg_channel_type(eeg, channel=channel, type=type)
    obj.header = obj_tmp.header

    return nothing
end

"""
    eeg_edit_electrode(eeg; <keyword arguments>)

Edit EEG electrode.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
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

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_edit_electrode(obj::NeuroAnalyzer.NEURO; channel::Union{String, Int64}, x::Union{Real, Nothing}=nothing, y::Union{Real, Nothing}=nothing, z::Union{Real, Nothing}=nothing, theta::Union{Real, Nothing}=nothing, radius::Union{Real, Nothing}=nothing, theta_sph::Union{Real, Nothing}=nothing, radius_sph::Union{Real, Nothing}=nothing, phi_sph::Union{Real, Nothing}=nothing, name::String="", type::String="")

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

    (x !== nothing || y !== nothing || z !== nothing || theta !== nothing || radius !== nothing || theta_sph !== nothing  || radius_sph !== nothing || phi_sph !== nothing) && (eeg_new.header.recording[:channel_locations] == true)

    eeg_reset_components!(eeg_new)
    push!(eeg_new.header.recording[:history], "eeg_edit_electrode(EEG; channel=$channel, x=$x, y=$y, z=$z, theta=$theta, radius=$radius, theta_sph=$theta_sph, radius_sph=$radius_sph, phi_sph=$phi_sph, name=$name, type=$type)")

    return eeg_new
end

"""
    eeg_edit_electrode!(eeg; <keyword arguments>)

Edit EEG electrode.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
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
function eeg_edit_electrode!(obj::NeuroAnalyzer.NEURO; channel::Union{String, Int64}, x::Union{Real, Nothing}=nothing, y::Union{Real, Nothing}=nothing, z::Union{Real, Nothing}=nothing, theta::Union{Real, Nothing}=nothing, radius::Union{Real, Nothing}=nothing, theta_sph::Union{Real, Nothing}=nothing, radius_sph::Union{Real, Nothing}=nothing, phi_sph::Union{Real, Nothing}=nothing, name::String="", type::String="")

    eeg_tmp = eeg_edit_electrode(eeg, channel=channel, x=x, y=y, z=z, theta=theta, radius=radius, theta_sph=theta_sph, radius_sph=radius_sph, phi_sph=phi_sph, name=name, type=type)
    obj.header = obj_tmp.header
    obj.locs = eeg_tmp.eeg_locs
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_electrode_loc(eeg; channel, output)

Return locations of EEG channel electrode.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
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
function eeg_electrode_loc(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, String}, output::Bool=true)

    obj.header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available, use eeg_load_electrodes() or eeg_add_electrodes() first."))

    channel = _get_channel_idx(eeg_labels(eeg), channel)

    x = obj.locs[!, :loc_x][channel]
    y = obj.locs[!, :loc_y][channel]
    z = obj.locs[!, :loc_z][channel]
    theta = obj.locs[!, :loc_theta][channel]
    radius = obj.locs[!, :loc_radius][channel]
    theta_sph = obj.locs[!, :loc_theta_sph][channel]
    radius_sph = obj.locs[!, :loc_radius_sph][channel]
    phi_sph = obj.locs[!, :loc_phi_sph][channel]

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
    eeg_view_marker(eeg)

Show markers.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_view_marker(obj::NeuroAnalyzer.NEURO)
    obj.header.markers == true || throw(ArgumentError("EEG has no markers."))
    for marker_idx in 1:size(obj.markers, 1)
        println("ID: $(rpad(("'" * obj.markers[marker_idx, :id] * "'"), 24, " ")) start [sample]: $(rpad(obj.markers[marker_idx, :start], 8, " ")) length [samples]: $(rpad(obj.markers[marker_idx, :length], 8, " ")) description: $(rpad(("'" * obj.markers[marker_idx, :description] * "'"), 24, " ")) channel: $(obj.markers[marker_idx, :channel])")
    end
end

"""
    eeg_delete_marker(eeg; n)

Delete marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `n::Int64`: marker number

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_delete_marker(obj::NeuroAnalyzer.NEURO; n::Int64)
    eeg_new = deepcopy(eeg)
    eeg_new.header.recording[:markers] == true || throw(ArgumentError("EEG has no markers."))
    nn = size(eeg_new.eeg_markers, 1)
    (n < 1 || n > nn) && throw(ArgumentError("n has to be ≥ 1 and ≤ $nn."))
    deleteat!(eeg_new.eeg_markers, n)
    size(eeg_new.eeg_markers, 1) == 0 && (eeg_new.header.recording[:markers] = false)
    eeg_reset_components!(eeg_new)
    push!(eeg_new.header.recording[:history], "eeg_delete_marker(EEG; n=$n)")
    
    return eeg_new
end

"""
    eeg_delete_marker!(eeg; n)

Delete marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `n::Int64`: marker number
"""
function eeg_delete_marker!(obj::NeuroAnalyzer.NEURO; n::Int64)

    eeg_tmp = eeg_delete_marker(eeg, n=n)
    obj.header.markers = obj_tmp.header[:markers]
    obj.markers = eeg_tmp.eeg_markers
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_add_marker(eeg; id, start, len, desc, channel)

Add marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `id::String`: marker ID
- `start::Int64`: marker time in samples
- `len::Int64=1`: marker length in samples
- `desc::String`: marker description
- `channel::Int64=0`: channel number, if 0 then marker is related to all channels

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_add_marker(obj::NeuroAnalyzer.NEURO; id::String, start::Int64, len::Int64=1, desc::String, channel::Int64=0)

    start < 1 && throw(ArgumentError("start must be > 0."))
    len < 1 && throw(ArgumentError("len must be > 0."))
    start >= eeg_signal_len(eeg) && throw(ArgumentError("start must be < $(eeg_signal_len(eeg) - 1)."))
    start + len > eeg_signal_len(eeg) && throw(ArgumentError("start + len must be ≤ $(eeg_signal_len(eeg))."))

    eeg_new = deepcopy(eeg)
    eeg_new.header.recording[:markers] = true
    append!(eeg_new.eeg_markers, DataFrame(:id => id, :start => start, :length => len, :description => desc, :channel => channel))
    sort!(eeg_new.eeg_markers)
    eeg_reset_components!(eeg_new)
    push!(eeg_new.header.recording[:history], "eeg_add_marker(EEG; id=$id, start=$start, len=$len, desc=$desc, channel=$channel)")

    return eeg_new
end

"""
    eeg_add_marker!(eeg; id, start, len, desc, channel)

Add marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `id::String`: marker ID
- `start::Int64`: marker time in samples
- `len::Int64=1`: marker length in samples
- `desc::String`: marker description
- `channel::Int64=0`: channel number, if 0 then marker is related to all channels
"""
function eeg_add_marker!(obj::NeuroAnalyzer.NEURO; id::String, start::Int64, len::Int64=1, desc::String, channel::Int64=0)

    eeg_tmp = eeg_add_marker(eeg, id=id, start=start, len=len, desc=desc, channel=channel)
    obj.header.markers = obj_tmp.header[:markers]
    obj.markers = eeg_tmp.eeg_markers
    eeg_reset_components!(eeg)

    return nothing
end


"""
    eeg_vch(eeg; f)

Calculate virtual channel using formula `f`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `f::String`: channel calculation formula, e.g. `"cz / mean(fp1 + fp2)"`; case of labels in the formula is ignored, all standard Julia math operators are available, channel labels must be the same as of the EEG object

# Returns
 
- `vc::Array{Float64, 3}`: single channel × time × epochs
"""
function eeg_vch(obj::NeuroAnalyzer.NEURO; f::String)

    ep_n = eeg_epoch_n(eeg)
    f = lowercase(f)
    labels = lowercase.(eeg_labels(eeg))
    vc = zeros(1, epoch_len(eeg), ep_n)
    Threads.@threads for epoch_idx in 1:ep_n
        f_tmp = f
        @inbounds for channel_idx in eachindex(labels)
            occursin(labels[channel_idx], f) == true && (f_tmp = replace(f_tmp, labels[channel_idx] => "$(obj.data[channel_idx, :, epoch_idx])"))
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

- `obj::NeuroAnalyzer.NEURO`
- `n::Int64`: marker number
- `id::String`: marker ID
- `start::Int64`: marker time in samples
- `len::Int64=1`: marker length in samples
- `desc::String`: marker description
- `channel::Int64`: channel number, if 0 then marker is related to all channels

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_edit_marker(obj::NeuroAnalyzer.NEURO; n::Int64, id::String, start::Int64, len::Int64=1, desc::String, channel::Int64)

    obj.header.markers == true || throw(ArgumentError("EEG has no markers."))
    start < 1 && throw(ArgumentError("start must be > 0."))
    len < 1 && throw(ArgumentError("len must be > 0."))
    start >= eeg_signal_len(eeg) && throw(ArgumentError("start must be < $(eeg_signal_len(eeg))."))
    start + len > eeg_signal_len(eeg) && throw(ArgumentError("start + len must be ≤ $(eeg_signal_len(eeg))."))

    nn = size(obj.markers, 1)
    n < 1 || n > nn && throw(ArgumentError("n has to be ≥ 1 and ≤ $nn."))
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_markers[n, :] = Dict(:id => id, :start => start, :length => len, :description => desc, :channel => channel)
     eeg_reset_components!(eeg_new)
    push!(eeg_new.header.recording[:history], "eeg_edit_marker(EEG; id=$id, start=$start, len=$len, desc=$desc, channel=$channel)")

    return eeg_new
end

"""
    eeg_edit_marker!(eeg; n, id, start, len, desc)

Edit EEG marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `n::Int64`: marker number
- `id::String`: marker ID
- `start::Int64`: marker time in samples
- `len::Int64=1`: marker length in samples
- `desc::String`: marker description
- `channel::Int64`: channel number, if 0 then marker is related to all channels
"""
function eeg_edit_marker!(obj::NeuroAnalyzer.NEURO; n::Int64, id::String, start::Int64, len::Int64=1, desc::String, channel::Int64)

    eeg_tmp = eeg_edit_marker(eeg, n=n, id=id, start=start, len=len, desc=desc, channel=channel)
    obj.header.markers = obj_tmp.header[:markers]
    obj.markers = eeg_tmp.eeg_markers
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_lrinterpolate_channel(eeg; channel, epoch)

Interpolate EEG channel using linear regression.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number to interpolate
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) within to interpolate

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_lrinterpolate_channel(obj::NeuroAnalyzer.NEURO; channel::Int64, epoch::Union{Int64, Vector{Int64}, AbstractRange})

    channels = eeg_get_channel_bytype(eeg, type=Symbol(obj.header.recording[:data_type]))
    channel in channels || throw(ArgumentError("channel must be EEG/MEG signal channel; cannot interpolate non-EEG/MEG channels."))

    bad_signal = obj.data[:, :, epoch]
    good_epochs = setdiff(1:eeg_epoch_n(eeg), epoch)
    good_channels = setdiff(channels, channel)
    good_signal = _make_epochs(obj.data[:, :, good_epochs], ep_n=1)

    # train
    df = @views DataFrame(hcat(good_signal[channel, :, 1], good_signal[good_channels, :, 1]'), :auto)
    train, test = _split(df, 0.80)
    fm = Term(:x1) ~ sum(Term.(Symbol.(names(df[!, Not(:x1)]))))
    linear_regressor = lm(fm, train)
    acc_r2 = r2(linear_regressor)
    prediction = GLM.predict(linear_regressor, test)
    accuracy_testdf = DataFrame(signal_actual = test[!, :x1], signal_predicted = prediction)
    accuracy_testdf.error = accuracy_testdf[!, :signal_actual] 
    accuracy_testdf[!, :signal_predicted]
    acc_mae = mean(abs.(accuracy_testdf.error))
    aic, bic = infcrit(linear_regressor)
    _info("R² for the linear regressor: $(round(acc_r2, digits=3))")
    _info("MAE (test dataset): $(round(acc_mae, digits=3))")
    _info("AIC: $(round(aic, digits=3))")
    _info("BIC: $(round(bic, digits=3))")

    # predict
    eeg_new = eeg_copy(eeg) 
    @inbounds @simd for epoch_idx in epoch
        df = @views DataFrame(hcat(bad_signal[channel, :, epoch_idx], bad_signal[good_channels, :, epoch_idx]'), :auto)
        eeg_new.eeg_signals[channel, :, epoch_idx] = GLM.predict(linear_regressor, df)
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.header.recording[:history], "eeg_lrinterpolate_channel(EEG, channel=$channel, epoch=$epoch)")

    return eeg_new
end

"""
    eeg_lrinterpolate_channel!(eeg; channel, epoch)

Interpolate EEG channel using linear regression.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number to interpolate
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) within to interpolate

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_lrinterpolate_channel!(obj::NeuroAnalyzer.NEURO; channel::Int64, epoch::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_tmp = eeg_lrinterpolate_channel(eeg, channel=channel, epoch=epoch)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_reflect(eeg; n)

Expand signal by adding reflected signal before the signal and after the signal, i.e. a signal 1234 becomes 432112344321. This may reduce edge artifacts, but will also affect amplitude of the filtered signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `n::Int64=eeg_sr(eeg)`: number of samples to add, default is 1 second

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_reflect(obj::NeuroAnalyzer.NEURO; n::Int64=eeg_sr(eeg))

    # add up to one epoch
    n > epoch_len(eeg) && (n = epoch_len(eeg))

    eeg_new = deepcopy(eeg)
    ch_n = channel_n(eeg)
    ep_n = eeg_epoch_n(eeg)
    s = zeros(ch_n, epoch_len(eeg) + 2 * n, ep_n)

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            s1 = eeg_new.eeg_signals[:, 1:n, epoch_idx]
            s2 = eeg_new.eeg_signals[:, end:-1:(end - n + 1), epoch_idx]
            @views s[channel_idx, :, epoch_idx] = _reflect(obj.data[channel_idx, :, epoch_idx], s1[channel_idx, :], s2[channel_idx, :])
        end
    end
    eeg_new.eeg_signals = s

    t = collect(0:(1 / eeg_sr(eeg)):(size(eeg_new.eeg_signals, 2) / eeg_sr(eeg)))[1:(end - 1)]
    eeg_new.eeg_time = t
    eeg_new.eeg_epoch_time = t .+ obj.epoch_time[1]
    eeg_new.header.recording[:eeg_duration_samples] = length(t) * ep_n
    eeg_new.header.recording[:eeg_duration_seconds] = length(t) * ep_n * (1 / eeg_sr(eeg))
    eeg_new.header.recording[:epoch_duration_samples] = size(eeg_new.eeg_signals, 2)
    eeg_new.header.recording[:epoch_duration_seconds] = size(eeg_new.eeg_signals, 2) * (1 / eeg_sr(eeg))

    push!(eeg_new.header.recording[:history], "eeg_reflect(EEG, n=$n)")

    return eeg_new
end

"""
    eeg_reflect!(eeg; n)

Expand signal by adding reflected signal before the signal and after the signal, i.e. a signal 1234 becomes 432112344321. This may reduce edge artifacts, but will also affect amplitude of the filtered signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `n::Int64=eeg_sr(eeg)`: number of samples to add, default is 1 second
"""
function eeg_reflect!(obj::NeuroAnalyzer.NEURO; n::Int64=eeg_sr(eeg))

    eeg_tmp = eeg_reflect(eeg, n=n)

    obj.header = obj_tmp.header
    obj.data = obj_tmp.data

    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_chop(eeg; n)

Reduce signal by removing reflected signal before the signal and after the signal, i.e. a signal 432112344321 becomes 1234.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `n::Int64=eeg_sr(eeg)`: number of samples to remove, default is 1 second

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_chop(obj::NeuroAnalyzer.NEURO; n::Int64=eeg_sr(eeg))

    # add up to one epoch
    n > epoch_len(eeg) && (n = epoch_len(eeg))

    eeg_new = deepcopy(eeg)
    ch_n = channel_n(eeg)
    ep_n = eeg_epoch_n(eeg)
    s = zeros(ch_n, epoch_len(eeg) - 2 * n, ep_n)

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            @views s[channel_idx, :, epoch_idx] = _chop(obj.data[channel_idx, :, epoch_idx], n)
        end
    end
    eeg_new.eeg_signals = s

    t = collect(0:(1 / eeg_sr(eeg)):(size(eeg_new.eeg_signals, 2) / eeg_sr(eeg)))[1:(end - 1)]
    eeg_new.eeg_time = t
    eeg_new.eeg_epoch_time = t .+ obj.epoch_time[1]
    eeg_new.header.recording[:eeg_duration_samples] = length(t) * ep_n
    eeg_new.header.recording[:eeg_duration_seconds] = length(t) * ep_n * (1 / eeg_sr(eeg))
    eeg_new.header.recording[:epoch_duration_samples] = size(eeg_new.eeg_signals, 2)
    eeg_new.header.recording[:epoch_duration_seconds] = size(eeg_new.eeg_signals, 2) * (1 / eeg_sr(eeg))

    push!(eeg_new.header.recording[:history], "eeg_chop(EEG, n=$n)")

    return eeg_new
end

"""
    eeg_chop!(eeg; c, v)

Reduce signal by removing reflected signal before the signal and after the signal, i.e. a signal 432112344321 becomes 1234.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `n::Int64=eeg_sr(eeg)`: number of samples to remove, default is 1 second
"""
function eeg_chop!(obj::NeuroAnalyzer.NEURO; n::Int64=eeg_sr(eeg))

    eeg_tmp = eeg_chop(eeg, n=n)

    obj.header = obj_tmp.header
    obj.data = obj_tmp.data

    eeg_reset_components!(eeg)

    return nothing
end

"""
    eeg_extract_data(eeg; channel)

Extract EEG data.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg)`: index of channels, default is all signal channels
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}=eeg_epoch_n(eeg)`: index of epochs, default is all epochs
- `time::Bool=false`: return time vector
- `etime::Bool=false`: return epoch time vector

# Returns

- `signal::Array{Float64, 3}`
- `time::Vector{Float64}`
- `etime::Vector{Float64}`
"""
function eeg_extract_data(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=eeg_signal_channels(eeg), epoch::Union{Int64, Vector{Int64}, AbstractRange}=eeg_epoch_n(eeg), time::Bool=false, etime::Bool=false)
    _check_channels(eeg, channel)
    _check_epochs(eeg, epoch)
    if time == false && etime == false
        return obj.data[channel, :, epoch][:, :, :]
    elseif time == true && etime == false
        return obj.data[channel, :, epoch][:, :, :], obj.time_pts
    elseif time == false && etime == true
        return obj.data[channel, :, epoch][:, :, :], obj.epoch_time
    else
        return obj.data[channel, :, epoch][:, :, :], obj.time_pts, obj.epoch_time
    end
end

"""
    eeg_extract_time(eeg)

Extract EEG time.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `time::Array{Float64, 3}`
"""
function eeg_extract_time(obj::NeuroAnalyzer.NEURO)
    return obj.time_pts
end

"""
    eeg_extract_etime(eeg)

Extract EEG epochs time.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `time::Array{Float64, 3}`
"""
function eeg_extract_etime(obj::NeuroAnalyzer.NEURO)
    return obj.epoch_time
end