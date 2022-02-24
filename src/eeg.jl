"""
    eeg_reset_components(eeg)

Resets `eeg` components.

# Arguments

- `eeg:EEG`

# Returns

- `eeg::EEG`
"""
function eeg_reset_components(eeg::EEG)
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_header[:components] = []
    eeg_new.eeg_components = []

    return eeg_new
end

"""
    eeg_reset_components!(eeg)

Resets `eeg` components.

# Arguments

- `eeg:EEG`
"""
function eeg_reset_components!(eeg::EEG)
    eeg.eeg_header[:components] = []
    eeg.eeg_components = []

    return
end

"""
    eeg_delete_channel(eeg; channel)

Removes `channel` from the `eeg`.

# Arguments

- `eeg::EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}` - channel index to be removed, vector of numbers or range

# Returns

- `eeg::EEG`
"""
function eeg_delete_channel(eeg::EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})
    typeof(channel) <: AbstractRange && (channel = collect(channel))
    channel_n = size(eeg.eeg_signals, 1)
    length(channel) == channel_n && throw(ArgumentError("You cannot delete all channels."))

    length(channel) > 1 && (channel = sort!(channel, rev=true))

    if channel[end] < 1 || channel[1] > size(eeg.eeg_signals, 1)
        throw(ArgumentError("Channel index does not match signal channels."))
    end

    eeg_header = deepcopy(eeg.eeg_header)
    eeg_time = deepcopy(eeg.eeg_time)
    eeg_signals = deepcopy(eeg.eeg_signals)
    eeg_components = []

    # update headers
    eeg_header[:channel_n] = channel_n - length(channel)
    for idx1 in 1:length(channel)
        for idx2 in 1:channel_n
            if idx2 == channel[idx1]
                deleteat!(eeg_header[:labels], idx2)
                eeg_header[:channel_locations] == true && (deleteat!(eeg_header[:xlocs], idx2))
                eeg_header[:channel_locations] == true && (deleteat!(eeg_header[:xlocs], idx2))
                deleteat!(eeg_header[:transducers], idx2)
                deleteat!(eeg_header[:physical_dimension], idx2)
                deleteat!(eeg_header[:physical_minimum], idx2)
                deleteat!(eeg_header[:physical_maximum], idx2)
                deleteat!(eeg_header[:digital_minimum], idx2)
                deleteat!(eeg_header[:digital_maximum], idx2)
                deleteat!(eeg_header[:prefiltering], idx2)
                deleteat!(eeg_header[:samples_per_datarecord], idx2)
                deleteat!(eeg_header[:sampling_rate], idx2)
                deleteat!(eeg_header[:gain], idx2)
            end
        end 
    end

    # remove channel
    eeg_signals = eeg_signals[setdiff(1:end, (channel)), :, :]

    # create new dataset
    eeg_new = EEG(eeg_header, eeg.eeg_time, eeg.eeg_signals, deepcopy(eeg_components))

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_delete_channel(EEG, $channel)")
    
    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_delete_channel!(eeg; channel)

Removes `channel` from the `eeg`.

# Arguments

- `eeg::EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}` - channel index to be removed
"""
function eeg_delete_channel!(eeg::EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})
    typeof(channel) <: AbstractRange && (channel = collect(channel))
    channel_n = size(eeg.eeg_signals, 1)
    length(channel) == channel_n && throw(ArgumentError("You cannot delete all channels."))

    length(channel) > 1 && (channel = sort!(channel, rev=true))

    if channel[end] < 1 || channel[1] > size(eeg.eeg_signals, 1)
        throw(ArgumentError("Channel index does not match signal channels."))
    end

    # update headers
    eeg.eeg_header[:channel_n] = channel_n - length(channel)
    for idx1 in 1:length(channel)
        for idx2 in 1:channel_n
            if idx2 == channel[idx1]
                deleteat!(eeg.eeg_header[:labels], idx2)
                eeg.eeg_header[:channel_locations] == true && (deleteat!(eeg.eeg_header[:xlocs], idx2))
                eeg.eeg_header[:channel_locations] == true && (deleteat!(eeg.eeg_header[:xlocs], idx2))
                deleteat!(eeg.eeg_header[:transducers], idx2)
                deleteat!(eeg.eeg_header[:physical_dimension], idx2)
                deleteat!(eeg.eeg_header[:physical_minimum], idx2)
                deleteat!(eeg.eeg_header[:physical_maximum], idx2)
                deleteat!(eeg.eeg_header[:digital_minimum], idx2)
                deleteat!(eeg.eeg_header[:digital_maximum], idx2)
                deleteat!(eeg.eeg_header[:prefiltering], idx2)
                deleteat!(eeg.eeg_header[:samples_per_datarecord], idx2)
                deleteat!(eeg.eeg_header[:sampling_rate], idx2)
                deleteat!(eeg.eeg_header[:gain], idx2)
            end
        end 
    end

    # remove channel
    eeg.eeg_signals = eeg.eeg_signals[setdiff(1:end, (channel)), :, :]

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_delete_channel!(EEG, $channel)")
    
    eeg_reset_components!(eeg)

    return
end

"""
    eeg_keep_channel(eeg; channel)

Keeps `channels` in the `eeg`.

# Arguments

- `eeg::EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}` - channel index to keep

# Returns

- `eeg::EEG`
"""
function eeg_keep_channel(eeg::EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})
    typeof(channel) <: AbstractRange && (channel = collect(channel))

    length(channel) > 1 && (channel = sort!(channel, rev=true))
    if channel[end] < 1 || channel[1] > size(eeg.eeg_signals, 1)
        throw(ArgumentError("Channel index does not match signal channels."))
    end

    channel_list = collect(1:size(eeg.eeg_signals, 1))
    channel_to_remove = setdiff(channel_list, channel)

    length(channel_to_remove) > 1 && (channel_to_remove = sort!(channel_to_remove, rev=true))

    eeg_header = deepcopy(eeg.eeg_header)
    eeg_time = deepcopy(eeg.eeg_time)
    eeg_signals = deepcopy(eeg.eeg_signals)

    channel_n = eeg_header[:channel_n]

    # update headers
    eeg_header[:channel_n] = channel_n - length(channel_to_remove)
    for idx1 in 1:length(channel_to_remove)
        for idx2 in 1:channel_n
            if idx2 == channel_to_remove[idx1]
                deleteat!(eeg_header[:labels], idx2)
                eeg_header[:channel_locations] == true && deleteat!(eeg_header[:xlocs], idx2)
                eeg_header[:channel_locations] == true && deleteat!(eeg_header[:ylocs], idx2)
                deleteat!(eeg_header[:transducers], idx2)
                deleteat!(eeg_header[:physical_dimension], idx2)
                deleteat!(eeg_header[:physical_minimum], idx2)
                deleteat!(eeg_header[:physical_maximum], idx2)
                deleteat!(eeg_header[:digital_minimum], idx2)
                deleteat!(eeg_header[:digital_maximum], idx2)
                deleteat!(eeg_header[:prefiltering], idx2)
                deleteat!(eeg_header[:samples_per_datarecord], idx2)
                deleteat!(eeg_header[:sampling_rate], idx2)
                deleteat!(eeg_header[:gain], idx2)
            end
        end
    end

    # remove channel
    eeg_signals = eeg_signals[setdiff(1:end, (channel_to_remove)), :, :]

    # create new dataset
    eeg_new = EEG(eeg_header, eeg_time, eeg_signals, deepcopy(eeg.eeg_components))

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_keep_channel(EEG, $channel)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_keep_channel!(eeg; channel)

Keeps `channels` in the `eeg`.

# Arguments

- `eeg::EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}` - channel index to keep
"""
function eeg_keep_channel!(eeg::EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})
    typeof(channel) <: AbstractRange && (channel = collect(channel))

    length(channel) > 1 && (channel = sort!(channel, rev=true))
    if channel[end] < 1 || channel[1] > size(eeg.eeg_signals, 1)
        throw(ArgumentError("Channel index does not match signal channels."))
    end

    channel_list = collect(1:size(eeg.eeg_signals, 1))
    channel_to_remove = setdiff(channel_list, channel)

    length(channel_to_remove) > 1 && (channel_to_remove = sort!(channel_to_remove, rev=true))

    channel_n = eeg.eeg_header[:channel_n]

    # update headers
    eeg.eeg_header[:channel_n] = channel_n - length(channel_to_remove)
    for idx1 in 1:length(channel_to_remove)
        for idx2 in 1:channel_n
            if idx2 == channel_to_remove[idx1]
                deleteat!(eeg.eeg_header[:labels], idx2)
                eeg.eeg_header[:channel_locations] == true && deleteat!(eeg.eeg_header[:xlocs], idx2)
                eeg.eeg_header[:channel_locations] == true && deleteat!(eeg.eeg_header[:ylocs], idx2)
                deleteat!(eeg.eeg_header[:transducers], idx2)
                deleteat!(eeg.eeg_header[:physical_dimension], idx2)
                deleteat!(eeg.eeg_header[:physical_minimum], idx2)
                deleteat!(eeg.eeg_header[:physical_maximum], idx2)
                deleteat!(eeg.eeg_header[:digital_minimum], idx2)
                deleteat!(eeg.eeg_header[:digital_maximum], idx2)
                deleteat!(eeg.eeg_header[:prefiltering], idx2)
                deleteat!(eeg.eeg_header[:samples_per_datarecord], idx2)
                deleteat!(eeg.eeg_header[:sampling_rate], idx2)
                deleteat!(eeg.eeg_header[:gain], idx2)
            end
        end
    end

    # remove channel
    eeg.eeg_signals = eeg.eeg_signals[setdiff(1:end, (channel_to_remove)), :, :]

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_keep_channel!(EEG, $channel)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_derivative(eeg)

Returns the derivative of the `eeg` with length same as the signal.

# Arguments

- `eeg::EEG` - EEG object

# Returns

- `eeg::EEG`
"""
function eeg_derivative(eeg::EEG)
    s_der = signal_derivative(eeg.eeg_signals)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_der, deepcopy(eeg.eeg_components))
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_derivative(EEG)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_derivative!(eeg)

Returns the derivative of the `eeg` with length same as the signal.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_derivative!(eeg::EEG)
    eeg.eeg_signals = signal_derivative(eeg.eeg_signals)

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_derivative!(EEG)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_total_power(eeg)

Calculates total power of the `eeg`.

# Arguments

- `eeg::EEG`

# Returns

- `stp::Vector{Float64}`
"""
function eeg_total_power(eeg::EEG)
    fs = eeg.eeg_header[:sampling_rate][1]
    stp = signal_total_power(eeg.eeg_signals, fs=fs)
    size(stp, 3) == 1 && (stp = reshape(stp, size(stp, 1), size(stp, 2)))

    return stp
end

"""
    eeg_total_power!(eeg)

Calculates total power of the `eeg`.

# Arguments

- `eeg::EEG`
"""
function eeg_total_power!(eeg::EEG)
    :total_power in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:total_power)
    fs = eeg.eeg_header[:sampling_rate][1]
    stp = signal_total_power(eeg.eeg_signals, fs=fs)
    push!(eeg.eeg_components, stp)
    push!(eeg.eeg_header[:components], :total_power)
    push!(eeg.eeg_header[:history], "eeg_total_power!(EEG)")

    return
end

"""
    eeg_band_power(eeg; f)

Calculates absolute band power between frequencies `f[1]` and `f[2]` of the `eeg`.

# Arguments

- `eeg::EEG`
- `f::Tuple(Union(Int64, Float64}, Union(Int64, Float64}}` - lower and upper frequency bounds

# Returns

- `sbp::Vector{Float64}`
"""
function eeg_band_power(eeg::EEG; f::Tuple)
    fs = eeg.eeg_header[:sampling_rate][1]
    sbp = signal_band_power(eeg.eeg_signals, fs=fs, f=f)
    (size(sbp, 3) == 1) && (sbp = reshape(sbp, size(sbp, 1), size(sbp, 2)))

    return sbp
end

"""
    eeg_detrend(eeg; type=:linear)

Removes linear trend from the `eeg`.

# Arguments

- `eeg::EEG`
- `type::Symbol[:linear, :constant]`, optional
    - `linear` - the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `constant` - the mean of `signal` is subtracted

# Returns

- `eeg::EEG`
"""
function eeg_detrend(eeg::EEG; type::Symbol=:linear)
    s_det = signal_detrend(eeg.eeg_signals, type=type)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_det, deepcopy(eeg.eeg_components))
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_detrend(EEG, type=$type)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_detrend!(eeg; type=:linear)

Removes linear trend from the `eeg`.

# Arguments

- `eeg::EEG`
- `type::Symbol[:linear, :constant]`, optional
    - `linear` - the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `constant` - the mean of `signal` is subtracted
"""
function eeg_detrend!(eeg::EEG; type::Symbol=:linear)
    eeg.eeg_signals = signal_detrend(eeg.eeg_signals, type=type)
    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_detrend!(EEG, type=$type)")

    return
end

"""
    eeg_reference_channel(eeg; channel)

References the `eeg` to specific channel `channel`.

# Arguments

- `eeg::EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}` - index of channels used as reference; if multiple channels are specified, their average is used as the reference

# Returns

- `eeg::EEG`
"""
function eeg_reference_channel(eeg::EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})
    typeof(channel) <: AbstractRange && (channel = collect(channel))

    s_referenced = signal_reference_channel(eeg.eeg_signals, channel=channel)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_referenced, deepcopy(eeg.eeg_components))
    eeg_new.eeg_header[:reference] = "channel: $channel"
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_reference_channel(EEG, channel=$channel)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_reference_channel!(eeg; channel)

References the `eeg` to specific channel `channel`.

# Arguments

- `eeg::EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}` - index of channels used as reference; if multiple channels are specified, their average is used as the reference
"""
function eeg_reference_channel!(eeg::EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})
    eeg.eeg_signals = signal_reference_channel(eeg.eeg_signals, channel=channel)
    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_reference_channel!(EEG, channel=$channel)")

    return
end

"""
    eeg_reference_car(eeg)

References the `eeg` to common average reference.

# Arguments

- `eeg::EEG`

# Returns

- `eeg::EEG`
"""
function eeg_reference_car(eeg::EEG)
    s_referenced = signal_reference_car(eeg.eeg_signals)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_referenced, deepcopy(eeg.eeg_components))
    eeg_new.eeg_header[:reference] = "CAR"
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_reference_car(EEG)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_reference_car!(eeg)

References the `eeg` to common average reference.

# Arguments

- `eeg::EEG`
"""
function eeg_reference_car!(eeg::EEG)
    eeg.eeg_signals = signal_reference_car(eeg.eeg_signals)

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_reference_car!(EEG)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_get_channel(eeg; channel)

Returns the `channel` index / name.

# Arguments

- `eeg::EEG`
- `channel::Union{Int64, String}` - channel name

# Returns

- `channel_idx::Int64`
"""
function eeg_get_channel(eeg::EEG; channel::Union{Int64, String})
    labels = eeg.eeg_header[:labels]
    if typeof(channel) == String
        channel_idx = nothing
        for idx in 1:length(labels)
            if channel == labels[idx]
                channel_idx = idx
            end
        end
        if channel_idx === nothing
            throw(ArgumentError("Channel name does not match signal labels."))
        end
        return channel_idx
    else
        if channel < 1 || channel > length(labels)
            throw(ArgumentError("Channel index does not match signal channels."))
        end
        return labels[channel]
    end
end

"""
    eeg_rename_channel(eeg; channel, new_name)

Renames the `eeg` `channel`.

# Arguments

- `eeg::EEG`
- `channel::Union{Int64, String}`
- `new_name::String`

# Returns

- `eeg::EEG`
"""
function eeg_rename_channel(eeg::EEG; channel::Union{Int64, String}, new_name::String)
    # create new dataset
    eeg_new = deepcopy(eeg)
    labels = eeg_new.eeg_header[:labels]
    
    if typeof(channel) == String
        channel_found = nothing
        for idx in 1:length(labels)
            if channel == labels[idx]
                labels[idx] = new_name
                channel_found = idx
            end
        end
        if channel_found === nothing
            throw(ArgumentError("Channel name does not match signal labels."))
        end
    else
        if channel < 1 || channel > length(labels)
            throw(ArgumentError("Channel index does not match signal channels."))
        else
            labels[channel] = new_name
        end
    end
    eeg_new.eeg_header[:labels] = labels
    
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_rename_channel(EEG, channel=$channel, new_name=$new_name)")

    return eeg_new
end

"""
    eeg_rename_channel!(eeg; channel, new_name)

Renames the `eeg` `channel`.

# Arguments

- `eeg::EEG`
- `channel::Union{Int64, String}`
- `new_name::String`
"""
function eeg_rename_channel!(eeg::EEG; channel::Union{Int64, String}, new_name::String)
    labels = eeg.eeg_header[:labels]
    
    if typeof(channel) == String
        channel_found = nothing
        for idx in 1:length(labels)
            if channel == labels[idx]
                labels[idx] = new_name
                channel_found = idx
            end
        end
        if channel_found === nothing
            throw(ArgumentError("Channel name does not match signal labels."))
        end
    else
        if channel < 1 || channel > length(labels)
            throw(ArgumentError("Channel index does not match signal channels."))
        else
            labels[channel] = new_name
        end
    end
    eeg.eeg_header[:labels] = labels
    
    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_rename_channel!(EEG, channel=$channel, new_name=$new_name)")

    return
end

"""
    eeg_taper(eeg; taper)

Taper `eeg` with `taper`.

# Arguments

- `eeg::EEG`
- `taper::Vector`

# Returns

- `eeg::EEG`
"""
function eeg_taper(eeg::EEG; taper::Vector)
    s_tapered = signal_taper(eeg.eeg_signals, taper=taper)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_tapered, deepcopy(eeg.eeg_components))
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_taper(EEG, taper=$taper)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_taper!(eeg; taper)

Taper `eeg` with `taper`.

# Arguments

- `eeg::EEG`
- `taper::Vector`
"""
function eeg_taper!(eeg::EEG; taper::Vector)
    eeg.eeg_signals = signal_taper(eeg.eeg_signals, taper=taper)

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_taper!(EEG, taper=$taper)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_demean(eeg)

Removes mean value (DC offset).

# Arguments

- `eeg::EEG`

# Returns

- `eeg::EEG`
"""
function eeg_demean(eeg::EEG)
    s_demeaned = signal_demean(eeg.eeg_signals)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_demeaned, deepcopy(eeg.eeg_components))
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_demean(EEG)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_demean!(eeg)

Removes mean value (DC offset).

# Arguments

- `eeg::EEG`
"""
function eeg_demean!(eeg::EEG)
    eeg.eeg_signals = signal_demean(eeg.eeg_signals)

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_demean!(EEG)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_normalize_zscore(eeg)

Normalize by z-score.

# Arguments

- `eeg::EEG` - EEG object

# Returns

- `eeg::EEG`
"""
function eeg_normalize_zscore(eeg::EEG)
    s_normalized = signal_normalize_zscore(eeg.eeg_signals)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_normalized, deepcopy(eeg.eeg_components))
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_normalize_zscore(EEG)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_normalize_zscore!(eeg)

Normalize by z-score.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_normalize_zscore!(eeg::EEG)
    eeg.eeg_signals = signal_normalize_zscore(eeg.eeg_signals)

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_normalize_zscore!(EEG)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_normalize_minmax(eeg)

Normalize to 0...1

# Arguments

- `eeg::EEG`

# Returns

- `eeg::EEG`
"""
function eeg_normalize_minmax(eeg::EEG)
    s_normalized = signal_normalize_minmax(eeg.eeg_signals)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_normalized, deepcopy(eeg.eeg_components))
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_normalize_minmax(EEG)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_normalize_minmax!(eeg)

Normalize to 0...1

# Arguments

- `eeg::EEG`
"""
function eeg_normalize_minmax!(eeg::EEG)
    eeg.eeg_signals = signal_normalize_minmax(eeg.eeg_signals)

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_normalize_minmax!(EEG)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_extract_channel(eeg; channel)

Extracts `channel` number or name.

# Arguments

- `eeg::EEG`
- `channel::Union{Int64, String}`

# Returns

- `channel::Vector{Float64}`
"""
function eeg_extract_channel(eeg::EEG; channel::Union{Int64, String})
    labels = eeg.eeg_header[:labels]
    if typeof(channel) == String
        channel_idx = nothing
        for idx in 1:length(labels)
            if channel == labels[idx]
                channel_idx = idx
            end
        end
        if channel_idx === nothing
            throw(ArgumentError("Channel name does not match signal labels."))
        end
    else
        if channel < 1 || channel > length(labels)
            throw(ArgumentError("Channel index does not match signal channels."))
        end
        channel_idx = channel
    end    
    eeg_channel = vec(eeg.eeg_signals[channel_idx, :, :])

    return eeg_channel
end

"""
    eeg_cov(eeg; norm=true)

Calculates covariance between all channels of `eeg`.

# Arguments

- `eeg::EEG`
- `norm::Bool` - normalize covariance

# Returns

- `cov_mat::Array{Float64, 3}`
"""
function eeg_cov(eeg::EEG; norm=true)
    cov_mat = signal_cov(eeg.eeg_signals, norm=norm)

    return cov_mat
end

"""
    eeg_cov!(eeg; norm=true)

Calculates covariance between all channels of `eeg`.

# Arguments

- `eeg::EEG`
- `norm::Bool` - normalize covariance
"""
function eeg_cov!(eeg::EEG; norm=true)
    :cov_mat in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:cov_mat)
    cov_mat = eeg_cov(eeg, norm=norm)
    push!(eeg.eeg_components, cov_mat)
    push!(eeg.eeg_header[:components], :cov_mat)
    push!(eeg.eeg_header[:history], "eeg_cov!(EEG)")

    return
end

"""
    eeg_cor(eeg)

Calculates correlation coefficients between all channels of `eeg`.

# Arguments

- `eeg::EEG`

# Returns

- `cov_mat::Array{Float64, 3}`
"""
function eeg_cor(eeg::EEG)
    cor_mat = signal_cor(eeg.eeg_signals)

    return cor_mat
end

"""
    eeg_cor!(eeg)

Calculates correlation coefficients between all channels of `eeg`.

# Arguments

- `eeg::EEG`
"""
function eeg_cor!(eeg::EEG)
    :cor_mat in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:cor_mat)
    cor_mat = eeg_cor(eeg)
    push!(eeg.eeg_components, cor_mat)
    push!(eeg.eeg_header[:components], :cor_mat)
    push!(eeg.eeg_header[:history], "eeg_cor!(EEG)")

    return
end

"""
    eeg_upsample(eeg; new_sr)

Upsamples all channels of `eeg` to `new_sr` sampling frequency.

# Arguments

- `eeg::EEG`
- `new_sr::Int64` - new sampling rate

# Returns

- `eeg::EEG`
"""
function eeg_upsample(eeg::EEG; new_sr::Int64)
    t = eeg.eeg_time[1]:(1 / eeg.eeg_header[:sampling_rate][1]):eeg.eeg_time[end]
    s_upsampled, t_upsampled = signal_upsample(eeg.eeg_signals, t=t, new_sr=new_sr)

    # create new dataset
    eeg_time = collect(t_upsampled)
    eeg_new = EEG(deepcopy(eeg.eeg_header), eeg_time, s_upsampled, deepcopy(eeg.eeg_components))
    eeg_new.eeg_header[:eeg_duration_samples] = size(s_upsampled, 2) * size(s_upsampled, 3)
    eeg_new.eeg_header[:eeg_duration_seconds] = (size(s_upsampled, 2) * size(s_upsampled, 3)) / new_sr
    eeg_new.eeg_header[:epoch_duration_samples] = size(s_upsampled, 2)
    eeg_new.eeg_header[:epoch_duration_seconds] = size(s_upsampled, 2) / new_sr
    eeg_new.eeg_header[:sampling_rate] = repeat([new_sr], eeg_new.eeg_header[:channel_n])

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_upsample(EEG, new_sr=$new_sr)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_upsample!(eeg; new_sr)

Upsamples all channels of `eeg` to `new_sr` sampling frequency.

# Arguments

- `eeg::EEG`
- `new_sr::Int64` - new sampling rate

# Returns

- `eeg::EEG`
"""
function eeg_upsample!(eeg::EEG; new_sr::Int64)
    t = eeg.eeg_time[1]:(1 / eeg.eeg_header[:sampling_rate][1]):eeg.eeg_time[end]
    eeg.eeg_signals, t_upsampled = signal_upsample(eeg.eeg_signals, t=t, new_sr=new_sr)

    eeg.eeg_time = collect(t_upsampled)
    eeg.eeg_header[:eeg_duration_samples] = size(eeg.eeg_signals, 2) * size(eeg.eeg_signals, 3)
    eeg.eeg_header[:eeg_duration_seconds] = (size(eeg.eeg_signals, 2) * size(eeg.eeg_signals, 3)) / new_sr
    eeg.eeg_header[:epoch_duration_samples] = size(eeg.eeg_signals, 2)
    eeg.eeg_header[:epoch_duration_seconds] = size(eeg.eeg_signals, 2) / new_sr
    eeg.eeg_header[:sampling_rate] = repeat([new_sr], eeg.eeg_header[:channel_n])

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_upsample!(EEG, new_sr=$new_sr)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_history(eeg)

Shows processing history.

# Arguments

- `eeg::EEG`

# Returns

- `eeg::EEG`
"""
function eeg_history(eeg::EEG)
    return eeg.eeg_header[:history]
end

"""
    eeg_labels(eeg)

Returns labels.

# Arguments

- `eeg::EEG`

# Returns

- `eeg::EEG`
"""
function eeg_labels(eeg::EEG)
    return eeg.eeg_header[:labels]
end

"""
    eeg_sr(eeg)

Returns sampling rate.

# Arguments

- `eeg::EEG`

# Returns

- `eeg::EEG`
"""
function eeg_sr(eeg::EEG)
    return eeg.eeg_header[:sampling_rate][1]
end

"""
    eeg_info(eeg)

Shows info.

# Arguments

- `eeg::EEG`
"""
function eeg_info(eeg::EEG)
    println("          EEG file name: $(eeg.eeg_header[:eeg_filename])")
    println("          EEG size [Mb]: $(eeg.eeg_header[:eeg_filesize_mb])")
    println("     Sampling rate (Hz): $(eeg.eeg_header[:sampling_rate][1])")
    println("Signal length (samples): $(eeg.eeg_header[:eeg_duration_samples])")
    println("Signal length (seconds): $(round(eeg.eeg_header[:eeg_duration_seconds], digits=2))")
    println("     Number of channels: $(eeg.eeg_header[:channel_n])")
    println("       Number of epochs: $(eeg.eeg_header[:epoch_n])")
    println(" Epoch length (samples): $(eeg.eeg_header[:epoch_duration_samples])")
    println(" Epoch length (seconds): $(round(eeg.eeg_header[:epoch_duration_seconds], digits=2))")
    if eeg.eeg_header[:reference] == ""
        println("         Reference type: unknown")
    else
        println("         Reference type: $(eeg.eeg_header[:reference])")
    end
    if length(eeg.eeg_header[:labels]) == 0
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
        print("             Components: no")
    end

end

"""
    eeg_epochs(eeg; epoch_n=nothing, epoch_len=nothing, average=false)

Splits `eeg` into epochs.

# Arguments

- `eeg::EEG`
- `epoch_n::Union{Int64, Nothing}` - number of epochs
- `epoch_len::Union{Int64, Nothing}` - epoch length in samples
- `average::Bool` - average all epochs, returns one averaged epoch; if false than returns array of epochs, each row is one epoch

# Returns

- `eeg::EEG`
"""
function eeg_epochs(eeg::EEG; epoch_n::Union{Int64, Nothing}=nothing, epoch_len::Union{Int64, Nothing}=nothing, average::Bool=false)
    # unsplit epochs
    s_merged = reshape(eeg.eeg_signals,
                       size(eeg.eeg_signals, 1),
                       size(eeg.eeg_signals, 2) * size(eeg.eeg_signals, 3))
    
    # split into epochs
    s_split = signal_epochs(s_merged, epoch_n=epoch_n, epoch_len=epoch_len, average=average)

    # convert into Array{Float64, 3}
    s_split = reshape(s_split, size(s_split, 1), size(s_split, 2), size(s_split, 3))

    # create new dataset
    epoch_n = size(s_split, 3)
    epoch_duration_samples = size(s_split, 2)
    epoch_duration_seconds = size(s_split, 2) / eeg.eeg_header[:sampling_rate][1]
    eeg_duration_samples = size(s_split, 2) * size(s_split, 3)
    eeg_duration_seconds = eeg_duration_samples / eeg.eeg_header[:sampling_rate][1]
    eeg_time = collect(0:(1 / eeg.eeg_header[:sampling_rate][1]):epoch_duration_seconds)
    eeg_time = eeg_time[1:(end - 1)]
    eeg_new = EEG(deepcopy(eeg.eeg_header), eeg_time, s_split, deepcopy(eeg.eeg_components))
    eeg_new.eeg_header[:eeg_duration_samples] = eeg_duration_samples
    eeg_new.eeg_header[:eeg_duration_seconds] = eeg_duration_seconds
    eeg_new.eeg_header[:epoch_n] = epoch_n
    eeg_new.eeg_header[:epoch_duration_samples] = epoch_duration_samples
    eeg_new.eeg_header[:epoch_duration_seconds] = epoch_duration_seconds

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_epochs(EEG, epoch_n=$epoch_n, epoch_len=$epoch_len, average=$average)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_epochs!(eeg; epoch_n=nothing, epoch_len=nothing, average=false)

Splits `eeg` into epochs.

# Arguments

- `eeg::EEG`
- `epoch_n::Union{Int64, Nothing}` - number of epochs
- `epoch_len::Union{Int64, Nothing}` - epoch length in samples
- `average::Bool` - average all epochs, returns one averaged epoch; if false than returns array of epochs, each row is one epoch
"""
function eeg_epochs!(eeg::EEG; epoch_n::Union{Int64, Nothing}=nothing, epoch_len::Union{Int64, Nothing}=nothing, average::Bool=false)
    # unsplit epochs
    s_merged = reshape(eeg.eeg_signals,
                       size(eeg.eeg_signals, 1),
                       size(eeg.eeg_signals, 2) * size(eeg.eeg_signals, 3))
    
    # split into epochs
    s_split = signal_epochs(s_merged, epoch_n=epoch_n, epoch_len=epoch_len, average=average)

    # convert into Array{Float64, 3}
    s_split = reshape(s_split, size(s_split, 1), size(s_split, 2), size(s_split, 3))

    # create new dataset
    epoch_n = size(s_split, 3)
    epoch_duration_samples = size(s_split, 2)
    epoch_duration_seconds = size(s_split, 2) / eeg.eeg_header[:sampling_rate][1]
    eeg_duration_samples = size(s_split, 2) * size(s_split, 3)
    eeg_duration_seconds = eeg_duration_samples / eeg.eeg_header[:sampling_rate][1]
    eeg_time = collect(0:(1 / eeg.eeg_header[:sampling_rate][1]):epoch_duration_seconds)
    eeg_time = eeg_time[1:(end - 1)]
    eeg.eeg_time = eeg_time
    eeg.eeg_signals = s_split
    eeg.eeg_header[:eeg_duration_samples] = eeg_duration_samples
    eeg.eeg_header[:eeg_duration_seconds] = eeg_duration_seconds
    eeg.eeg_header[:epoch_n] = epoch_n
    eeg.eeg_header[:epoch_duration_samples] = epoch_duration_samples
    eeg.eeg_header[:epoch_duration_seconds] = epoch_duration_seconds

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_epochs!(EEG, epoch_n=$epoch_n, epoch_len=$epoch_len, average=$average)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_extract_epoch(eeg; epoch)

Extracts the `epoch` epoch.

# Arguments

- `eeg::EEG`
- `epoch::Int64` - epoch index

# Returns

- `eeg::EEG`
"""
function eeg_extract_epoch(eeg::EEG; epoch::Int64)
    if epoch < 1 || epoch > eeg.eeg_header[:epoch_n]
        throw(ArgumentError("Epoch index out of range."))
    end

    s_new = reshape(eeg.eeg_signals[:, :, epoch], size(eeg.eeg_signals, 1), size(eeg.eeg_signals, 2), 1)
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_new, deepcopy(eeg.eeg_components))
    eeg_new.eeg_header[:epoch_n] = 1
    eeg_new.eeg_header[:eeg_duration_samples] = eeg_new.eeg_header[:epoch_duration_samples]
    eeg_new.eeg_header[:eeg_duration_seconds] = eeg_new.eeg_header[:epoch_duration_seconds]

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_get_epoch(EEG, epoch=$epoch)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_tconv(eeg; kernel)

Performs convolution in the time domain.

# Arguments

- `eeg::EEG`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `eeg::EEG`
"""
function eeg_tconv(eeg::EEG; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})
    s_convoluted = signal_tconv(eeg.eeg_signals, kernel=kernel)

    ## EEG signal can only store Float64
    typeof(kernel) == Vector{ComplexF64} && (s_convoluted = abs.(s_convoluted))

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_convoluted, deepcopy(eeg.eeg_components))
    eeg_new.eeg_header[:eeg_duration_samples] = size(s_convoluted, 2) * size(s_convoluted, 3)
    eeg_new.eeg_header[:eeg_duration_seconds] = (size(s_convoluted, 2) * size(s_convoluted, 3)) / eeg_sr(eeg_new)
    eeg_new.eeg_header[:epoch_duration_samples] = size(s_convoluted, 2)
    eeg_new.eeg_header[:epoch_duration_seconds] = size(s_convoluted, 2) / eeg_sr(eeg_new)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_tconv(EEG, kernel=$kernel)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_tconv!(eeg; kernel)

Performs convolution in the time domain.

# Arguments

- `eeg::EEG`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`
"""
function eeg_tconv!(eeg::EEG; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

    # EEG signal can only store Float64
    if typeof(kernel) == Vector{ComplexF64}
        eeg.eeg_signals = abs.(signal_tconv(eeg.eeg_signals, kernel=kernel))
    else
        eeg.eeg_signals = signal_tconv(eeg.eeg_signals, kernel=kernel)
    end

    # create new dataset
    eeg.eeg_header[:eeg_duration_samples] = size(s_convoluted, 2) * size(s_convoluted, 3)
    eeg.eeg_header[:eeg_duration_seconds] = (size(s_convoluted, 2) * size(s_convoluted, 3)) / eeg_sr(eeg)
    eeg.eeg_header[:epoch_duration_samples] = size(s_convoluted, 2)
    eeg.eeg_header[:epoch_duration_seconds] = size(s_convoluted, 2) / eeg_sr(eeg)
    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_tconv!(EEG, kernel=$kernel)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_filter(eeg; fprototype, ftype=nothing, cutoff=0, order=-1, rp=-1, rs=-1, dir=:twopass, d=1, t=0, window=nothing)

Filters `eeg` using zero phase distortion filter.

# Arguments

- `eeg::EEG`
- `fprototype::Symbol[:mavg, :mmed, :poly, :butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir]` - filter prototype:
    - `:mavg` - moving average (with threshold and/or weight window)
    - `:mmed` - moving median (with threshold and/or weight window)
    - `:poly` - polynomial of `order` order
- `ftype::Symbol[:lp, :hp, :bp, :bs]` - filter type
- `cutoff::Union{Int64, Float64, Tuple}` - filter cutoff in Hz (vector for `:bp` and `:bs`)
- `order::Int64` - filter order
- `rp::Union{Int64, Float64}` - dB ripple in the passband
- `rs::Union{Int64, Float64}` - dB attentuation in the stopband
- `dir:Symbol[:onepass, :onepass_reverse, :twopass]` - filter direction
- `d::Int64` - window length for mean average and median average filter
- `t::Union{Int64, Float64}` - threshold for :mavg and :mmed filters; threshold = threshold * std(signal) + mean(signal) for :mavg or threshold = threshold * std(signal) + median(signal) for :mmed filter
- `window::Union{Vector{Float64}, Nothing} - window, required for FIR filter

# Returns

- `eeg::EEG`
"""
function eeg_filter(eeg::EEG; fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Int64, Float64, Tuple}=0, order::Int64=0, rp::Union{Int64, Float64}=-1, rs::Union{Int64, Float64}=-1, dir::Symbol=:twopass, d::Int64=1, t::Union{Int64, Float64}=0, window::Union{Vector{Float64}, Nothing}=nothing)

    s_filtered = signal_filter(eeg.eeg_signals,
                                    fprototype=fprototype,
                                    ftype=ftype,
                                    cutoff=cutoff,
                                    fs=eeg_sr(eeg),
                                    order=order,
                                    rp=rp,
                                    rs=rs,
                                    dir=dir,
                                    d=d,
                                    t=t,
                                    window=window)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_filtered, deepcopy(eeg.eeg_components))
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_filter(EEG, fprototype=$fprototype, ftype=$ftype, cutoff=$cutoff, order=$order, rp=$rp, rs=$rs, dir=$dir, d=$d, window=$window)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_filter!(eeg; fprototype, ftype=nothing, cutoff, fs, order, rp, rs, dir=:twopass, d=1, window)

Filters `eeg` using zero phase distortion filter.

# Arguments

- `eeg::EEG`
- `fprototype::Symbol[:mavg, :mmed, :poly, :butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir]` - filter prototype:
    - `:mavg` - moving average (with threshold and/or weight window)
    - `:mmed` - moving median (with threshold and/or weight window)
    - `:poly` - polynomial of `order` order
- `ftype::Union{Symbol[:lp, :hp, :bp, :bs], Nothing}` - filter type
- `cutoff::Union{Int64, Float64, Tuple}` - filter cutoff in Hz (vector for `:bp` and `:bs`)
- `order::Int64` - filter order
- `rp::Union{Int64, Float64}` - dB ripple in the passband
- `rs::Union{Int64, Float64}` - dB attentuation in the stopband
- `dir:Symbol[:onepass, :onepass_reverse, :twopass]` - filter direction
- `d::Int64` - window length for mean average and median average filter
- `t::Union{Int64, Float64}` - threshold for :mavg and :mmed filters; threshold = threshold * std(signal) + mean(signal) for :mavg or threshold = threshold * std(signal) + median(signal) for :mmed filter
- `window::Union{Vector{Float64}, Nothing} - window, required for FIR filter
"""
function eeg_filter!(eeg::EEG; fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Int64, Float64, Tuple}=0, order::Int64=0, rp::Union{Int64, Float64}=-1, rs::Union{Int64, Float64}=-1, dir::Symbol=:twopass, d::Int64=1, t::Union{Int64, Float64}=0, window::Union{Vector{Float64}, Nothing}=nothing)

    eeg.eeg_signals = signal_filter(eeg.eeg_signals,
                                    fprototype=fprototype,
                                    ftype=ftype,
                                    cutoff=cutoff,
                                    fs=eeg_sr(eeg),
                                    order=order,
                                    rp=rp,
                                    rs=rs,
                                    dir=dir,
                                    d=d,
                                    t=t,
                                    window=window)

    # create new dataset
    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_filter!(EEG, fprototype=$fprototype, ftype=$ftype, cutoff=$cutoff, order=$order, rp=$rp, rs=$rs, dir=$dir, d=$d, window=$window)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_downsample(eeg; new_sr)

Downsamples all channels of `eeg` to `new_sr` sampling frequency.

# Arguments

- `eeg::EEG`
- `new_sr::Int64` - new sampling rate

# Returns

- `eeg::EEG`
"""
function eeg_downsample(eeg::EEG; new_sr::Int64)
    t = eeg.eeg_time[1]:(1 / eeg.eeg_header[:sampling_rate][1]):eeg.eeg_time[end]
    s_downsampled, t_downsampled = signal_downsample(eeg.eeg_signals, t=t, new_sr=new_sr)

    # create new dataset
    eeg_time = collect(t_downsampled)
    eeg_new = EEG(deepcopy(eeg.eeg_header), eeg_time, s_downsampled, deepcopy(eeg.eeg_components))
    eeg_new.eeg_header[:eeg_duration_samples] = size(s_downsampled, 2) * size(s_downsampled, 3)
    eeg_new.eeg_header[:eeg_duration_seconds] = (size(s_downsampled, 2) * size(s_downsampled, 3)) / new_sr
    eeg_new.eeg_header[:epoch_duration_samples] = size(s_downsampled, 2)
    eeg_new.eeg_header[:epoch_duration_seconds] = size(s_downsampled, 2) / new_sr
    eeg_new.eeg_header[:sampling_rate] = repeat([new_sr], eeg_new.eeg_header[:channel_n])

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_downsample(EEG, new_sr=$new_sr)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_downsample!(eeg; new_sr)

Downsamples all channels of `eeg` to `new_sr` sampling frequency.

# Arguments

- `eeg::EEG`
- `new_sr::Int64` - new sampling rate
"""
function eeg_downsample!(eeg::EEG; new_sr::Int64)
    t = eeg.eeg_time[1]:(1 / eeg.eeg_header[:sampling_rate][1]):eeg.eeg_time[end]
    eeg.eeg_signals, t_downsampled = signal_downsample(eeg.eeg_signals, t=t, new_sr=new_sr)

    # create new dataset
    eeg.eeg_time = collect(t_downsampled)
    eeg.eeg_header[:eeg_duration_samples] = size(eeg.eeg_signals, 2) * size(eeg.eeg_signals, 3)
    eeg.eeg_header[:eeg_duration_seconds] = (size(eeg.eeg_signals, 2) * size(eeg.eeg_signals, 3)) / new_sr
    eeg.eeg_header[:epoch_duration_samples] = size(eeg.eeg_signals, 2)
    eeg.eeg_header[:epoch_duration_seconds] = size(eeg.eeg_signals, 2) / new_sr
    eeg.eeg_header[:sampling_rate] = repeat([new_sr], eeg.eeg_header[:channel_n])

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_downsample!(EEG, new_sr=$new_sr)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_autocov(eeg; lag=1, demean=false, norm=false)

Calculates autocovariance of each the `eeg` channels.

# Arguments

- `eeg::EEG`
- `lag::Int64` - lags range is `-lag:lag`
- `demean::Bool` - demean signal prior to analysis
- `norm::Bool` - normalize autocovariance

# Returns

- `acov::Matrix{Float64}`
- `lags::Vector{Float64}
"""
function eeg_autocov(eeg::EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)
    acov, lags = signal_autocov(eeg.eeg_signals, lag=lag, demean=demean, norm=norm)
    size(acov, 3) == 1 && (acov = reshape(acov, size(acov, 1), size(acov, 2)))
    lags = (eeg.eeg_time[2] - eeg.eeg_time[1]) .* collect(-lag:lag)

    return acov, lags
end

"""
    eeg_autocov!(eeg; lag=1, demean=false, norm=false)

Calculates autocovariance of each the `eeg` channels.

# Arguments

- `eeg::EEG`
- `lag::Int64` - lags range is `-lag:lag`
- `demean::Bool` - demean signal prior to analysis
- `norm::Bool` - normalize autocovariance
"""
function eeg_autocov!(eeg::EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)
    :acov in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:acov)
    :acov_lags in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:acov_lags)
    acov, lags = eeg_autocov(eeg, lag=lag, demean=demean, norm=norm)
    push!(eeg.eeg_components, acov)
    push!(eeg.eeg_components, lags)
    push!(eeg.eeg_header[:components], :acov)
    push!(eeg.eeg_header[:components], :acov_lags)
    push!(eeg.eeg_header[:history], "eeg_autocov!(EEG, lag=$lag, demean=$demean, norm=$norm)")

    return
end

"""
    eeg_crosscov(eeg; lag=1, demean=false, norm=false)

Calculates cross-covariance of each the `eeg` channels.

# Arguments

- `eeg::EEG`
- `lag::Int64` - lags range is `-lag:lag`
- `demean::Bool` - demean signal prior to analysis
- `norm::Bool` - normalize cross-covariance

# Returns

- `ccov::Matrix{Float64}`
- `lags::Vector{Float64}
"""
function eeg_crosscov(eeg::EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)
    ccov, lags = signal_crosscov(eeg.eeg_signals, lag=lag, demean=demean, norm=norm)
    size(ccov, 3) == 1 && (ccov = reshape(ccov, size(ccov, 1), size(ccov, 2)))
    lags = (eeg.eeg_time[2] - eeg.eeg_time[1]) .* collect(-lag:lag)

    return ccov, lags
end

"""
    eeg_crosscov!(eeg; lag=1, demean=false, norm=false)

Calculates cross-covariance of each the `eeg` channels.

# Arguments

- `eeg::EEG`
- `lag::Int64` - lags range is `-lag:lag`
- `demean::Bool` - demean signal prior to analysis
- `norm::Bool` - normalize cross-covariance
"""
function eeg_crosscov!(eeg::EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)
    :ccov in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:ccov)
    :ccov_lags in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:ccov_lags)
    ccov, lags = eeg_crosscov(eeg, lag=lag, demean=demean, norm=norm)
    push!(eeg.eeg_components, ccov)
    push!(eeg.eeg_components, lags)
    push!(eeg.eeg_header[:components], :ccov)
    push!(eeg.eeg_header[:components], :ccov_lags)
    push!(eeg.eeg_header[:history], "eeg_crosscov!(EEG, lag=$lag, demean=$demean, norm=$norm)")

    return
end

"""
    eeg_crosscov(eeg1, eeg1; lag=1, demean=false, norm=false)

Calculates cross-covariance between same channels in `eeg1` and `eeg2`.

# Arguments

- `eeg1::EEG`
- `eeg2::EEG`
- `lag::Int64` - lags range is `-lag:lag`
- `demean::Bool` - demean signal prior to analysis
- `norm::Bool` - normalize cross-covariance

# Returns

- `ccov::Matrix{Float64}`
- `lags::Vector{Float64}
"""
function eeg_crosscov(eeg1::EEG, eeg2::EEG; lag::Int64=1, demean::Bool=false, norm::Bool=false)
    ccov, lags = signal_crosscov(eeg1.eeg_signals, eeg2.eeg_signals, lag=lag, demean=demean, norm=norm)
    size(ccov, 3) == 1 && (ccov = reshape(ccov, size(ccov, 1), size(ccov, 2)))
    lags = (eeg1.eeg_time[2] - eeg1.eeg_time[1]) .* collect(-lag:lag)

    return ccov, lags
end

"""
    eeg_psd(eeg; norm=false)

Calculates total power for each the `eeg` channels.

# Arguments

- `eeg::EEG`
- `norm::Bool` - normalize do dB

# Returns

- `powers::Array{Float64, 3}`
- `frequencies::Array{Float64, 3}`
"""
function eeg_psd(eeg::EEG; norm::Bool=false)
    s_psd_powers, s_psd_frequencies = signal_psd(eeg.eeg_signals, fs=eeg_sr(eeg), norm=norm)
    size(s_psd_powers, 3) == 1 && (s_psd_powers = reshape(s_psd_powers, size(s_psd_powers, 1), size(s_psd_powers, 2)))
    size(s_psd_frequencies, 3) == 1 && (s_psd_frequencies = reshape(s_psd_frequencies, size(s_psd_frequencies, 1), size(s_psd_frequencies, 2)))

    return s_psd_powers, s_psd_frequencies
end

"""
    eeg_psd!(eeg; norm=false)

Calculates total power for each the `eeg` channels.

# Arguments

- `eeg::EEG`
- `norm::Bool` - normalize do dB
"""
function eeg_psd!(eeg::EEG; norm::Bool=false)
    :psd_powers in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:psd_powers)
    :psd_frequencies in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:psd_frequencies)
    s_psd_powers, s_psd_frequencies = eeg_psd(eeg, norm=norm)
    push!(eeg.eeg_components, s_psd_powers)
    push!(eeg.eeg_components, s_psd_frequencies)
    push!(eeg.eeg_header[:components], :psd_powers)
    push!(eeg.eeg_header[:components], :psd_frequencies)
    push!(eeg.eeg_header[:history], "eeg_psd!(EEG, norm=$norm)")

    return
end

"""
    eeg_stationarity(eeg:EEG; window=10, method=:euclid)

Calculates stationarity.

# Arguments

- `eeg:EEG`
- `window::Int64` - time window in samples
- `method::Symbol[:mean, :var, :euclid, :hilbert]

# Returns

- `stationarity::Union{Matrix{Float64}, Array{Float64, 3}}`

"""
function eeg_stationarity(eeg::EEG; window::Int64=10, method::Symbol=:hilbert)
    s_stationarity = signal_stationarity(eeg.eeg_signals, window=window, method=method)

    return s_stationarity
end

"""
    eeg_stationarity!(eeg:EEG; window=10, method=:euclid)

Calculates stationarity.

# Arguments

- `eeg:EEG`
- `window::Int64` - time window in samples
- `method::Symbol[:mean, :var, :euclid, :hilbert]
"""
function eeg_stationarity!(eeg::EEG; window::Int64=10, method::Symbol=:hilbert)
    :stationarity in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:stationarity)
    s_stationarity = eeg_stationarity(eeg, window=window, method=method)
    push!(eeg.eeg_components, s_stationarity)
    push!(eeg.eeg_header[:components], :stationarity)
    push!(eeg.eeg_header[:history], "eeg_stationarity!(EEG, window=$window, method=$method)")

    return
end

"""
    eeg_trim(eeg:EEG; len, offset=0, from=:start)

Removes `len` samples from the beginning (`from` = :start, default) or end (`from` = :end) of the `signal`.

# Arguments

- `eeg:EEG`
- `len::Int64` - number of samples to remove
- `offset::Int64` - offset from which trimming starts, only works for `from` = :start
- `from::Symbol[:start, :end]`

# Returns

- `eeg:EEG`
"""
function eeg_trim(eeg::EEG; len::Int64, offset::Int64=0, from::Symbol=:start)
    # create new dataset
    eeg_signals = deepcopy(eeg.eeg_signals)
    eeg_time = deepcopy(eeg.eeg_time)
    eeg_signals = signal_trim(eeg_signals, len=len, from=from)
    eeg_time = signal_trim(eeg_time, len=len, from=from)

    eeg_trimmed = EEG(deepcopy(eeg.eeg_header), eeg_time, eeg_signals, deepcopy(eeg.eeg_components))
    eeg_trimmed.eeg_header[:eeg_duration_samples] -= len
    eeg_trimmed.eeg_header[:eeg_duration_seconds] -= len * (1 / eeg_sr(eeg))
    eeg_trimmed.eeg_header[:epoch_duration_samples] -= len
    eeg_trimmed.eeg_header[:epoch_duration_seconds] -= len * (1 / eeg_sr(eeg))

    # add entry to :history field
    push!(eeg_trimmed.eeg_header[:history], "eeg_trim(EEG, len=$len, from=$from)")

    eeg_reset_components!(eeg_trimmed)

    return eeg_trimmed
end

"""
    eeg_trim!(eeg:EEG; len, offset=0, from=:start)

Removes `len` samples from the beginning (`from` = :start, default) or end (`from` = :end) of the `signal`.

# Arguments

- `eeg:EEG`
- `len::Int64` - number of samples to remove
- `offset::Int64` - offset from which trimming starts, only works for `from` = :start
- `from::Symbol[:start, :end]`

"""
function eeg_trim!(eeg::EEG; len::Int64, offset::Int64=0, from::Symbol=:start)
    eeg.eeg_signals = signal_trim(eeg.eeg_signals, len=len, from=from)
    eeg.eeg_time = signal_trim(eeg.eeg_time, len=len, from=from)

    eeg.eeg_header[:eeg_duration_samples] -= len
    eeg.eeg_header[:eeg_duration_seconds] -= len * (1 / eeg_sr(eeg))
    eeg.eeg_header[:epoch_duration_samples] -= len
    eeg.eeg_header[:epoch_duration_seconds] -= len * (1 / eeg_sr(eeg))

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_trim!(EEG, len=$len, from=$from)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_mi(eeg)

Calculates mutual information between all channels of `eeg`.

# Arguments

- `eeg::EEG`

# Returns

- `mi::Array{Float64, 3}`
"""
function eeg_mi(eeg::EEG)
    mi = signal_mi(eeg.eeg_signals)
    size(mi, 3) == 1 && (mi = reshape(mi, size(mi, 1), size(mi, 2)))

    return mi
end

"""
    eeg_mi!(eeg)

Calculates mutual information between all channels of `eeg`.

# Arguments

- `eeg::EEG`
"""
function eeg_mi!(eeg::EEG)
    :mi in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:mi)
    mi = signal_mi(eeg.eeg_signals)
    size(mi, 3) == 1 && (mi = reshape(mi, size(mi, 1), size(mi, 2)))
    push!(eeg.eeg_components, mi)
    push!(eeg.eeg_header[:components], :mi)
    push!(eeg.eeg_header[:history], "eeg_mi!(EEG)")

    return
end

"""
    eeg_mi(eeg1, eeg2)

Calculates mutual information between all channels of `eeg1` and `eeg2`.

# Arguments

- `eeg1::EEG`
- `eeg2::EEG`

# Returns

- `mi::Array{Float64, 3}`
"""
function eeg_mi(eeg1::EEG, eeg2::EEG)
    mi = signal_mi(eeg1.eeg_signals, eeg2.eeg_signals)
    size(mi, 3) == 1 && (mi = reshape(mi, size(mi, 1), size(mi, 2)))

    return mi
end

"""
    eeg_entropy(eeg)

Calculates entropy of all channels of `eeg`.

# Arguments

- `eeg::EEG`

# Returns

- `entropy::Matrix{Float64}`
"""
function eeg_entropy(eeg::EEG)
    ent = signal_entropy(eeg.eeg_signals)
    size(ent, 3) == 1 && (ent = reshape(ent, size(ent, 1), size(ent, 2)))

    return ent
end

"""
    eeg_entropy!(eeg)

Calculates entropy of all channels of `eeg1`.

# Arguments

- `eeg::EEG`
"""
function eeg_entropy!(eeg::EEG)
    :entropy in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:entropy)
    ent = eeg_entropy(eeg)
    push!(eeg.eeg_components, ent)
    push!(eeg.eeg_header[:components], :entropy)
    push!(eeg.eeg_header[:history], "eeg_entropy!(EEG)")

    return
end

"""
    eeg_band(band)

Return EEG band frequency limits.

# Arguments

- `band::Symbol`

# Returns

- `band_frequency::Tuple{Float64, Float64}`

"""
function eeg_band(band::Symbol)
    band in [:delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher] || throw(ArgumentError("Available bands: :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher."))

    band === :delta && (band_frequency = (0.5, 4.0))
    band === :theta && (band_frequency = (4.0, 8.0))
    band === :alpha && (band_frequency = (8.0, 13.0))
    band === :beta && (band_frequency = (14.0, 30.0))
    band === :beta_high && (band_frequency = (25.0, 30.0))
    band === :gamma && (band_frequency = (30.0, 150.0))
    band === :gamma_1 && (band_frequency = (31.0, 40.0))
    band === :gamma_2 && (band_frequency = (41.0, 50.0))
    band === :gamma_lower && (band_frequency = (30.0, 80.0))
    band === :gamma_higher && (band_frequency = (80.0, 150.0))
    
    return band_frequency
end

"""
    eeg_coherence(eeg1, eeg2)

Calculates coherence between all channels of `eeg1` and `eeg2`.

# Arguments

- `eeg1::EEG`
- `eeg2::EEG`

# Returns

- `coherence::Union{Matrix{Float64}, Array{ComplexF64, 3}}`
"""
function eeg_coherence(eeg1::EEG, eeg2::EEG)
    coherence = signal_coherence(eeg1.eeg_signals, eeg2.eeg_signals)
    size(coherence, 3) == 1 && (coherence = reshape(coherence, size(coherence, 1), size(coherence, 2)))

    return coherence
end

"""
    eeg_coherence(eeg; channel1, channel2, epoch1, epoch2)

Calculates coherence between `channel1`/`epoch1` and `channel2` of `epoch2` of `eeg`.

# Arguments

- `eeg::EEG`
- `channel1::Int64`
- `channel2::Int64`
- `epoch1::Int64`
- `epoch2::Int64`

# Returns

- `coherence::Vector{ComplexF64}`
"""
function eeg_coherence(eeg::EEG; channel1::Int64, channel2::Int64, epoch1::Int64, epoch2::Int64)
    (channel1 < 0 || channel2 < 0 || epoch1 < 0 || epoch2 < 0) && throw(ArgumentError("Channels/epochs indices must be > 0."))
    channel_n = eeg.eeg_header[:channel_n]
    epoch_n = eeg.eeg_header[:epoch_n]
    (channel1 > channel_n || channel2 > channel_n || epoch1 > epoch_n || epoch2 > epoch_n) && throw(ArgumentError("Channels/epochs indices out of range."))

    coherence = signal_coherence(eeg.eeg_signals[channel1, :, epoch1], eeg.eeg_signals[channel2, :, epoch2])

    return coherence
end

"""
    eeg_freqs(eeg)

Returns vector of frequencies and Nyquist frequency for `eeg`.

# Arguments

- `eeg::EEG`

# Returns

- `hz::Vector{Float64}`
- `nyquist::Float64`
"""
function eeg_freqs(eeg::EEG)
    hz, nyq = freqs(eeg.eeg_signals[1, :, 1], eeg_sr(eeg))

    return hz, nyq
end

"""
    eeg_freqs!(eeg)

Returns vector of frequencies and Nyquist frequency for `eeg`.

# Arguments

- `eeg::EEG`
"""
function eeg_freqs!(eeg::EEG)
    :hz in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:hz)
    :nyq in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:nyq)
    hz, nyq = eeg_freqs(eeg)
    push!(eeg.eeg_components, hz)
    push!(eeg.eeg_components, nyq)
    push!(eeg.eeg_header[:components], :hz)
    push!(eeg.eeg_header[:components], :nyq)
    push!(eeg.eeg_header[:history], "eeg_freqs!(EEG)")

    return
end

"""
    eeg_pca(eeg; n)

Calculates `n` first PCs for `eeg`.

# Arguments

- `eeg::EEG`
- `n::Int64` - number of PCs

# Returns

- `pc::Array{Float64, 3}:` - PC(1)..PC(n) × epoch
- `pc_var::Matrix{Float64}` - PC_VAR(1)..PC_VAR(n) × epoch
"""
function eeg_pca(eeg::EEG; n::Int64)
    pc, pc_var = signal_pca(eeg.eeg_signals, n=n)

    return pc, pc_var
end

"""
    eeg_pca!(eeg; n)

Calculates `n` first PCs for `eeg`.

# Arguments

- `eeg::EEG`
- `n::Int64` - number of PCs
"""
function eeg_pca!(eeg::EEG; n::Int64)
    :pc in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:pc)
    :pc_var in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:pc_var)
    pc, pc_var = eeg_pca(eeg, n=n)
    push!(eeg.eeg_components, pc)
    push!(eeg.eeg_components, pc_var)
    push!(eeg.eeg_header[:components], :pc)
    push!(eeg.eeg_header[:components], :pc_var)
    push!(eeg.eeg_header[:history], "eeg_pca!(EEG, n=$n)")

    return
end

"""
    eeg_difference(eeg1, eeg2; n=3, method=:absdiff)

Calculates mean difference and 95% confidence interval for `eeg1` and `eeg2`.

# Arguments

- `eeg1::EEG`
- `eeg2::EEG`
- `n::Int64` - number of bootstraps
- `method::Symbol[:absdiff, :diff2int]`
    - `:absdiff` - maximum difference
    - `:diff2int` - integrated area of the squared difference

# Returns

- `signals_statistic::Matrix{Float64}`
- `signals_statistic_single::Vector{Float64}`
- `p::Vector{Float64}`
"""
function eeg_difference(eeg1::EEG, eeg2::EEG; n::Int64=3, method::Symbol=:absdiff)
    epoch_n = size(eeg1.eeg_signals, 3)
    signals_statistic = zeros(epoch_n, size(eeg1.eeg_signals, 1) * n)
    signals_statistic_single = zeros(epoch_n)
    p = zeros(epoch_n)

    Threads.@threads for epoch in 1:epoch_n
        signals_statistic[epoch, :], signals_statistic_single[epoch], p[epoch] = signal_difference(eeg1.eeg_signals[:, :, epoch], eeg2.eeg_signals[:, :, epoch], n=n, method=method)
    end

    return signals_statistic, signals_statistic_single, p
end

"""
    eeg_fconv(eeg, kernel)

Performs convolution in the time domain.

# Arguments

- `eeg::EEG`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `eeg::EEG`
"""
function eeg_fconv(eeg::EEG; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})
    s_convoluted = signal_fconv(eeg.eeg_signals, kernel=kernel)

    ## EEG signal can only store Float64
    s_convoluted = abs.(s_convoluted)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), s_convoluted, deepcopy(eeg.eeg_components))
    eeg_new.eeg_header[:eeg_duration_samples] = size(s_convoluted, 2) * size(s_convoluted, 3)
    eeg_new.eeg_header[:eeg_duration_seconds] = (size(s_convoluted, 2) * size(s_convoluted, 3)) / eeg_sr(eeg_new)
    eeg_new.eeg_header[:epoch_duration_samples] = size(s_convoluted, 2)
    eeg_new.eeg_header[:epoch_duration_seconds] = size(s_convoluted, 2) / eeg_sr(eeg_new)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_fconv(EEG, kernel=$kernel)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_fconv!(eeg, kernel)

Performs convolution in the time domain.

# Arguments

- `eeg::EEG`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`
"""
function eeg_fconv!(eeg::EEG; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})
    # EEG signal can only store Float64
    eeg.eeg_signals = abs.(signal_fconv(eeg.eeg_signals, kernel=kernel))

    eeg.eeg_header[:eeg_duration_samples] = size(s_convoluted, 2) * size(s_convoluted, 3)
    eeg.eeg_header[:eeg_duration_seconds] = (size(s_convoluted, 2) * size(s_convoluted, 3)) / eeg_sr(eeg)
    eeg.eeg_header[:epoch_duration_samples] = size(s_convoluted, 2)
    eeg.eeg_header[:epoch_duration_seconds] = size(s_convoluted, 2) / eeg_sr(eeg)
    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_fconv!(EEG, kernel=$kernel)")

    eeg_reset_components!(eeg)

    return
end

"""
    eeg_edit_header(eeg; field, value)

Changes value of `eeg` `field` to `value`.

# Arguments

- `eeg::EEG`
- `field::Symbol`
- `value::Any`

# Returns

- `eeg:EEG`
"""
function eeg_edit_header(eeg::EEG; field::Symbol, value::Any)
  value === nothing && throw(ArgumentError("Value cannot be empty."))
  
  eeg_new = deepcopy(eeg)
  fields = keys(eeg_new.eeg_header)
  field in fields || throw(ArgumentError("Field does not exist."))
  typeof(eeg_new.eeg_header[field]) == typeof(value) || throw(ArgumentError("Field type does not mach value type."))
  eeg_new.eeg_header[field] = value
  # add entry to :history field
  push!(eeg_new.eeg_header[:history], "eeg_edit(EEG, field=$field, value=$value)")    

  return eeg_new
end

"""
    eeg_edit_header!(eeg; field, value)

Changes value of `eeg` `field` to `value`.

# Arguments

- `eeg::EEG`
- `field::Symbol`
- `value::Any`

# Returns

- `eeg:EEG`
"""
function eeg_edit_header!(eeg::EEG; field::Symbol, value::Any)
  value === nothing && throw(ArgumentError("Value cannot be empty."))
  
  fields = keys(eeg.eeg_header)
  field in fields || throw(ArgumentError("Field does not exist."))
  typeof(eeg.eeg_header[field]) == typeof(value) || throw(ArgumentError("Field type does not mach value type."))
  eeg.eeg_header[field] = value
  # add entry to :history field
  push!(eeg.eeg_header[:history], "eeg_edit!(EEG, field=$field, value=$value)")    

  return
end

"""
    eeg_show_header(eeg)

Shows keys and values of `eeg` header.

# Arguments

- `eeg::EEG`
"""
function eeg_show_header(eeg::EEG)
    for (key, value) in eeg.eeg_header
        println("$key: $value")
    end
end

"""
    eeg_delete_epoch(eeg; epoch)

Removes `epoch` from the `eeg`.

# Arguments

- `eeg::EEG`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}` - epoch index to be removed, vector of numbers or range

# Returns

- `eeg::EEG`
"""
function eeg_delete_epoch(eeg::EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange})
    size(eeg.eeg_signals, 3) == 1 && throw(ArgumentError("You cannot delete the last epoch."))

    if typeof(epoch) <: AbstractRange
        epoch = collect(epoch)
    end

    length(epoch) == size(eeg.eeg_signals, 3) && throw(ArgumentError("You cannot delete all epochs."))

    length(epoch) > 1 && (epoch = sort!(epoch, rev=true))

    if epoch[end] < 1 || epoch[1] > size(eeg.eeg_signals, 3)
        throw(ArgumentError("Epoch index does not match signal epochs."))
    end

    eeg_header = deepcopy(eeg.eeg_header)
    eeg_time = deepcopy(eeg.eeg_time)
    eeg_signals = deepcopy(eeg.eeg_signals)

    # remove epoch
    eeg_signals = eeg_signals[:, :, setdiff(1:end, (epoch))]

    # update headers
    eeg_header[:epoch_n] = eeg_header[:epoch_n] - length(epoch)
    epoch_n = eeg_header[:epoch_n]
    eeg_header[:eeg_duration_samples] = epoch_n * size(eeg.eeg_signals, 2)
    eeg_header[:eeg_duration_seconds] = round(epoch_n * (size(eeg.eeg_signals, 2) / eeg_sr(eeg)), digits=2)
    eeg_header[:epoch_duration_samples] = size(eeg.eeg_signals, 2)
    eeg_header[:epoch_duration_seconds] = round(size(eeg.eeg_signals, 2) / eeg_sr(eeg), digits=2)

    # create new dataset
    eeg_new = EEG(eeg_header, eeg_time, eeg_signals, deepcopy(eeg.eeg_components))
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_delete_epoch(EEG, $epoch)")
    
    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_delete_epoch!(eeg; epoch)

Removes `epoch` from the `eeg`.

# Arguments

- `eeg::EEG`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}` - epoch index to be removed, vector of numbers or range
"""
function eeg_delete_epoch!(eeg::EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange})
    size(eeg.eeg_signals, 3) == 1 && throw(ArgumentError("You cannot delete the last epoch."))

    if typeof(epoch) <: AbstractRange
        epoch = collect(epoch)
    end

    length(epoch) == size(eeg.eeg_signals, 3) && throw(ArgumentError("You cannot delete all epochs."))

    length(epoch) > 1 && (epoch = sort!(epoch, rev=true))

    if epoch[end] < 1 || epoch[1] > size(eeg.eeg_signals, 3)
        throw(ArgumentError("Epoch index does not match signal epochs."))
    end

      # remove epoch
    eeg.eeg_signals = eeg.eeg_signals[:, :, setdiff(1:end, (epoch))]

    # update headers
    eeg.eeg_header[:epoch_n] = eeg.eeg_header[:epoch_n] - length(epoch)
    epoch_n = eeg.eeg_header[:epoch_n]
    eeg.eeg_header[:eeg_duration_samples] = epoch_n * size(eeg.eeg_signals, 2)
    eeg.eeg_header[:eeg_duration_seconds] = round(epoch_n * (size(eeg.eeg_signals, 2) / eeg_sr(eeg)), digits=2)
    eeg.eeg_header[:epoch_duration_samples] = size(eeg.eeg_signals, 2)
    eeg.eeg_header[:epoch_duration_seconds] = round(size(eeg.eeg_signals, 2) / eeg_sr(eeg), digits=2)

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_delete_epoch!(EEG, $epoch)")
    
    eeg_reset_components!(eeg)

    return
end

"""
    eeg_keep_epoch(eeg; epoch)

Keeps `epoch` in the `eeg`.

# Arguments

- `eeg::EEG`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}` - epoch index to keep, vector of numbers or range

# Returns

- `eeg::EEG`
"""
function eeg_keep_epoch(eeg::EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange})
    size(eeg.eeg_signals, 3) == 1 && throw(ArgumentError("EEG contains only one epoch."))

    if typeof(epoch) <: AbstractRange
        epoch = collect(epoch)
    end

    length(epoch) > 1 && (epoch = sort!(epoch, rev=true))
    if epoch[end] < 1 || epoch[1] > size(eeg.eeg_signals, 1)
        throw(ArgumentError("Epoch index does not match signal epochs."))
    end

    epoch_list = collect(1:size(eeg.eeg_signals, 3))
    epoch_to_remove = setdiff(epoch_list, epoch)

    length(epoch_to_remove) > 1 && (epoch_to_remove = sort!(epoch_to_remove, rev=true))

    eeg_header = deepcopy(eeg.eeg_header)
    eeg_time = deepcopy(eeg.eeg_time)
    eeg_signals = deepcopy(eeg.eeg_signals)

    # remove epoch
    eeg_signals = eeg_signals[:, :, setdiff(1:end, (epoch_to_remove))]

    # update headers
    eeg_header[:epoch_n] = eeg_header[:epoch_n] - length(epoch_to_remove)
    epoch_n = eeg_header[:epoch_n]
    eeg_header[:eeg_duration_samples] = epoch_n * size(eeg.eeg_signals, 2)
    eeg_header[:eeg_duration_seconds] = round(epoch_n * (size(eeg.eeg_signals, 2) / eeg_sr(eeg)), digits=2)
    eeg_header[:epoch_duration_samples] = size(eeg.eeg_signals, 2)
    eeg_header[:epoch_duration_seconds] = round(size(eeg.eeg_signals, 2) / eeg_sr(eeg), digits=2)

    # create new dataset
    eeg_new = EEG(eeg_header, eeg_time, eeg_signals, deepcopy(eeg.eeg_components))
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_keep_epoch(EEG, $epoch)")
    
    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_keep_epoch!(eeg; epoch)

Keeps `epoch` in the `eeg`.

# Arguments

- `eeg::EEG`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}` - epoch index to keep, vector of numbers or range
"""
function eeg_keep_epoch!(eeg::EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange})
    size(eeg.eeg_signals, 3) == 1 && throw(ArgumentError("EEG contains only one epoch."))

    if typeof(epoch) <: AbstractRange
        epoch = collect(epoch)
    end

    length(epoch) > 1 && (epoch = sort!(epoch, rev=true))
    if epoch[end] < 1 || epoch[1] > size(eeg.eeg_signals, 1)
        throw(ArgumentError("Epoch index does not match signal epochs."))
    end

    epoch_list = collect(1:size(eeg.eeg_signals, 3))
    epoch_to_remove = setdiff(epoch_list, epoch)

    length(epoch_to_remove) > 1 && (epoch_to_remove = sort!(epoch_to_remove, rev=true))

    # remove epoch
    eeg.eeg_signals = eeg.eeg_signals[:, :, setdiff(1:end, (epoch_to_remove))]

    # update headers
    eeg.eeg_header[:epoch_n] = eeg.eeg_header[:epoch_n] - length(epoch_to_remove)
    epoch_n = eeg.eeg_header[:epoch_n]
    eeg.eeg_header[:eeg_duration_samples] = epoch_n * size(eeg.eeg_signals, 2)
    eeg.eeg_header[:eeg_duration_seconds] = round(epoch_n * (size(eeg.eeg_signals, 2) / eeg_sr(eeg)), digits=2)
    eeg.eeg_header[:epoch_duration_samples] = size(eeg.eeg_signals, 2)
    eeg.eeg_header[:epoch_duration_seconds] = round(size(eeg.eeg_signals, 2) / eeg_sr(eeg), digits=2)

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_keep_epoch(EEG, $epoch)")
    
    eeg_reset_components!(eeg)

    return
end

"""
    eeg_picks(eeg; pick)

Return `pick` of electrodes for `eeg` electrodes.

# Arguments

- `pick::Vector{Symbol}`

# Returns

- `channels::Vector{Int64}`
"""
function eeg_pick(eeg::EEG; pick::Union{Symbol, Vector{Symbol}})
    length(eeg.eeg_header[:labels]) == 0 && throw(ArgumentError("EEG does not contain channel labels."))

    if typeof(pick) == Vector{Symbol}
        for idx in 1:length(pick)
            pick[idx] in [:list, :central, :c, :left, :l, :right, :r, :frontal, :f, :temporal, :t, :parietal, :p, :occipital, :o] || throw(ArgumentError("Available picks: :central, :c, :left, :l, :right, :r, :frontal, :f, :temporal, :t, :parietal, :p, :occipital, :o"))
        end

        c = Vector{Char}()
        for idx in 1:length(pick)
            (pick[idx] === :central || pick[idx] === :c) && push!(c, 'z')
            (pick[idx] === :frontal || pick[idx] === :f) && push!(c, 'F')
            (pick[idx] === :temporal || pick[idx] === :t) && push!(c, 'T')
            (pick[idx] === :parietal || pick[idx] === :p) && push!(c, 'P')
            (pick[idx] === :occipital || pick[idx] === :o) && push!(c, 'O')
        end
        
        labels = deepcopy(eeg.eeg_header[:labels])
        channels = Vector{Int64}()
        for idx1 in 1:length(labels)
            for idx2 in 1:length(c)
                in(c[idx2], labels[idx1]) && push!(channels, idx1)
            end
        end

        # check for both :l and :r
        for idx1 in 1:length(pick)
            if (pick[idx1] === :left || pick[idx1] === :l)
                for idx2 in 1:length(pick)
                    if (pick[idx2] === :right || pick[idx2] === :r)
                        return channels
                    end
                end
            end
            if (pick[idx1] === :right || pick[idx1] === :r)
                for idx2 in 1:length(pick)
                    if (pick[idx2] === :left || pick[idx2] === :l)
                        return channels
                    end
                end
            end
        end

        labels = deepcopy(eeg.eeg_header[:labels])
        labels = labels[channels]
        pat = nothing
        for idx in 1:length(pick)
            # for :right remove lefts
            (pick[idx] === :right || pick[idx] === :r) && (pat = r"[z13579]$")
            # for :left remove rights
            (pick[idx] === :left || pick[idx] === :l) && (pat = r"[z02468]$")
        end
        if typeof(pat) == Regex
            for idx in length(labels):-1:1
                typeof(match(pat, labels[idx])) == RegexMatch && deleteat!(channels, idx)
            end
        end

        return channels
    else
        pick in [:central, :c, :left, :l, :right, :r, :frontal, :f, :temporal, :t, :parietal, :p, :occipital, :o] || throw(ArgumentError("Available picks: :central, :c, :left, :l, :right, :r, :frontal, :f, :temporal, :t, :parietal, :p, :occipital, :o"))

        c = Vector{Char}()
        (pick === :central || pick === :c) && (c = ['z'])
        (pick === :left || pick === :l) && (c = ['1', '3', '5', '7', '9'])
        (pick === :right || pick === :r) && (c = ['2', '4', '6', '8'])
        (pick === :frontal || pick === :f) && (c = ['F'])
        (pick === :temporal || pick === :t) && (c = ['T'])
        (pick === :parietal || pick === :p) && (c = ['P'])
        (pick === :occipital || pick === :o) && (c = ['O'])

        labels = deepcopy(eeg.eeg_header[:labels])
        channels = Vector{Int64}()
        for idx1 in 1:length(c)
            for idx2 in 1:length(labels)
                in(c[idx1], labels[idx2]) && push!(channels, idx2)
            end
        end

        return channels
    end
end

"""
    eeg_ica(eeg; n)

Calculates `n` first ICs for `eeg`.

# Arguments

- `eeg::EEG`
- `n::Int64` - number of ICs
- `tol::Float64` - tolerance for ICA
- `iter::Int64` - maximum number of iterations
- `f::Symbol` - neg-entropy functor
# Returns

- `ic::Array{Float64, 3}` - IC(1)..IC(n) × epoch
"""
function eeg_ica(eeg::EEG; n::Int64, tol::Float64=1.0e-6, iter::Int64=100, f::Symbol=:tanh)
    ic = signal_ica(eeg.eeg_signals, n=n, tol=tol, iter=iter, f=f)

    return ic
end

"""
    eeg_ica!(eeg; n, tol)

Calculates `n` first ICs for `eeg`.

# Arguments

- `eeg::EEG`
- `n::Int64` - number of ICs
- `tol::Float64` - tolerance for ICA
- `iter::Int64` - maximum number of iterations
- `f::Symbol` - neg-entropy functor
# Returns

- `eeg::EEG`
"""
function eeg_ica!(eeg::EEG; n::Int64, tol::Float64=1.0e-6, iter::Int64=100, f::Symbol=:tanh)
    :ica in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:ica)
    ic = signal_ica(eeg.eeg_signals, n=n, tol=tol, iter=iter, f=f)
    push!(eeg.eeg_components, ic)
    push!(eeg.eeg_header[:components], :ica)
    push!(eeg.eeg_header[:history], "eeg_ica!(EEG, n=$n, tol=$tol, iter=$iter, f=$f))")

    return
end

"""
    eeg_epochs_stats(eeg)

Calculates mean, sd and variance of `eeg` epochs.

# Arguments

- `eeg::EEG`

# Returns

- `mean::Vector{Float64}`
- `sd::Vector{Float64}`
- `var::Vector{Float64}`
"""
function eeg_epochs_stats(eeg::EEG)
    e_mean, e_sd, e_var = signal_epochs_stats(eeg.eeg_signals)

    return e_mean, e_sd, e_var
end

"""
    eeg_epochs_stats!(eeg)

Calculates mean, sd and variance of `eeg` epochs and stores in `eeg`.

# Arguments

- `eeg::EEG`
"""
function eeg_epochs_stats!(eeg::EEG)
    :epochs_mean in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:epochs_mean)
    :epochs_sd in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:epochs_sd)
    :epochs_var in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:epochs_var)
    e_mean, e_sd, e_var = signal_epochs_stats(eeg.eeg_signals)
    push!(eeg.eeg_components, e_mean)
    push!(eeg.eeg_components, e_sd)
    push!(eeg.eeg_components, e_var)
    push!(eeg.eeg_header[:components], :epochs_mean)
    push!(eeg.eeg_header[:components], :epochs_sd)
    push!(eeg.eeg_header[:components], :epochs_var)
    push!(eeg.eeg_header[:history], "eeg_epochs_stats!(EEG)")

    return
end

"""
    eeg_list_components(eeg)

Lists `eeg` components.

# Arguments

- `eeg::EEG`

# Returns

- `components::Vector{Symbol}`
"""
function eeg_list_components(eeg::EEG)
    return eeg.eeg_header[:components]
end

"""
    eeg_extract_component(eeg, c)

Extracts `component` of `eeg`.

# Arguments

- `eeg::EEG`
- `component::Symbol`

# Returns

- `component::Any`
"""
function eeg_extract_component(eeg::EEG; c::Symbol)
    c in eeg.eeg_header[:components] || throw(ArgumentError("Component $c does not exist. Use eeg_list_component() to view available components."))
    for idx in 1:length(eeg.eeg_header[:components])
        if c == eeg.eeg_header[:components][idx]
            return eeg.eeg_components[idx]
        end
    end
end

"""
    eeg_delete_component(eeg, c)

Deletes `component` of `eeg`.

# Arguments

- `eeg::EEG`
- `component::Symbol`

# Returns

- `eeg::EEG`
"""
function eeg_delete_component(eeg::EEG; c::Symbol)
    eeg_new = deepcopy(eeg)
    c in eeg_new.eeg_header[:components] || throw(ArgumentError("Component $c does not exist. Use eeg_list_component() to view available components."))
    for idx in 1:length(eeg.eeg_header[:components])
        if c == eeg_new.eeg_header[:components][idx]
            deleteat!(eeg_new.eeg_components, idx)
            deleteat!(eeg_new.eeg_header[:components], idx)
            push!(eeg_new.eeg_header[:history], "eeg_delete_component(EEG, c=$c)")
            return eeg_new
        end
    end
end

"""
    eeg_delete_component!(eeg, c)

Deletes `component` of `eeg`.

# Arguments

- `eeg::EEG`
- `component::Symbol`
"""
function eeg_delete_component!(eeg::EEG; c::Symbol)
    c in eeg.eeg_header[:components] || throw(ArgumentError("Component $c does not exist. Use eeg_list_component() to view available components."))
    for idx in length(eeg.eeg_header[:components]):-1:1
        if c == eeg.eeg_header[:components][idx]
            deleteat!(eeg.eeg_components, idx)
            deleteat!(eeg.eeg_header[:components], idx)
            push!(eeg.eeg_header[:history], "eeg_delete_component(EEG, c=$c)")
        end
    end

    return
end

"""
    eeg_spectrogram(eeg; norm=true, demean=true)

Calculates spectrogram of `eeg`.

# Arguments

- `eeg::EEG`
- `norm::Bool` - normalize powers to dB
- `demean::Bool` - demean signal prior to analysis

# Returns

- `spec.power::Array{Float64, 3}`
- `spec.freq::Matrix{Float64}`
- `spec.time::Matrix{Float64}`
"""
function eeg_spectrogram(eeg::EEG; norm::Bool=true, demean::Bool=true)
    s_pow, s_frq, s_t = signal_spectrogram(eeg.eeg_signals, fs=eeg_sr(eeg), norm=norm, demean=demean)

    return s_pow, s_frq, s_t
end


"""
    eeg_spectrogram!(eeg, norm=true, demean=true)

Calculates spectrogram of `eeg`.

# Arguments

- `eeg::EEG`
- `norm::Bool` - normalize powers to dB
- `demean::Bool` - demean signal prior to analysis
"""
function eeg_spectrogram!(eeg::EEG; norm::Bool=true, demean::Bool=true)
    :spectrogram_pow in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:spec_pow)
    :spectrogram_frq in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:spec_frq)
    :spectrogram_t in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:spec_t)
    s_pow, s_frq, s_t = signal_spectrogram(eeg.eeg_signals, fs=eeg_sr(eeg), norm=norm, demean=demean)
    push!(eeg.eeg_components, s_pow)
    push!(eeg.eeg_components, s_frq)
    push!(eeg.eeg_components, s_t)
    push!(eeg.eeg_header[:components], :spec_pow)
    push!(eeg.eeg_header[:components], :spec_frq)
    push!(eeg.eeg_header[:components], :spec_t)
    push!(eeg.eeg_header[:history], "eeg_spectrogram!(EEG, norm=$norm, demean=$demean)")

    return
end

"""
    eeg_average(eeg)

Returns the average signal of all `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object

# Returns

- `eeg::EEG`
"""
function eeg_average(eeg::EEG)
    # create new dataset
    eeg_new = deepcopy(eeg)
    eeg_keep_channel!(eeg_new, channel=1)
    eeg_new.eeg_signals = signal_average(eeg.eeg_signals)

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_average(EEG)")

    eeg_reset_components!(eeg_new)

    return eeg_new
end

"""
    eeg_average!(eeg)

Returns the average signal of all `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_average!(eeg::EEG)
    s_avg = signal_average(eeg.eeg_signals)
    eeg_delete_channel!(eeg, channel=2:size(eeg.eeg_signals, 1))
    eeg.eeg_signals = s_avg

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_average!(EEG)")

    eeg_reset_components!(eeg)

    return
end


"""
    eeg_spectrum(eeg; pad=0)

Calculates FFT, amplitudes, powers and phases for each channel of the `eeg`.

# Arguments

- `eeg::EEG` - the signal
- `pad::Int64` - pad channels `pad` zeros

# Returns

- `fft::Array{ComplexF64, 3}`
- `amplitudes::Array{Float64, 3}`
- `powers::Array{Float64, 3}`
- `phases::Array{Float64, 3}
"""
function eeg_spectrum(eeg::EEG; pad::Int64=0)
    s_fft, s_amp, s_pow, s_pha = signal_spectrum(eeg.eeg_signals, pad=pad)

    return s_fft, s_amp, s_pow, s_pha
end

"""
    eeg_spectrum!(eeg; pad=0)

Calculates FFT, amplitudes, powers and phases for each channel of the `eeg`.

# Arguments

- `eeg::EEG` - the signal
- `pad::Int64` - pad channels `pad` zeros
"""
function eeg_spectrum!(eeg::EEG; pad::Int64=0)
    :spectrum_fft in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:spectrum_fft)
    :spectrum_amp in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:spectrum_amp)
    :spectrum_pow in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:spectrum_pow)
    :spectrum_phase in eeg.eeg_header[:components] && eeg_delete_component!(eeg, c=:spectrum_phase)
    s_fft, s_amplitudes, s_powers, s_phases = signal_spectrum(eeg.eeg_signals)
    push!(eeg.eeg_components, s_fft)
    push!(eeg.eeg_components, s_amplitudes)
    push!(eeg.eeg_components, s_powers)
    push!(eeg.eeg_components, s_phases)
    push!(eeg.eeg_header[:components], :spectrum_fft)
    push!(eeg.eeg_header[:components], :spectrum_amp)
    push!(eeg.eeg_header[:components], :spectrum_pow)
    push!(eeg.eeg_header[:components], :spectrum_phase)
    push!(eeg.eeg_header[:history], "eeg_spectrum!(EEG, pad=$pad)")

    return
end

