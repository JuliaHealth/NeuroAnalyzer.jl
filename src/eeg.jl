"""
    eeg_drop_channel(eeg, channels)

Removes `channels` from the `eeg`.

# Arguments

- `eeg::EEG`
- `channels::Float64` - channels to be removed, vector of numbers or range
"""
function eeg_drop_channel(eeg::EEG, channels::Union{Int64, Vector{Int64}, UnitRange{Int64}})
    if typeof(channels) == UnitRange{Int64}
        channels = collect(channels)
    end

    length(channels) > 1 && (channels = sort!(channels, rev=true))

    if channels[end] < 1 || channels[1] > length(eeg.eeg_header[:labels])
        throw(ArgumentError("Channel index does not match signal channels."))
    end

    eeg_header = deepcopy(eeg.eeg_header)
    eeg_time = deepcopy(eeg.eeg_time)
    eeg_signals = deepcopy(eeg.eeg_signals)

    channels_no = eeg_header[:channels_no]

    # update headers
    eeg_header[:channels_no] = channels_no - length(channels)
    for idx1 in 1:length(channels)
        for idx2 in 1:channels_no
            if idx2 == channels[idx1]
                deleteat!(eeg_header[:labels], idx2)
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

    # remove channels
    eeg_signals = eeg_signals[setdiff(1:end, (channels)), :, :]

    # create new dataset
    eeg_new = EEG(eeg_header, eeg_time, eeg_signals)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_drop_channel(EEG, $channels)")
    
    return eeg_new
end

"""
    eeg_derivative(eeg)

Returns the derivative of the `eeg` with length same as the signal.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_derivative(eeg)
    signal_der = signal_derivative(eeg.eeg_signals)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), signal_der)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_derivative(EEG)")

    return eeg_new
end

"""
    eeg_total_power(eeg)

Calculates total power of the `eeg`.

# Arguments

- `eeg::EEG`
"""
function eeg_total_power(eeg)
    fs = eeg.eeg_header[:sampling_rate][1]
    stp = signal_total_power(eeg.eeg_signals, fs=fs)
    size(stp, 3) == 1 && (stp = reshape(stp, size(stp, 1), size(stp, 2)))

    return stp
end

"""
    eeg_band_power(eeg; f1, f2)

Calculates absolute band power between frequencies `f1` and `f2` of the `eeg`.

# Arguments

- `eeg::EEG`
- `f1::Float64` - lower frequency bound
- `f2::Float64` - upper frequency bound
"""
function eeg_band_power(eeg; f1, f2)
    fs = eeg.eeg_header[:sampling_rate][1]
    sbp = signal_band_power(eeg.eeg_signals, fs=fs, f1=f1, f2=f2)
    size(sbp, 3) == 1 && (sbp = reshape(sbp, size(sbp, 1), size(sbp, 2)))

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
"""
function eeg_detrend(eeg::EEG; type::Symbol=:linear)
    signal_det = signal_detrend(eeg.eeg_signals, type=type)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), signal_det)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_detrend(EEG, type=$type)")

    return eeg_new
end

"""
    eeg_reference_channel(eeg, reference_idx)

References the `eeg` to specific channel `reference_idx`.

# Arguments

- `eeg::EEG`
- `reference_idx::Union{Int64, Vector{Int64}}` - index of channels used as reference; if multiple channels are specified, their average is used as the reference
"""
function eeg_reference_channel(eeg::EEG, reference_idx::Union{Int64, Vector{Int64}, UnitRange{Int64}})
    if typeof(reference_idx) == UnitRange{Int64}
        reference_idx = collect(reference_idx)
    end
    signal_referenced = signal_reference_channel(eeg.eeg_signals, reference_idx)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), signal_referenced)
    eeg_new.eeg_header[:reference] = "channel: $reference_idx"
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_reference_channel(EEG, reference_idx=$reference_idx)")

    return eeg_new
end

"""
    eeg_reference_car(eeg)

References the `eeg` to common average reference.

# Arguments

- `eeg::EEG`
"""
function eeg_reference_car(eeg::EEG)
    signal_referenced = signal_reference_car(eeg.eeg_signals)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), signal_referenced)
    eeg_new.eeg_header[:reference] = "CAR"
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_reference_car(EEG)")

    return eeg_new
end

"""
    eeg_save(eeg, file_name; overwrite=false)

Saves the `eeg` to `file_name` file (HDF5-based).

# Arguments

- `eeg::EEG`
- `file_name::String` - file name
- `overwrite::Bool`
"""
function eeg_save(eeg::EEG, file_name::String; overwrite::Bool=false)
    if isfile(file_name) & overwrite == false
        throw(ArgumentError("""File $file_name already exists. To overwrite, add "overwrite=true" argument."""))
    end
    save_object(file_name, eeg)
end

"""
    eeg_load(file_name)

Loads the `eeg` from `file_name` file (HDF5-based).

# Arguments

- `file_name::String` - file name
"""
function eeg_load(file_name::String)
    eeg = load_object(file_name)
    return eeg
end

"""
    eeg_get_channel_idx(eeg, channel_name)

Returns the `channel_name` index.

# Arguments

- `eeg::EEG`
- `channel_name::String` - channel name
"""
function eeg_get_channel_idx(eeg::EEG, channel_name::String)
    labels = eeg.eeg_header[:labels]
    channel_idx = nothing
    for idx in 1:length(labels)
        if channel_name == labels[idx]
            channel_idx = idx
        end
    end
    if channel_idx == nothing
        throw(ArgumentError("Channel name does not match signal labels."))
    end
    return channel_idx
end

"""
    eeg_get_channel_idx(eeg, channel_idx)

Returns the `channel_idx` name.

# Arguments

- `eeg::EEG`
- `channel_idx::Int64` - channel index
"""
function eeg_get_channel_name(eeg::EEG, channel_idx::Int64)
    labels = eeg.eeg_header[:labels]
    if channel_idx < 1 || channel_idx > length(labels)
        throw(ArgumentError("Channel index does not match signal channels."))
    else
        channel_name = labels[channel_idx]
    end
    return channel_name
end

"""
    eeg_rename_channel(eeg, old_channel_name, new_channel_name)

Renames the `eeg` channel.

# Arguments

- `eeg::EEG`
- `old_channel_name::String`
- `new_name::String`
"""
function eeg_rename_channel(eeg::EEG, old_channel_name::String, new_channel_name::String)
    labels = eeg.eeg_header[:labels]
    channel_idx = nothing
    for idx in 1:length(labels)
        if old_channel_name == labels[idx]
            labels[idx] = new_channel_name
            channel_idx = idx
        end
    end
    if channel_idx == nothing
        throw(ArgumentError("Channel name does not match signal labels."))
    end

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), deepcopy(eeg.eeg_signals))
    eeg_new.eeg_header[:labels] = labels
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_rename_channel(EEG, old_channel_name=$old_channel_name, new_channel_name=$new_channel_name)")

    return eeg_new
end

"""
    eeg_rename_channel(eeg, channel_idx, new_channel_name)

Rename the `eeg` channel.

# Arguments

- `eeg::EEG`
- `channel_idx::Int64`
- `new_name::String`
"""
function eeg_rename_channel(eeg::EEG, channel_idx::Int64, new_channel_name::String)
    labels = eeg.eeg_header[:labels]
    if channel_idx < 1 || channel_idx > length(labels)
        throw(ArgumentError("Channel index does not match signal channels."))
    else
        labels[channel_idx] = new_channel_name
    end

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), deepcopy(eeg.eeg_signals))
    eeg_new.eeg_header[:labels] = labels
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_rename_channel(EEG, channel_idx=$channel_idx, new_channel_name=$new_channel_name)")

    return eeg_new
end

"""
    eeg_taper(eeg, taper)

Taper `eeg` with `taper`.

# Arguments

- `eeg::EEG`
- `taper::Vector`
"""
function eeg_taper(eeg::EEG, taper::Vector)
    signal_tapered = signal_taper(eeg.eeg_signals, taper)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), signal_tapered)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_taper(EEG, taper=$taper)")

    return eeg_new
end

"""
    eeg_demean(eeg)

Removes mean value (DC offset).

# Arguments

- `eeg::EEG`
"""
function eeg_demean(eeg::EEG)
    signal_demeaned = signal_demean(eeg.eeg_signals)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), signal_demeaned)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_demean(EEG)")

    return eeg_new
end

"""
    eeg_normalize_mean(eeg)

Normalize (scales around the mean).

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_normalize_mean(eeg::EEG)
    signal_normalized = signal_normalize_mean(eeg.eeg_signals)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), signal_normalized)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_normalize_mean(EEG)")

    return eeg_new
end

"""
    eeg_normalize_minmax(eeg)

Normalize to 0...1

# Arguments

- `eeg::EEG`
"""
function eeg_normalize_minmax(eeg::EEG)
    signal_normalized = signal_normalize_minmax(eeg.eeg_signals)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), signal_normalized)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_normalize_minmax(EEG)")

    return eeg_new
end

"""
    eeg_get_channel(eeg, channel_name)

Get the `eeg` channel by `channel_name`.

# Arguments

- `eeg::EEG`
- `channel_name::String`
"""
function eeg_get_channel(eeg::EEG, channel_name::String)
    labels = eeg.eeg_header[:labels]
    channel_idx = nothing
    for idx in 1:length(labels)
        if channel_name == labels[idx]
            channel_idx = idx
        end
    end
    if channel_idx == nothing
        throw(ArgumentError("Channel name does not match signal labels."))
    end
    channel = eeg.eeg_signals[channel_idx, :, :]

    return channel
end

"""
    eeg_get_channel(eeg, channel_idx)

Get the `eeg` channel by `channel_idx`.

# Arguments

- `eeg::EEG`
- `channel_idx::Int64`
"""
function eeg_get_channel(eeg::EEG, channel_idx::Int64)
    labels = eeg.eeg_header[:labels]
    if channel_idx < 1 || channel_idx > length(labels)
        throw(ArgumentError("Channel index does not match signal channels."))
    end
    channel = eeg.eeg_signals[channel_idx, :, :]

    return channel
end

"""
    eeg_cov(eeg; normalize=true)

Calculates covariance between all channels of `eeg`.

# Arguments

- `eeg::EEG`
"""
function eeg_cov(eeg::EEG; normalize=true)
    cov_mat = signal_cov(eeg.eeg_signals, normalize=normalize)
    size(cov_mat, 3) == 1 && (cov_mat = reshape(cov_mat, size(cov_mat, 1), size(cov_mat, 2)))

    return cov_mat
end

"""
    eeg_cor(eeg; normalize=true)

Calculates correlation coefficients between all channels of `eeg`.

# Arguments

- `eeg::EEG`
"""
function eeg_cor(eeg::EEG)
    cor_mat = signal_cor(eeg.eeg_signals)
    size(cor_mat, 3) == 1 && (cor_mat = reshape(cor_mat, size(cor_mat, 1), size(cor_mat, 2)))

    return cor_mat
end

"""
    eeg_upsample(eeg; new_sr)

Upsamples all channels of `eeg` to `new_sr` sampling frequency.

# Arguments

- `eeg::EEG`
- `new_sr::Int64` - new sampling rate
"""
function eeg_upsample(eeg::EEG; new_sr::Int64)
    t = eeg.eeg_time[1]:(1 / eeg.eeg_header[:sampling_rate][1]):eeg.eeg_time[end]
    signal_upsampled, t_upsampled = signal_upsample(eeg.eeg_signals, t=t, new_sr=new_sr)

    # create new dataset
    eeg_duration_samples = size(signal_upsampled, 2)
    eeg_duration_seconds = size(signal_upsampled, 2) / new_sr
    eeg_time = collect(t_upsampled)
    eeg_new = EEG(deepcopy(eeg.eeg_header), eeg_time, signal_upsampled)
    eeg_new.eeg_header[:eeg_duration_samples] = eeg_duration_samples
    eeg_new.eeg_header[:eeg_duration_seconds] = eeg_duration_seconds
    eeg_new.eeg_header[:sampling_rate] = repeat([new_sr], eeg_new.eeg_header[:channels_no])

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_upsample(EEG, new_sr=$new_sr)")

    return eeg_new
end

"""
    eeg_history(eeg)

Shows processing history.

# Arguments

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
"""
function eeg_labels(eeg::EEG)
    return eeg.eeg_header[:labels]
end

"""
    eeg_samplingrate(eeg)

Returns sampling rate.

# Arguments

- `eeg::EEG`
"""
function eeg_samplingrate(eeg::EEG)
    return eeg.eeg_header[:sampling_rate][1]
end

"""
    eeg_info(eeg)

Shows info.

# Arguments

- `eeg::EEG`
"""
function eeg_info(eeg::EEG)
    println("         EEG file name: $(eeg.eeg_header[:eeg_filename])")
    println("         EEG size [Mb]: $(eeg.eeg_header[:eeg_filesize_mb])")
    println("    Number of channels: $(eeg.eeg_header[:channels_no])")
    println("    Sampling rate (Hz): $(eeg.eeg_header[:sampling_rate][1])")
    println("      Number of epochs: $(eeg.eeg_header[:epochs_no])")
    println("Epoch length (samples): $(eeg.eeg_header[:epoch_duration_samples])")
    println("Epoch length (seconds): $(round(eeg.eeg_header[:epoch_duration_seconds], digits=2))")
    if eeg.eeg_header[:reference] == ""
        println("        Reference type: unknown")
    else
        println("        Reference type: $(eeg.eeg_header[:reference])")
    end
    if length(eeg.eeg_header[:labels]) == 0
        println("                Labels: no")
    else
        println("                Labels: yes")
    end
    if length(eeg.eeg_header[:channel_locations]) == false
        println("     Channel locations: no")
    else
        println("     Channel locations: yes")
    end
end

"""
    eeg_epochs(eeg; epochs_no=nothing, epochs_len=nothing, average=true)

Splits `eeg` into epochs.

# Arguments

- `eeg::EEG`
- `epochs_no::Union{Int64, Nothing}` - number of epochs
- `epochs_len::Union{Int64, Nothing}` - epoch length in samples
- `average::Bool` - average all epochs, returns one averaged epoch; if false than returns array of epochs, each row is one epoch
"""
function eeg_epochs(eeg::EEG; epochs_no::Union{Int64, Nothing}=nothing, epochs_len::Union{Int64, Nothing}=nothing, average::Bool=false)
    # unsplit epochs
    signal_merged = reshape(eeg.eeg_signals,
                            size(eeg.eeg_signals, 1), size(eeg.eeg_signals, 2) * size(eeg.eeg_signals, 3))
    
    # split into epochs
    signal_split = signal_epochs(signal_merged, epochs_no=epochs_no, epochs_len=epochs_len, average=average)

    # convert into Array{Float64, 3}
    signal_split = reshape(signal_split, size(signal_split, 1), size(signal_split, 2), size(signal_split, 3))

    # create new dataset
    epochs_no = size(signal_split, 3)
    epoch_duration_samples = size(signal_split, 2)
    epoch_duration_seconds = size(signal_split, 2) / eeg.eeg_header[:sampling_rate][1]
    eeg_duration_samples = size(signal_split, 2) * size(signal_split, 3)
    eeg_duration_seconds = eeg_duration_samples / eeg.eeg_header[:sampling_rate][1]
    eeg_time = collect(0:(1 / eeg.eeg_header[:sampling_rate][1]):epoch_duration_seconds)
    eeg_time = eeg_time[1:(end - 1)]
    eeg_new = EEG(deepcopy(eeg.eeg_header), eeg_time, signal_split)
    eeg_new.eeg_header[:eeg_duration_samples] = eeg_duration_samples
    eeg_new.eeg_header[:eeg_duration_seconds] = eeg_duration_seconds
    eeg_new.eeg_header[:epochs_no] = epochs_no
    eeg_new.eeg_header[:epoch_duration_samples] = epoch_duration_samples
    eeg_new.eeg_header[:epoch_duration_seconds] = epoch_duration_seconds

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_epochs(EEG, epochs_no=$epochs_no, epochs_len=$epochs_len, average=$average)")

    return eeg_new
end

"""
    eeg_get_epoch(eeg, epoch_idx)

Returns the `epoch_idx` epoch.

# Arguments

- `eeg::EEG`
- `epoch_idx::Int64` - epoch index
"""
function eeg_get_epoch(eeg::EEG, epoch_idx::Int64)
    if epoch_idx < 1 || epoch_idx > eeg.eeg_header[:epochs_no]
        throw(ArgumentError("Epoch index out of range."))
    end

    signal_new = eeg.eeg_signals[:, :, epoch_idx]
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), signal_new)
    eeg_new.eeg_header[:epochs_no] = 1
    eeg_new.eeg_header[:eeg_duration_samples] = eeg_new.eeg_header[:epoch_duration_samples]
    eeg_new.eeg_header[:eeg_duration_seconds] = eeg_new.eeg_header[:epoch_duration_seconds]

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_get_epoch(EEG, epoch_idx=$epoch_idx)")

    return eeg_new
end

"""
    eeg_tconv(eeg, kernel)

Performs convolution in the time domain.

# Arguments

- `eeg::EEG`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`
"""
function eeg_tconv(eeg::EEG; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})
    signal_convoluted = signal_tconv(eeg.eeg_signals, kernel)

    ## EEG signal can only store Float64
    typeof(kernel) == Vector{ComplexF64} && (signal_convoluted = abs.(signal_convoluted))

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), signal_convoluted)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_tconv(EEG, kernel=$kernel)")

    return eeg_new
end

"""
    eeg_filter(signal; prototype, filter_type, cutoff, fs, order=8, window=hanning(64))

Filters `signal` using zero phase distortion filter.

# Arguments

- `eeg::EEG`
- `fprototype::Symbol[:butterworth, :fir]
- `ftype::Symbol[:lp, :hp, :bp, :bs]` - filter type
- `cutoff::Union{Int64, Float64, Vector{Int64}, Vector{Float64}}` - filter cutoff in Hz (vector for `:bp` and `:bs`)
- `fs::Int64` - sampling rate
- `order::Int64` - filter order
- `window::Vector{Float64} - window, required for FIR filter
"""
function eeg_filter(eeg::EEG; fprototype::Symbol, ftype::Symbol, cutoff::Union{Int64, Float64, Vector{Int64}, Vector{Float64}}, order::Int64=8, window::Vector{Float64}=hanning(64))
    fs = eeg.eeg_header[:sampling_rate][1]

    signal_filtered = signal_filter(eeg.eeg_signals,
                                    fprototype=fprototype,
                                    ftype=ftype,
                                    cutoff=cutoff,
                                    fs=fs,
                                    order=order,
                                    window=window)

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), signal_filtered)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_filter(EEG, fprototype=$fprototype,ftype=$ftype, cutoff=$cutoff, order=$order, window=$window)")

    return eeg_new
end