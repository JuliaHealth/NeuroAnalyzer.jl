"""
    eeg_plot(eeg; t=nothing, epoch=1, offset=0, labels=[], normalize=false, xlabel="Time [s]", ylabel="Channels", figure=nothing)

Plots `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
- `t::Union{Vector{Float64}, UnitRange{Int64}, Nothing}` - the time vector
- `epoch::Int64` - epoch number to display
- `offset::Int64` - displayed segment offset in samples
- `len::Float64` - length in seconds
- `labels::Vector{String}` - channel labels vector
- `normalize::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis lable
- `figure::String` - name of the output figure file
"""
function eeg_plot(eeg::EEG; t::Union{Vector{Float64}, UnitRange{Int64}, Nothing}=nothing, epoch::Int64=1, offset::Int64=1, len::Float64=10.0, labels::Vector{String}=[], normalize::Bool=true, xlabel::String="Time [s]", ylabel::String="Channels", figure::String="")

    if epoch < 1 || epoch > eeg.eeg_header[:epochs_no]
        throw(ArgumentError("Epoch index out of range."))
    end

    if typeof(t) == UnitRange{Int64}
        t = collect(t)
    end

    signal = eeg.eeg_signals[:, :, epoch]
    labels = eeg.eeg_header[:labels]
    fs = eeg.eeg_header[:sampling_rate][1]

    # default time is 10 seconds or epoch_duration_seconds
    len > eeg.eeg_header[:epoch_duration_seconds] && (len = eeg.eeg_header[:epoch_duration_seconds])
    t === nothing && (t = collect(0:1/fs:len))

    p = signal_plot(t, signal, offset=offset, labels=labels, xlabel=xlabel, ylabel=ylabel, normalize=normalize)

    plot(p)

    # TO DO: catching error while saving
    figure !== "" && (savefig(p, figure))
    
    return p
end

"""
    eeg_drop_channel(eeg, channels)

Removes `channels` from the `eeg` object.

# Arguments

- `eeg::EEG` - EEG object
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

    eeg_header = eeg.eeg_header
    eeg_time = eeg.eeg_time
    eeg_signals = eeg.eeg_signals

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
    eeg_filter_butter(eeg; filter_type, cutoff, fs, poles=8)

Filters `eeg` channels using Butterworth filter.

# Arguments

- `eeg::EEG` - EEG object
- `filter_type::Symbol[:lp, :hp, :bp, :bs]` - filter type
- `cutoff::Float64` - filter cutoff in Hz (tuple or vector for `:bp` and `:bs`)
- `fs::Float64` - sampling rate
- `poles::Int` - filter pole
"""
function eeg_filter_butter(eeg::EEG; filter_type, cutoff, poles=8)
    fs = eeg.eeg_header[:sampling_rate][1]
    signal_filtered = signal_filter_butter(eeg.eeg_signals, filter_type=filter_type, cutoff=cutoff, fs=fs, poles=poles)

    # create new dataset
    eeg_new = EEG(eeg.eeg_header, eeg.eeg_time, signal_filtered)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_filter_butter(EEG, filter_type=$filter_type, cutoff=$cutoff, poles=$poles)")

    return eeg_new
end

"""
    eeg_derivative(eeg)

Returns the derivative of each the `eeg` channels with length same as the signal.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_derivative(eeg)
    signal_der = signal_derivative(eeg.eeg_signals)

    # create new dataset
    eeg_new = EEG(eeg.eeg_header, eeg.eeg_time, signal_der)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_derivative(EEG)")

    return eeg_new
end

"""
    eeg_total_power(eeg)

Calculates total power for each the `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_total_power(eeg)
    fs = eeg.eeg_header[:sampling_rate][1]
    stp = signal_total_power(eeg.eeg_signals, fs=fs)

    return stp
end

"""
    eeg_band_power(eeg; f1, f2)

Calculates absolute band power between frequencies `f1` and `f2` for each the `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
- `f1::Float64` - Lower frequency bound
- `f2::Float64` - Upper frequency bound
"""
function eeg_band_power(eeg; f1, f2)
    fs = eeg.eeg_header[:sampling_rate][1]
    sbp = signal_band_power(eeg.eeg_signals, fs=fs, f1=f1, f2=f2)

    return sbp
end

"""
    eeg_make_spectrum(eeg)

Returns FFT and DFT sample frequencies for a DFT for each the `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_make_spectrum(eeg)
    fs = eeg.eeg_header[:sampling_rate][1]
    signal_fft, signal_sf = signal_make_spectrum(eeg.eeg_signals, fs)

    return signal_fft, signal_sf
end

"""
    eeg_detrend(eeg, type=:linear)

Removes linear trend for each the `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
- `type::Symbol[:linear, :constant]`, optional
    - `linear` - the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `constant` - the mean of `signal` is subtracted
"""
function eeg_detrend(eeg, type=:linear)
    signal_det = signal_detrend(eeg.eeg_signals, type=type)

    # create new dataset
    eeg_new = EEG(eeg.eeg_header, eeg.eeg_time, signal_det)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_detrend(EEG, type=$type)")

    return eeg_new
end

"""
    eeg_draw_head(p, loc_x, loc_y, add_labels=true)

Draws head over a topographical plot `p`.

# Arguments

- `p::Plot` - toppgraphical plot
- `loc_x::Vector{Float64}` - vector of x electrode position
- `loc_y::Vector{Float64}` - vector of y electrode position
- `add_labels::Bool` - add text labels to the plot
"""
function eeg_draw_head(p, loc_x::Vector{Float64}, loc_y::Vector{Float64}, add_labels=true)
    pts = Plots.partialcircle(0, 2π, 100, maximum(loc_x))
    x, y = Plots.unzip(pts)
    x = x .* 1.1
    y = y .* 1.1
    head = Shape(x, y)
    nose = Shape([(-0.1, maximum(y)), (0, maximum(y) + 0.1 * maximum(y)), (0.1, maximum(y))])
    ear_l = Shape([(minimum(x), -0.1), (minimum(x) + 0.1 * minimum(x), -0.1), (minimum(x) + 0.1 * minimum(x), 0.1), (minimum(x), 0.1)])
    ear_r = Shape([(maximum(x), -0.1), (maximum(x) + 0.1 * maximum(x), -0.1), (maximum(x) + 0.1 * maximum(x), 0.1), (maximum(x), 0.1)])
    plot!(p, head, fill=nothing, label="")
    plot!(p, nose, fill=nothing, label="")
    plot!(p, ear_l, fill=nothing, label="")
    plot!(p, ear_r, fill=nothing, label="")
    if add_labels == true
        plot!(p, annotation=(0, 1 - maximum(y) / 5, text("Inion", pointsize=12, halign=:center, valign=:center)))
        plot!(p, annotation=(0, -1 - minimum(y) / 5, text("Nasion", pointsize=12, halign=:center, valign=:center)))
        plot!(p, annotation=(-1 - minimum(x) / 5, 0, text("Left", pointsize=12, halign=:center, valign=:center, rotation=90)))
        plot!(p, annotation=(1 - maximum(x) / 5, 0, text("Right", pointsize=12, halign=:center, valign=:center, rotation=-90)))
    end
end

"""
    eeg_reference_channel(eeg, reference)

References the `eeg` channels to specific signal channel.

# Arguments

- `eeg::EEG` - EEG object
- `reference::Float64` - index of channels used as reference; if multiple channels are specified, their average is used as the reference
"""
function eeg_reference_channel(eeg::EEG, reference_idx)
    signal_referenced = signal_reference_channel(eeg.eeg_signals, reference_idx)
    eeg.eeg_header[:reference_type] = "common reference"
    eeg.eeg_header[:reference_channel] = reference_idx

    # create new dataset
    eeg_new = EEG(eeg.eeg_header, eeg.eeg_time, signal_referenced)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_reference_channel(EEG, reference_idx=$reference_idx)")

    return eeg_new
end

"""
    eeg_reference_car(eeg)

References the `eeg` channels to common average reference.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_reference_car(eeg::EEG)
    signal_referenced = signal_reference_car(eeg.eeg_signals)
    eeg.eeg_header[:reference_type] = "CAR"
    eeg.eeg_header[:reference_channel] = []

    # create new dataset
    eeg_new = EEG(eeg.eeg_header, eeg.eeg_time, signal_referenced)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_reference_car(EEG)")

    return eeg_new
end

"""
    eeg_save(eeg, file_name; overwrite=false)

Saves the `eeg` object to `file_name` file (HDF5-based).

# Arguments

- `eeg::EEG` - EEG object
- `file_name::String` - file name
"""
function eeg_save(eeg::EEG, file_name; overwrite=false)
    if isfile(file_name) & overwrite == false
        throw(ArgumentError("""File $file_name already exists. To overwrite, add "overwrite=true" argument."""))
    end
    save_object(file_name, eeg)
end

"""
    eeg_load(file_name)

Loads the `eeg` object from `file_name` file (HDF5-based).

# Arguments

- `file_name::String` - file name
"""
function eeg_load(file_name)
    eeg = load_object(file_name)
    return eeg
end

"""
    eeg_get_channel_idx(eeg, channel_name)

Returns the `channel_name` index.

# Arguments

- `eeg::EEG` - EEG object
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

- `eeg::EEG` - EEG object
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

Rename the `eeg` channel.

# Arguments

- `eeg::EEG` - EEG object
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
    eeg.eeg_header[:labels] = labels

    # create new dataset
    eeg_new = EEG(eeg.eeg_header, eeg.eeg_time, eeg.eeg_signals)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_rename_channel(EEG, old_channel_name=$old_channel_name, new_channel_name=$new_channel_name)")

    return eeg_new
end

"""
    eeg_rename_channel(eeg, channel_idx, new_channel_name)

Rename the `eeg` channel.

# Arguments

- `eeg::EEG` - EEG object
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
    eeg.eeg_header[:labels] = labels

    # create new dataset
    eeg_new = EEG(eeg.eeg_header, eeg.eeg_time, eeg.eeg_signals)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_rename_channel(EEG, channel_idx=$channel_idx, new_channel_name=$new_channel_name)")

    return eeg_new
end

"""
    eeg_taper(eeg, taper)

Taper each the `eeg` channels with `taper`.

# Arguments

- `eeg::EEG` - EEG object
- `taper::Vector`
"""
function eeg_taper(eeg::EEG, taper::Vector)
    signal_tapered = signal_taper(eeg.eeg_signals, taper)

    # create new dataset
    eeg_new = EEG(eeg.eeg_header, eeg.eeg_time, signal_tapered)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_taper(EEG, taper=$taper)")

    return eeg_new
end

"""
    eeg_demean(eeg)

Removes mean value (DC offset) for each the `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_demean(eeg::EEG)
    signal_demeaned = signal_demean(eeg.eeg_signals)

    # create new dataset
    eeg_new = EEG(eeg.eeg_header, eeg.eeg_time, signal_demeaned)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_demean(EEG)")

    return eeg_new
end

"""
    eeg_normalize_mean(eeg)

Normalize (scales around the mean) each the `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_normalize_mean(eeg::EEG)
    signal_normalized = signal_normalize_mean(eeg.eeg_signals)

    # create new dataset
    eeg_new = EEG(eeg.eeg_header, eeg.eeg_time, signal_normalized)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_normalize_mean(EEG)")

    return eeg_new
end

"""
    eeg_normalize_minmax(eeg)

Normalize (to 0…1) each the `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_normalize_minmax(eeg::EEG)
    signal_normalized = signal_normalize_minmax(eeg.eeg_signals)

    # create new dataset
    eeg_new = EEG(eeg.eeg_header, eeg.eeg_time, signal_normalized)
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_normalize_minmax(EEG)")

    return eeg_new
end

"""
    eeg_get_channel(eeg, channel_name)

Get the `eeg` channel by `channel_name`.

# Arguments

- `eeg::EEG` - EEG object
- `channel_name::String`
"""
function eeg_get_channel(eeg::EEG, channel_name::String)
    labels = eeg.eeg_header[:labels]
    channel_idx = nothing
    for idx in 1:length(labels)
        if channel_name == labels[idx]
            labels[idx] = new_channel_name
            channel_idx = idx
        end
    end
    if channel_idx == nothing
        throw(ArgumentError("Channel name does not match signal labels."))
    end
    channel = eeg.eeg_signals[channel_idx, :]

    return channel
end

"""
    eeg_get_channel(eeg, channel_idx)

Get the `eeg` channel by `channel_idx`.

# Arguments

- `eeg::EEG` - EEG object
- `channel_idx::Int64`
"""
function eeg_get_channel(eeg::EEG, channel_idx::Int64)
    labels = eeg.eeg_header[:labels]
    if channel_idx < 1 || channel_idx > length(labels)
        throw(ArgumentError("Channel index does not match signal channels."))
    else
        labels[channel_idx] = new_channel_name
    end
    eeg_header[:labels] = labels
    channel = eeg.eeg_signals[channel_idx, :]

    return channel
end

"""
    eeg_cov(eeg; normalize=true)

Calculates covariance between all channels of the `eeg` object.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_cov(eeg::EEG; normalize=true)
    result = signal_cov(eeg.eeg_signals, normalize=normalize)
    return result
end

"""
    eeg_cor(eeg; normalize=true)

Calculates correlation coefficients between all channels of the `eeg` object.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_cor(eeg::EEG)
    result = signal_cor(eeg.eeg_signals)
    return result
end

"""
    eeg_upsample(eeg; new_sr)

Upsamples all channels of the `eeg` object to `new_sr` sampling frequency.

# Arguments

- `eeg::EEG` - EEG object
- `new_sr::Int64` - new sampling rate
"""
function eeg_upsample(eeg::EEG; new_sr::Int64)
    t = eeg.eeg_time[1]:(1 / eeg.eeg_header[:sampling_rate][1]):eeg.eeg_time[end]
    signal_upsampled, t_upsampled = signal_upsample(eeg.eeg_signals, t=t, new_sr=new_sr)

    # create new dataset
    eeg_duration_samples = size(signal_upsampled, 2)
    eeg_duration_seconds = size(signal_upsampled, 2) / new_sr
    eeg_time = collect(t_upsampled)
    eeg_new = EEG(eeg.eeg_header, eeg_time, signal_upsampled)
    eeg_new.eeg_header[:eeg_duration_samples] = eeg_duration_samples
    eeg_new.eeg_header[:eeg_duration_seconds] = eeg_duration_seconds
    eeg_new.eeg_header[:sampling_rate] = repeat([new_sr], eeg_new.eeg_header[:channels_no])

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_upsample(EEG, new_sr=$new_sr)")

    return eeg_new
end

"""
    eeg_show_processing_history(eeg)

Shows processing history of the `eeg` object.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_show_processing_history(eeg::EEG)
    return eeg.eeg_header[:history]
end

"""
    eeg_info(eeg)

Shows info of the `eeg` object.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_info(eeg::EEG)
    println("         EEG file name: $(eeg.eeg_header[:eeg_filename])")
    println("         EEG size [Mb]: $(eeg.eeg_header[:eeg_filesize_mb])")
    println("    Number of channels: $(eeg.eeg_header[:channels_no])")
    println("    Sampling rate (Hz): $(eeg.eeg_header[:sampling_rate][1])")
    println("      Number of epochs: $(eeg.eeg_header[:epochs_no])")
    println("Epoch length (samples): $(eeg.eeg_header[:epoch_duration_samples])")
    println("Epoch length (seconds): $(round(eeg.eeg_header[:epoch_duration_seconds], digits=2))")
    if eeg.eeg_header[:reference_type] == ""
        println("        Reference type: unknown")
    elseif eeg.eeg_header[:reference_type] == "common reference"
        println("        Reference type: $(eeg.eeg_header[:reference_type]), channel: $(eeg.eeg_header[:reference_channel])")
    else
        println("        Reference type: $(eeg.eeg_header[:reference_type])")
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

Splits `eeg` signals into epochs.

# Arguments

- `eeg::EEG` - EEG object
- `epochs_no::Union{Int64, Nothing}=nothing` - number of epochs
- `epochs_len::Union{Int64, Nothing}` - epoch length in samples
- `average::Bool` - average all epochs, returns one averaged epoch; if false than returns array of epochs, each row is one epoch
"""
function eeg_epochs(eeg::EEG; epochs_no::Union{Int64, Nothing}=nothing, epochs_len::Union{Int64, Nothing}=nothing, average::Bool=false)
    # unsplit epochs
    signal_merged = reshape(eeg.eeg_signals,
                            size(eeg.eeg_signals, 1), size(eeg.eeg_signals, 2) * size(eeg.eeg_signals, 3))
    
    # split into epochs
    signal_split = signal_epochs(signal_merged, epochs_no=epochs_no, epochs_len=epochs_len, average=average)

    # create new dataset
    epochs_no = size(signal_split, 3)
    epoch_duration_samples = size(signal_split, 2)
    epoch_duration_seconds = size(signal_split, 2) / eeg.eeg_header[:sampling_rate][1]
    eeg_duration_samples = size(signal_split)[2] * size(signal_split)[3]
    eeg_duration_seconds = eeg_duration_samples / eeg.eeg_header[:sampling_rate][1]
    eeg_time = collect(1:1/eeg.eeg_header[:sampling_rate][1]:epoch_duration_samples)
    eeg_new = EEG(eeg.eeg_header, eeg_time, signal_split)
    eeg_new.eeg_header[:eeg_duration_samples] = eeg_duration_samples
    eeg_new.eeg_header[:eeg_duration_seconds] = eeg_duration_seconds
    eeg_new.eeg_header[:epochs_no] = epochs_no
    eeg_new.eeg_header[:epoch_duration_samples] = epoch_duration_samples
    eeg_new.eeg_header[:epoch_duration_seconds] = epoch_duration_seconds

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_epochs(EEG, epochs_no=$epochs_no, epochs_len=$epochs_len, average=$average")

    return eeg_new
end

"""
    eeg_get_epoch(eeg, epoch_idx)

Returns the `epoch_idx` epoch.

# Arguments

- `eeg::EEG` - EEG object
- `epoch_idx::Int64` - epoch index
"""
function eeg_get_epoch(eeg::EEG, epoch_idx::Int64)
    if epoch_idx < 1 || epoch_idx > eeg.eeg_header[:epochs_no]
        throw(ArgumentError("Epoch index out of range."))
    end

    signal_new = eeg.eeg_signals[:, :, epoch_idx]
    eeg_new = EEG(eeg.eeg_header, eeg_time, signal_new)
    eeg_new.eeg_header[:epochs_no] = 1

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_get_epoch(EEG, epoch_idx=$epoch_idx")

    return eeg_new
end