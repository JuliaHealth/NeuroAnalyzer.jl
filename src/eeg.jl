"""
    eeg_plot(eeg; t=nothing, offset=0, labels=[], normalize=false, xlabel="Time [s]", ylabel="Channels", figure=nothing)

Plots `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
- `t::Vector{Float64}` - the time vector
- `offset::Float64` - displayed segment offset in samples
- `labels::Vector{String}` - channel labels vector
- `normalize::Bool` - normalize the `signal` prior to calculations
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis lable
- `figure::String` - name of the output figure file
"""
function eeg_plot(eeg::EEG; t=nothing, offset=1, labels=[], normalize=true, xlabel="Time [s]", ylabel="Channels", figure::String="")
    
    if typeof(t) == UnitRange{Int64}
        t = collect(t)
    end

    signal = eeg.eeg_signals
    labels = eeg.eeg_signal_header[:labels]
    fs = eeg.eeg_signal_header[:sampling_rate][1]

    # default time is 10 seconds
    t === nothing && (t = collect(0:1/fs:10))

    p, signal_new = signal_plot(t, signal, offset=offset, labels=labels, xlabel=xlabel, ylabel=ylabel, normalize=normalize, figure=figure)

    plot(p)

    return p
end

"""
    eeg_drop_channel(eeg, channels)

Removes `channels` from the `eeg` object.

# Arguments

- `eeg::EEG` - EEG object
- `channels::Float64` - channels to be removed, vector of numbers or range
"""
function eeg_drop_channel(eeg::EEG, channels)
    if typeof(channels) == UnitRange{Int64}
        channels = collect(channels)
    end

    length(channels) > 1 && (channels = sort!(channels, rev=true))

    eeg_object_header = eeg.eeg_object_header
    eeg_signal_header = eeg.eeg_signal_header
    eeg_signals = eeg.eeg_signals

    channels_no = eeg_object_header[:channels_no]

    # update headers
    eeg_object_header[:channels_no] = channels_no - length(channels)
    for idx1 in 1:length(channels)
        for idx2 in 1:channels_no
            if idx2 == channels[idx1]
                deleteat!(eeg_signal_header[:labels], idx2)
                deleteat!(eeg_signal_header[:transducers], idx2)
                deleteat!(eeg_signal_header[:physical_dimension], idx2)
                deleteat!(eeg_signal_header[:physical_minimum], idx2)
                deleteat!(eeg_signal_header[:physical_maximum], idx2)
                deleteat!(eeg_signal_header[:digital_minimum], idx2)
                deleteat!(eeg_signal_header[:digital_maximum], idx2)
                deleteat!(eeg_signal_header[:prefiltering], idx2)
                deleteat!(eeg_signal_header[:samples_per_datarecord], idx2)
                deleteat!(eeg_signal_header[:sampling_rate], idx2)
                deleteat!(eeg_signal_header[:gain], idx2)
            end
        end 
    end

    # remove channels
    eeg_signals = eeg_signals[setdiff(1:end, (channels)), :]

    # add entry to :history field
    push!(eeg_object_header[:history], "eeg_drop_channel(EEG, $channels)")

    # create new dataset    
    eeg_new = EEG(eeg_object_header, eeg_signal_header, eeg_signals)
    
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
    fs = eeg.eeg_signal_header[:sampling_rate][1]
    signal_filtered = signal_filter_butter(eeg.eeg_signals, filter_type=filter_type, cutoff=cutoff, fs=fs, poles=poles)
    eeg_new = EEG(eeg.eeg_object_header, eeg.eeg_signal_header, signal_filtered)

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
    eeg_new = EEG(eeg.eeg_object_header, eeg.eeg_signal_header, signal_der)

    return eeg_new
end

"""
    eeg_total_power(eeg)

Calculates total power for each the `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_total_power(eeg)
    fs = eeg.eeg_signal_header[:sampling_rate][1]
    stp = signal_total_power(eeg.eeg_signals, fs)

    return stp
end

"""
    eeg_band_power(eeg, f1, f2)

Calculates absolute band power between frequencies `f1` and `f2` for each the `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
- `f1::Float64` - Lower frequency bound
- `f2::Float64` - Upper frequency bound
"""
function eeg_total_power(eeg, f1, f2)
    fs = eeg.eeg_signal_header[:sampling_rate][1]
    sbp = signal_band_power(eeg.eeg_signals, fs, f1, f2)

    return sbp
end

"""
    eeg_make_spectrum(eeg)

Returns FFT and DFT sample frequencies for a DFT for each the `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_make_spectrum(eeg)
    fs = eeg.eeg_signal_header[:sampling_rate][1]
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
    signal_det = signal_detrend(eeg.eeg_signals, type)
    eeg_new = EEG(eeg.eeg_object_header, eeg.eeg_signal_header, signal_det)

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
    eeg.eeg_object_header[:reference_type] = "common reference"
    eeg.eeg_object_header[:reference_channel] = reference_idx
    eeg_new = EEG(eeg.eeg_object_header, eeg.eeg_signal_header, signal_referenced)

    return eeg_new
end

"""
    eeg_reference_car(eeg)

References the `eeg` channels to common average reference.

# Arguments

- `eeg::EEG` - EEG object
- `reference::Float64` - index of channels used as reference; if multiple channels are specified, their average is
"""
function eeg_reference_car(eeg::EEG, reference_idx)
    signal_referenced = eeg_reference_car(eeg.eeg_signals, reference_idx)
    eeg.eeg_object_header[:reference_type] = "CAR"
    eeg.eeg_object_header[:reference_channel] = []
    eeg_new = EEG(eeg.eeg_object_header, eeg.eeg_signal_header, signal_referenced)

    return eeg_new
end

"""
    eeg_save(eeg, file_name)

Saves the `eeg` object to `file_name` file (HDF5-based).

# Arguments

- `eeg::EEG` - EEG object
- `file_name::String` - file name
"""
function eeg_save(eeg::EEG, file_name)
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
    labels = eeg.eeg_signal_header[:labels]
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
- `channel_idx::Int` - channel index
"""
function eeg_get_channel_name(eeg::EEG, channel_idx::Int)
    labels = eeg.eeg_signal_header[:labels]
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
    labels = eeg.eeg_signal_header[:labels]
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
    eeg_signal_header[:labels] = labels
    eeg_new = EEG(eeg.eeg_object_header, eeg_signal_header, eeg.eeg_signals)

    return eeg_new
end

"""
    eeg_rename_channel(eeg, channel_idx, new_channel_name)

Rename the `eeg` channel.

# Arguments

- `eeg::EEG` - EEG object
- `channel_idx::Int`
- `new_name::String`
"""
function eeg_rename_channel(eeg::EEG, channel_idx::Int, new_channel_name::String)
    labels = eeg.eeg_signal_header[:labels]
    if channel_idx < 1 || channel_idx > length(labels)
        throw(ArgumentError("Channel index does not match signal channels."))
    else
        labels[channel_idx] = new_channel_name
    end
    eeg_signal_header[:labels] = labels
    eeg_new = EEG(eeg.eeg_object_header, eeg_signal_header, eeg.eeg_signals)

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
    eeg_new = EEG(eeg.eeg_object_header, eeg.eeg_signal_header, signal_tapered)

    return eeg_new
end

"""
    eeg_demean(eeg)

Removes mean value (DC offset) for each the `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_demean(eeg::EEG)
    signal_demeaned = signal_demean(eeg.eeg_signals, taper)
    eeg_new = EEG(eeg.eeg_object_header, eeg.eeg_signal_header, signal_demeaned)

    return eeg_new
end

"""
    eeg_normalize_mean(eeg)

Normalize (scales around the mean) each the `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_normalize_mean(eeg::EEG)
    signal_normalized = signal_normalize_mean(eeg.eeg_signals, taper)
    eeg_new = EEG(eeg.eeg_object_header, eeg.eeg_signal_header, signal_tapered)

    return eeg_new
end

"""
    eeg_normalize_minmax(eeg)

Normalize (to 0…1) each the `eeg` channels.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_normalize_minmax(eeg::EEG)
    signal_normalized = signal_normalize_minmax(eeg.eeg_signals, taper)
    eeg_new = EEG(eeg.eeg_object_header, eeg.eeg_signal_header, signal_tapered)

    return eeg_new
end
