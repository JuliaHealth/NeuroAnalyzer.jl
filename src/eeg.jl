"""
    eeg_plot(eeg; t=nothing, offset=0, labels=[], normalize=false, xlabel="Time [s]", ylabel="Channels", figure=nothing)

Plots `eeg` signals.

# Arguments

- `eeg::EEG` - EEG object
- `t::Vector{Float64} - the time vector
- `offset::Float64` - displayed segment offset in samples
- `labels::Vector{String}` - channel labels vector
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis lable
- `normalize::Bool` - normalize the `signal` prior to calculations
- `remove_dc::Bool` - demean the `signal` prior to calculations
- `detrend::Bool` - detrend the `signal` prior to calculations
- `derivative::Bool` - derivate `signal` prior to calculations
- `taper::Bool` - taper the `signal` with `taper`-window prior to calculations
- `figure::String` - name of the output figure file
"""
function eeg_plot(eeg::EEG; t=nothing, offset=1, labels=[], xlabel="Time [s]", ylabel="Channels", normalize=true, remove_dc=false, detrend=false, derivative=false, taper=nothing, figure::String="")
    
    if typeof(t) == UnitRange{Int64}
        t = collect(t)
    end

    signal = eeg.eeg_signals
    labels = eeg.eeg_signal_header[:labels]
    fs = eeg.eeg_signal_header[:sampling_rate][1]

    # default time is 5 seconds
    t === nothing && (t = collect(0:1/fs:5))

    p, signal_new = signal_plot(t, signal, offset=offset, labels=labels, xlabel=xlabel, ylabel=ylabel, normalize=normalize, remove_dc=remove_dc, detrend=detrend, derivative=derivative, taper=taper, figure=figure)

    plot(p)
    
    # create new dataset    
    eeg_new = EEG(eeg.eeg_file_header, eeg.eeg_signal_header, signal_new)

    return p, eeg_new
end

"""
    eeg_drop_channel(eeg, channels)

Removes `channels` from the `eeg` set.

# Arguments

- `eeg::EEG` - EEG object
- `channels::Float64` - channels to be removed, vector of numbers or range
"""
function eeg_drop_channel(eeg::EEG, channels)
    if typeof(channels) == UnitRange{Int64}
        channels = collect(channels)
    end

    channels = sort!(channels, rev=true)

    eeg_file_header = eeg.eeg_file_header
    eeg_signal_header = eeg.eeg_signal_header
    eeg_signals = eeg.eeg_signals

    channels_no = eeg_file_header[:channels_no]

    # update headers
    eeg_file_header[:channels_no] = channels_no - length(channels)
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

    # create new dataset    
    eeg_new = EEG(eeg_file_header, eeg_signal_header, eeg_signals)
    
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
    filter_type in [:lp, :hp, :bp, :bs] || throw(ArgumentError("""Filter type must be ":bp", ":hp", ":bp" or ":bs"."""))

    signal = eeg.eeg_signals
    fs = eeg.eeg_signal_header[:sampling_rate][1]

    signal_filtered = zeros(size(signal))

    for idx in 1:channels_no
        signal_filtered[idx, :] = signal_filter_butter(signal[idx, :], filter_type=filter_type, cutoff=cutoff, fs=fs, poles=poles)
    end

    eeg_new = EEG(eeg.eeg_file_header, eeg.eeg_signal_header, signal_filtered)

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
    eeg_new = EEG(eeg.eeg_file_header, eeg.eeg_signal_header, signal_der)

    return eeg_new
end

"""
    eeg_total_power(eeg)

Calculates total power for each the `eeg` signal channels.

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

Calculates absolute band power between frequencies `f1` and `f2` for each the `eeg` signal channels.

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

Returns FFT and DFT sample frequencies for a DFT for each the `eeg` signal channels.

# Arguments

- `eeg::EEG` - EEG object
"""
function eeg_make_spectrum(eeg)
    fs = eeg.eeg_signal_header[:sampling_rate][1]
    signal_fft, signal_sf = signal_make_spectrum(eeg.eeg_signals, fs)

    return signal_fft, signal_sf
end

"""
    eeg_detrend(eeg, type=:linear; channels=[])

Removes linear trend for each the `eeg` signal channels.

# Arguments

- `eeg::EEG` - EEG object
- `type::Symbol[:linear, :constant]`, optional
    - `linear` - the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `constant` - the mean of `signal` is subtracted
"""
function eeg_detrend(eeg, type=:linear)
    trend in [:linear, :constant] || throw(ArgumentError("""Trend type must be ":linear" or ":constant"."""))

    signal_det = signal_detrend(eeg.eeg_signals)
    eeg_new = EEG(eeg.eeg_file_header, eeg.eeg_signal_header, signal_det)

    return eeg_new
end

"""
    eeg_draw_head(p, loc_x, loc_y)

Draws head over a topographical plot `p`.

# Arguments

- `p::Plot` - toppgraphical plot
- `loc_x::Vector{Float64` - vector of x electrode position
- `loc_y::Vector{Float64` - vector of y electrode position
"""
function eeg_draw_head(p, loc_x::Vector{Float64}, loc_y::Vector{Float64})
    pts = Plots.partialcircle(0, 2π, 100, maximum(loc_x))
    x, y = Plots.unzip(pts)
    x = x .* 1.1
    y = y .* 0.91
    head = Shape(x, y)
    nose = Shape([(-0.1, maximum(y)), (0, maximum(y) + 0.1 * maximum(y)), (0.1, maximum(y))])
    ear_l = Shape([(minimum(x), -0.1), (minimum(x) + 0.1 * minimum(x), -0.1), (minimum(x) + 0.1 * minimum(x), 0.1), (minimum(x), 0.1)])
    ear_r = Shape([(maximum(x), -0.1), (maximum(x) + 0.1 * maximum(x), -0.1), (maximum(x) + 0.1 * maximum(x), 0.1), (maximum(x), 0.1)])
    plot!(p, head, fill=nothing, label="")
    plot!(p, nose, fill=nothing, label="")
    plot!(p, ear_l, fill=nothing, label="")
    plot!(p, ear_r, fill=nothing, label="")
end

"""
    eeg_rereference_channel(eeg, reference)

Re-references the `eeg` signal channels to specific signal channel.

# Arguments

- `eeg::EEG` - EEG object
- `reference::Float64` - index of channels used as reference; if multiple channels are specififed, their average is used as the reference
"""
function eeg_rereference_channel(eeg::EEG, reference_idx)
    signal_rereferenced = signal_rereference_channel(eeg.eeg_signals, reference_idx)
    eeg.eeg_file_header[:reference_type] = "channel"
    eeg.eeg_file_header[:reference_channel] = reference_idx
    eeg_new = EEG(eeg.eeg_file_header, eeg.eeg_signal_header, signal_rereferenced)

    return eeg_new
end

"""
    eeg_save(eeg, file_name)

Saves the `eeg` object to `file_name` file.

# Arguments

- `eeg::EEG` - EEG object
- `file_name::String` - file name
"""
function eeg_save(eeg::EEG, file_name)
    save_object(file_name, eeg)
end

"""
    eeg_load(file_name)

Loads the `eeg` object from `file_name` file.

# Arguments

- `file_name::String` - file name
"""
function eeg_load(file_name)
    eeg = load_object(file_name)
    return eeg
end