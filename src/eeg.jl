"""
    eeg_plot(eeg; channels=[], labels=[], normalize=false, xlabel="Time [s]", ylabel="Channels")

Plots `signal` matrix.

# Arguments
- `eee::EEG` - EEG object
- `t::Vector{Float64} - the time vector
- `offset::Float64` - displayed segment offset in samples
- `channels::Float64` - channels to be plotted (all if empty), vector or range
- `labels::Vector{String}` - channel labels vector
- `xlabel::String` - x-axis label
- `ylabel::String` - y-axis lable
- `normalize::Bool` - normalize the `signal` prior to calculations
- `remove_dc::Bool` - demean the `signal` prior to calculations
- `detrend::Bool` - detrend the `signal` prior to calculations
- `derivative::Bool` - derivate `signal` prior to calculations
- `taper::Bool` - taper the `signal` with `taper`-window prior to calculations
"""
function eeg_plot(eeg::EEG; t=nothing, offset=0, channels=[], labels=[], xlabel="Time [s]", ylabel="Channels", normalize=true, remove_dc=false, detrend=false, derivative=false, taper=nothing)
    
    if typeof(channels) == UnitRange{Int64}
        channels = collect(channels)
    end

    if typeof(t) == UnitRange{Int64}
        t = collect(t)
    end

    channels_no = eeg.eeg_file_header[:channels_no]

    # drop channels not in the list
    channels_to_drop = collect(1:channels_no)
    if length(channels) > 1
        for idx in length(channels):-1:1
            channels_to_drop = deleteat!(channels_to_drop, channels[idx])
        end
        eeg = eeg_drop_channel(eeg, channels_to_drop)
    end

    signal = eeg.eeg_signals
    labels = eeg.eeg_signal_header[:labels]
    fs = eeg.eeg_signal_header[:sampling_rate][1]

    # default time is 5 seconds
    t === nothing && (t = collect(0:1/fs:5))

    p = signal_plot(t, signal, offset=offset, channels=[], labels=labels, xlabel=xlabel, ylabel=ylabel, normalize=normalize, remove_dc=remove_dc, detrend=detrend, derivative=derivative, taper=taper)

    return p
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

    length(channels) > 1 && (channels = sort!(channels, rev=true))

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
    eeg_filter_butter(eeg; channels=[], filter_type, cutoff, fs, poles=8)

Filters `eeg` channels using Butterworth filter.

# Arguments
- `eeg::EEG` - EEG object
- `channels::Float64` - channels to filter, vector of numbers or range
- `filter_type::Symbol[:lp, :hp, :bp, :bs]` - filter type
- `cutoff::Float64` - filter cutoff in Hz (tuple or vector for `:bp` and `:bs`)
- `fs::Float64` - sampling rate
- `poles::Int` - filter pole
"""
function eeg_filter_butter(eeg::EEG; channels=[], filter_type, cutoff, poles=8)
    filter_type in [:lp, :hp, :bp, :bs] || throw(ArgumentError("""Filter type must be ":bp", ":hp", ":bp" or ":bs"."""))

    if typeof(channels) == UnitRange{Int64}
        channels = collect(channels)
    end

    channels_no = eeg.eeg_file_header[:channels_no]

    # drop channels not in the list
    channels_to_drop = collect(1:channels_no)
    if length(channels) > 1
        for idx in length(channels):-1:1
            channels_to_drop = deleteat!(channels_to_drop, channels[idx])
        end
        eeg = eeg_drop_channel(eeg, channels_to_drop)
    end

    signal = eeg.eeg_signals
    fs = eeg.eeg_signal_header[:sampling_rate][1]

    signal_filtered = zeros(size(signal))

    for idx in 1:channels_no
        signal_filtered[idx, :] = signal_filter_butter(signal[idx, :], filter_type=filter_type, cutoff=cutoff, fs=fs, poles=poles)
    end

    eeg_new = EEG(eeg.eeg_file_header, eeg.eeg_signal_header, signal_filtered)

    return eeg_new
end