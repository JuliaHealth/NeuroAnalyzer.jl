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
