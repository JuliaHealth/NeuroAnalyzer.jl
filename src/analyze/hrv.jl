export hrv

function hrv(obj::NeuroAnalyzer.NEURO)

    "ecg" in obj.header.recording[:channel_type] || throw(ArgumentError("OBJ does not contain ECG channel."))
    ch = findfirst(obj.header.recording[:channel_type] .== "ecg")
    NeuroAnalyzer._info("ECG channel found: $ch")
    ecg = eeg.data[ch, :, :][:]
    nn_idx, _ = findpeaks1d(ecg, height=mean(ecg) + 2*std(ecg))

    # convert to ms
    nn_diff = diff(nn_idx) ./ sr(eeg) * 1000
    
    NeuroAnalyzer._info("Detected NN intervals: $(length(nn_diff))")

    return nn_diff

end

mean(rr_diff)
std(rr_diff)

var(rr_diff)

median(rr_diff)

SDNN, the standard deviation of NN intervals. Often calculated over a 24-hour period. SDANN, the standard deviation of the average NN intervals calculated over short periods, usually 5 minutes. SDANN is therefore a measure of changes in heart rate due to cycles longer than 5 minutes. SDNN reflects all the cyclic components responsible for variability in the period of recording, therefore it represents total variability.
RMSSD ("root mean square of successive differences"), the square root of the mean of the squares of the successive differences between adjacent NNs.[39]
SDSD ("standard deviation of successive differences"), the standard deviation of the successive differences between adjacent NNs.[39]
NN50, the number of pairs of successive NNs that differ by more than 50 ms.
pNN50, the proportion of NN50 divided by total number of NNs.
NN20, the number of pairs of successive NNs that differ by more than 20 ms.[40]
pNN20, the proportion of NN20 divided by total number of NNs.