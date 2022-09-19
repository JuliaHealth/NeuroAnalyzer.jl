"""
    test_pipeline(eeg::NeuroAnalyzer.EEG)

Basic NeuroAnalyzer EEG processing pipeline.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `dc::Real=50`: DC frequency
- `lp::Real=45`: low-pass filtering cutoff
- `hp::Real=0.1`: low-pass filtering cutoff

# Returns

- `eeg_processed::NeuroAnalyzer.EEG`
"""
function test_pipeline(eeg::NeuroAnalyzer.EEG, dc::Real=50, lp::Real=45, hp::Real=0.1)

    # re-reference to common average
    eeg_processed = eeg_reference_car(eeg)

    # DC filtering
    eeg_filter!(eeg_processed, fprototype=:iirnotch, cutoff=dc, bw=2)
    # HP filtering
    eeg_filter!(eeg_processed, fprototype=:butterworth, ftype=:hp, cutoff=hp, order=8)
    # LP filtering
    eeg_filter!(eeg_processed, fprototype=:butterworth, ftype=:lp, cutoff=lp, order=8)
    
    return eeg_processed
end