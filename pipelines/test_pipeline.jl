"""
    test_pipeline(eeg::NeuroAnalyzer.NEURO)

Basic NeuroAnalyzer EEG processing pipeline.

# Arguments

- `eeg::NeuroAnalyzer.NEURO`
- `dc::Real=50`: DC frequency
- `lp::Real=45`: low-pass filtering cutoff
- `hp::Real=0.1`: low-pass filtering cutoff

# Returns

- `eeg_processed::NeuroAnalyzer.NEURO`
"""
function test_pipeline(eeg::NeuroAnalyzer.NEURO, dc::Real=50, lp::Real=45, hp::Real=0.1)

    # re-reference to common average
    eeg_processed = reference_car(eeg)

    # DC filtering
    filter!(eeg_processed, fprototype=:iirnotch, cutoff=dc, bw=2)
    # HP filtering
    filter!(eeg_processed, fprototype=:butterworth, ftype=:hp, cutoff=hp, order=8)
    # LP filtering
    filter!(eeg_processed, fprototype=:butterworth, ftype=:lp, cutoff=lp, order=8)
    
    return eeg_processed
end