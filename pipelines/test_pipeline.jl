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
function test_pipeline(eeg::NeuroAnalyzer.NEURO, pl_frq::Real=50, lp::Real=45, hp::Real=0.1)

    # DC filtering
    remove_powerline!(eeg_processed, pl_frq=pl_frq)

    # HP filtering
    NeuroAnalyzer.filter!(eeg_processed, fprototype=:fir, ftype=:hp, cutoff=hp, order=91)

    # LP filtering
    NeuroAnalyzer.filter!(eeg_processed, fprototype=:fir, ftype=:lp, cutoff=lp, order=91)

    # re-reference to common average
    eeg_processed = reference_avg(eeg)

    return eeg_processed

end