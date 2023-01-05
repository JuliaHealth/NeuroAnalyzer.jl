"""
    eeg_study_create(eeg, group)

Create EEG study object.

# Arguments

- `eeg::Vector{NeuroAnalyzer.EEG}`
- `group::Vector{Symbol}`

# Returns

- `study::NeuroAnalyzer.STUDY`
"""
function eeg_study_create(eeg::Vector{NeuroAnalyzer.EEG}, group::Vector{Symbol})
    length(eeg) == length(group) || throw(ArgumentError("Length of EEGs and groups must be equal."))

    channel_n = eeg_channel_n(eeg[1])
    epoch_n = eeg_epoch_n(eeg[1])
    epoch_len = eeg_epoch_len(eeg[1])
    sr = eeg_sr(eeg[1])

    for eeg_idx in 1:length(eeg)
        eeg_epoch_n(eeg[eeg_idx]) == epoch_n || throw(ArgumentError("All EEG object in the study must have the same number of epochs."))
        eeg_epoch_len(eeg[eeg_idx]) == epoch_len || throw(ArgumentError("All EEG object in the study must have the same length of epochs."))
        eeg_sr(eeg[eeg_idx]) == sr || throw(ArgumentError("All EEG object in the study must have the same sampling rate."))
        eeg_channel_n(eeg[eeg_idx]) == channel_n || throw(ArgumentError("All EEG object in the study must have the same number of channels."))
    end

    study = STUDY(Dict{Symbol, Any}(), eeg, group)
    
    return study
end

"""
    eeg_study_n(study)

Return number of EEG object in the study.

# Arguments

- `study::NeuroAnalyzer.STUDY`

# Returns

- `n::Int64`
"""
function eeg_study_n(study::NeuroAnalyzer.STUDY)
    return length(study.study_eeg)
end

"""
    eeg_study_channel_n(study)

Return number of channels per EEG object in the study.

# Arguments

- `study::NeuroAnalyzer.STUDY`

# Returns

- `n::Int64`
"""
function eeg_study_channel_n(study::NeuroAnalyzer.STUDY)
    return eeg_channel_n(study.study_eeg[1])
end

"""
    eeg_study_epoch_n(study)

Return number of epochs per EEG object in the study.

# Arguments

- `study::NeuroAnalyzer.STUDY`

# Returns

- `n::Int64`
"""
function eeg_study_epoch_n(study::NeuroAnalyzer.STUDY)
    return eeg_epoch_n(study.study_eeg[1])
end

"""
    eeg_study_epoch_len(study)

Return length of epochs per EEG object in the study.

# Arguments

- `study::NeuroAnalyzer.STUDY`

# Returns

- `len::Int64`
"""
function eeg_study_epoch_len(study::NeuroAnalyzer.STUDY)
    return eeg_epoch_len(study.study_eeg[1])
end

"""
    eeg_study_sr(study)

Return sampling rate of EEG objects in the study.

# Arguments

- `study::NeuroAnalyzer.STUDY`

# Returns

- `sr::Int64`
"""
function eeg_study_sr(study::NeuroAnalyzer.STUDY)
    return eeg_sr(study.study_eeg[1])
end
