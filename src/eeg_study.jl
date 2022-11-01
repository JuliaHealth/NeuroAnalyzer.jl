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
    study = STUDY(Dict{Symbol, Any}(), eeg, group)
    
    return study
end