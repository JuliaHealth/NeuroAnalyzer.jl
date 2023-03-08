export study_create

"""
    study_create(obj, group)

Create NeuroAnalyzer STUDY object.

# Arguments

- `obj::Vector{NeuroAnalyzer.NEURO}`
- `group::Vector{Symbol}`

# Returns

- `study::NeuroAnalyzer.STUDY`
"""
function study_create(obj::Vector{NeuroAnalyzer.NEURO}, group::Vector{Symbol})
    length(obj) == length(group) || throw(ArgumentError("Length of OBJs and groups must be equal."))

    ch_n = channel_n(obj[1])
    ep_n = epoch_n(obj[1])
    epoch_len = epoch_len(obj[1])
    sr = sr(obj[1])

    for idx in 1:length(obj)
        epoch_n(obj[idx]) == ep_n || throw(ArgumentError("All OBJs in the study must have the same number of epochs."))
        epoch_len(obj[idx]) == epoch_len || throw(ArgumentError("All OBJs in the study must have the same length of epochs."))
        sr(obj[idx]) == sr || throw(ArgumentError("All OBJs in the study must have the same sampling rate."))
        channel_n(obj[idx]) == ch_n || throw(ArgumentError("All OBJs in the study must have the same number of channels."))
    end

    study = STUDY(Dict{Symbol, Any}(), obj, group)
    
    return study
end