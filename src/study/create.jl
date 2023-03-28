export create_study

"""
    study(obj, group)

Create NeuroAnalyzer STUDY object.

# Arguments

- `obj::Vector{NeuroAnalyzer.NEURO}`
- `group::Vector{Symbol}`

# Returns

- `study::NeuroAnalyzer.STUDY`
"""
function create_study(obj::Vector{NeuroAnalyzer.NEURO}, group::Vector{Symbol})
    length(obj) == length(group) || throw(ArgumentError("Length of OBJs and groups must be equal."))

    ch_n = channel_n(obj[1])
    ep_n = epoch_n(obj[1])
    ep_len = epoch_len(obj[1])
    fs = sr(obj[1])

    for idx in 1:length(obj)
        epoch_n(obj[idx]) == ep_n || throw(ArgumentError("All OBJs in the study must have the same number of epochs."))
        epoch_len(obj[idx]) == ep_len || throw(ArgumentError("All OBJs in the study must have the same length of epochs."))
        sr(obj[idx]) == fs || throw(ArgumentError("All OBJs in the study must have the same sampling rate."))
        channel_n(obj[idx]) == ch_n || throw(ArgumentError("All OBJs in the study must have the same number of channels."))
    end

    return NeuroAnalyzer.STUDY(Dict{Symbol, Any}(), obj, group)
    
end