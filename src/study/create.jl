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
    @assert length(obj) == length(group) "Length of OBJs and groups must be equal."

    ch_n = nchannels(obj[1])
    ep_n = nepochs(obj[1])
    ep_len = epoch_len(obj[1])
    fs = sr(obj[1])

    for idx in eachindex(obj)
        @assert nepochs(obj[idx]) == ep_n "All OBJs in the study must have the same number of epochs."
        @assert epoch_len(obj[idx]) == ep_len "All OBJs in the study must have the same length of epochs."
        @assert sr(obj[idx]) == fs "All OBJs in the study must have the same sampling rate."
        @assert nchannels(obj[idx]) == ch_n "All OBJs in the study must have the same number of channels."
    end

    return NeuroAnalyzer.STUDY(Dict{Symbol, Any}(), obj, group)
    
end