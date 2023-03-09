export study_n
export study_channel_n
export study_epoch_n
export study_epoch_len
export study_sr

"""
    study_n(study)

Return number of NeuroAnalyzer NEURO objects in the study.

# Arguments

- `study::NeuroAnalyzer.STUDY`

# Returns

- `n::Int64`
"""
function study_n(study::NeuroAnalyzer.STUDY)
    return length(study.objects)
end

"""
    study_channel_n(study)

Return number of channels per NeuroAnalyzer NEURO object in the study.

# Arguments

- `study::NeuroAnalyzer.STUDY`

# Returns

- `n::Int64`
"""
function study_channel_n(study::NeuroAnalyzer.STUDY)
    return channel_n(study.objects[1])
end

"""
    study_epoch_n(study)

Return number of epochs per NeuroAnalyzer NEURO object in the study.

# Arguments

- `study::NeuroAnalyzer.STUDY`

# Returns

- `n::Int64`
"""
function study_epoch_n(study::NeuroAnalyzer.STUDY)
    return epoch_n(study.objects[1])
end

"""
    study_epoch_len(study)

Return length of epochs per NeuroAnalyzer NEURO object in the study.

# Arguments

- `study::NeuroAnalyzer.STUDY`

# Returns

- `len::Int64`
"""
function study_epoch_len(study::NeuroAnalyzer.STUDY)
    return epoch_len(study.objects[1])
end

"""
    study_sr(study)

Return sampling rate of NeuroAnalyzer NEURO objects in the study.

# Arguments

- `study::NeuroAnalyzer.STUDY`

# Returns

- `sr::Int64`
"""
function study_sr(study::NeuroAnalyzer.STUDY)
    return sr(study.objects[1])
end
