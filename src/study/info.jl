export size
export channel_n
export epoch_n
export epoch_len
export sr

"""
    size(study)

Return number of NeuroAnalyzer NEURO objects in the study.

# Arguments

- `study::NeuroAnalyzer.STUDY`

# Returns

- `n::Int64`
"""
function size(study::NeuroAnalyzer.STUDY)
    return length(study.objects)
end

"""
    channel_n(study)

Return number of channels per NeuroAnalyzer NEURO object in the study.

# Arguments

- `study::NeuroAnalyzer.STUDY`

# Returns

- `n::Int64`
"""
function channel_n(study::NeuroAnalyzer.STUDY)
    return channel_n(study.objects[1])
end

"""
    epoch_n(study)

Return number of epochs per NeuroAnalyzer NEURO object in the study.

# Arguments

- `study::NeuroAnalyzer.STUDY`

# Returns

- `n::Int64`
"""
function epoch_n(study::NeuroAnalyzer.STUDY)
    return epoch_n(study.objects[1])
end

"""
    epoch_len(study)

Return length of epochs per NeuroAnalyzer NEURO object in the study.

# Arguments

- `study::NeuroAnalyzer.STUDY`

# Returns

- `len::Int64`
"""
function epoch_len(study::NeuroAnalyzer.STUDY)
    return epoch_len(study.objects[1])
end

"""
    sr(study)

Return sampling rate of NeuroAnalyzer NEURO objects in the study.

# Arguments

- `study::NeuroAnalyzer.STUDY`

# Returns

- `sr::Int64`
"""
function sr(study::NeuroAnalyzer.STUDY)
    return sr(study.objects[1])
end
