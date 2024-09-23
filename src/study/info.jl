export obj_n
export nchannels
export nepochs
export epoch_len
export sr

"""
    obj_n(study)

Return number of NeuroAnalyzer NEURO objects in the study.

# Arguments

- `study::NeuroAnalyzer.STUDY`

# Returns

- `n::Int64`
"""
function obj_n(study::NeuroAnalyzer.STUDY)::Int64

    return length(study.objects)

end

"""
    nchannels(study)

Return number of channels per NeuroAnalyzer NEURO object in the study.

# Arguments

- `study::NeuroAnalyzer.STUDY`

# Returns

- `n::Int64`
"""
function nchannels(study::NeuroAnalyzer.STUDY)::Int64

    return nchannels(study.objects[1])

end

"""
    nepochs(study)

Return number of epochs per NeuroAnalyzer NEURO object in the study.

# Arguments

- `study::NeuroAnalyzer.STUDY`

# Returns

- `n::Int64`
"""
function nepochs(study::NeuroAnalyzer.STUDY)::Int64

    return nepochs(study.objects[1])

end

"""
    epoch_len(study)

Return length of epochs per NeuroAnalyzer NEURO object in the study.

# Arguments

- `study::NeuroAnalyzer.STUDY`

# Returns

- `len::Int64`
"""
function epoch_len(study::NeuroAnalyzer.STUDY)::Int64

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
function sr(study::NeuroAnalyzer.STUDY)::Int64

    return sr(study.objects[1])

end
