function _clean_labels(clabels::Vector{String})
    l = deepcopy(clabels)
    l = replace.(l, "eeg " => "")
    l = replace.(l, "EEG " => "")
    l = replace.(l, "eeg" => "")
    l = replace.(l, "EEG" => "")
    l = replace.(l, "eog eog" => "EOG")
    l = replace.(l, "EOG EOG" => "EOG")
    l = replace.(l, "ecg ekg" => "ECG")
    l = replace.(l, "ECG EKG" => "ECG")
    l = replace.(l, "edf" => "")
    l = replace.(l, "edf " => "")
    l = replace.(l, "EDF" => "")
    l = replace.(l, "EDF " => "")
    l = replace.(l, "bdf" => "")
    l = replace.(l, "BDF" => "")
    l = replace.(l, "bdf " => "")
    l = replace.(l, "BDF " => "")
    l = replace.(l, ".." => "")
    l = replace.(l, "." => "")
    return l
end

function _gen_clabels(obj::NeuroAnalyzer.NEURO, c::Symbol)
    c = _get_component(obj, c)
    clabels = Vector{String}()
    for idx in 1:size(c, 1)
        push!(clabels, lpad(string(idx), length(string(size(c, 1))), "0"))
    end
    return clabels
end

function _gen_clabels(c::Union{AbstractVector, AbstractArray})
    clabels = Vector{String}()
    if ndims(c) == 1
        push!(clabels, "1")
    else
        for idx in 1:size(c, 1)
            push!(clabels, lpad(string(idx), length(string(size(c, 1))), "0"))
        end
    end
    return clabels
end
