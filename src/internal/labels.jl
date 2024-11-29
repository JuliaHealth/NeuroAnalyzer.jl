function _clean_labels(clabels::Vector{String})::Vector{String}
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
    l = replace.(l, "MEG" => "MEG ")
    l = replace.(l, "MEG 0" => "MEG ")
    l = replace.(l, "  " => " ")
    return l
end

function _clean_meg_labels(clabels::Vector{String})::Vector{String}
    l = deepcopy(clabels)
    l = replace.(l, "MEG" => "MEG ")
    l = replace.(l, "EEG" => "EEG ")
    l = replace.(l, "EOG" => "EOG ")
    l = replace.(l, "EMG" => "EMG ")
    l = replace.(l, "  " => " ")
    l = replace.(l, "MEG 0" => "MEG ")
    l = replace.(l, "EEG 0" => "EEG ")
    l = replace.(l, "EOG 0" => "EOG ")
    l = replace.(l, "EMG 0" => "EMG ")
    l = replace.(l, "MEG 0" => "MEG ")
    l = replace.(l, "EEG 0" => "EEG ")
    l = replace.(l, "EOG 0" => "EOG ")
    l = replace.(l, "EMG 0" => "EMG ")
    return l
end

function _clean_eeg_labels(clabels::Vector{String})::Vector{String}
    l = deepcopy(clabels)
    l = replace.(l, "EEG " => "")
    l = replace.(l, "  " => " ")
    l = replace.(l, "EOG EOG" => "EOG")
    l = replace.(l, "ECG EKG" => "ECG")
    l = replace.(l, "EEG 0" => "EEG ")
    l = replace.(l, "EMG 0" => "EMG ")
    l = replace.(l, "EEG 0" => "EEG ")
    l = replace.(l, "EOG 0" => "EOG ")
    l = replace.(l, "EEG 0" => "EEG ")
    l = replace.(l, "EMG 0" => "EMG ")
    l = replace.(l, "EEG 0" => "EEG ")
    l = replace.(l, "EOG 0" => "EOG ")
    return l
end

function _gen_clabels(obj::NeuroAnalyzer.NEURO, c::Symbol)::Vector{String}
    c = _get_component(obj, c)
    clabels = Vector{String}()
    for idx in axes(c, 1)
        push!(clabels, lpad(string(idx), length(string(size(c, 1))), "0"))
    end
    return clabels
end

function _gen_clabels(c::Union{AbstractVector, AbstractArray})::Vector{String}
    clabels = Vector{String}()
    if ndims(c) == 1
        push!(clabels, "1")
    else
        for idx in axes(c, 1)
            push!(clabels, lpad(string(idx), length(string(size(c, 1))), "0"))
        end
    end
    return clabels
end
