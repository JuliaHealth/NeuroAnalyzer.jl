function _clean_labels(clabels::Vector{String})
    clabels = replace.(lowercase.(clabels), "eeg " => "")
    clabels = replace.(lowercase.(clabels), "eog eog" => "EOG")
    clabels = replace.(lowercase.(clabels), "ecg ekg" => "ECG")
    clabels = replace.(lowercase.(clabels), "edf " => "")
    clabels = replace.(lowercase.(clabels), "bdf " => "")
    clabels = replace.(clabels, ".." => "")
    clabels = replace.(clabels, "." => "")
    return clabels
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
