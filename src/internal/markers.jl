function _delete_markers(markers::DataFrame, segment::Tuple{Int64, Int64})
    for marker_idx in nrow(markers):-1:1
        markers[marker_idx, :start] in segment[1]:segment[2] && deleteat!(markers, marker_idx)
    end
    return markers
end

function _shift_markers(m::DataFrame, pos::Int64, offset::Int64)
    markers = deepcopy(m)
    for marker_idx in 1:nrow(markers)
        markers[marker_idx, :start] > pos && (markers[marker_idx, :start] -= offset)
    end
    return markers
end

function _get_epoch_markers(obj::NeuroAnalyzer.NEURO)
    return round.(s2t.(collect(1:epoch_len(obj):epoch_len(obj) * epoch_n(obj)), sr(obj)), digits=2)
end
