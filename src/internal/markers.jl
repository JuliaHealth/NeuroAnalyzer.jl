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

function _has_markers(channel_types::Vector{String})
    markers = false
    markers_channel = 0
    if "mrk" in channel_types
        markers = true
        markers_channel = nothing
        for ch_idx in eachindex(channel_types)
            channel_types[ch_idx] == "mrk" && (markers_channel = ch_idx)
        end
    end
    return markers, markers_channel
end

function _has_markers(obj::NeuroAnalyzer.NEURO)
    return nrow(obj.markers) > 0 ? true : false
end

function _a2df(annotations::Vector{String})
    # convert EDF/BDF annotations to markers DataFrame
    annotations = replace.(annotations, "\x14\x14\0" => "|")
    annotations = replace.(annotations, "\x14\x14" => "|")
    annotations = replace.(annotations, "\x14" => "|")
    annotations = replace.(annotations, "\0" => "")
    a_start = Vector{Float64}()
    a_event = Vector{String}()
    # what about annotations containing event duration?
    for idx in length(annotations):-1:1
        length(annotations[idx]) == 0 && deleteat!(annotations, idx)
    end
    for idx in length(annotations):-1:1
        length(split(annotations[idx], "|")) < 3 && deleteat!(annotations, idx)
    end
    for idx in 1:length(annotations)
        s = split(annotations[idx], "|")
        if length(s) > 3
            push!(a_start, parse(Float64, strip(s[2])))
            push!(a_event, strip(s[3]))
        elseif length(s) < 4
            push!(a_start, parse(Float64, strip(s[1])))
            push!(a_event, strip(s[2]))
        end
    end
    return DataFrame(:id=>repeat([""], length(a_event)), :start=>a_start, :length=>zeros(Int64, length(a_event)), :description=>a_event, :channel=>zeros(Int64, length(a_event)))
end
