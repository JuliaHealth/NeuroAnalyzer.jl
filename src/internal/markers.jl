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
    mrk = replace.(annotations, "\x14\x14\0" => "|")
    mrk = replace.(mrk, "\x14\x14" => "|")
    mrk = replace.(mrk, "\x14" => "|")
    mrk = replace.(mrk, "\x15" => "|")
    mrk = replace.(mrk, "\0" => "")
    mrk = replace.(mrk, r"\|$" => "")
    a_start = Vector{Float64}()
    a_length = Vector{Float64}()
    a_event = Vector{String}()

    # remove empty
    for idx in length(mrk):-1:1
        (length(mrk[idx]) == 0 || occursin('|', mrk[idx]) == false) && deleteat!(mrk, idx)
    end
    # for idx in length(mrk):-1:1
    #     length(split(mrk[idx], "|")) < 3 && deleteat!(mrk, idx)
    # end

    if length(mrk) == 1
        s = split(mrk[1], "|")
        for idx in length(s):-1:1
            s[idx] == "" && deleteat!(s, idx)
        end
        if length(s) % 3 == 0
            for idx in 1:3:length(s) ÷ 3
                push!(a_start, parse(Float64, strip(s[idx])))
                push!(a_length, parse(Float64, strip(s[idx + 1])))
                push!(a_event, strip(s[idx + 2]))
            end
        else
            offset = parse(Float64, strip(s[1]))
            deleteat!(s, 1)
            for idx in 1:3:length(s) ÷ 3
                push!(a_start, parse(Float64, strip(s[idx])))
                push!(a_length, parse(Float64, strip(s[idx + 1])))
                push!(a_event, strip(s[idx + 2]))
            end
        end
        all(isascii.(a_event)) == false && _info("Unicode labels were not converted.")
        return DataFrame(:id=>string.(collect(1:length(a_event))), :start=>a_start, :length=>a_length, :description=>a_event, :channel=>zeros(Int64, length(a_event)))
    else
        for idx in 1:length(mrk)
            s = split(mrk[idx], "|")
            # drop the first column if it contains annotation number
            s[1][1] == '+' && s[2][1] == '+' && (s = s[2:end])
            if length(s) == 3
                push!(a_start, parse(Float64, strip(s[1])))
                push!(a_length, parse(Float64, strip(s[2])))
                push!(a_event, strip(s[3]))
            elseif length(s) == 2
                push!(a_start, parse(Float64, strip(s[1])))
                push!(a_length, 0.0)
                push!(a_event, strip(s[2]))
            end
        end
        all(isascii.(a_event)) == false && _info("Unicode labels were not converted.")
        return DataFrame(:id=>string.(collect(1:length(a_event))), :start=>a_start, :length=>a_length, :description=>a_event, :channel=>zeros(Int64, length(a_event)))
    end
end
