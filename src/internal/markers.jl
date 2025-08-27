function _delete_markers(markers::DataFrame, segment::Tuple{Real, Real}, fs::Int64)::DataFrame
    for marker_idx in nrow(markers):-1:1
        round(Int64, fs * markers[marker_idx, :start]) in segment[1]:segment[2] && deleteat!(markers, marker_idx)
    end
    return markers
end

function _shift_markers(m::DataFrame, pos::Real, offset::Real, fs::Int64)::DataFrame
    markers = deepcopy(m)
    for marker_idx in 1:nrow(markers)
        round(Int64, fs * markers[marker_idx, :start]) > pos && (markers[marker_idx, :start] -= (offset / fs))
    end
    return markers
end

function _get_epoch_markers(obj::NeuroAnalyzer.NEURO)::Vector{Float64}
    return round.(s2t.(collect(1:epoch_len(obj):epoch_len(obj) * nepochs(obj)), sr(obj)), digits=3)
end

function _has_markers(channel_types::Vector{String})::Tuple{Bool, Int64}
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

_has_markers(obj::NeuroAnalyzer.NEURO)::Bool = nrow(obj.markers) > 0 ? true : false

function _a2df(annotations::Vector{String})::DataFrame
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
        !(length(mrk[idx]) == 0 || occursin('|', mrk[idx])) && deleteat!(mrk, idx)
    end

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
            # TO DO: use offset if provided
            offset = parse(Float64, strip(s[1]))
            deleteat!(s, 1)
            for idx in 1:3:length(s) ÷ 3
                push!(a_start, parse(Float64, strip(s[idx])))
                push!(a_length, parse(Float64, strip(s[idx + 1])))
                push!(a_event, strip(s[idx + 2]))
            end
        end
        !all(isascii.(a_event)) && _warn("Unicode labels were not converted.")
        return DataFrame(:id=>string.(collect(eachindex(a_event))), :start=>a_start, :length=>a_length, :value=>a_event, :channel=>zeros(Int64, length(a_event)))
    else
        for idx in eachindex(mrk)
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
        evnts = unique(a_event)
        id = zeros(Int64, length(a_event))
        for idx1 in eachindex(a_event)
            for idx2 in eachindex(unique(a_event))
                id[idx1] = findfirst(a_event[idx1] .== unique(a_event) )
            end
        end
        !all(isascii.(a_event)) && _warn("Unicode labels were not converted.")
        return DataFrame(:id=>string.(id), :start=>a_start, :length=>a_length, :value=>a_event, :channel=>zeros(Int64, length(a_event)))
    end
end
