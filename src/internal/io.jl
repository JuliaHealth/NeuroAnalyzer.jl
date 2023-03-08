function _has_markers(channel_types::Vector{String})
    markers = false
    markers_channel = 0
    if "mrk" in channel_types
        markers = true
        markers_channel = nothing
        for channel_idx in eachindex(channel_types)
            channel_types[channel_idx] == "mrk" && (markers_channel = channel_idx)
        end
    end
    return markers, markers_channel
end

function _set_channel_types(clabels::Vector{String})
    channel_names = ["af3", "af4", "af7", "af8", "afz", "c1", "c2", "c3", "c4", "c5", "c6", "cp1", "cp2", "cp3", "cp4", "cp5", "cp6", "cpz", "cz", "f1", "f10", "f2", "f3", "f4", "f5", "f6", "f7", "f8", "f9", "fc1", "fc2", "fc3", "fc4", "fc5", "fc6", "fcz", "fp1", "fp2", "fpz", "ft10", "ft7", "ft8", "ft9", "fz", "nz", "o1", "o2", "oz", "p1", "p10", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "po3", "po4", "po7", "po8", "poz", "pz", "t10", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "tp10", "tp7", "tp8", "tp9"]
    ref_channels = ["a1", "a2", "m1", "m2", "pg1", "pg2"]
    eog_channel = ["e", "e1", "e2"]
    channel_type = repeat(["???"], length(clabels))
    for idx in eachindex(clabels)
        occursin("meg", lowercase(clabels[idx])) && (channel_type[idx] = "meg")
        occursin("ecg", lowercase(clabels[idx])) && (channel_type[idx] = "ecg")
        occursin("ekg", lowercase(clabels[idx])) && (channel_type[idx] = "ecg")
        occursin("eog", lowercase(clabels[idx])) && (channel_type[idx] = "eog")
        occursin("rr", lowercase(clabels[idx])) && (channel_type[idx] = "misc")
        occursin("mic", lowercase(clabels[idx])) && (channel_type[idx] = "misc")
        occursin("flw", lowercase(clabels[idx])) && (channel_type[idx] = "misc")
        occursin("tho", lowercase(clabels[idx])) && (channel_type[idx] = "misc")
        occursin("abd", lowercase(clabels[idx])) && (channel_type[idx] = "misc")
        occursin("sao2", lowercase(clabels[idx])) && (channel_type[idx] = "misc")
        occursin("sa02", lowercase(clabels[idx])) && (channel_type[idx] = "misc")
        occursin("plr", lowercase(clabels[idx])) && (channel_type[idx] = "misc")
        occursin("body", lowercase(clabels[idx])) && (channel_type[idx] = "misc")
        occursin("ux", lowercase(clabels[idx])) && (channel_type[idx] = "misc")
        for idx2 in eachindex(eog_channel)
            occursin(channel_names[idx2], lowercase(clabels[idx])) && (channel_type[idx] = "eog")
        end
        occursin("emg", lowercase(clabels[idx])) && (channel_type[idx] = "emg")
        in(lowercase(clabels[idx]), ref_channels) && (channel_type[idx] = "ref")
        for idx2 in eachindex(ref_channels)
            occursin(ref_channels[idx2], lowercase(clabels[idx])) && (channel_type[idx] = "ref")
        end
        occursin("mark", lowercase(clabels[idx])) && (channel_type[idx] = "mrk")
        occursin("marker", lowercase(clabels[idx])) && (channel_type[idx] = "mrk")
        occursin("markers", lowercase(clabels[idx])) && (channel_type[idx] = "mrk")
        occursin("event", lowercase(clabels[idx])) && (channel_type[idx] = "mrk")
        occursin("events", lowercase(clabels[idx])) && (channel_type[idx] = "mrk")
        occursin("annotation", lowercase(clabels[idx])) && (channel_type[idx] = "mrk")
        occursin("annotations", lowercase(clabels[idx])) && (channel_type[idx] = "mrk")
        occursin("status", lowercase(clabels[idx])) && (channel_type[idx] = "mrk")
        # eeg channels should have priority, e.g. C3A1 (C3 referenced to A1 should be of eeg type, not ref)
        in(lowercase(clabels[idx]), channel_names) && (channel_type[idx] = "eeg")
        for idx2 in eachindex(channel_names)
            occursin(channel_names[idx2], lowercase(clabels[idx])) && (channel_type[idx] = "eeg")
        end
        lowercase(clabels[idx])[1] == 'c' && (channel_type[idx] = "eeg")
        lowercase(clabels[idx])[1] == 'f' && (channel_type[idx] = "eeg")
        lowercase(clabels[idx])[1] == 'n' && (channel_type[idx] = "eeg")
        lowercase(clabels[idx])[1] == 'o' && (channel_type[idx] = "eeg")
        lowercase(clabels[idx])[1] == 'p' && (channel_type[idx] = "eeg")
        lowercase(clabels[idx])[1] == 't' && (channel_type[idx] = "eeg")
        lowercase(clabels[idx])[1] == 'i' && (channel_type[idx] = "eeg")
        (length(clabels[idx]) > 1 && lowercase(clabels[idx])[1:2] == "af") && (channel_type[idx] = "eeg")
    end
    return channel_type
end

function _m2df(markers::Vector{String})
    # convert EDF/BDF markers to DataFrame
    markers = replace.(markers, "\x14\x14\0" => "|")
    markers = replace.(markers, "\x14\x14" => "|")
    markers = replace.(markers, "\x14" => "|")
    markers = replace.(markers, "\0" => "")
    a_start = Vector{Float64}()
    a_event = Vector{String}()
    # what about markers containing event duration?
    for idx in eachindex(markers)
        s = split(markers[idx], "|")
        if length(s) > 2
            push!(a_start, parse(Float64, strip(s[2])))
            push!(a_event, strip(s[3]))
        end
    end
    return DataFrame(:id => repeat([""], length(a_event)), :start => a_start, :length => zeros(Int64, length(a_event)), :description => a_event, :channel => zeros(Int64, length(a_event)))
end

function _sort_channels(ch_t::Vector{String})
    replace!(ch_t, "eeg" => "1")
    replace!(ch_t, "meg" => "2")
    replace!(ch_t, "ref" => "3")
    replace!(ch_t, "eog" => "4")
    replace!(ch_t, "ecg" => "5")
    replace!(ch_t, "emg" => "6")
    replace!(ch_t, "other" => "7")
    replace!(ch_t, "markers" => "7")
    return sortperm(ch_t)
end
