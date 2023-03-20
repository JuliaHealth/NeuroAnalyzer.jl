function _set_units(obj::NeuroAnalyzer.NEURO, ch::Int64)
    units = ""
    obj.header.recording[:channel_type][ch] == "eeg" && (units = "μV")
    obj.header.recording[:channel_type][ch] == "csd" && (units = "μV/m²")
    obj.header.recording[:channel_type][ch] == "meg" && (units = "fT")
    obj.header.recording[:channel_type][ch] == "grad" && (units = "fT/cm")
    obj.header.recording[:channel_type][ch] == "ecg" && (units = "mV")
    return units
end

function _channel2channel_name(channel::Union{Int64, Vector{Int64}, <:AbstractRange})
    if typeof(channel) == Int64
        return channel
    else
        if collect(channel[1]:channel[end]) == channel
            channel_name = string(channel[1]) * ":" * string(channel[end])
        else
            channel_name = ""
            for idx in 1:(length(channel) - 1)
                channel_name *= string(channel[idx])
                channel_name *= ", "
            end
            channel_name *= string(channel[end])
        end
    end
    return channel_name
end

function _map_channels(channel::Union{Int64, Vector{Int64}, <:AbstractRange}, channels=Vector{Int64})
    channel_orig = channel
    if typeof(channel) == Int64
        channel = vsearch(channel, channels)
    else
        for idx in eachindex(channel)
            channel[idx] = vsearch(channel[idx], channels)
        end
    end
    return channel, channel_orig
end

function _get_ch_idx(clabels::Vector{String}, ch::Union{String, Int64})
    if typeof(ch) == String
        ch_found = nothing
        for idx in eachindex(clabels)
            if ch == clabels[idx]
                ch_found = idx
            end
        end
        if ch_found === nothing
            throw(ArgumentError("ch name does not match signal labels."))
        end
    else
        ch < 1 || ch > length(clabels) && throw(ArgumentError("channel index does not match signal channels."))
        ch_found = ch
    end

    return ch_found
end

function _set_channel_types(clabels::Vector{String})
    channel_names = ["af3", "af4", "af7", "af8", "afz", "c1", "c2", "c3", "c4", "c5", "c6", "cp1", "cp2", "cp3", "cp4", "cp5", "cp6", "cpz", "cz", "f1", "f10", "f2", "f3", "f4", "f5", "f6", "f7", "f8", "f9", "fc1", "fc2", "fc3", "fc4", "fc5", "fc6", "fcz", "fp1", "fp2", "fpz", "ft10", "ft7", "ft8", "ft9", "fz", "nz", "o1", "o2", "oz", "p1", "p10", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "po3", "po4", "po7", "po8", "poz", "pz", "t10", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "tp10", "tp7", "tp8", "tp9"]
    ref_channels = ["a1", "a2", "m1", "m2", "pg1", "pg2"]
    eog_channel = ["e", "e1", "e2"]
    channel_type = repeat(["???"], length(clabels))
    for idx in eachindex(clabels)
        occursin("meg", lowercase(clabels[idx])) && (channel_type[idx] = "meg")
        occursin("mag", lowercase(clabels[idx])) && (channel_type[idx] = "mag")
        occursin("grad", lowercase(clabels[idx])) && (channel_type[idx] = "grad")
        occursin("ecg", lowercase(clabels[idx])) && (channel_type[idx] = "ecg")
        occursin("ekg", lowercase(clabels[idx])) && (channel_type[idx] = "ecg")
        occursin("eog", lowercase(clabels[idx])) && (channel_type[idx] = "eog")
        occursin("rr", lowercase(clabels[idx])) && (channel_type[idx] = "other")
        occursin("mic", lowercase(clabels[idx])) && (channel_type[idx] = "other")
        occursin("flw", lowercase(clabels[idx])) && (channel_type[idx] = "other")
        occursin("tho", lowercase(clabels[idx])) && (channel_type[idx] = "other")
        occursin("abd", lowercase(clabels[idx])) && (channel_type[idx] = "other")
        occursin("sao2", lowercase(clabels[idx])) && (channel_type[idx] = "other")
        occursin("sa02", lowercase(clabels[idx])) && (channel_type[idx] = "other")
        occursin("plr", lowercase(clabels[idx])) && (channel_type[idx] = "other")
        occursin("body", lowercase(clabels[idx])) && (channel_type[idx] = "other")
        occursin("ux", lowercase(clabels[idx])) && (channel_type[idx] = "other")
        for idx2 in eachindex(eog_channel)
            occursin(channel_names[idx2], lowercase(clabels[idx])) && (channel_type[idx] = "eog")
        end
        occursin("emg", lowercase(clabels[idx])) && (channel_type[idx] = "emg")
        in(lowercase(clabels[idx]), ref_channels) && (channel_type[idx] = "ref")
        for idx2 in eachindex(ref_channels)
            occursin(ref_channels[idx2], lowercase(clabels[idx])) && (channel_type[idx] = "ref")
        end
        occursin("stim", lowercase(clabels[idx])) && (channel_type[idx] = "mrk")
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

function _sort_channels(ch_t::Vector{String})
    replace!(ch_t, "eeg" => "1")
    replace!(ch_t, "meg" => "2")
    replace!(ch_t, "mag" => "2")
    replace!(ch_t, "grad" => "3")
    replace!(ch_t, "csd" => "1")
    replace!(ch_t, "ref" => "4")
    replace!(ch_t, "eog" => "5")
    replace!(ch_t, "ecg" => "6")
    replace!(ch_t, "emg" => "7")
    replace!(ch_t, "other" => "8")
    replace!(ch_t, "mrk" => "9")
    return sortperm(ch_t)
end
