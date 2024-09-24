function _v2r(v::Vector{Int64})
    if sort(v) == sort(v)[1]:sort(v)[end]
        return sort(v)[1]:sort(v)[end]
    else
        return v
    end
end
_v2r(v::AbstractRange) = return v
_v2r(v::Int64) = return v

function _ch_rename(ch_type::String)::String
    lowercase(ch_type) == "eeg" && return "EEG"
    lowercase(ch_type) == "ecog" && return "ECoG"
    lowercase(ch_type) == "ecog" && return "ECoG"
    lowercase(ch_type) == "seeg" && return "sEEG"
    lowercase(ch_type) == "grad" && return "MEG gradiometer"
    lowercase(ch_type) == "mag" && return "MEG magnetometer"
    lowercase(ch_type) == "nirs_int" && return "NIRS intensity"
    lowercase(ch_type) == "nirs_od" && return "NIRS optical density"
    lowercase(ch_type) == "nirs_hbo" && return "NIRS HbO concentration"
    lowercase(ch_type) == "nirs_hbr" && return "NIRS HbR concentration"
    lowercase(ch_type) == "nirs_hbt" && return "NIRS HbT concentration"
    lowercase(ch_type) == "ref" && return "reference"
    lowercase(ch_type) == "other" && return "other"
    lowercase(ch_type) == "mrk" && return "marker"
    lowercase(ch_type) == "accel" && return "acceleration"
    lowercase(ch_type) == "magfld" && return "magnetic field"
    lowercase(ch_type) == "orient" && return "orientation"
    lowercase(ch_type) == "angvel" && return "angular velocity"
    return ch_type
end

function _def_ylabel(ch_type::String, u::String)::String
    yl = "Amplitude [$u]"
    lowercase(ch_type) == "nirs_int" && return "Intensity [$u]"
    lowercase(ch_type) == "nirs_od" && return "OD [$u]"
    lowercase(ch_type) == "nirs_hbo" && return "HbO concentration [$u]"
    lowercase(ch_type) == "nirs_hbr" && return "HbR concentration [$u]"
    lowercase(ch_type) == "nirs_hbt" && return "HbT concentration [$u]"
    return yl
end

function _ch_units(ch_type::String)::String
    u = ""
    lowercase(ch_type) == "eeg" && return "μV"
    lowercase(ch_type) == "ieeg" && return "μV"
    lowercase(ch_type) == "ecog" && return "μV"
    lowercase(ch_type) == "seeg" && return "μV"
    lowercase(ch_type) == "csd" && return "μV/m²"
    lowercase(ch_type) == "mag" && return "fT"
    lowercase(ch_type) == "grad" && return "fT/cm"
    lowercase(ch_type) == "emg" && return "μV"
    lowercase(ch_type) == "eog" && return "μV"
    lowercase(ch_type) == "ref" && return "μV"
    lowercase(ch_type) == "ecg" && return "mV"
    lowercase(ch_type) == "nirs_int" && return "V"
    lowercase(ch_type) == "nirs_od" && return ""
    lowercase(ch_type) == "nirs_hbo" && return "μM/mm"
    lowercase(ch_type) == "nirs_hbr" && return "μM/mm"
    lowercase(ch_type) == "nirs_hbt" && return "μM/mm"
    lowercase(ch_type) == "nirs_aux" && return ""
    lowercase(ch_type) == "mrk" && return ""
    lowercase(ch_type) == "accel" && return "m/s²"
    lowercase(ch_type) == "magfld" && return "µT"
    lowercase(ch_type) == "orient" && return "°"
    lowercase(ch_type) == "angvel" && return "rad/s"
    lowercase(ch_type) == "mep" && return "μV"
    lowercase(ch_type) == "eda" && return "μS"
    lowercase(ch_type) == "other" && return ""
    return u
end

_ch_units(obj::NeuroAnalyzer.NEURO, ch::String) = _ch_units(obj.header.recording[:channel_type][_ch_idx(obj, ch)[1]])

function _ch_idx(cl::Union{String, Vector{String}}, l::Union{String, Vector{String}})
    l == "" && return Int64[]
    isa(l, String) && (l = [l])
    isa(cl, String) && (cl = [l])
    any(occursin.("all", l)) && (l = cl)
    _check_channels(cl, l)
    ch = Int64[]
    for ch_idx in eachindex(l)
        @assert l[ch_idx] in cl "$(l[ch_idx]) does not match signal labels."
        push!(ch, findfirst(isequal(l[ch_idx]), cl))
    end
    return unique(ch)
end

function _ch_idx(obj::NeuroAnalyzer.NEURO, l::Union{String, Vector{String}})
    l == "" && return Int64[]
    cl = labels(obj)
    isa(l, String) && (l = [l])
    isa(cl, String) && (cl = [l])
    any(occursin.("all", l)) && (l = cl)
    if any(occursin.("meg", l))
        idx = findfirst(isequal("meg"), l)
        l[idx] = "mag"
        if idx < length(l)
            l = [l[1:idx]; "grad"; l[idx+1:end]]
        else
            l = append!(l, ["grad"])
        end
    end
    if any(occursin.("nirs", l))
        idx = findfirst(isequal("nirs"), l)
        l[idx] = "nirs_int"
        if idx < length(l)
            l = [l[1:idx]; "nirs_od"; "nirs_dmean"; "nirs_dvar"; "nirs_dskew"; "nirs_mua"; "nirs_musp"; "nirs_hbo"; "nirs_hbr"; "nirs_hbt"; "nirs_h2o"; "nirs_lipid"; "nirs_bfi"; "nirs_hrf_dod"; "nirs_hrf_dmean"; "nirs_hrf_dvar"; "nirs_hrf_dskew"; "nirs_hrf_hbo"; "nirs_hrf_hbr"; "nirs_hrf_hbt"; "nirs_hrf_bfi"; "nirs_aux"; l[idx+1:end]]
        else
            append!(l, ["nirs_od", "nirs_dmean", "nirs_dvar", "nirs_dskew", "nirs_mua", "nirs_musp", "nirs_hbo", "nirs_hbr", "nirs_hbt", "nirs_h2o", "nirs_lipid", "nirs_bfi", "nirs_hrf_dod", "nirs_hrf_dmean", "nirs_hrf_dvar", "nirs_hrf_dskew", "nirs_hrf_hbo", "nirs_hrf_hbr", "nirs_hrf_hbt", "nirs_hrf_bfi", "nirs_aux"])
        end
    end
    if any(occursin.("sensors", l))
        idx = findfirst(isequal("sensors"), l)
        l[idx] = "accel"
        if idx < length(l)
            l = [l[1:idx]; "magfld"; "orient"; "angvel"; l[idx+1:end]]
        else
            append!(l, ["magfld", "orient", "angvel"])
        end
    end
    if any(occursin.("bad", l))
        idx = findfirst(isequal("bad"), l)
        bads = vec(sum(obj.header.recording[:bad_channel], dims=2))
        bad_ch = labels(obj)[bads.>0]
        deleteat!(l, idx)
        if length(bad_ch) > 0
            if idx < length(l)
                l = [l[1:idx]; bad_ch l[idx+1:end]]
            else
                append!(l, bad_ch)
            end
        end
    end
    length(l) == 0 && return Int64[]
    l_tmp = String[]
    for idx1 in eachindex(l)
        if l[idx1] in NeuroAnalyzer.channel_types
            for idx2 in NeuroAnalyzer.channel_types
                l[idx1] == idx2 && append!(l_tmp, get_channel(obj, type=idx2))
            end
        else
            push!(l_tmp, l[idx1])
        end
    end
    l = l_tmp
    NeuroAnalyzer._check_channels(cl, l)
    ch = Int64[]
    for ch_idx in eachindex(l)
        @assert l[ch_idx] in cl "$(l[ch_idx]) does not match signal labels."
        push!(ch, findfirst(isequal(l[ch_idx]), cl))
    end
    return unique(ch)
end

function _set_channel_types(clabels::Vector{String}, default::String="other")
    channel_names = ["af3", "af4", "af7", "af8", "afz", "c1", "c2", "c3", "c4", "c5", "c6", "cp1", "cp2", "cp3", "cp4", "cp5", "cp6", "cpz", "cz", "f1", "f10", "f2", "f3", "f4", "f5", "f6", "f7", "f8", "f9", "fc1", "fc2", "fc3", "fc4", "fc5", "fc6", "fcz", "fp1", "fp2", "fpz", "ft10", "ft7", "ft8", "ft9", "fz", "nz", "o1", "o2", "oz", "p1", "p10", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "po3", "po4", "po7", "po8", "poz", "pz", "t10", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "tp10", "tp7", "tp8", "tp9"]
    ref_channels = ["a1", "a2", "m1", "m2", "pg1", "pg2"]
    eog_channel = ["e", "e1", "e2"]
    channel_type = repeat([default], length(clabels))
    for idx in eachindex(clabels)
        occursin("ecg", lowercase(clabels[idx])) && (channel_type[idx] = "ecg")
        occursin("ekg", lowercase(clabels[idx])) && (channel_type[idx] = "ecg")

        occursin("eog", lowercase(clabels[idx])) && (channel_type[idx] = "eog")
        for idx2 in eachindex(eog_channel)
            occursin(channel_names[idx2], lowercase(clabels[idx])) && (channel_type[idx] = "eog")
        end

        occursin("emg", lowercase(clabels[idx])) && (channel_type[idx] = "emg")

        in(lowercase(clabels[idx]), ref_channels) && (channel_type[idx] = "ref")

        occursin("meg", lowercase(clabels[idx])) && (channel_type[idx] = "meg")
        occursin("mag", lowercase(clabels[idx])) && (channel_type[idx] = "mag")
        occursin("grad", lowercase(clabels[idx])) && (channel_type[idx] = "grad")

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
        (length(clabels[idx]) > 2 && lowercase(clabels[idx])[1] == "a") && (channel_type[idx] = "eeg")

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
        occursin("acc", lowercase(clabels[idx])) && (channel_type[idx] = "accel")

        occursin("stim", lowercase(clabels[idx])) && (channel_type[idx] = "mrk")
        occursin("mark", lowercase(clabels[idx])) && (channel_type[idx] = "mrk")
        occursin("marker", lowercase(clabels[idx])) && (channel_type[idx] = "mrk")
        occursin("markers", lowercase(clabels[idx])) && (channel_type[idx] = "mrk")
        occursin("event", lowercase(clabels[idx])) && (channel_type[idx] = "mrk")
        occursin("event", lowercase(clabels[idx])) && (channel_type[idx] = "mrk")
        occursin("trigger", lowercase(clabels[idx])) && (channel_type[idx] = "mrk")
        lowercase(clabels[idx]) == "e" && (channel_type[idx] = "mrk")
        lowercase(clabels[idx]) == "stim" && (channel_type[idx] = "mrk")
        occursin("annotation", lowercase(clabels[idx])) && (channel_type[idx] = "mrk")
        occursin("annotations", lowercase(clabels[idx])) && (channel_type[idx] = "mrk")
        occursin("status", lowercase(clabels[idx])) && (channel_type[idx] = "mrk")
    end
    return channel_type
end

function _sort_channels(ch_t::Vector{String})
    ch_order = deepcopy(ch_t)
    replace!(ch_order, "eeg" => "3")
    replace!(ch_order, "seeg" => "3")
    replace!(ch_order, "ecog" => "3")
    replace!(ch_order, "ieeg" => "3")
    replace!(ch_order, "meg" => "1")
    replace!(ch_order, "grad" => "1")
    replace!(ch_order, "mag" => "2")
    replace!(ch_order, "csd" => "4")
    replace!(ch_order, "ref" => "4")
    replace!(ch_order, "eog" => "5")
    replace!(ch_order, "ecg" => "6")
    replace!(ch_order, "emg" => "7")
    replace!(ch_order, "other" => "8")
    replace!(ch_order, "nirs_int" => "1")
    replace!(ch_order, "nirs_od" => "1")
    replace!(ch_order, "nirs_dmean" => "3")
    replace!(ch_order, "nirs_dvar" => "3")
    replace!(ch_order, "nirs_dskew" => "3")
    replace!(ch_order, "nirs_mua" => "3")
    replace!(ch_order, "nirs_musp" => "3")
    replace!(ch_order, "nirs_hbo" => "2")
    replace!(ch_order, "nirs_hbr" => "2")
    replace!(ch_order, "nirs_hbt" => "2")
    replace!(ch_order, "nirs_h2o" => "2")
    replace!(ch_order, "nirs_lipid" => "3")
    replace!(ch_order, "nirs_bfi" => "3")
    replace!(ch_order, "nirs_hrf_dod" => "3")
    replace!(ch_order, "nirs_hrf_dmean" => "3")
    replace!(ch_order, "nirs_hrf_dvar" => "3")
    replace!(ch_order, "nirs_hrf_dskew" => "3")
    replace!(ch_order, "nirs_hrf_hbo" => "3")
    replace!(ch_order, "nirs_hrf_hbr" => "3")
    replace!(ch_order, "nirs_hrf_hbt" => "3")
    replace!(ch_order, "nirs_hrf_bfi" => "3")
    replace!(ch_order, "nirs_aux" => "4")
    replace!(ch_order, "accel" => "1")
    replace!(ch_order, "magfld" => "2")
    replace!(ch_order, "orient" => "3")
    replace!(ch_order, "angvel" => "4")
    replace!(ch_order, "mrk" => "9")
    replace!(ch_order, "stim" => "9")
    return sortperm(ch_order)
end
