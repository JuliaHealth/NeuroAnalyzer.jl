function _v2r(v::Vector{Int64})::Union{AbstractRange, Vector{Int64}}
    sv = sort(v)
    # return a contiguous range only if the sorted vector equals its own range
    if sv == sv[1]:sv[end]
        return sv[1]:sv[end]
    else
        return sv
    end
end
_v2r(v::AbstractRange)::AbstractRange = v
_v2r(v::Int64)::Int64 = v

function _ch_rename(ch_type::String)::String
    ct = lowercase(ch_type)
    ct == "eeg"       && return "EEG"
    ct == "ecog"      && return "ECoG"
    ct == "seeg"      && return "sEEG"
    ct == "grad"      && return "MEG gradiometer"
    ct == "mag"       && return "MEG magnetometer"
    ct == "nirs_int"  && return "NIRS intensity"
    ct == "nirs_od"   && return "NIRS optical density"
    ct == "nirs_hbo"  && return "NIRS HbO concentration"
    ct == "nirs_hbr"  && return "NIRS HbR concentration"
    ct == "nirs_hbt"  && return "NIRS HbT concentration"
    ct == "ref"       && return "reference"
    ct == "other"     && return "other"
    ct == "mrk"       && return "marker"
    ct == "accel"     && return "acceleration"
    ct == "magfld"    && return "magnetic field"
    ct == "orient"    && return "orientation"
    ct == "angvel"    && return "angular velocity"
    return ch_type
end

function _def_ylabel(ch_type::String, u::String)::String
    ct = lowercase(ch_type)
    ct == "nirs_int" && return "Intensity [$u]"
    ct == "nirs_od"  && return "OD [$u]"
    ct == "nirs_hbo" && return "HbO concentration [$u]"
    ct == "nirs_hbr" && return "HbR concentration [$u]"
    ct == "nirs_hbt" && return "HbT concentration [$u]"
    return "Amplitude [$u]"
end

function _ch_units(ch_type::String)::String
    ct = lowercase(ch_type)
    ct == "eeg"       && return "μV"
    ct == "ieeg"      && return "μV"
    ct == "ecog"      && return "μV"
    ct == "seeg"      && return "μV"
    ct == "csd"       && return "μV/m²"
    ct == "mag"       && return "fT"
    ct == "grad"      && return "fT/cm"
    ct == "emg"       && return "μV"
    ct == "eog"       && return "μV"
    ct == "ref"       && return "μV"
    ct == "ecg"       && return "mV"
    ct == "nirs_int"  && return "V"
    ct == "nirs_od"   && return ""
    ct == "nirs_hbo"  && return "μM/mm"
    ct == "nirs_hbr"  && return "μM/mm"
    ct == "nirs_hbt"  && return "μM/mm"
    ct == "nirs_aux"  && return ""
    ct == "mrk"       && return ""
    ct == "accel"     && return "m/s²"
    ct == "magfld"    && return "µT"
    ct == "orient"    && return "°"
    ct == "angvel"    && return "rad/s"
    ct == "mep"       && return "μV"
    ct == "eda"       && return "μS"
    ct == "other"     && return ""
    return ""
end

_ch_units(obj::NeuroAnalyzer.NEURO, ch::String)::String =
    _ch_units(obj.header.recording[:channel_type][_ch_idx(obj, ch)[1]])

function _ch_idx(
    cl::Union{String, Vector{String}},
    l::Union{String, Vector{String}, Regex}
)::Vector{Int64}
    if isa(l, Regex)
        matches = Base.filter(!isnothing, match.(l, cl))
        l = [m.match for m in matches]
        length(l) == 0 && return Int64[]
    end
    l == "" && return Int64[]
    isa(l, String) && (l = [l])
    isa(cl, String) && (cl = [cl])
    any(isequal("all"), l) && (l = cl)
    _check_channels(cl, l)
    ch = Int64[]
    for label in l
        !(label in cl) && throw(ArgumentError("$label does not match signal labels."))
        push!(ch, findfirst(isequal(label), cl))
    end
    return unique(ch)
end

# ---------------------------------------------------------------------------
# NIRS channel type list — defined once to avoid repetition
# ---------------------------------------------------------------------------
const _NIRS_TYPES = [
    "nirs_od", "nirs_dmean", "nirs_dvar", "nirs_dskew", "nirs_mua", "nirs_musp",
    "nirs_hbo", "nirs_hbr", "nirs_hbt", "nirs_h2o", "nirs_lipid", "nirs_bfi",
    "nirs_hrf_dod", "nirs_hrf_dmean", "nirs_hrf_dvar", "nirs_hrf_dskew",
    "nirs_hrf_hbo", "nirs_hrf_hbr", "nirs_hrf_hbt", "nirs_hrf_bfi", "nirs_aux",
]

function _ch_idx(
    obj::NeuroAnalyzer.NEURO,
    l::Union{String, Vector{String}, Regex}
)::Vector{Int64}

    cl = labels(obj)

    if isa(l, Regex)
        matches = Base.filter(!isnothing, match.(l, cl))
        l = [m.match for m in matches]
        length(l) == 0 && return Int64[]
    end

    l == "" && return Int64[]
    isa(l, String) && (l = [l])
    isa(cl, String) && (cl = [cl])

    any(isequal("all"), l) && (l = cl)

    # expand "meg" → ["mag", "grad"]
    if any(isequal("meg"), l)
        idx = findfirst(isequal("meg"), l)
        l = [l[1:(idx - 1)]; "mag"; "grad"; l[(idx + 1):end]]
    end

    # expand "nirs" → all NIRS sub-types
    if any(isequal("nirs"), l)
        idx = findfirst(isequal("nirs"), l)
        l = [l[1:(idx - 1)]; "nirs_int"; _NIRS_TYPES; l[(idx + 1):end]]
    end

    # expand "sensors" → ["accel", "magfld", "orient", "angvel"]
    if any(isequal("sensors"), l)
        idx = findfirst(isequal("sensors"), l)
        l = [l[1:(idx - 1)]; "accel"; "magfld"; "orient"; "angvel"; l[(idx + 1):end]]
    end

    # expand "bad" → labels of channels flagged as bad
    if any(isequal("bad"), l)
        idx  = findfirst(isequal("bad"), l)
        bads = labels(obj)[obj.header.recording[:bad_channel] .== true]
        l    = [l[1:(idx - 1)]; bads; l[(idx + 1):end]]
    end

    length(l) == 0 && return Int64[]

    # replace channel-type tokens with the matching channel labels
    l_tmp = String[]
    for token in l
        if token in NeuroAnalyzer.channel_types
            append!(l_tmp, get_channel(obj; type=token))
        else
            push!(l_tmp, token)
        end
    end
    l = l_tmp

    _check_channels(cl, l)

    ch = Int64[]
    for label in l
        !(label in cl) && throw(ArgumentError("$label does not match signal labels."))
        push!(ch, findfirst(isequal(label), cl))
    end
    return unique(ch)
end

function _set_channel_types(
    clabels::Vector{String},
    default::String = "other"
)::Vector{String}
    channel_names = [
        "af3","af4","af7","af8","afz",
        "c1","c2","c3","c4","c5","c6",
        "cp1","cp2","cp3","cp4","cp5","cp6","cpz","cz",
        "f1","f10","f2","f3","f4","f5","f6","f7","f8","f9",
        "fc1","fc2","fc3","fc4","fc5","fc6","fcz",
        "fp1","fp2","fpz",
        "ft10","ft7","ft8","ft9","fz","nz",
        "o1","o2","oz",
        "p1","p10","p2","p3","p4","p5","p6","p7","p8","p9",
        "po3","po4","po7","po8","poz","pz",
        "t10","t3","t4","t5","t6","t7","t8","t9",
        "tp10","tp7","tp8","tp9",
    ]
    ref_channels = ["a1", "a2", "m1", "m2", "pg1", "pg2"]
    eog_channels = ["e", "e1", "e2"]

    channel_type = fill(default, length(clabels))

    for idx in eachindex(clabels)
        lbl = lowercase(clabels[idx])

        occursin("ecg",  lbl) && (channel_type[idx] = "ecg")
        occursin("ekg",  lbl) && (channel_type[idx] = "ecg")
        occursin("eog",  lbl) && (channel_type[idx] = "eog")

        for eog in eog_channels
            lbl == eog && (channel_type[idx] = "eog")
        end

        occursin("emg",  lbl) && (channel_type[idx] = "emg")
        lbl in ref_channels   && (channel_type[idx] = "ref")
        occursin("mag",  lbl) && (channel_type[idx] = "mag")
        occursin("grad", lbl) && (channel_type[idx] = "grad")
        occursin("meg",  lbl) && (channel_type[idx] = "meg")

        # EEG channels take priority (e.g. "C3A1" → eeg, not ref)
        lbl in channel_names  && (channel_type[idx] = "eeg")
        for ch_name in channel_names
            occursin(ch_name, lbl) && (channel_type[idx] = "eeg")
        end

        length(lbl) >= 1 && lbl[1] in ('c','f','n','o','p','t','i') && (channel_type[idx] = "eeg")
        length(lbl) >= 2 && lbl[1:2] == "af" && (channel_type[idx] = "eeg")
        length(lbl) >= 3 && lbl[1]   == 'a'  && (channel_type[idx] = "eeg")

        # non-neural / auxiliary channels
        for pattern in ("rr","mic","flw","tho","abd","sao2","sa02","plr","body","ux","ias","sys","aux")
            occursin(pattern, lbl) && (channel_type[idx] = "other")
        end

        occursin("acc", lbl) && (channel_type[idx] = "accel")

        # marker / event channels
        for pattern in ("sti","stim","mark","marker","markers","event","trigger","annotation","annotations","status")
            occursin(pattern, lbl) && (channel_type[idx] = "mrk")
        end
        lbl == "e" && (channel_type[idx] = "mrk")
        lbl == "stim" && (channel_type[idx] = "mrk") # already covered by occursin above; kept for explicitness

        occursin("nm", lbl) && (channel_type[idx] = "nirs_od")
    end

    return channel_type
end

function _sort_channels(ch_t::Vector{String})::Vector{Int64}
    # map each channel type to a sort-priority string (lower string = sorted first)
    priority = Dict(
        "meg"           => "1", "grad"          => "1", "mag"           => "2",
        "eeg"           => "3", "seeg"          => "3", "ecog"          => "3",
        "ieeg"          => "3", "csd"           => "4", "ref"           => "4",
        "eog"           => "5", "ecg"           => "6", "emg"           => "7",
        "other"         => "8", "mrk"           => "9", "stim"          => "9",
        "nirs_int"      => "1", "nirs_od"       => "1", "nirs_hbo"      => "2",
        "nirs_hbr"      => "2", "nirs_hbt"      => "2", "nirs_h2o"      => "2",
        "nirs_dmean"    => "3", "nirs_dvar"     => "3", "nirs_dskew"    => "3",
        "nirs_mua"      => "3", "nirs_musp"     => "3", "nirs_lipid"    => "3",
        "nirs_bfi"      => "3", "nirs_hrf_dod"  => "3", "nirs_hrf_dmean"=> "3",
        "nirs_hrf_dvar" => "3", "nirs_hrf_dskew"=> "3", "nirs_hrf_hbo"  => "3",
        "nirs_hrf_hbr"  => "3", "nirs_hrf_hbt"  => "3", "nirs_hrf_bfi"  => "3",
        "nirs_aux"      => "4", "accel"         => "1", "magfld"        => "2",
        "orient"        => "3", "angvel"        => "4",
    )
    ch_order = [get(priority, t, "8") for t in ch_t]   # default to "8" for unknowns
    return sortperm(ch_order)
end
