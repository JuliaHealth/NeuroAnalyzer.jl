"""Convert a vector to range."""
function _v2r(v::Vector{Int64})::Union{AbstractUnitRange{Int64}, Vector{Int64}}
    isempty(v) && return v
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

"""Get channel names from their types."""
function _ch_rename(ch_type::String)::String
    ct = lowercase(ch_type)
    ct == "eeg" && return "EEG"
    ct == "ecog" && return "ECoG"
    ct == "seeg" && return "sEEG"
    ct == "grad" && return "MEG gradiometer"
    ct == "mag" && return "MEG magnetometer"
    ct == "nirs_int" && return "NIRS intensity"
    ct == "nirs_od" && return "NIRS optical density"
    ct == "nirs_hbo" && return "NIRS HbO concentration"
    ct == "nirs_hbr" && return "NIRS HbR concentration"
    ct == "nirs_hbt" && return "NIRS HbT concentration"
    ct == "ref" && return "reference"
    ct == "other" && return "other"
    ct == "mrk" && return "marker"
    ct == "accel" && return "acceleration"
    ct == "magfld" && return "magnetic field"
    ct == "orient" && return "orientation"
    ct == "angvel" && return "angular velocity"
    return ch_type
end

"""Set default y axis label based on channel type and units."""
function _def_ylabel(ch_type::String, u::String)::String
    ct = lowercase(ch_type)
    ct == "nirs_int" && return "Intensity [$u]"
    ct == "nirs_od" && return "OD [$u]"
    ct == "nirs_hbo" && return "HbO concentration [$u]"
    ct == "nirs_hbr" && return "HbR concentration [$u]"
    ct == "nirs_hbt" && return "HbT concentration [$u]"
    return "Amplitude [$u]"
end

"""Get channel units based on channel type."""
function _ch_units(ch_type::String)::String
    ct = lowercase(ch_type)
    ct == "eeg" && return "μV"
    ct == "ieeg" && return "μV"
    ct == "ecog" && return "μV"
    ct == "seeg" && return "μV"
    ct == "csd" && return "μV/m²"
    ct == "mag" && return "fT"
    ct == "grad" && return "fT/cm"
    ct == "emg" && return "μV"
    ct == "eog" && return "μV"
    ct == "ref" && return "μV"
    ct == "ecg" && return "mV"
    ct == "nirs_int" && return "V"
    ct == "nirs_od" && return ""
    ct == "nirs_hbo" && return "μM/mm"
    ct == "nirs_hbr" && return "μM/mm"
    ct == "nirs_hbt" && return "μM/mm"
    ct == "nirs_aux" && return ""
    ct == "mrk" && return ""
    ct == "accel" && return "m/s²"
    ct == "magfld" && return "μT"
    ct == "orient" && return "°"
    ct == "angvel" && return "rad/s"
    ct == "mep" && return "μV"
    ct == "eda" && return "μS"
    ct == "other" && return ""
    return ""
end

"""Get channel units from the object header."""
_ch_units(obj::NeuroAnalyzer.NEURO, ch::String)::String =
    _ch_units(obj.header.recording[:channel_type][_ch_idx(obj, ch)[1]])

# ------------------------------------------------------------------------------
# data types sharing the EEG-like / MEG-like recording header schema
# (per-channel :transducers, :gain / :coil_type, …; :lowpass, :highpass for all data types)
# ------------------------------------------------------------------------------
const _EEG_LIKE_DATA = ("eeg", "seeg", "ecog", "ieeg", "csd", "erp")
const _MEG_LIKE_DATA = ("meg", "erf")

# ------------------------------------------------------------------------------
# NIRS channel type list - defined once to avoid repetition
# ------------------------------------------------------------------------------
const _NIRS_TYPES = [
    "nirs_od", "nirs_dmean", "nirs_dvar", "nirs_dskew", "nirs_mua", "nirs_musp",
    "nirs_hbo", "nirs_hbr", "nirs_hbt", "nirs_h2o", "nirs_lipid", "nirs_bfi",
    "nirs_hrf_dod", "nirs_hrf_dmean", "nirs_hrf_dvar", "nirs_hrf_dskew",
    "nirs_hrf_hbo", "nirs_hrf_hbr", "nirs_hrf_hbt", "nirs_hrf_bfi", "nirs_aux",
]

# channel type groups expanded to their member types
const _CH_TYPE_GROUPS = Dict(
    "meg" => ["mag", "grad"],
    "nirs" => ["nirs_int"; _NIRS_TYPES],
    "sensors" => ["accel", "magfld", "orient", "angvel"],
)

# ------------------------------------------------------------------------------
# channel types present in an object
#
# An object may contain channels of several modalities (e.g. EEG recorded with NIRS):
# functions for a modality check that the object contains channels of its type(s), not the
# data type of the object (`datatype()` describes the recording; "erp", "erf", "mep", ...
# also describe its state)
# ------------------------------------------------------------------------------

# EEG-like channel types (per-channel :gain, :transducers; :reference)
const _EEG_LIKE_TYPES = ["eeg", "ecog", "seeg", "ieeg", "csd"]
# channel types of brain signals (not auxiliary channels)
const _SIGNAL_CH_TYPES =
    [_EEG_LIKE_TYPES; "mag"; "grad"; "nirs_int"; Base.filter(!=("nirs_aux"), _NIRS_TYPES)]

# channel types and type groups ("meg", "nirs", "sensors") expanded to channel types
_ch_types_expand(type::Union{String, AbstractVector{String}})::Vector{String} =
    unique(reduce(vcat, [get(_CH_TYPE_GROUPS, t, [t]) for t in (type isa String ? [type] : type)]))

# true if the object contains channels of the type(s)
_has_chtype(obj::NeuroAnalyzer.NEURO, type::Union{String, AbstractVector{String}})::Bool =
    any(in(_ch_types_expand(type)), obj.header.recording[:channel_type])

# error unless the object contains channels of the type(s)
function _check_chtype(
    obj::NeuroAnalyzer.NEURO,
    type::Union{String, AbstractVector{String}},
)::Nothing
    _has_chtype(obj, type) && return nothing
    t = type isa String ? [type] : type
    ct = unique(obj.header.recording[:channel_type])
    throw(
        ArgumentError(
            "This function requires $(join(uppercase.(t), "/")) channels; OBJ has none (channel types: $(join(ct, ", "))).",
        ),
    )
end

# ------------------------------------------------------------------------------
# per-channel fields of the recording header
#
# Fields are handled by their presence (not by the data type of the object), so that objects
# with channels of several modalities keep all of them consistent: a field present in the
# header has one value per channel; channels of other types get the default value.
# ------------------------------------------------------------------------------

# per-channel vector fields (besides :label, :channel_type, :unit, :bad_channel, :lowpass,
# :highpass, :notch) and their default values
const _REC_CH_FIELDS = Dict{Symbol, Any}(
    :transducers => "",             # EEG-like
    :gain => 1.0,                   # EEG-like
    :coil_type => "",               # MEG (sensor coil), MEP (stimulation coil)
    :ssp_channels => false,         # MEG: channels of the SSP projectors
    :coil_id => 0,                  # MEG: FIFF coil type code
    :stimulation_intensity => 0,    # MEP
)

# per-channel fields needed by channels of the type
function _rec_fields_for_type(t::String)::Vector{Symbol}
    t in _EEG_LIKE_TYPES && return [:transducers, :gain]
    t in ("mag", "grad") && return [:coil_type, :ssp_channels]
    t == "mep" && return [:stimulation_intensity, :coil_type]
    return Symbol[]
end

# true if the per-channel field `f` is present in the header (one value per channel; empty
# optional fields such as :coil_id or :coil_loc are absent)
function _rec_has_ch_field(rec::Dict, f::Symbol, ch_n::Int64)::Bool
    haskey(rec, f) || return false
    v = rec[f]
    v isa AbstractMatrix && return size(v, 1) == ch_n && ch_n > 0
    return v isa AbstractVector && length(v) == ch_n && ch_n > 0
end

# add the per-channel field `f` (default values for the `ch_n` channels) if it is missing
function _rec_ensure_ch_field!(rec::Dict, f::Symbol, ch_n::Int64)::Nothing
    (haskey(rec, f) && rec[f] isa AbstractVector && length(rec[f]) == ch_n) && return nothing
    d = _REC_CH_FIELDS[f]
    rec[f] = fill(d, ch_n)
    return nothing
end

# recording-level fields needed by channels of the type (added with these values if missing)
function _rec_scalar_defaults(t::String)::Dict{Symbol, Any}
    if t in _EEG_LIKE_TYPES || t in ("mag", "grad")
        d = Dict{Symbol, Any}(:reference => "", :line_frequency => 50)
        if t in ("mag", "grad")
            merge!(d, Dict{Symbol, Any}(
                :magnetometers => Int64[], :gradiometers => Int64[],
                :ssp_labels => String[], :ssp_data => Matrix{Float64}(undef, 0, 0),
                :maxshield => false, :sss => false, :tsss => false,
            ))
        end
        return merge!(d, _locs_coord_fields())
    end
    return Dict{Symbol, Any}()
end

# signal channels (labels) of the type of channel `ch` (e.g. the EEG channels for an EEG
# channel of an object with EEG and NIRS channels); empty if `ch` is not a signal channel
function _signal_channels_like(obj::NeuroAnalyzer.NEURO, ch::String)::Vector{String}
    i = findfirst(==(ch), labels(obj))
    isnothing(i) && return String[]
    t = obj.header.recording[:channel_type][i]
    t in _SIGNAL_CH_TYPES || return String[]
    return labels(obj)[obj.header.recording[:channel_type] .== t]
end

# true for empty header values ("", empty vectors, matrices and dictionaries)
_rec_empty(v)::Bool = (v isa AbstractString || v isa AbstractArray || v isa AbstractDict) && isempty(v)

# number of NIRS optical channels described by optode pairs (the first channels)
_nirs_n_opt(rec::Dict)::Int64 =
    haskey(rec, :wavelength_index) ? length(rec[:wavelength_index]) : 0

"""Labels from `cl` matched by a regex anywhere in the label (use `^...\$` for exact matches)."""
_ch_regex(cl::Vector{String}, l::Regex)::Vector{String} = cl[occursin.(l, cl)]

"""Get integer indices from channel names."""
function _ch_idx(
    cl::Union{String, Vector{String}},
    l::Union{String, Vector{String}, Regex},
)::Vector{Int64}
    isa(cl, String) && (cl = [cl])
    if isa(l, Regex)
        l = _ch_regex(cl, l)
    elseif isa(l, String)
        l == "" && return Int64[]
        l = [l]
    end
    isempty(l) && return Int64[]
    any(isequal("all"), l) && return collect(1:length(cl))
    _check_channels(cl, l)
    ch = Int64[]
    for label in l
        label in cl || throw(ArgumentError("$label does not match signal labels."))
        push!(ch, findfirst(isequal(label), cl))
    end
    return unique(ch)
end

"""Get integer indices from channel names."""
function _ch_idx(
    obj::NeuroAnalyzer.NEURO,
    l::Union{String, Vector{String}, Regex},
)::Vector{Int64}
    cl = labels(obj)

    if isa(l, Regex)
        l = _ch_regex(cl, l)
    elseif isa(l, String)
        l == "" && return Int64[]
        l = [l]
    end
    isempty(l) && return Int64[]

    any(isequal("all"), l) && return collect(1:length(cl))

    # expand every keyword token (group, "bad", channel type) in place, preserving order
    l_tmp = String[]
    for token in l
        if token == "bad"
            append!(l_tmp, cl[obj.header.recording[:bad_channel] .== true])
        elseif haskey(_CH_TYPE_GROUPS, token)
            for type in _CH_TYPE_GROUPS[token]
                append!(l_tmp, get_channel(obj; type = type))
            end
        elseif token in NeuroAnalyzer.channel_types
            append!(l_tmp, get_channel(obj; type = token))
        else
            push!(l_tmp, token)
        end
    end
    l = l_tmp

    isempty(l) && return Int64[]
    _check_channels(cl, l)

    ch = Int64[]
    for label in l
        label in cl || throw(ArgumentError("$label does not match signal labels."))
        push!(ch, findfirst(isequal(label), cl))
    end
    return unique(ch)
end

"""Set channel types based on their names."""
function _set_channel_types(
    clabels::Vector{String},
    default::String = "other",
)::Vector{String}
    channel_names = [
        "af3", "af4", "af7", "af8", "afz",
        "c1", "c2", "c3", "c4", "c5", "c6",
        "cp1", "cp2", "cp3", "cp4", "cp5", "cp6", "cpz", "cz",
        "f1", "f10", "f2", "f3", "f4", "f5", "f6", "f7", "f8", "f9",
        "fc1", "fc2", "fc3", "fc4", "fc5", "fc6", "fcz",
        "fp1", "fp2", "fpz",
        "ft10", "ft7", "ft8", "ft9", "fz", "nz",
        "o1", "o2", "oz",
        "p1", "p10", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9",
        "po3", "po4", "po7", "po8", "poz", "pz",
        "t10", "t3", "t4", "t5", "t6", "t7", "t8", "t9",
        "tp10", "tp7", "tp8", "tp9",
    ]
    ref_channels = ["a1", "a2", "m1", "m2", "pg1", "pg2"]
    eog_channels = ["e", "e1", "e2"]

    channel_type = fill(default, length(clabels))

    for idx in eachindex(clabels)
        lbl = lowercase(clabels[idx])

        occursin("ecg", lbl) && (channel_type[idx] = "ecg")
        occursin("ekg", lbl) && (channel_type[idx] = "ecg")
        occursin("eog", lbl) && (channel_type[idx] = "eog")

        for eog in eog_channels
            lbl == eog && (channel_type[idx] = "eog")
        end

        occursin("emg", lbl) && (channel_type[idx] = "emg")
        lbl in ref_channels && (channel_type[idx] = "ref")
        occursin("mag", lbl) && (channel_type[idx] = "mag")
        occursin("grad", lbl) && (channel_type[idx] = "grad")
        occursin("meg", lbl) && (channel_type[idx] = "meg")

        # EEG channels take priority (e.g. "C3A1" → eeg, not ref)
        lbl in channel_names && (channel_type[idx] = "eeg")
        for ch_name in channel_names
            occursin(ch_name, lbl) && (channel_type[idx] = "eeg")
        end

        length(lbl) >= 1 && lbl[1] in ('c', 'f', 'n', 'o', 'p', 't', 'i') &&
            (channel_type[idx] = "eeg")
        length(lbl) >= 2 && lbl[1:2] == "af" && (channel_type[idx] = "eeg")
        length(lbl) >= 3 && lbl[1] == 'a' && (channel_type[idx] = "eeg")

        # non-neural / auxiliary channels
        for pattern in (
            "rr",
            "mic",
            "flw",
            "tho",
            "abd",
            "sao2",
            "sa02",
            "plr",
            "body",
            "ux",
            "ias",
            "sys",
            "aux",
        )
            occursin(pattern, lbl) && (channel_type[idx] = "other")
        end

        occursin("acc", lbl) && (channel_type[idx] = "accel")

        # marker / event channels
        for pattern in (
            "sti",
            "stim",
            "mark",
            "marker",
            "markers",
            "event",
            "trigger",
            "annotation",
            "annotations",
            "status",
        )
            occursin(pattern, lbl) && (channel_type[idx] = "mrk")
        end
        lbl == "stim" && (channel_type[idx] = "mrk") # already covered by occursin above; kept for explicitness

        occursin("nm", lbl) && (channel_type[idx] = "nirs_od")
    end

    return channel_type
end

"""Sort channel names based on their types."""
function _sort_channels(ch_t::Vector{String})::Vector{Int64}
    # map each channel type to a sort-priority string (lower string = sorted first)
    priority = Dict(
        "meg" => "1", "grad" => "1", "mag" => "2",
        "eeg" => "3", "seeg" => "3", "ecog" => "3",
        "ieeg" => "3", "csd" => "4", "ref" => "4",
        "eog" => "5", "ecg" => "6", "emg" => "7",
        "other" => "8", "mrk" => "9", "stim" => "9",
        # NIRS: each data type forms its own block (the channels of the pairs are
        # interleaved by od2conc(): S1_D1 HbO, S1_D1 HbR, S1_D1 HbT, S1_D2 HbO, …)
        "nirs_int" => "1", "nirs_od" => "2", "nirs_dmean" => "3",
        "nirs_dvar" => "3", "nirs_dskew" => "3", "nirs_mua" => "3",
        "nirs_musp" => "3", "nirs_hbo" => "4", "nirs_hbr" => "5", "nirs_hbt" => "6",
        "nirs_h2o" => "7", "nirs_lipid" => "7", "nirs_bfi" => "7",
        "nirs_hrf_dod" => "8", "nirs_hrf_dmean" => "8", "nirs_hrf_dvar" => "8",
        "nirs_hrf_dskew" => "8", "nirs_hrf_hbo" => "8", "nirs_hrf_hbr" => "8",
        "nirs_hrf_hbt" => "8", "nirs_hrf_bfi" => "8",
        "nirs_aux" => "9", "accel" => "1", "magfld" => "2",
        "orient" => "3", "angvel" => "4",
    )
    ch_order = [get(priority, t, "8") for t in ch_t] # default to "8" for unknowns
    return sortperm(ch_order)
end

"""
    _ch_pairs(obj1, ch1, obj2, ch2, exclude_bad_channels)

Resolve the channel pairs of two-object functions (`ch1[i]` of `obj1` with `ch2[i]` of `obj2`).

Unlike `get_channel()`, the order of the channels given by the user is kept and channels may repeat (e.g. one seed channel paired with several others). Keywords (channel types, groups) and regular expressions are expanded in the order of the channels in the object. If `exclude_bad_channels` is set, pairs in which either channel is marked as bad are removed as whole pairs.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}, Regex}`: channel name(s) in `obj1`
- `obj2::NeuroAnalyzer.NEURO`
- `ch2::Union{String, Vector{String}, Regex}`: channel name(s) in `obj2`
- `exclude_bad_channels::Bool=NeuroAnalyzer.exclude_bad_channels`: if `true`, pairs with a bad channel are removed

# Returns

- `Tuple{Vector{Int64}, Vector{Int64}}`: channel numbers in `obj1` and `obj2`
"""
function _ch_pairs(
    obj1::NeuroAnalyzer.NEURO,
    ch1::Union{String, Vector{String}, Regex},
    obj2::NeuroAnalyzer.NEURO,
    ch2::Union{String, Vector{String}, Regex},
    exclude_bad_channels::Bool = NeuroAnalyzer.exclude_bad_channels,
)::Tuple{Vector{Int64}, Vector{Int64}}
    # resolve each token separately: keeps the order of the tokens and repeated channels;
    # a token expanding to several channels (keyword, regex) is in the object order
    function _resolve(obj::NeuroAnalyzer.NEURO, ch)::Vector{Int64}
        ch isa Vector{String} || return sort(_ch_idx(obj, ch))
        idx = Int64[]
        for token in ch
            append!(idx, sort(_ch_idx(obj, token)))
        end
        return idx
    end

    c1 = _resolve(obj1, ch1)
    c2 = _resolve(obj2, ch2)
    isempty(c1) && throw(ArgumentError("No channels selected in ch1."))
    isempty(c2) && throw(ArgumentError("No channels selected in ch2."))
    length(c1) == length(c2) || throw(
        ArgumentError(
            "Lengths of ch1 ($(length(c1))) and ch2 ($(length(c2))) must be equal.",
        ),
    )

    if exclude_bad_channels
        bad1 = obj1.header.recording[:bad_channel]
        bad2 = obj2.header.recording[:bad_channel]
        keep = [!bad1[i1] && !bad2[i2] for (i1, i2) in zip(c1, c2)]
        if !all(keep)
            dropped = join(
                [labels(obj1)[i1] * "-" * labels(obj2)[i2] for (i1, i2) in
                 zip(c1[.!keep], c2[.!keep])],
                ", ",
            )
            _warn("Pairs with bad channels removed: $dropped")
            c1 = c1[keep]
            c2 = c2[keep]
        end
        isempty(c1) && throw(ArgumentError("No channel pairs without bad channels."))
    end

    return c1, c2
end

"""
    _get_ch(obj, ch, exclude_bad_channels)

Channel numbers of the channel name(s) `ch` (see `get_channel()`); if `exclude_bad_channels` is `true`, the channels marked as bad are not included.
"""
_get_ch(
    obj::NeuroAnalyzer.NEURO,
    ch::Union{String, Vector{String}, Regex},
    exclude_bad_channels::Bool,
)::Vector{Int64} = get_channel(obj; ch = ch, exclude = exclude_bad_channels ? "bad" : "")

"""
    _get_single_ch(obj, ch, exclude_bad_channels)

Channel number of the channel name `ch`, which must resolve to exactly one channel; if `exclude_bad_channels` is `true` and the channel is marked as bad, an error is thrown.
"""
function _get_single_ch(
    obj::NeuroAnalyzer.NEURO,
    ch::Union{String, Vector{String}, Regex},
    exclude_bad_channels::Bool,
)::Int64
    c = get_channel(obj; ch = ch)
    isempty(c) && throw(ArgumentError("No channels selected."))
    length(c) == 1 || throw(ArgumentError("ch must resolve to exactly one channel."))
    exclude_bad_channels && obj.header.recording[:bad_channel][c[1]] && throw(
        ArgumentError(
            "Channel $(labels(obj)[c[1]]) is marked as bad (exclude_bad_channels = true).",
        ),
    )
    return c[1]
end
