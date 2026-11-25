"""
    _delete_markers(markers, seg)

Return a copy of `markers` with all entries whose start time falls within `[seg[1], seg[2]]` removed.

# Arguments

- `markers::DataFrame`: markers table
- `seg::Tuple{Real, Real}`: time range `(start, stop)` in seconds (inclusive)

# Returns

- `DataFrame`: filtered copy
"""
function _delete_markers(markers::DataFrame, seg::Tuple{Real, Real})::DataFrame
    markers_new = copy(markers)
    for idx = DataFrames.nrow(markers_new):-1:1
        if markers_new[idx, :start] >= seg[1] && markers_new[idx, :start] <= seg[2]
            deleteat!(markers_new, idx)
        end
    end
    return markers_new
end

"""
    _shift_markers(markers, seg)

Return a copy of `markers` with all entries whose start time is after `seg[2]` shifted back by `seg[2] - seg[1]` seconds (i.e. the duration of the removed segment).

# Arguments

- `markers::DataFrame`: marker table
- `seg::Tuple{Real, Real}`: removed time range `(start, stop)` in seconds

# Returns

- `DataFrame`: adjusted copy
"""
function _shift_markers(markers::DataFrame, seg::Tuple{Real, Real})::DataFrame
    markers_new = copy(markers)
    duration = seg[2] - seg[1]
    for idx = 1:DataFrames.nrow(markers_new)
        if markers_new[idx, :start] > seg[2]
            markers_new[idx, :start] -= duration
        end
    end
    return markers_new
end

"""
    _get_epoch_markers(obj)

Return a vector of time points (in seconds) marking the start of each epoch in `obj`.
"""
function _get_epoch_markers(obj::NeuroAnalyzer.NEURO)::Vector{Float64}
    return s2t.(collect(1:epoch_len(obj):(epoch_len(obj) * nepochs(obj))), sr(obj))
end

"""
    _has_markers(channel_types)

Check whether a `"mrk"` channel exists in `channel_types`.

# Arguments

- `channel_types::Vector{String}`: vector of channel type strings

# Returns

- `Tuple{Bool, Int64}`: `(has_markers, marker_channel_index)`; `marker_channel_index` is `0` if no marker channel is found
"""
function _has_markers(channel_types::Vector{String})::Tuple{Bool, Int64}
    idx = findfirst(==("mrk"), channel_types)
    isnothing(idx) && return false, 0
    return true, idx
end

"""
    _has_markers(obj)

Return `true` if `obj.markers` is non-empty.
"""
_has_markers(obj::NeuroAnalyzer.NEURO)::Bool = !isempty(obj.markers)

"""
    _a2df(annotations)

Convert a vector of raw EDF/BDF annotation strings into a standardized markers DataFrame.

Each string (the annotation channel of one data record) is parsed into time-stamped annotation lists (TALs, EDF+ specification): TALs are separated by 0x00; a TAL is `onset[0x15 duration]0x14 annotation 0x14 [annotation 0x14 …]`. Every non-empty annotation of a TAL is one event (several annotations may share an onset, and a data record may contain several TALs); empty annotations (the time-keeping TAL of each data record) are skipped.

# Arguments

- `annotations::Vector{String}`: raw annotation strings from an EDF/BDF file

# Returns

- `DataFrame` with columns `id`, `start`, `length`, `value`, `channel`

# Notes

- the marker type is `"annotation"`, or a standard type recognized in the text (e.g. `"BAD_ACQ_SKIP"`: `"BAD"`, see `_marker_type()`)
- Non-ASCII labels trigger a warning but are not dropped
"""
function _a2df(annotations::Vector{String})::DataFrame
    a_start  = Float64[]
    a_length = Float64[]
    a_event  = String[]

    # TAL: onset (+/- number), optional 0x15 duration, 0x14, annotations terminated by 0x14
    tal_re = r"^\s*([+-]\d+(?:\.\d*)?)(?:\x15(\d*(?:\.\d*)?))?\x14(.*)$"s
    for rec in annotations, tal in split(rec, '\0')
        m = match(tal_re, tal)
        isnothing(m) && continue
        onset = parse(Float64, m.captures[1])
        dur = m.captures[2]
        dur = (isnothing(dur) || isempty(dur)) ? 0.0 : something(tryparse(Float64, dur), 0.0)
        for txt in split(m.captures[3], '\x14')
            txt = strip(txt)
            isempty(txt) && continue
            push!(a_start, onset)
            push!(a_length, dur)
            push!(a_event, txt)
        end
    end

    isempty(a_event) && return DataFrame(
        :type    => String[],
        :start   => Float64[],
        :length  => Float64[],
        :value   => String[],
        :channel => Vector{String}[],
    )

    all(isascii.(a_event)) || _warn("Unicode labels were not converted.")

    return DataFrame(
        :type    => [_marker_type(String(ev); default = "annotation") for ev in a_event],
        :start   => a_start,
        :length  => a_length,
        :value   => a_event,
        :channel => [String[] for _ = 1:length(a_event)],
    )
end

"""
    _remap_markers(markers, fs, f)

Map markers through a transformation of sample positions.

`f(s)` receives the 1-based sample number of a marker start (in the concatenated signal, i.e. `t2s(start, fs)`) and returns its new sample number, or `nothing` to drop the marker. Marker lengths and other columns are kept.

# Returns

- `DataFrame`: remapped markers
"""
function _remap_markers(markers::DataFrame, fs::Int64, f::Function)::DataFrame
    m = copy(markers)
    keep = trues(DataFrames.nrow(m))
    for idx = 1:DataFrames.nrow(m)
        s_new = f(t2s(m[idx, :start], fs))
        if isnothing(s_new)
            keep[idx] = false
        else
            m[idx, :start] = s2t(s_new, fs)
        end
    end
    return m[keep, :]
end

# ---------------------------------------------------------------------------
# marker channels: every marker holds a vector of channel labels (`channel` column,
# `Vector{String}`); an empty vector means that the marker applies to all channels

"""
    _markers_df()

Return an empty markers table (columns `type`, `start`, `length`, `value`, `channel`).
"""
_markers_df()::DataFrame = DataFrame(
    :type    => String[],
    :start   => Float64[],
    :length  => Float64[],
    :value   => String[],
    :channel => Vector{String}[],
)

"""
    _mrk_ch(idx, clabels)

Channel labels of a marker given as a channel number (`0`, or a number outside `clabels`: all channels).

# Arguments

- `idx::Integer`: channel number
- `clabels::Vector{String}`: channel labels

# Returns

- `Vector{String}`: `[clabels[idx]]`, or `String[]` (all channels)
"""
_mrk_ch(idx::Integer, clabels::Vector{String})::Vector{String} =
    idx in eachindex(clabels) ? [clabels[idx]] : String[]

"""
    _mrk_ch(ch, clabels)

Validate the channel(s) of a marker given as labels.

# Arguments

- `ch::Union{String, Vector{String}}`: channel label(s); `""` or `String[]`: all channels
- `clabels::Vector{String}`: channel labels of the object

# Returns

- `Vector{String}`: unique channel labels (`String[]`: all channels)
"""
function _mrk_ch(ch::Union{String, Vector{String}}, clabels::Vector{String})::Vector{String}
    ch = ch isa String ? (ch == "" ? String[] : [ch]) : unique(ch)
    for c in ch
        c in clabels || throw(ArgumentError("Channel $c does not exist."))
    end
    return ch
end

# string of the marker channels (for printing and CSV export)
_mrk_ch_str(ch::Vector{String})::String = isempty(ch) ? "" : join(ch, ";")

"""
    _markers_relabel!(markers, old, new)

Rename the marker channels after channel labels were changed (`old[i]` → `new[i]`).
"""
function _markers_relabel!(
    markers::DataFrame,
    old::AbstractVector{<:AbstractString},
    new::AbstractVector{<:AbstractString},
)::Nothing
    DataFrames.nrow(markers) == 0 && return nothing
    length(old) == length(new) ||
        throw(ArgumentError("old and new labels must have the same length."))
    m = Dict(String(o) => String(n) for (o, n) in zip(old, new) if o != n)
    isempty(m) && return nothing
    # new vectors are assigned: rows may share one vector object
    markers[!, :channel] = Vector{String}[[get(m, c, c) for c in ch] for ch in markers[!, :channel]]
    return nothing
end

"""
    _markers_delete_channels!(markers, del)

Remove deleted channels from the marker channels. A marker whose channels were all deleted is removed (otherwise its empty channel list would mean all channels); markers of all channels are kept.

# Returns

- `Int64`: number of removed markers
"""
function _markers_delete_channels!(
    markers::DataFrame,
    del::AbstractVector{<:AbstractString},
)::Int64
    (DataFrames.nrow(markers) == 0 || isempty(del)) && return 0
    del = Set(String.(del))
    drop = falses(DataFrames.nrow(markers))
    # new vectors are assigned: rows may share one vector object
    new_ch = Vector{String}[]
    for (idx, ch) in enumerate(markers[!, :channel])
        nch = String[c for c in ch if !(c in del)]
        drop[idx] = !isempty(ch) && isempty(nch)
        push!(new_ch, nch)
    end
    markers[!, :channel] = new_ch
    any(drop) && deleteat!(markers, findall(drop))
    return count(drop)
end

"""
    _markers_upgrade!(obj)

Convert the markers of an object saved by an older version of NeuroAnalyzer: the `id` column is renamed to `type` and its values are mapped to the standard marker types (see `_marker_type()`), channel numbers (`0` = all channels) are converted to channel labels.
"""
function _markers_upgrade!(obj::NeuroAnalyzer.NEURO)::Nothing
    # marker `id` (before v0.26.12): free-form → standard marker type
    if "id" in names(obj.markers) && !("type" in names(obj.markers))
        DataFrames.rename!(obj.markers, :id => :type)
        obj.markers[!, :type] = String[
            _marker_type_old(string(t), string(v)) for
            (t, v) in zip(obj.markers[!, :type], obj.markers[!, :value])
        ]
    end
    "channel" in names(obj.markers) || return nothing
    eltype(obj.markers[!, :channel]) <: AbstractVector{<:AbstractString} && return nothing
    clabels = labels(obj)
    obj.markers[!, :channel] = Vector{String}[
        _mrk_ch(Int64(c), clabels) for c in obj.markers[!, :channel]
    ]
    return nothing
end

"""
    _markers_map_channels!(markers, m)

Replace the marker channels using the map `m` (label → labels of the channels derived from it, e.g. montage channels); labels not in `m` are kept. A marker whose channels map to no channel is removed; markers of all channels are kept.

# Returns

- `Int64`: number of removed markers
"""
function _markers_map_channels!(markers::DataFrame, m::Dict{String, Vector{String}})::Int64
    DataFrames.nrow(markers) == 0 && return 0
    drop = falses(DataFrames.nrow(markers))
    new_ch = Vector{String}[]
    for (idx, ch) in enumerate(markers[!, :channel])
        nch = unique(reduce(vcat, [get(m, c, [c]) for c in ch]; init = String[]))
        drop[idx] = !isempty(ch) && isempty(nch)
        push!(new_ch, nch)
    end
    markers[!, :channel] = new_ch
    any(drop) && deleteat!(markers, findall(drop))
    return count(drop)
end

# ---------------------------------------------------------------------------
# marker types: `type` is the kind of marker, `value` its code or condition (used by `epoch()`)

"""Standard marker types (the `type` column of `obj.markers`)."""
const _MARKER_TYPES = ["stim", "response", "annotation", "BAD", "boundary", "epoch", "other"]

"""
    _marker_type(t; <keyword arguments>)

Map a marker type name used by a file format or another program to a standard marker type (`_MARKER_TYPES`).

| standard type | recognized names (case-insensitive) |
| ------------- | ----------------------------------- |
| `stim`        | stim, stimulus, trigger, event, mrk, marker |
| `response`    | response, resp, button, keypress |
| `annotation`  | annotation, comment, note |
| `BAD`         | BAD, BAD_*, BAD … (e.g. MNE `BAD_ACQ_SKIP`) |
| `boundary`    | boundary, new segment, EDGE, EDGE boundary, BAD boundary, deleted |
| `epoch`       | epoch |

# Arguments

- `t::AbstractString`: type name
- `default::String="other"`: type of names that are not recognized

# Returns

- `String`
"""
function _marker_type(t::AbstractString; default::String = "other")::String
    tl = lowercase(strip(t))
    tl in ("stim", "stimulus", "trigger", "event", "mrk", "marker") && return "stim"
    tl in ("response", "resp", "button", "keypress") && return "response"
    tl in ("annotation", "comment", "note") && return "annotation"
    tl in ("boundary", "new segment", "edge", "edge boundary", "bad boundary", "deleted") &&
        return "boundary"
    tl == "epoch" && return "epoch"
    occursin(r"^bad([_\s]|$)", tl) && return "BAD"
    return default
end

# marker `id` of objects saved before v0.26.12 → marker type
function _marker_type_old(id::String, value::String)::String
    # boundary markers placed by trim() (id "NA", value "DELETED")
    value == "DELETED" && return "boundary"
    t = _marker_type(id; default = "")
    isempty(t) || return t
    # numeric IDs (one per annotation text: EDF+, DigiTrack, SNIRF, NIRS) and other names
    return occursin(r"^\d+$", id) ? "stim" : "other"
end

"""
    _gdf_marker_type(code)

Standard marker type of a GDF event code (BioSig `eventcodes.txt`): artifacts (`0x0101`–`0x01FF`): `BAD`, start of a new segment (`0x7FFE`): `boundary`, user-specific events and trial/cue codes (`0x0001`–`0x00FF`, `0x0300`–`0x03FF`): `stim`, other codes (e.g. sleep stages, ECG waves): `annotation`.
"""
function _gdf_marker_type(code::UInt16)::String
    c = code & 0x7fff
    0x0101 <= c <= 0x01ff && return "BAD"
    c == 0x7ffe && return "boundary"
    (0x0001 <= c <= 0x00ff || 0x0300 <= c <= 0x03ff) && return "stim"
    return "annotation"
end

"""
    _fiff_an_desc(type, value)

Description of the MNE annotation of a marker (`export_fiff()`): the value (or the type, if the value is empty); `BAD` markers: `BAD_<value>` (unless the value starts with `BAD`), `boundary` markers: `BAD boundary` (MNE-Python treats both as bad segments).
"""
function _fiff_an_desc(type::String, value::String)::String
    type == "boundary" && return "BAD boundary"
    if type == "BAD"
        startswith(uppercase(value), "BAD") && return value
        return isempty(value) ? "BAD" : "BAD_" * value
    end
    return isempty(value) ? type : value
end

# ---------------------------------------------------------------------------
# extra marker columns: a markers table has the 5 fixed columns (`type`, `start`, `length`,
# `value`, `channel`) followed by optional extra columns (e.g. `trial_type`, `response_time`
# from BIDS events.tsv); rows added by NeuroAnalyzer get `missing` in the extra columns

# names of the fixed marker columns
const _MARKER_FIXED = ["type", "start", "length", "value", "channel"]

# names of the extra columns of a markers table
_marker_extra_cols(markers::DataFrame)::Vector{String} = setdiff(names(markers), _MARKER_FIXED)

"""
    _markers_push!(markers, row)

Append a marker (`NamedTuple` with the fixed columns) to `markers`; extra columns get `missing`.
"""
function _markers_push!(markers::DataFrame, row::NamedTuple)::Nothing
    push!(markers, row; cols = :subset, promote = true)
    return nothing
end

"""
    _markers_append!(markers, new)

Append the rows of `new` to `markers`; columns missing in either table get `missing`.
"""
function _markers_append!(markers::DataFrame, new::DataFrame)::Nothing
    append!(markers, new; cols = :union, promote = true)
    return nothing
end

"""
    _bids_events(events, fs)

Markers from a BIDS `events.tsv` table.

- `start`, `length`: `onset` and `duration` (seconds), or `sample` (0-based) / `fs`
- `type`, `value`: a `trial_type` that names a marker type (e.g. `stimulus`, `response`, `BAD_…`, see `_marker_type()`) gives the type and the value is the event code (`value`); a condition (e.g. `face`) gives type `stim` and the value is the condition (used by `epoch()`); without `trial_type` the type is `stim` and the value is the event code
- all other columns are kept as extra columns (`value` as `bids_value`)

# Arguments

- `events::DataFrame`: events table (`n/a` read as `missing`)
- `fs::Real`: sampling rate in Hz

# Returns

- `DataFrame`: markers
"""
function _bids_events(events::DataFrame, fs::Real)::DataFrame
    n = DataFrames.nrow(events)
    cols = names(events)
    col(c, default) = c in cols ? [ismissing(x) ? default : x for x in events[!, c]] : fill(default, n)
    start = if "onset" in cols
        Float64.(col("onset", 0.0))
    else
        Float64.(col("sample", 0)) ./ fs
    end
    len = Float64.(col("duration", 0.0))
    tt = string.(col("trial_type", ""))
    ev = string.(col("value", ""))
    mtype = [_marker_type(t; default = "stim") for t in tt]
    mvalue = [
        (isempty(t) || _marker_type(t; default = "") != "") ? (isempty(v) ? t : v) : t for
        (t, v) in zip(tt, ev)
    ]
    markers = DataFrame(
        :type => mtype,
        :start => start,
        :length => len,
        :value => mvalue,
        :channel => [String[] for _ = 1:n],
    )
    for c in cols
        c in ("onset", "duration") && continue
        name = c in _MARKER_FIXED ? "bids_" * c : c
        markers[!, name] = events[!, c]
    end
    return markers
end

"""
    _marker_select(markers, value, type)

Markers with the value(s) `value` (a value, a list of values or a regular expression; `""`: all markers) and the type `type` (`""`: any type).

# Returns

- `BitVector`
"""
function _marker_select(
    markers::DataFrame,
    value::Union{String, Vector{String}, Regex},
    type::String,
)::BitVector
    v = string.(markers[!, :value])
    sel = if value isa Regex
        occursin.(value, v)
    elseif value isa Vector
        in(Set(value)).(v)
    else
        value == "" ? trues(length(v)) : v .== value
    end
    type == "" || (sel .&= markers[!, :type] .== type)
    return BitVector(sel)
end

"""
    _marker_idx(markers, n, value, type; required=false)

Numbers of the markers selected by number (`n`) or by value and type (`nothing`: any; a string, a vector of strings or a regex). Without any selection all markers are selected, or an error is thrown if `required`.
"""
function _marker_idx(
    markers::DataFrame,
    n::Union{Nothing, Int64, AbstractVector{Int64}},
    value::Union{Nothing, String, Vector{String}, Regex},
    type::Union{Nothing, String, Vector{String}, Regex};
    required::Bool = false,
)::Vector{Int64}
    nn = DataFrames.nrow(markers)
    if !isnothing(n)
        isnothing(value) && isnothing(type) ||
            throw(ArgumentError("n cannot be combined with value or type."))
        idx = sort(unique(n isa Integer ? [n] : collect(n)))
        isempty(idx) && throw(ArgumentError("n must not be empty."))
        for i in idx
            _in(i, (1, nn), "n")
        end
        return idx
    end
    required && isnothing(value) && isnothing(type) &&
        throw(ArgumentError("n, value or type must be given."))
    m(x, s) = isnothing(s) ? true : s isa Regex ? occursin(s, x) : s isa Vector ? x in s : x == s
    return Int64[
        i for i in 1:nn if
        m(string(markers[i, :value]), value) && m(string(markers[i, :type]), type)
    ]
end

"""
    _marker_conditions(markers, condition)

Condition label of each marker: the values of the marker column `condition` (`missing`: `""`), or the marker values mapped by the dictionary `condition` (values not in the dictionary are their own labels).

# Returns

- `Vector{String}`
"""
function _marker_conditions(
    markers::DataFrame,
    condition::Union{String, AbstractDict},
)::Vector{String}
    v = String[string(x) for x in markers[!, :value]]
    condition isa AbstractDict && return String[string(get(condition, x, x)) for x in v]
    condition in names(markers) || throw(
        ArgumentError(
            "Marker column $condition does not exist; available: $(join(names(markers), ", ")).",
        ),
    )
    condition == "channel" && throw(ArgumentError("condition cannot be the channel column."))
    return String[ismissing(x) ? "" : string(x) for x in markers[!, condition]]
end

"""
    _mark_bad!(obj, t1, t2, ch)

Mark the segment `t1`–`t2` (seconds) as bad: add a marker of type `"BAD"` (value `"manual"`) for the channels `ch` (labels; empty: all channels). Used by the interactive signal plot (`b`).
"""
function _mark_bad!(
    obj::NeuroAnalyzer.NEURO,
    t1::Real,
    t2::Real,
    ch::Vector{String},
)::Nothing
    t2 > t1 || throw(ArgumentError("The end of the segment must be after its start."))
    add_marker!(obj; type = "BAD", start = t1, len = t2 - t1, value = "manual", ch = ch)
    return nothing
end

"""
    _unmark_bad!(obj, t1, t2)

Delete the markers of type `"BAD"` overlapping the segment `t1`–`t2` (seconds). Used by the interactive signal plot (`u`).

# Returns

- `Int64`: number of deleted markers
"""
function _unmark_bad!(obj::NeuroAnalyzer.NEURO, t1::Real, t2::Real)::Int64
    m = obj.markers
    del = findall(
        i -> m[i, :type] == "BAD" && m[i, :start] < t2 && m[i, :start] + m[i, :length] > t1,
        1:DataFrames.nrow(m),
    )
    isempty(del) || delete_marker!(obj; n = del)
    return length(del)
end

"""
    _bad_marker_epochs(obj, ch)

Epochs that contain a bad segment: a marker of type `"BAD"` of all channels, or of at least one of the channels `ch` (labels), overlapping the epoch (marker times are on the timeline of the concatenated epochs).

# Returns

- `Vector{Int64}`: sorted epoch numbers
"""
function _bad_marker_epochs(obj::NeuroAnalyzer.NEURO, ch::Vector{String})::Vector{Int64}
    fs = sr(obj)
    n = epoch_len(obj)
    ep_n = nepochs(obj)
    ep = Int64[]
    for m in eachrow(obj.markers)
        m[:type] == "BAD" || continue
        (isempty(m[:channel]) || any(in(ch), m[:channel])) || continue
        # first and last sample of the marker (1-based), its epochs
        s1 = t2s(m[:start], fs)
        s2 = max(s1, round(Int64, (m[:start] + m[:length]) * fs))
        e1 = clamp((s1 - 1) ÷ n + 1, 1, ep_n)
        e2 = clamp((s2 - 1) ÷ n + 1, 1, ep_n)
        append!(ep, e1:e2)
    end
    return sort(unique(ep))
end

"""
    _bad_samples(obj, ch)

Samples of a continuous object covered by bad segments: markers of type `"BAD"` of all channels or of at least one of the channels `ch` (channel numbers).

# Returns

- `BitVector`: `true` for the samples in bad segments
"""
function _bad_samples(obj::NeuroAnalyzer.NEURO, ch::AbstractVector{<:Integer})::BitVector
    n = epoch_len(obj)
    fs = sr(obj)
    lbl = Set(labels(obj)[ch])
    bad = falses(n)
    for m in eachrow(obj.markers)
        m[:type] == "BAD" || continue
        (isempty(m[:channel]) || any(in(lbl), m[:channel])) || continue
        s1 = clamp(t2s(m[:start], fs), 1, n)
        s2 = clamp(max(s1, round(Int64, (m[:start] + m[:length]) * fs)), 1, n)
        bad[s1:s2] .= true
    end
    return bad
end

"""
    _seg_data(obj, ch, exclude_bad_segments)

Signal of the channels `ch` (channel numbers), shape (channels, samples, epochs). For a continuous object with `exclude_bad_segments = true`, the samples of bad segments (markers of type `"BAD"` of all channels or of any of the channels `ch`) are removed and the remaining parts are joined (as `reject_by_annotation = "omit"` of MNE-Python); epoched objects are returned unchanged (see the `reject` column of the epoch metadata).
"""
function _seg_data(
    obj::NeuroAnalyzer.NEURO,
    ch::AbstractVector{<:Integer},
    exclude_bad_segments::Bool,
)::AbstractArray{Float64, 3}
    (exclude_bad_segments && nepochs(obj) == 1) || return @view(obj.data[ch, :, :])
    bad = _bad_samples(obj, ch)
    any(bad) || return @view(obj.data[ch, :, :])
    all(bad) && throw(ArgumentError("All samples are in bad segments."))
    _info("$(round(count(bad) / sr(obj); digits = 2)) s of bad segments excluded")
    return obj.data[ch, .!bad, :]
end
