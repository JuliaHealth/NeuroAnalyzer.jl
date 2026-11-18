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

- IDs are assigned per unique event label (matching events share the same numeric ID)
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
        :id      => String[],
        :start   => Float64[],
        :length  => Float64[],
        :value   => String[],
        :channel => Int64[],
    )

    all(isascii.(a_event)) || _warn("Unicode labels were not converted.")

    # build numeric IDs: each unique event label gets one integer ID, shared by all occurrences of that label
    unique_events = unique(a_event)
    event_id_map = Dict(ev => i for (i, ev) in enumerate(unique_events))
    id = [event_id_map[ev] for ev in a_event]

    return DataFrame(
        :id      => string.(id),
        :start   => a_start,
        :length  => a_length,
        :value   => a_event,
        :channel => zeros(Int64, length(a_event)),
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
