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
    for idx in DataFrames.nrow(markers_new):-1:1
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
    for idx in 1:DataFrames.nrow(markers_new)
        if markers_new[idx, :start] > seg[2]
            markers_new[idx, :start] -= duration
        end
    end
    return markers_new
end

"""
    _get_epoch_markers(obj)

Return a vector of time points (in seconds) marking the start of each epoch in `obj`, rounded to 4 decimal places.
"""
function _get_epoch_markers(obj::NeuroAnalyzer.NEURO)::Vector{Float64}
    return round.(
        s2t.(collect(1:epoch_len(obj):(epoch_len(obj) * nepochs(obj))), sr(obj));
        digits = 4,
    )
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

Each annotation is pre-processed to normalize control characters (0x14, 0x15, 0x00) to pipe `|` separators or removed entirely. The resulting token stream is parsed into `(start, length, event)` triples.

# Arguments

- `annotations::Vector{String}`: raw annotation strings from an EDF/BDF file

# Returns

- `DataFrame` with columns `id`, `start`, `length`, `value`, `channel`

# Notes

- IDs are assigned per unique event label (matching events share the same numeric ID)
- Non-ASCII labels trigger a warning but are not dropped
"""
function _a2df(annotations::Vector{String})::DataFrame
    # normalize control characters to "|" separators
    mrk = replace.(annotations,
        "\x14\x14\0" => "|",
        "\x14\x14"   => "|",
        "\x14"       => "|",
        "\x15"       => "|",
        "\0"         => "",
        r"\|$"       => "",
    )
 
    # remove entries that are empty or contain no pipe separators
    for idx in length(mrk):-1:1
        if isempty(mrk[idx]) || !occursin('|', mrk[idx])
            deleteat!(mrk, idx)
        end
    end
 
    a_start  = Float64[]
    a_length = Float64[]
    a_event  = String[]
 
    if length(mrk) == 1
        s = split(mrk[1], "|")
 
        for idx in length(s):-1:1
            s[idx] == "" && deleteat!(s, idx)
        end
 
        if length(s) % 3 == 0
            for idx in 1:3:(length(s) - 2)
                push!(a_start,  parse(Float64, strip(s[idx])))
                push!(a_length, parse(Float64, strip(s[idx + 1])))
                push!(a_event,  strip(s[idx + 2]))
            end
        else

            # TODO: use offset if provided

            _offset = parse(Float64, strip(s[1]))
            for idx in 2:3:(length(s) - 1)
                push!(a_start,  parse(Float64, strip(s[idx])))
                push!(a_length, parse(Float64, strip(s[idx + 1])))
                push!(a_event,  strip(s[idx + 2]))
            end
        end
 
    else
        for idx in eachindex(mrk)
            s = split(mrk[idx], "|")
            # drop leading annotation-number column when both first two tokens start with '+'
            length(s) >= 2 && s[1][1] == '+' && s[2][1] == '+' && (s = s[2:end])
 
            if length(s) == 3
                push!(a_start,  parse(Float64, strip(s[1])))
                push!(a_length, parse(Float64, strip(s[2])))
                push!(a_event,  strip(s[3]))
            elseif length(s) == 2
                push!(a_start,  parse(Float64, strip(s[1])))
                push!(a_length, 0.0)
                push!(a_event,  strip(s[2]))
            end
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