export view_marker
export delete_marker
export delete_marker!
export add_marker
export add_marker!
export edit_marker
export edit_marker!
export channel2marker
export channel2marker!
export add_markers
export add_markers!


"""
    _check_marker_cols(markers)

Throw `ArgumentError` if `markers` does not start with the fixed marker columns (`type`, `start`, `length`, `value`, `channel`; `id` instead of `type` before v0.26.12); extra columns may follow.
"""
function _check_marker_cols(markers::DataFrame)::Nothing
    n = names(markers)
    k = length(_MARKER_FIXED)
    (length(n) >= k && (n[1:k] == _MARKER_FIXED || n[1:k] == ["id"; _MARKER_FIXED[2:end]])) ||
        throw(
            ArgumentError(
                "Markers DataFrame must start with columns: $(_MARKER_FIXED) (extra columns may follow); got $n.",
            ),
        )
    return nothing
end

"""
    _marker_channels(markers, clabels)

Return a copy of `markers` with the `channel` column as vectors of channel labels (`String[]`: all channels). Accepted input: vectors of labels, strings (labels separated by `;`, `""`: all channels, as written by `export_markers()`) or channel numbers (`0`: all channels). Channel labels must exist in `clabels`.
"""
function _marker_channels(markers::DataFrame, clabels::Vector{String})::DataFrame
    m = copy(markers)
    # `id` column of NeuroAnalyzer before v0.26.12
    if "id" in names(m)
        DataFrames.rename!(m, :id => :type)
        m[!, :type] = String[_marker_type_old(string(t), string(v)) for (t, v) in zip(m[!, :type], m[!, :value])]
    end
    m[!, :channel] = Vector{String}[
        if c isa Integer
            0 <= c <= length(clabels) ||
                throw(ArgumentError("Marker channel $c must be in [0, $(length(clabels))]."))
            _mrk_ch(c, clabels)
        elseif c isa AbstractString
            _mrk_ch(String[strip(x) for x in split(c, ';') if strip(x) != ""], clabels)
        elseif c isa AbstractVector
            _mrk_ch(String.(c), clabels)
        elseif ismissing(c)
            String[]
        else
            throw(ArgumentError("Marker channel must be a label, a vector of labels or a channel number; got $(typeof(c))."))
        end for c in markers[!, :channel]
    ]
    return m
end

"""
    view_marker(obj)

Print a formatted table of all markers (extra columns as `name=value`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Nothing`
"""
function view_marker(obj::NeuroAnalyzer.NEURO)::Nothing
    _has_markers(obj) || throw(ArgumentError("OBJ has no markers."))

    println(
        rpad("n", 5) *
        rpad("type", 14) *
        rpad("start [s]", 12) *
        rpad("length [s]", 12) *
        rpad("value", 24) *
        rpad("channel", 1),
    )
    ex = _marker_extra_cols(obj.markers)
    for i = 1:DataFrames.nrow(obj.markers)
        println(
            rpad(string(i), 5) *
            rpad(obj.markers[i, :type], 14) *
            rpad(string(round(obj.markers[i, :start]; digits = 3)), 12) *
            rpad(string(round(obj.markers[i, :length]; digits = 3)), 12) *
            rpad("'" * obj.markers[i, :value] * "'", 24) *
            (isempty(obj.markers[i, :channel]) ? "all" : _mrk_ch_str(obj.markers[i, :channel])) *
            # extra columns: name=value
            join(["  $c=$(obj.markers[i, c])" for c in ex if !ismissing(obj.markers[i, c])]),
        )
    end
    return nothing
end

"""
    delete_marker(obj; <keyword arguments>)

Delete a marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `n::Int64`: marker number

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function delete_marker(obj::NeuroAnalyzer.NEURO; n::Int64)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    delete_marker!(obj_new; n = n)

    return obj_new
end

"""
    delete_marker!(obj; <keyword arguments>)

Delete a marker in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `n::Int64`: marker number

# Returns

- `Nothing`
"""
function delete_marker!(obj::NeuroAnalyzer.NEURO; n::Int64)::Nothing
    _has_markers(obj) || throw(ArgumentError("OBJ has no markers."))

    nn = DataFrames.nrow(obj.markers)
    _in(n, (1, nn), "n")

    deleteat!(obj.markers, n)

    _history!(obj, "delete_marker", (; n))

    return nothing
end

"""
    add_marker(obj; <keyword arguments>)

Append a new marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `type::String="stim"`: marker type: `"stim"`, `"response"`, `"annotation"`, `"BAD"` (bad segment), `"boundary"` (discontinuity), `"epoch"` or `"other"`
- `start::Real`: marker start time in seconds (must be ≥ 0 and within signal)
- `len::Real=0.0`: marker duration in seconds (must be ≥ 0)
- `value::String`: marker value
- `ch::Union{String, Vector{String}}=""`: channel label(s); `""` (or `String[]`) means the marker applies to all channels

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function add_marker(
    obj::NeuroAnalyzer.NEURO;
    type::String = "stim",
    start::Real,
    len::Real = 0.0,
    value::String,
    ch::Union{String, Vector{String}} = "",
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    add_marker!(obj_new; type = type, start = start, len = len, value = value, ch = ch)

    return obj_new
end

"""
    add_marker!(obj; <keyword arguments>)

Add a marker in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `type::String="stim"`: marker type: `"stim"`, `"response"`, `"annotation"`, `"BAD"` (bad segment), `"boundary"` (discontinuity), `"epoch"` or `"other"`
- `start::Real`: marker start time in seconds (must be ≥ 0 and within signal)
- `len::Real=0.0`: marker duration in seconds (must be ≥ 0)
- `value::String`: marker value
- `ch::Union{String, Vector{String}}=""`: channel label(s); `""` (or `String[]`) means the marker applies to all channels

# Returns

- `Nothing`
"""
function add_marker!(
    obj::NeuroAnalyzer.NEURO;
    type::String = "stim",
    start::Real,
    len::Real = 0.0,
    value::String,
    ch::Union{String, Vector{String}} = "",
)::Nothing
    # validate
    start >= 0 || throw(ArgumentError("start must be ≥ 0."))
    len >= 0 || throw(ArgumentError("len must be ≥ 0."))
    start <= obj.time_pts[end] ||
        throw(ArgumentError("start must be ≤ $(obj.time_pts[end])."))
    start + len <= obj.time_pts[end] ||
        throw(ArgumentError("start + len must be ≤ $(obj.time_pts[end])."))
    # empty = marker for all channels
    mch = _mrk_ch(ch, labels(obj))

    # extra marker columns get missing
    _markers_push!(
        obj.markers,
        (type = type, start = Float64(start), length = Float64(len), value = value, channel = mch),
    )
    sort!(obj.markers, :start)

    _history!(obj, "add_marker", (; type, start, len, value, ch))

    return nothing
end

"""
    edit_marker(obj; <keyword arguments>)

Edit a marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `n::Int64`: marker number to edit
- `type::String="stim"`: marker type: `"stim"`, `"response"`, `"annotation"`, `"BAD"` (bad segment), `"boundary"` (discontinuity), `"epoch"` or `"other"`
- `start::Real`: marker start time in seconds (must be ≥ 0 and within signal)
- `len::Real=0.0`: marker duration in seconds (must be ≥ 0)
- `value::String`: marker value
- `ch::Union{String, Vector{String}}=""`: channel label(s); `""` (or `String[]`) means the marker applies to all channels

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function edit_marker(
    obj::NeuroAnalyzer.NEURO;
    n::Int64,
    type::String = "stim",
    start::Real,
    len::Real = 0.0,
    value::String,
    ch::Union{String, Vector{String}} = "",
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    edit_marker!(
        obj_new; n = n, type = type, start = start, len = len, value = value, ch = ch,
    )

    return obj_new
end

"""
    edit_marker!(obj; <keyword arguments>)

Edit a marker in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `n::Int64`: marker number to edit
- `type::String="stim"`: marker type: `"stim"`, `"response"`, `"annotation"`, `"BAD"` (bad segment), `"boundary"` (discontinuity), `"epoch"` or `"other"`
- `start::Real`: marker start time in seconds (must be ≥ 0 and within signal)
- `len::Real=0.0`: marker duration in seconds (must be ≥ 0)
- `value::String`: marker value
- `ch::Union{String, Vector{String}}=""`: channel label(s); `""` (or `String[]`) means the marker applies to all channels

# Returns

- `Nothing`
"""
function edit_marker!(
    obj::NeuroAnalyzer.NEURO;
    n::Int64,
    type::String = "stim",
    start::Real,
    len::Real = 0.0,
    value::String,
    ch::Union{String, Vector{String}} = "",
)::Nothing
    _has_markers(obj) || throw(ArgumentError("OBJ has no markers."))
    start >= 0 || throw(ArgumentError("start must be > 0."))
    len >= 0 || throw(ArgumentError("len must be ≥ 0."))
    start < signal_len(obj) / sr(obj) ||
        throw(ArgumentError("start must be < $(signal_len(obj) / sr(obj))."))
    start + len <= signal_len(obj) / sr(obj) ||
        throw(ArgumentError("start + len must be ≤ $(signal_len(obj) / sr(obj))."))

    nn = size(obj.markers, 1)
    1 <= n <= nn || throw(ArgumentError("n must be in [1, $nn]."))
    mch = _mrk_ch(ch, labels(obj))

    obj.markers[n, :type] = type
    obj.markers[n, :start] = start
    obj.markers[n, :length] = len
    obj.markers[n, :value] = value
    obj.markers[n, :channel] = mch
    sort!(obj.markers, :start)
    _history!(obj, "edit_marker", (; n, type, start, len, value, ch))

    return nothing
end

"""
    channel2marker(obj; <keyword arguments>)

Return a copy of `obj` with events detected in an event channel added as markers.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `v::Real=1.0`: signal value interpreted as an event
- `type::String="stim"`: marker type
- `value::String=""`: marker value string (default: channel label)

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function channel2marker(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    v::Real = 1.0,
    type::String = "stim",
    value::String = "",
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    channel2marker!(obj_new; ch = ch, v = v, type = type, value = value)

    return obj_new
end

"""
    channel2marker!(obj; <keyword arguments>)

Convert an event channel to markers in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `v::Real=1.0`: signal value interpreted as an event
- `type::String="stim"`: marker type
- `value::String=""`: marker value string (default: channel label)

# Returns

- `Nothing`
"""
function channel2marker!(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    v::Real = 1.0,
    type::String = "stim",
    value::String = "",
)::Nothing
    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    length(ch) == 1 || throw(ArgumentError("ch must resolve to exactly one channel."))
    ch = ch[1]

    # check if the event channel contain events
    ev_ch = obj.data[ch, :, :][:]
    length(unique(ev_ch)) > 1 ||
        throw(ArgumentError("Channel $ch does not contain events."))

    # extract events
    ev_v = unique(ev_ch)
    v in ev_v || throw(ArgumentError("Event channel does not contain value $v."))
    _info("Event channel contains values: $ev_v")

    # edges of the samples equal to v (a trigger channel may contain several event codes,
    # so jumps of the channel value itself are not reliable: e.g. 1 → 32 is the end of
    # event 1 and the start of event 32)
    ev_segs = diff(Int8.(ev_ch .== v))

    ev_start = findall(ev_segs .== 1)
    ev_end   = findall(ev_segs .== -1)

    # an edge at position i lies between samples i and i + 1, so an event detected at
    # edge i starts at sample i + 1, i.e. at i / fs seconds
    # does the signal start with an active event? (edge 0: starts at sample 1, t = 0)
    !isempty(ev_end) && (isempty(ev_start) || ev_end[1] < ev_start[1]) &&
        pushfirst!(ev_start, 0)

    # does the signal end with an active event? (it lasts through the last sample)
    !isempty(ev_start) && (isempty(ev_end) || ev_end[end] < ev_start[end]) &&
        push!(ev_end, length(ev_ch))

    length(ev_start) == length(ev_end) || throw(
        ArgumentError(
            "Mismatched event start/end edges in channel $(labels(obj)[ch]).",
        ),
    )

    ev_len = ev_end .- ev_start

    ch_label = labels(obj)[ch]
    ev_desc = fill(value == "" ? ch_label : value, length(ev_start))
    ev_type = fill(type, length(ev_start))
    # the events apply to all channels
    ev_ch_v = [String[] for _ in eachindex(ev_start)]

    _info("$(length(ev_start)) events found and added as markers")

    _markers_append!(
        obj.markers,
        DataFrame(
            :type    => ev_type,
            :start   => ev_start ./ sr(obj),
            :length  => ev_len ./ sr(obj),
            :value   => ev_desc,
            :channel => ev_ch_v,
        ),
    )
    sort!(obj.markers, :start)

    _history!(obj, "channel2marker", (; ch, v, type, value))

    return nothing
end

"""
    add_markers(obj; <keyword arguments>)

Return a copy of `obj` with its marker table replaced by `markers`.

The provided DataFrame must start with the columns `type`, `start`, `length`, `value`, `channel` (a table with an `id` column instead of `type`, as used by NeuroAnalyzer before v0.26.12, is converted); extra columns (e.g. `trial_type`, `response_time`) may follow and are kept by all marker functions. The `channel` column holds vectors of channel labels (`String[]`: all channels); strings with labels separated by `;` (`""`: all channels, as written by `export_markers()`) and channel numbers (`0`: all channels) are converted.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `markers::DataFrame`: replacement marker table

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function add_markers(obj::NeuroAnalyzer.NEURO; markers::DataFrame)::NeuroAnalyzer.NEURO
    # validate
    _check_marker_cols(markers)

    # create new dataset
    obj_tmp         = deepcopy(obj)
    obj_tmp.markers = _marker_channels(markers, labels(obj))

    _history!(obj_tmp, "add_markers", (; markers))

    return obj_tmp
end

"""
    add_markers!(obj; <keyword arguments>)

Replace `obj`'s marker table in-place.

The provided DataFrame must start with the columns `type`, `start`, `length`, `value`, `channel` (a table with an `id` column instead of `type`, as used by NeuroAnalyzer before v0.26.12, is converted); extra columns (e.g. `trial_type`, `response_time`) may follow and are kept by all marker functions. The `channel` column holds vectors of channel labels (`String[]`: all channels); strings with labels separated by `;` (`""`: all channels, as written by `export_markers()`) and channel numbers (`0`: all channels) are converted.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `markers::DataFrame`: replacement marker table

# Returns

- `Nothing`
"""
function add_markers!(obj::NeuroAnalyzer.NEURO; markers::DataFrame)::Nothing
    _check_marker_cols(markers)
    obj.markers = _marker_channels(markers, labels(obj))
    _history!(obj, "add_markers", (; markers))

    return nothing
end
