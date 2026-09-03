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

# Expected column schema for a markers DataFrame
const _MARKER_COLS = ["id", "start", "length", "value", "channel"]

"""
    _check_marker_cols(markers)

Throw `ArgumentError` if `markers` does not have the expected column schema.
"""
function _check_marker_cols(markers::DataFrame)::Nothing
    names(markers) == _MARKER_COLS ||
        throw(
            ArgumentError(
                "Markers DataFrame must have columns: $(_MARKER_COLS); got $(names(markers)).",
            ),
        )
    return nothing
end

"""
    view_marker(obj)

Print a formatted table of all markers.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Nothing`
"""
function view_marker(obj::NeuroAnalyzer.NEURO)::Nothing
    _has_markers(obj) || throw(ArgumentError("OBJ has no markers."))

    println(
        rpad("n", 5) *
        rpad("ID", 24) *
        rpad("start [s]", 12) *
        rpad("length [s]", 12) *
        rpad("value", 24) *
        rpad("channel", 1),
    )
    for i = 1:DataFrames.nrow(obj.markers)
        println(
            rpad(string(i), 5) *
            rpad("'" * obj.markers[i, :id] * "'", 24) *
            rpad(string(round(obj.markers[i, :start]; digits = 3)), 12) *
            rpad(string(round(obj.markers[i, :length]; digits = 3)), 12) *
            rpad("'" * obj.markers[i, :value] * "'", 24) *
            rpad(string(obj.markers[i, :channel]), 1),
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
    _has_markers(obj) || throw(ArgumentError("OBJ has no markers."))

    nn = DataFrames.nrow(obj.markers)
    _in(n, (1, nn), "n")

    # create new dataset
    obj_tmp = deepcopy(obj)

    deleteat!(obj_tmp.markers, n)

    push!(obj_tmp.history, "delete_marker(obj; n=$n)")

    return obj_tmp
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
    obj_tmp = delete_marker(obj; n = n)
    obj.history = obj_tmp.history
    obj.markers = obj_tmp.markers
    obj_tmp = nothing

    return nothing
end

"""
    add_marker(obj; <keyword arguments>)

Append a new marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `id::String`: marker ID
- `start::Real`: marker start time in seconds (must be ≥ 0 and within signal)
- `len::Real=0.0`: marker duration in seconds (must be ≥ 0)
- `value::String`: marker value
- `ch::Int64=0`: channel number; `0` means the marker applies to all channels

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function add_marker(
    obj::NeuroAnalyzer.NEURO;
    id::String,
    start::Real,
    len::Real = 0.0,
    value::String,
    ch::Int64 = 0,
)::NeuroAnalyzer.NEURO
    # validate
    start >= 0 || throw(ArgumentError("start must be ≥ 0."))
    len >= 0 || throw(ArgumentError("len must be ≥ 0."))
    start <= obj.time_pts[end] ||
        throw(ArgumentError("start must be ≤ $(obj.time_pts[end])."))
    start + len <= obj.time_pts[end] ||
        throw(ArgumentError("start + len must be ≤ $(obj.time_pts[end])."))

    # create new dataset
    obj_tmp = deepcopy(obj)
    append!(
        obj_tmp.markers,
        DataFrame(
            :id => id,
            :start => start,
            :length => len,
            :value => value,
            :channel => ch,
        ),
    )
    sort!(obj_tmp.markers, :start)

    push!(
        obj_tmp.history,
        "add_marker(obj; id=$id, start=$start, len=$len, value=$value, ch=$ch)",
    )

    return obj_tmp
end

"""
    add_marker!(obj; <keyword arguments>)

Add a marker in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `id::String`: marker ID
- `start::Real`: marker start time in seconds (must be ≥ 0 and within signal)
- `len::Real=0.0`: marker duration in seconds (must be ≥ 0)
- `value::String`: marker value
- `ch::Int64=0`: channel number; `0` means the marker applies to all channels

# Returns

- `Nothing`
"""
function add_marker!(
    obj::NeuroAnalyzer.NEURO;
    id::String,
    start::Real,
    len::Real = 0.0,
    value::String,
    ch::Int64 = 0,
)::Nothing
    obj_tmp = add_marker(obj; id = id, start = start, len = len, value = value, ch = ch)
    obj.history = obj_tmp.history
    obj.markers = obj_tmp.markers
    obj_tmp = nothing

    return nothing
end

"""
    edit_marker(obj; <keyword arguments>)

Edit a marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `n::Int64`: marker number to edit
- `id::String`: marker ID
- `start::Real`: marker start time in seconds (must be ≥ 0 and within signal)
- `len::Real=0.0`: marker duration in seconds (must be ≥ 0)
- `value::String`: marker value
- `ch::Int64=0`: channel number; `0` means the marker applies to all channels

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function edit_marker(
    obj::NeuroAnalyzer.NEURO;
    n::Int64,
    id::String,
    start::Real,
    len::Real = 0.0,
    value::String,
    ch::Int64 = 0,
)::NeuroAnalyzer.NEURO
    _has_markers(obj) || throw(ArgumentError("OBJ has no markers."))
    start >= 0 || throw(ArgumentError("start must be > 0."))
    len >= 0 || throw(ArgumentError("len must be ≥ 0."))
    start < signal_len(obj) / sr(obj) ||
        throw(ArgumentError("start must be < $(signal_len(obj) / sr(obj))."))
    start + len <= signal_len(obj) / sr(obj) ||
        throw(ArgumentError("start + len must be ≤ $(signal_len(obj) / sr(obj))."))

    nn = size(obj.markers, 1)
    n < 1 || n > nn && throw(ArgumentError("n must be in [1, $nn]."))

    # create new dataset
    obj_tmp = deepcopy(obj)

    obj_tmp.markers[n, :] = Dict(
        :id => id, :start => start, :length => len, :value => value, :channel => ch,
    )
    sort!(obj_tmp.markers, :start)
    push!(
        obj_tmp.history,
        "edit_marker(obj, id=$id, start=$start, len=$len, value=$value, ch=$ch)",
    )

    return obj_tmp
end

"""
    edit_marker!(obj; <keyword arguments>)

Edit a marker in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `n::Int64`: marker number to edit
- `id::String`: marker ID
- `start::Real`: marker start time in seconds (must be ≥ 0 and within signal)
- `len::Real=0.0`: marker duration in seconds (must be ≥ 0)
- `value::String`: marker value
- `ch::Int64=0`: channel number; `0` means the marker applies to all channels

# Returns

- `Nothing`
"""
function edit_marker!(
    obj::NeuroAnalyzer.NEURO;
    n::Int64,
    id::String,
    start::Real,
    len::Real = 0.0,
    value::String,
    ch::Int64 = 0,
)::Nothing
    obj_tmp = edit_marker(
        obj; n = n, id = id, start = start, len = len, value = value, ch = ch,
    )
    obj.history = obj_tmp.history
    obj.markers = obj_tmp.markers
    obj_tmp = nothing

    return nothing
end

"""
    channel2marker(obj; <keyword arguments>)

Return a copy of `obj` with events detected in an event channel added as markers.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `v::Real=1.0`: signal value interpreted as an event
- `id::String`: prefix for generated marker IDs (default: channel name + `"_"`)
- `value::String=""`: marker value string (default: channel label)

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function channel2marker(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    v::Real = 1.0,
    id::String = "",
    value::String = "",
)::NeuroAnalyzer.NEURO

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    length(ch) == 1 || throw(ArgumentError("ch must resolve to exactly one channel."))
    ch = ch[1]

    # resolve markers channel
    stim_ch = get_channel(obj; type = "mrk")
    isempty(stim_ch) && throw(ArgumentError("No MRK channel."))

    # check if the event channel contain events
    ev_ch = obj.data[ch, :, :][:]
    length(unique(ev_ch)) > 1 ||
        throw(ArgumentError("Channel $ch does not contain events."))

    # extract events
    ev_v = unique(ev_ch)
    v in ev_v || throw(ArgumentError("Event channel does not contain value $v."))
    _info("Event channel contains values: $ev_v")

    ev_segs = diff(ev_ch)

    ev_start = findall(ev_segs .== v)
    ev_end   = findall(ev_segs .== -v)

    # does the signal start with an active event?
    !isempty(ev_end) && !isempty(ev_start) && ev_end[1] < ev_start[1] &&
        pushfirst!(ev_start, 1)

    # does the signal end with an active event?
    !isempty(ev_start) && !isempty(ev_end) && ev_end[end] < ev_start[end] &&
        push!(ev_end, length(ev_ch))

    length(ev_start) == length(ev_end) || throw(
        ArgumentError(
            "Mismatched event start/end edges in channel $(labels(obj)[ch]).",
        ),
    )

    ev_len = ev_end .- ev_start

    ch_label = labels(obj)[ch]
    ev_desc = fill(value == "" ? ch_label : value, length(ev_start))
    id_prefix = id == "" ? ch_label * "_" : id
    ev_id = ["$id_prefix$i" for i in eachindex(ev_start)]
    ev_ch_v = zeros(Int64, length(ev_start))

    _info("$(length(ev_start)) events found and added as markers.")

    obj_tmp = deepcopy(obj)
    append!(
        obj_tmp.markers,
        DataFrame(
            :id      => ev_id,
            :start   => ev_start ./ sr(obj),
            :length  => ev_len ./ sr(obj),
            :value   => ev_desc,
            :channel => ev_ch_v,
        ),
    )
    sort!(obj_tmp.markers, :start)

    push!(obj_tmp.history, "channel2marker(obj; ch=$ch, v=$v, id=$id, value=$value)")

    return obj_tmp
end

"""
    channel2marker!(obj; <keyword arguments>)

Convert an event channel to markers in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `v::Real=1.0`: signal value interpreted as an event
- `id::String`: prefix for generated marker IDs (default: channel name + `"_"`)
- `value::String=""`: marker value string (default: channel label)

# Returns

- `Nothing`
"""
function channel2marker!(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    v::Real = 1.0,
    id::String = "",
    value::String = "",
)::Nothing
    obj_tmp = channel2marker(obj; ch = ch, v = v, id = id, value = value)
    obj.history = obj_tmp.history
    obj.markers = obj_tmp.markers
    obj_tmp = nothing

    return nothing
end

"""
    add_markers(obj; <keyword arguments>)

Return a copy of `obj` with its marker table replaced by `markers`.

The provided DataFrame must have exactly the columns: `id`, `start`, `length`, `value`, `channel`.

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
    obj_tmp.markers = markers

    push!(obj_tmp.history, "add_markers(obj; markers)")

    return obj_tmp
end

"""
    add_markers!(obj; <keyword arguments>)

Replace `obj`'s marker table in-place.

The provided DataFrame must have exactly the columns: `id`, `start`, `length`, `value`, `channel`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `markers::DataFrame`: replacement marker table

# Returns

- `Nothing`
"""
function add_markers!(obj::NeuroAnalyzer.NEURO; markers::DataFrame)::Nothing
    _check_marker_cols(markers)
    obj.markers = markers
    push!(obj.history, "add_markers!(obj; markers)")

    return nothing
end
