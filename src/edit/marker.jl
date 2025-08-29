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
    view_marker(obj)

Show markers.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

Nothing
"""
function view_marker(obj::NeuroAnalyzer.NEURO)::Nothing

    @assert _has_markers(obj) "OBJ has no markers."

    println(rpad("n", 5) *
            rpad("ID", 24) *
            rpad("start [s]", 12) *
            rpad("length [s]", 12) *
            rpad("value", 24) *
            rpad("channel", 1))

    for mrk_idx in 1:DataFrames.nrow(obj.markers)
        println(rpad(string(mrk_idx), 5) *
                rpad("'" * obj.markers[mrk_idx, :id] * "'", 24) *
                rpad(string(round(obj.markers[mrk_idx, :start], digits=3)), 12) *
                rpad(string(round(obj.markers[mrk_idx, :length], digits=3)), 12) *
                rpad("'" * obj.markers[mrk_idx, :value] * "'", 24) *
                rpad(string(obj.markers[mrk_idx, :channel]), 1))
    end

end

"""
    delete_marker(obj; <keyword arguments>)

Delete marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `n::Int64`: marker number

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function delete_marker(obj::NeuroAnalyzer.NEURO; n::Int64)::NeuroAnalyzer.NEURO

    @assert _has_markers(obj) "OBJ has no markers."

    obj_new = deepcopy(obj)
    nn = DataFrames.nrow(obj_new.markers)
    @assert !(n < 1 || n > nn) "n must be in [1, $nn]."
    deleteat!(obj_new.markers, n)
    reset_components!(obj_new)
    push!(obj_new.history, "delete_marker(OBJ; n=$n)")

    return obj_new

end

"""
    delete_marker!(obj; n)

Delete marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `n::Int64`: marker number

# Returns

Nothing
"""
function delete_marker!(obj::NeuroAnalyzer.NEURO; n::Int64)::Nothing

    obj_new = delete_marker(obj, n=n)
    obj.history = obj_new.history
    obj.markers = obj_new.markers

    return nothing

end

"""
    add_marker(obj; <keyword arguments>)

Add marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `id::String`: marker ID
- `start::Real`: marker time in seconds
- `len::Real=1.0`: marker length in seconds
- `value::String`: marker value
- `ch::Int64=0`: channel number, if 0 then marker is related to all channels

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function add_marker(obj::NeuroAnalyzer.NEURO; id::String, start::Real, len::Real=1.0, value::String, ch::Int64=0)::NeuroAnalyzer.NEURO

    @assert start > 0 "start must be > 0."
    @assert len > 0 "len must be > 0."
    @assert start < signal_len(obj) "start must be < $(signal_len(obj) - 1)."
    @assert start + len <= signal_len(obj) "start + len must be ≤ $(signal_len(obj))."

    obj_new = deepcopy(obj)
    append!(obj_new.markers, DataFrame(:id=>id, :start=>start, :length=>len, :value=>value, :channel=>ch))
    sort!(obj_new.markers, :start)
    reset_components!(obj_new)
    push!(obj_new.history, "add_marker(OBJ; id=$id, start=$start, len=$len, value=$value, ch=$ch)")

    return obj_new

end

"""
    add_marker!(obj; id, start, len, value, ch)

Add marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `id::String`: marker ID
- `start::Real`: marker time in seconds
- `len::Real=1.0`: marker length in seconds
- `value::String`: marker value
- `ch::Int64=0`: channel number, if 0 then marker is related to all channels

# Returns

Nothing
"""
function add_marker!(obj::NeuroAnalyzer.NEURO; id::String, start::Real, len::Real=1.0, value::String, ch::Int64=0)::Nothing

    obj_new = add_marker(obj, id=id, start=start, len=len, value=value, ch=ch)
    obj.history = obj_new.history
    obj.markers = obj_new.markers

    return nothing

end

"""
    edit_marker(obj; <keyword arguments>)

Edit marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `n::Int64`: marker number
- `id::String`: marker ID
- `start::Real`: marker time in seconds
- `len::Real=1.0`: marker length in seconds
- `value::String`: marker value
- `ch::Int64=0`: channel number, if 0 then marker is related to all channels

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function edit_marker(obj::NeuroAnalyzer.NEURO; n::Int64, id::String, start::Real, len::Real=1.0, value::String, ch::Int64=0)::NeuroAnalyzer.NEURO

    @assert _has_markers(obj) "OBJ has no markers."
    @assert start > 0 "start must be > 0."
    @assert len > 0 "len must be > 0."
    @assert start < signal_len(obj) / sr(obj) "start must be < $(signal_len(obj) / sr(obj))."
    @assert start + len <= signal_len(obj) / sr(obj) "start + len must be ≤ $(signal_len(obj) / sr(obj))."

    nn = size(obj.markers, 1)
    @assert !(n < 1 || n > nn) "n must be in [1, $nn]."
    obj_new = deepcopy(obj)
    obj_new.markers[n, :] = Dict(:id=>id, :start=>start, :length=>len, :value=>value, :channel=>ch)
     reset_components!(obj_new)
    sort!(obj_new.markers, :start)
    push!(obj_new.history, "edit_marker(OBJ; id=$id, start=$start, len=$len, value=$value, ch=$ch)")

    return obj_new

end

"""
    edit_marker!(obj; n, id, start, len, value, ch)

Edit marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `n::Int64`: marker number
- `id::String`: marker ID
- `start::Real`: marker time in seconds
- `len::Real=1`: marker length in seconds
- `value::String`: marker value
- `ch::Int64=0`: channel number, if 0 then marker is related to all channels

# Returns

Nothing
"""
function edit_marker!(obj::NeuroAnalyzer.NEURO; n::Int64, id::String, start::Real, len::Real=1.0, value::String, ch::Int64=0)::Nothing

    obj_new = edit_marker(obj, n=n, id=id, start=start, len=len, value=value, ch=ch)
    obj.history = obj_new.history
    obj.markers = obj_new.markers

    return nothing

end

"""
    channel2marker(obj; <keyword arguments>)

Convert event channel to markers.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `v::Real=1.0`: event channel value interpreted as an event
- `id::String`: prefix for marker ID; default is based on event channel name (e.g. "stim1_")
- `value::String=""`: marker value; default is based on event channel name (e.g. "stim1")

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function channel2marker(obj::NeuroAnalyzer.NEURO; ch::String, v::Real=1.0, id::String="", value::String="")::NeuroAnalyzer.NEURO

    stim_ch = get_channel(obj, type="mrk")
    ch = get_channel(obj, ch=ch)[1]

    # check if the event channel contain events
    ev_ch = obj.data[ch, :, :][:]
    @assert length(unique(ev_ch)) > 1 "Channel $ch does not contain events."

    # extract events
    ev_v = unique(ev_ch)
    @assert v in ev_v "Event channel does not contain value $v."
    _info("Event channel contains values: $ev_v")

    ev_start = getindex.(findall(ev_ch .== v), 1)

    ev_segs = diff(ev_ch)
    ev_segments = getindex.(findall(abs.(ev_segs) .== v), 1)
    ev_segments[2:end] .-= 1

    # segment starts
    ev_start = getindex.(findall(ev_segs .== v))

    # segment ends
    ev_end = getindex.(findall(ev_segs .== -v))

    # does signal start with a marker?
    ev_end[1] < ev_start[1] && pushfirst!(ev_start, 1)

    # does signal end with a marker?
    ev_end[end] < ev_start[end] && push!(ev_end, length(ev_ch))

    # segment lengths
    ev_len = ev_end .- ev_start

    # generate descriptors and IDs
    ev_desc = String[]
    if value == ""
        ev_desc = repeat([labels(obj)[ch]], length(ev_start))
    else
        ev_desc = repeat([value], length(ev_start))
    end

    ev_id = String[]
    id == "" && (id = labels(obj)[ch] * "_")
    for idx in eachindex(ev_start)
        push!(ev_id, "$id$idx")
    end

    ev_ch = zeros(Int64, length(ev_start))

    _info("$(length(ev_start)) events added")

    obj_new = deepcopy(obj)
    append!(obj_new.markers, DataFrame(:id=>ev_id, :start=>(ev_start ./ sr(obj)), :length=>(ev_len ./ sr(obj)), :value=>ev_desc, :channel=>ev_ch))
    sort!(obj_new.markers, :start)
    reset_components!(obj_new)
    push!(obj_new.history, "channel2marker(OBJ, ch=$ch, v=$v, id=$id, value=$value")

    return obj_new

end

"""
    channel2marker!(obj; <keyword arguments>)

Convert event channel to markers.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `v::Real=1.0`: event channel value interpreted as an event
- `id::String`: prefix for marker ID; default is "mrk_"
- `value::String=""`: prefix for marker value; default is based on event channel name (e.g. "stim1_")

# Returns

Nothing
"""
function channel2marker!(obj::NeuroAnalyzer.NEURO; ch::String, v::Real=1.0, id::String="", value::String="")::Nothing

    obj_new = channel2marker(obj, ch=ch, v=v, id=id, value=value)
    obj.history = obj_new.history
    obj.markers = obj_new.markers

    return nothing

end

"""
    add_markers(obj; markers)

Add markers.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `markers::DataFrame`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function add_markers(obj::NeuroAnalyzer.NEURO; markers::DataFrame)::NeuroAnalyzer.NEURO

    obj_new = deepcopy(obj)
    @assert names(markers) == ["id", "start", "length", "value", "channel"] "Markers column names are incorrect."
    obj_new.markers = markers

    return nothing

end

"""
    add_markers!(obj; markers)

Add markers.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `markers::DataFrame`

# Returns

Nothing
"""
function add_markers!(obj::NeuroAnalyzer.NEURO; markers::DataFrame)::Nothing

    @assert names(markers) == ["id", "start", "length", "value", "channel"] "Markers column names are incorrect."
    obj.markers = markers

    return nothing

end