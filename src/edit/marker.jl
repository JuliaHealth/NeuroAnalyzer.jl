export view_marker
export delete_marker
export delete_marker!
export add_marker
export add_marker!
export edit_marker
export edit_marker!
export channel2marker
export channel2marker!

"""
    view_marker(obj)

Show markers.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function view_marker(obj::NeuroAnalyzer.NEURO)

    _has_markers(obj) == true || throw(ArgumentError("OBJ has no markers."))
    
    for marker_idx in 1:nrow(obj.markers)
        println("n: $(rpad(string(marker_idx), 6)) ID: $(rpad(("'" * obj.markers[marker_idx, :id] * "'"), 24, " ")) start [sample]: $(rpad(obj.markers[marker_idx, :start], 8, " ")) length [samples]: $(rpad(obj.markers[marker_idx, :length], 8, " ")) description: $(rpad(("'" * obj.markers[marker_idx, :description] * "'"), 24, " ")) channel: $(obj.markers[marker_idx, :channel])")
    end

end

"""
    delete_marker(obj; n)

Delete marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `n::Int64`: marker number

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function delete_marker(obj::NeuroAnalyzer.NEURO; n::Int64)

    _has_markers(obj) == true || throw(ArgumentError("OBJ has no markers."))

    obj_new = deepcopy(obj)
    nn = nrow(obj_new.markers)
    (n < 1 || n > nn) && throw(ArgumentError("n has to be ≥ 1 and ≤ $nn."))
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
"""
function delete_marker!(obj::NeuroAnalyzer.NEURO; n::Int64)

    obj_new = delete_marker(obj, n=n)
    obj.history = obj_new.history
    obj.markers = obj_new.markers

    return nothing

end

"""
    add_marker(obj; id, start, len, desc, ch)

Add marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `id::String`: marker ID
- `start::Int64`: marker time in samples
- `len::Int64=1`: marker length in samples
- `desc::String`: marker description
- `ch::Int64=0`: channel number, if 0 then marker is related to all channels

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function add_marker(obj::NeuroAnalyzer.NEURO; id::String, start::Int64, len::Int64=1, desc::String, ch::Int64=0)

    start < 1 && throw(ArgumentError("start must be > 0."))
    len < 1 && throw(ArgumentError("len must be > 0."))
    start >= signal_len(obj) && throw(ArgumentError("start must be < $(signal_len(obj) - 1)."))
    start + len > signal_len(obj) && throw(ArgumentError("start + len must be ≤ $(signal_len(obj))."))

    obj_new = deepcopy(obj)
    append!(obj_new.markers, DataFrame(:id=>id, :start=>start, :length=>len, :description=>desc, :channel=>ch))
    sort!(obj_new.markers, :start)
    reset_components!(obj_new)
    push!(obj_new.history, "add_marker(OBJ; id=$id, start=$start, len=$len, desc=$desc, ch=$ch)")

    return obj_new

end

"""
    add_marker!(obj; id, start, len, desc, ch)

Add marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `id::String`: marker ID
- `start::Int64`: marker time in samples
- `len::Int64=1`: marker length in samples
- `desc::String`: marker description
- `ch::Int64=0`: channel number, if 0 then marker is related to all channels
"""
function add_marker!(obj::NeuroAnalyzer.NEURO; id::String, start::Int64, len::Int64=1, desc::String, ch::Int64=0)

    obj_new = add_marker(obj, id=id, start=start, len=len, desc=desc, ch=ch)
    obj.history = obj_new.history
    obj.markers = obj_new.markers

    return nothing

end

"""
    edit_marker(obj; n, id, start, len, desc)

Edit marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `n::Int64`: marker number
- `id::String`: marker ID
- `start::Int64`: marker time in samples
- `len::Int64=1`: marker length in samples
- `desc::String`: marker description
- `ch::Int64`: channel number, if 0 then marker is related to all channels

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function edit_marker(obj::NeuroAnalyzer.NEURO; n::Int64, id::String, start::Int64, len::Int64=1, desc::String, ch::Int64)

    _has_markers(obj) == true || throw(ArgumentError("OBJ has no markers."))
    start < 1 && throw(ArgumentError("start must be > 0."))
    len < 1 && throw(ArgumentError("len must be > 0."))
    start >= signal_len(obj) && throw(ArgumentError("start must be < $(signal_len(obj))."))
    start + len > signal_len(obj) && throw(ArgumentError("start + len must be ≤ $(signal_len(obj))."))

    nn = size(obj.markers, 1)
    n < 1 || n > nn && throw(ArgumentError("n has to be ≥ 1 and ≤ $nn."))
    obj_new = deepcopy(obj)
    obj_new.markers[n, :] = Dict(:id=>id, :start=>start, :length=>len, :description=>desc, :channel=>ch)
     reset_components!(obj_new)
    push!(obj_new.history, "edit_marker(OBJ; id=$id, start=$start, len=$len, desc=$desc, ch=$ch)")

    return obj_new

end

"""
    edit_marker!(obj; n, id, start, len, desc)

Edit marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `n::Int64`: marker number
- `id::String`: marker ID
- `start::Int64`: marker time in samples
- `len::Int64=1`: marker length in samples
- `desc::String`: marker description
- `ch::Int64`: channel number, if 0 then marker is related to all channels
"""
function edit_marker!(obj::NeuroAnalyzer.NEURO; n::Int64, id::String, start::Int64, len::Int64=1, desc::String, ch::Int64)

    obj_new = edit_marker(obj, n=n, id=id, start=start, len=len, desc=desc, ch=ch)
    obj.history = obj_new.history
    obj.markers = obj_new.markers

    return nothing

end

"""
    channel2marker(obj; ch, v, id, desc)

Convert event channel to markers.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Int64`: event channel number
- `v::Real=1.0`: event channel value interpreted as an event
- `id::String`: prefix for marker ID; default is "mrk_"
- `desc::String=""`: prefix for marker description; default is based on event channel name (e.g. "stim1_")

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function channel2marker(obj::NeuroAnalyzer.NEURO; ch::Int64, v::Real=1.0, id::String="", desc::String="")

    stim_ch = get_channel_bytype(obj, type=:mrk)
    _check_channels(stim_ch, ch)

    # check if the event channel contain events
    ev_ch = obj.data[ch, :, :][:]
    length(unique(ev_ch)) == 1 && throw(ArgumentError("Channel $ch does not contain events."))

    # extract events
    ev_v = unique(ev_ch)
    v in ev_v || throw(ArgumentError("Event channel does not contain value $v."))
    _info("Event channel contain values: $ev_v")

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
    desc == "" && (desc = labels(obj)[ch] * "_")
    id == "" && (id = "mrk_")
    ev_id = String[]
    ev_desc = String[]
    ev_ch = zeros(Int64, length(ev_start))
    for idx in 1:length(ev_start)
        push!(ev_id, "$id$idx")
        push!(ev_desc, "$desc$idx")
    end
    
    _info("$(length(ev_start)) events added.")

    obj_new = deepcopy(obj)
    append!(obj_new.markers, DataFrame(:id=>ev_id, :start=>ev_start, :length=>ev_len, :description=>ev_desc, :channel=>ev_ch))
    sort!(obj_new.markers, :start)
    reset_components!(obj_new)
    push!(obj_new.history, "channel2marker(OBJ, ch=$ch, v=$v, id=$id, desc=$desc")

    return obj_new

end

"""
    channel2marker!(obj; ch, v, id, desc)

Convert event channel to markers.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Int64`: event channel number
- `v::Real=1.0`: event channel value interpreted as an event
- `id::String`: prefix for marker ID; default is "mrk_"
- `desc::String=""`: prefix for marker description; default is based on event channel name (e.g. "stim1_")
"""
function channel2marker!(obj::NeuroAnalyzer.NEURO; ch::Int64, v::Real=1.0, id::String="", desc::String="")

    obj_new = channel2marker(obj, ch=ch, v=v, id=id, desc=desc)
    obj.history = obj_new.history
    obj.markers = obj_new.markers

    return nothing

end
