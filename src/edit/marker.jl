export view_marker
export delete_marker
export delete_marker!
export add_marker
export add_marker!
export edit_marker
export edit_marker!

"""
    view_marker(obj)

Show markers.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function view_marker(obj::NeuroAnalyzer.NEURO)
    obj.header.markers == true || throw(ArgumentError("OBJ has no markers."))
    for marker_idx in 1:size(obj.markers, 1)
        println("ID: $(rpad(("'" * obj.markers[marker_idx, :id] * "'"), 24, " ")) start [sample]: $(rpad(obj.markers[marker_idx, :start], 8, " ")) length [samples]: $(rpad(obj.markers[marker_idx, :length], 8, " ")) description: $(rpad(("'" * obj.markers[marker_idx, :description] * "'"), 24, " ")) channel: $(obj.markers[marker_idx, :channel])")
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
    obj_new = deepcopy(obj)
    obj_new.header.recording[:markers] == true || throw(ArgumentError("OBJ has no markers."))
    nn = size(obj_new.markers, 1)
    (n < 1 || n > nn) && throw(ArgumentError("n has to be ≥ 1 and ≤ $nn."))
    deleteat!(obj_new.markers, n)
    size(obj_new.markers, 1) == 0 && (obj_new.header.recording[:markers] = false)
    reset_components!(obj_new)
    push!(obj_new.header.recording.history, "delete_marker(OBJ; n=$n)")
    
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

    obj_tmp = delete_marker(obj, n=n)
    obj.header.markers = obj_tmp.header.markers
    obj.markers = obj_tmp.markers
    reset_components!(obj)

    return nothing
end

"""
    add_marker(obj; id, start, len, desc, channel)

Add marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `id::String`: marker ID
- `start::Int64`: marker time in samples
- `len::Int64=1`: marker length in samples
- `desc::String`: marker description
- `channel::Int64=0`: channel number, if 0 then marker is related to all channels

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function add_marker(obj::NeuroAnalyzer.NEURO; id::String, start::Int64, len::Int64=1, desc::String, channel::Int64=0)

    start < 1 && throw(ArgumentError("start must be > 0."))
    len < 1 && throw(ArgumentError("len must be > 0."))
    start >= signal_len(obj) && throw(ArgumentError("start must be < $(signal_len(obj) - 1)."))
    start + len > signal_len(obj) && throw(ArgumentError("start + len must be ≤ $(signal_len(obj))."))

    obj_new = deepcopy(obj)
    obj_new.header.recording[:markers] = true
    append!(obj_new.markers, DataFrame(:id=>id, :start=>start, :length=>len, :description=>desc, :channel=>channel))
    sort!(obj_new.markers)
    reset_components!(obj_new)
    push!(obj_new.header.recording.history, "add_marker(OBJ; id=$id, start=$start, len=$len, desc=$desc, channel=$channel)")

    return obj_new
end

"""
    add_marker!(obj; id, start, len, desc, channel)

Add marker.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `id::String`: marker ID
- `start::Int64`: marker time in samples
- `len::Int64=1`: marker length in samples
- `desc::String`: marker description
- `channel::Int64=0`: channel number, if 0 then marker is related to all channels
"""
function add_marker!(obj::NeuroAnalyzer.NEURO; id::String, start::Int64, len::Int64=1, desc::String, channel::Int64=0)

    obj_tmp = add_marker(obj, id=id, start=start, len=len, desc=desc, channel=channel)
    obj.header.markers = obj_tmp.header.markers
    obj.markers = obj_tmp.markers
    reset_components!(obj)

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
- `channel::Int64`: channel number, if 0 then marker is related to all channels

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function edit_marker(obj::NeuroAnalyzer.NEURO; n::Int64, id::String, start::Int64, len::Int64=1, desc::String, channel::Int64)

    obj.header.markers == true || throw(ArgumentError("OBJ has no markers."))
    start < 1 && throw(ArgumentError("start must be > 0."))
    len < 1 && throw(ArgumentError("len must be > 0."))
    start >= signal_len(obj) && throw(ArgumentError("start must be < $(signal_len(obj))."))
    start + len > signal_len(obj) && throw(ArgumentError("start + len must be ≤ $(signal_len(obj))."))

    nn = size(obj.markers, 1)
    n < 1 || n > nn && throw(ArgumentError("n has to be ≥ 1 and ≤ $nn."))
    obj_new = deepcopy(obj)
    obj_new.markers[n, :] = Dict(:id=>id, :start=>start, :length=>len, :description=>desc, :channel=>channel)
     reset_components!(obj_new)
    push!(obj_new.header.recording.history, "edit_marker(OBJ; id=$id, start=$start, len=$len, desc=$desc, channel=$channel)")

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
- `channel::Int64`: channel number, if 0 then marker is related to all channels
"""
function edit_marker!(obj::NeuroAnalyzer.NEURO; n::Int64, id::String, start::Int64, len::Int64=1, desc::String, channel::Int64)

    obj_tmp = edit_marker(obj, n=n, id=id, start=start, len=len, desc=desc, channel=channel)
    obj.header.markers = obj_tmp.header.markers
    obj.markers = obj_tmp.markers
    reset_components!(obj)

    return nothing
end
