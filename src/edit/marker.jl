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

    _has_markers(obj) == true || throw(ArgumentError("OBJ has no markers."))
    
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

    _has_markers(obj) == true || throw(ArgumentError("OBJ has no markers."))

    obj_new = deepcopy(obj)
    nn = size(obj_new.markers, 1)
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
    sort!(obj_new.markers)
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
