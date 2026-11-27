export view_note
export add_note
export add_note!
export delete_note
export delete_note!

"""
    view_note(obj)

Return the header notes (`obj.header.notes`) stored in the object header.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `String`: header notes (empty string if no note has been set)
"""
function view_note(obj::NeuroAnalyzer.NEURO)::String
    return obj.header.notes
end

"""
    add_note(obj; <keyword arguments>)

Return a copy of `obj` with the header notes set to `note`. The original object is not modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `note::String`: note text to store

# Returns

- `NeuroAnalyzer.NEURO`: new NEURO object with the updated header notes
"""
function add_note(obj::NeuroAnalyzer.NEURO; note::String)::NeuroAnalyzer.NEURO
    # create new dataset
    obj_tmp = deepcopy(obj)
    add_note!(obj_tmp; note = note)

    return obj_tmp
end

"""
    add_note!(obj; <keyword arguments>)

Set the header notes in `obj` in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `note::String`: note text to store

# Returns

- `Nothing`
"""
function add_note!(obj::NeuroAnalyzer.NEURO; note::String)::Nothing
    obj.header.notes = note
    _history!(obj, "add_note", (; note))

    return nothing
end

"""
    delete_note(obj)

Return a copy of `obj` with the header notes cleared. The original object is not modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `NeuroAnalyzer.NEURO`: new NEURO object with an empty header notes
"""
function delete_note(obj::NeuroAnalyzer.NEURO)::NeuroAnalyzer.NEURO
    # create new dataset
    obj_tmp = deepcopy(obj)
    delete_note!(obj_tmp)

    return obj_tmp
end

"""
    delete_note!(obj)

Clear the header notes in `obj` in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Nothing`
"""
function delete_note!(obj::NeuroAnalyzer.NEURO)::Nothing
    obj.header.notes = ""
    _history!(obj, "delete_note", NamedTuple())

    return nothing
end
