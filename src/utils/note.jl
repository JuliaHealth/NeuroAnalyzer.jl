export view_note
export add_note
export add_note!
export delete_note
export delete_note!

"""
    view_note(obj)

Return the object recording note.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `note::String`
"""
function view_note(obj::NeuroAnalyzer.NEURO)::String

    return obj.header.recording[:recording_notes]

end

"""
    add_note(obj; <keyword arguments>)

Add recording note to the object header.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `note::String`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function add_note(obj::NeuroAnalyzer.NEURO; note::String)::NeuroAnalyzer.NEURO

    obj_new = deepcopy(obj)
    obj_new.header.recording[:recording_notes] = note

    return obj_new

end

"""
    add_note!(obj; <keyword arguments>)

Add recording note to the object header.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `note::String`

# Returns

Nothing
"""
function add_note!(obj::NeuroAnalyzer.NEURO; note::String)::Nothing

    obj.header.recording[:recording_notes] = note

    return nothing

end

"""
    delete_note(obj)

Delete recording note from the object header.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function delete_note(obj::NeuroAnalyzer.NEURO)::NeuroAnalyzer.NEURO

    obj_new = deepcopy(obj)
    obj_new.header.recording[:recording_notes] = ""

    return obj_new

end

"""
    delete_note!(obj)

Delete recording note from the object header.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

Nothing
"""
function delete_note!(obj::NeuroAnalyzer.NEURO)::Nothing

    obj.header.recording[:recording_notes] = ""

    return nothing

end
