export view_note
export add_note
export add_note!
export delete_note
export delete_note!

"""
    view_note(obj)

Return recording note.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function view_note(obj::NeuroAnalyzer.NEURO)
    return obj.header.recording[:recording_notes]
end

"""
    add_note(obj; note)

Add recording note.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `note::String`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function add_note(obj::NeuroAnalyzer.NEURO; note::String)
    obj_new = deepcopy(obj)
    obj_new.header.recording[:recording_notes] = note
    return obj_new
end

"""
    add_note!(obj; note)

Add recording note.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `note::String`
"""
function add_note!(obj::NeuroAnalyzer.NEURO; note::String)
    obj.header.recording[:recording_notes] = note
    return nothing
end

"""
    delete_note(obj)

Delete recording note.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function delete_note(obj::NeuroAnalyzer.NEURO)
    obj_new = deepcopy(obj)
    obj_new.header.recording[:recording_notes] = ""
    return obj_new
end

"""
    delete_note!(obj)

Delete recording note.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function delete_note!(obj::NeuroAnalyzer.NEURO)
    obj.header.recording[:recording_notes] = ""
    return nothing
end
