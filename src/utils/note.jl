export view_note
export add_note
export add_note!
export delete_note
export delete_note!

"""
    view_note(obj)

Return the recording note stored in the object header.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `String`: recording note (empty string if no note has been set)

# See also

[`add_note!`](@ref), [`delete_note!`](@ref)
"""
function view_note(obj::NeuroAnalyzer.NEURO)::String

    return obj.header.recording[:recording_notes]

end

"""
    add_note(obj; <keyword arguments>)

Return a copy of `obj` with the recording note set to `note`. The original object is not modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `note::String`: note text to store

# Returns

- `NeuroAnalyzer.NEURO`: new NEURO object with the updated recording note

# See also

[`add_note!`](@ref), [`view_note`](@ref), [`delete_note`](@ref)
"""
function add_note(obj::NeuroAnalyzer.NEURO; note::String)::NeuroAnalyzer.NEURO

    obj_new = deepcopy(obj)
    obj_new.header.recording[:recording_notes] = note

    return obj_new

end

"""
    add_note!(obj; <keyword arguments>)

Set the recording note in `obj` in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `note::String`: note text to store

# Returns

- `Nothing`

# See also

[`add_note`](@ref), [`view_note`](@ref), [`delete_note!`](@ref)
"""
function add_note!(obj::NeuroAnalyzer.NEURO; note::String)::Nothing

    obj.header.recording[:recording_notes] = note

    return nothing

end

"""
    delete_note(obj)

Return a copy of `obj` with the recording note cleared. The original object is not modified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `NeuroAnalyzer.NEURO`: new NEURO object with an empty recording note

# See also

[`delete_note!`](@ref), [`view_note`](@ref), [`add_note`](@ref)
"""
function delete_note(obj::NeuroAnalyzer.NEURO)::NeuroAnalyzer.NEURO

    obj_new = deepcopy(obj)
    obj_new.header.recording[:recording_notes] = ""

    return obj_new

end

"""
    delete_note!(obj)

Clear the recording note in `obj` in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Nothing`

# See also

[`delete_note`](@ref), [`view_note`](@ref), [`add_note!`](@ref)
"""
function delete_note!(obj::NeuroAnalyzer.NEURO)::Nothing

    obj.header.recording[:recording_notes] = ""

    return nothing

end
