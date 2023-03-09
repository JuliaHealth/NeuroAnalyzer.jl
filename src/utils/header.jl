export edit_header
export edit_header!
export show_header

"""
    edit_header(obj; field, value)

Change value of OBJ header.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `field::Symbol`
- `value::Any`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function edit_header(obj::NeuroAnalyzer.NEURO; field::Symbol, value::Any)

    value === nothing && throw(ArgumentError("value cannot be empty."))

    obj_new = deepcopy(obj)
    field in keys(obj_new.header.recording) || throw(ArgumentError("$field does not exist."))
    typeof(obj_new.header.recording[field]) == typeof(value) || throw(ArgumentError("field type ($(typeof(obj_new.header.recording[field]))) does not mach value type ($(typeof(value)))."))
    obj_new.header.recording[field] = value
    push!(obj_new.header.history, "edit(OBJ, field=$field, value=$value)")    

    return obj_new
end

"""
    edit_header!(obj; field, value)

Change value of OBJ header.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `field::Symbol`
- `value::Any`
"""
function edit_header!(obj::NeuroAnalyzer.NEURO; field::Symbol, value::Any)
    obj.header.recording = edit_header(obj, field=field, value=value).header.recording
    return nothing
end

"""
    show_header(obj)

Show keys and values of OBJ header.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function show_header(obj::NeuroAnalyzer.NEURO)
    for (key, value) in obj.header
        println("$key: $value")
    end
end
