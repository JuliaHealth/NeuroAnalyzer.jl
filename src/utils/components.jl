export component_type
export add_component
export add_component!
export list_component
export extract_component
export delete_component
export delete_component!
export reset_components
export reset_components!
export rename_component
export rename_component!

"""
    component_type(obj, c)

Return component data type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `c::Symbol`: component name

# Return

- `c_type::DataType`
"""
function component_type(obj::NeuroAnalyzer.NEURO; c::Symbol)

    c in keys(obj.components) || throw(ArgumentError("Component :$c does not exist. Use list_component() to view existing components."))

    return typeof(obj.components[c])

end

"""
    rename_component(obj, c_old, c_new)

Rename component.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `c_old::Symbol`: old component name
- `c_new::Symbol`: new component name

# Return

- `obj_new::NEURO`
"""
function rename_component(obj::NeuroAnalyzer.NEURO; c_old::Symbol, c_new::Symbol)

    c_old in keys(obj.components) || throw(ArgumentError("Component $c_old does not exist. Use list_component() to view existing components."))
    c_new in keys(obj.components) && throw(ArgumentError("Component $c_new already exists. Use list_component() to view existing components."))

    obj_new = deepcopy(obj)
    c = pop!(obj.components, c_old)
    push!(obj_new.components, c_new=>c)

    push!(obj_new.header.history, "rename_component(OBJ, c_old=$c_old, c_new=$c_new)")

    return obj_new
    
end

"""
    rename_component!(obj, c_old, c_new)

Rename component.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `c_old::Symbol`: old component name
- `c_new::Symbol`: new component name
"""
function rename_component!(obj::NeuroAnalyzer.NEURO; c_old::Symbol, c_new::Symbol)

    obj_new = rename_component(obj, c_old=c_old, c_new=c_new)
    obj.components = obj_new.components

    return nothing

end

"""
    add_component(obj; c, v)

Add component.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `c::Symbol`: component name
- `v::Any`: component value

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function add_component(obj::NeuroAnalyzer.NEURO; c::Symbol, v::Any)

    obj_new = deepcopy(obj)
    c in keys(obj.components) && throw(ArgumentError("Component $c already exists. Use different component name or delete_component() to remove it prior the operation."))

    # add component
    push!(obj_new.components, c=>v)
    
    # update history
    push!(obj_new.header.history, "add_component(OBJ, c=$c, v=$v)")

    return obj_new

end

"""
    add_component!(obj; c, v)

Add component.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `c::Symbol`: component name
- `v::Any`: component value
"""
function add_component!(obj::NeuroAnalyzer.NEURO; c::Symbol, v::Any)

    obj_new = add_component(obj, c=c, v=v)
    obj.components = obj_new.components

    return nothing

end

"""
    list_component(obj)

List component names.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `components::Vector{Symbol}`
"""
function list_component(obj::NeuroAnalyzer.NEURO)

    return keys(obj.components)

end

"""
    extract_component(obj, c)

Extract component values.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `c::Symbol`: component name

# Returns

- `c::Any`
"""
function extract_component(obj::NeuroAnalyzer.NEURO; c::Symbol)

    c in keys(obj.components) || throw(ArgumentError("Component $c does not exist. Use list_component() to view existing components."))
    c = obj.components[c]

    return c

end

"""
    delete_component(obj; c)

Delete component. 

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `c::Symbol`: component name

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function delete_component(obj::NeuroAnalyzer.NEURO; c::Symbol)

    c in keys(obj.components) || throw(ArgumentError("Component $c does not exist. Use list_component() to view existing components."))
    
    # delete component values
    obj_new = deepcopy(obj)
    pop!(obj_new.components, c)

    push!(obj_new.header.history, "delete_component(OBJ, c=$c)")
    
    return obj_new

end

"""
    delete_component!(obj; c)

Delete component.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `c::Symbol`: component name
"""
function delete_component!(obj::NeuroAnalyzer.NEURO; c::Symbol)

    obj_new = delete_component(obj, c=c)
    obj.components = obj_new.components

    return nothing

end

"""
    reset_components(obj)

Remove all components.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function reset_components(obj::NeuroAnalyzer.NEURO)

    obj_new = deepcopy(obj)
    obj_new.components = Dict()

    return obj_new

end

"""
    reset_components!(obj)

Remove all components.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function reset_components!(obj::NeuroAnalyzer.NEURO)

    obj.components = Dict()

    return nothing

end
