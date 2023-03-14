export component_idx
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
    component_idx(obj, c)

Return component index.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `c::Symbol`: component name

# Return

- `c_idx::Int64`
"""
function component_idx(obj::NeuroAnalyzer.NEURO; c::Symbol)

    c in obj.header.component_names || throw(ArgumentError("Component $c does not exist. Use list_component() to view existing components."))
    
    return findfirst(isequal(c), obj.header.component_names)

end

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

    c in obj.header.component_names || throw(ArgumentError("Component $c does not exist. Use list_component() to view existing components."))
    c_idx = component_idx(obj; c=c)

    return typeof(obj.components[c_idx])

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

    c_old in obj.header.component_names || throw(ArgumentError("Component $c_old does not exist. Use list_component() to view existing components."))
    c_new in obj.header.component_names && throw(ArgumentError("Component $c_new already exists. Use list_component() to view existing components."))

    obj_new = deepcopy(obj)
    c_idx = component_idx(obj, c=c_old)
    obj_new.header.component_names[c_idx] = c_new

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
    obj.header.component_names = obj_new.header.component_names
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
    c in obj_new.header.component_names && throw(ArgumentError("Component $c already exists. Use delete_component() to remove it prior the operation."))

    # add component name
    push!(obj_new.header.component_names, c)
    
    # add component values
    push!(obj_new.components, v)
    
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
    obj.header.component_names = obj_new.header.component_names
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

    return obj.header.component_names

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

    c in obj.header.component_names || throw(ArgumentError("Component $c does not exist. Use list_component() to view existing components."))
    c_idx = component_idx(obj, c=c)
    
    c = obj.components[c_idx]

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

    c in obj.header.component_names || throw(ArgumentError("Component $c does not exist. Use list_component() to view existing components."))
    
    obj_new = deepcopy(obj)
    c_idx = component_idx(obj, c=c)

    # delete component values
    deleteat!(obj_new.components, c_idx)

    # delete component name
    deleteat!(obj_new.header.component_names, c_idx)

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
    obj.header.component_names = obj_new.header.component_names
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
    obj_new.header.component_names = Symbol[]
    obj_new.components = Any[]

    return obj_new

end

"""
    reset_components!(obj)

Remove all components.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function reset_components!(obj::NeuroAnalyzer.NEURO)

    obj.header.component_names = Symbol[]
    obj.components = Any[]

    return nothing

end
