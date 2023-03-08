export component_idx
export component_type
export add_component
export add_component!
export list_component
export extract_component
export delete_component
export delete_component!
export reset_component
export reset_component!
export rename_component
export rename_component!

function _get_component(obj::NeuroAnalyzer.NEURO, c::Symbol)
    # get component values
    c in obj.header.components || throw(ArgumentError("Component $c not found."))
    c_idx = findfirst(isequal(c), obj.header.components)
    c = obj.components[c_idx]
    return (c=c, c_idx=c_idx)
end

function _select_cidx(obj::NeuroAnalyzer.NEURO, c::Symbol, c_idx::Union{Int64, Vector{Int64}, AbstractRange}, def_cidx::Int64=0)
    # select component channels, default is all or def_cidx
    c, _ = _get_component(obj, c)
    def_cidx > size(c, 1) && (def_cidx = size(c, 1))
    def_cidx == 0 && (def_cidx = size(c, 1))
    c_idx == 0 && (c_idx = 1:def_cidx)
    typeof(c_idx) <: AbstractRange && (c_idx = collect(c_idx))
    length(c_idx) > 1 && sort!(c_idx)
    for idx in c_idx
        (idx < 1 || idx > size(c, 1)) && throw(ArgumentError("c_idx must be ≥ 1 and ≤ $(size(c, 1))."))
    end
    return c_idx
end

function _select_cidx(c::AbstractArray, c_idx::Union{Int64, Vector{Int64}, AbstractRange}, def_cidx::Int64=0)
    # select component channels, default is all or def_cidx
    def_cidx > size(c, 1) && (def_cidx = size(c, 1))
    def_cidx == 0 && (def_cidx = size(c, 1))
    c_idx == 0 && (c_idx = 1:def_cidx)
    typeof(c_idx) <: AbstractRange && (c_idx = collect(c_idx))
    length(c_idx) > 1 && sort!(c_idx)
    for idx in c_idx
        (idx < 1 || idx > size(c, 1)) && throw(ArgumentError("c_idx must be ≥ 1 and ≤ $(size(c, 1))."))
    end
    return c_idx
end

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

    c in obj.header.components || throw(ArgumentError("Component $c does not exist. Use list_component() to view existing components."))
    return findfirst(isequal(c), obj.header.components)

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

    c in obj.header.components || throw(ArgumentError("Component $c does not exist. Use list_component() to view existing components."))
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

    c_old in obj.header.components || throw(ArgumentError("Component $c_old does not exist. Use list_component() to view existing components."))
    c_new in obj.header.components && throw(ArgumentError("Component $c_new already exists. Use list_component() to view existing components."))

    obj_new = deepcopy(obj)
    c_idx = component_idx(obj, c=c_old)
    obj_new.header.components[c_idx] = c_new

    push!(obj_new.header.history, "rename_component(NEURO, c_old=$c_old, c_new=$c_new)")

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
    obj.header.components = obj_new.header.components
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
    c in obj_new.header.components && throw(ArgumentError("Component $c already exists. Use delete_component() to remove it prior the operation."))

    # add component name
    push!(obj_new.header.components, c)
    
    # add component values
    push!(obj_new.components, v)
    
    # update history
    push!(obj_new.header.history, "add_component(NEURO, c=$c, v=$v)")

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
    obj.header.components = obj_new.header.components
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
    return obj.header.components
end

"""
    extract_component(obj, c)

Extract component values.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `c::Symbol`: component name

# Returns

- `component::Any`
"""
function extract_component(obj::NeuroAnalyzer.NEURO; c::Symbol)

    c in obj.header.components || throw(ArgumentError("Component $c does not exist. Use list_component() to view existing components."))
    c_idx = component_idx(obj, c=c)
    
    return obj.components[c_idx]
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

    c in obj.header.components || throw(ArgumentError("Component $c does not exist. Use list_component() to view existing components."))
    
    obj_new = deepcopy(obj)
    c_idx = component_idx(obj, c=c)
    # delete component values
    deleteat!(obj_new.components, c_idx)
    # delete component name
    deleteat!(obj_new.header.components, c_idx)
    push!(obj_new.header.history, "delete_component(NEURO, c=$c)")
    
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
    obj.header.components = obj_new.header.components
    obj.components = obj_new.components

    return nothing
end

"""
    reset_component(obj)

Remove all components.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function reset_component(obj::NeuroAnalyzer.NEURO)

    obj_new = deepcopy(obj)
    obj_new.header.components = Symbol[]
    obj_new.components = Any[]

    return obj_new
end

"""
    reset_component!(obj)

Remove all components.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function reset_component!(obj::NeuroAnalyzer.NEURO)

    obj.header.components = Symbol[]
    obj.components = Any[]

    return nothing
end

