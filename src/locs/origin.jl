export locs_origin
export locs_origin!

"""
    locs_origin(locs; <keyword arguments>)

Move locs origin ([0, 0, 0]) along the axes.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `x::Real=0`: move origin along X axis
- `y::Real=0`: move origin along Y axis
- `z::Real=0`: move origin along Z axis

# Returns

- `locs_new::DataFrame`
"""
function locs_origin(locs::DataFrame; x::Real=0, y::Real=0, z::Real=0)::DataFrame

    locs_new = deepcopy(locs)
    locs_new[:, :loc_x] .+= x
    locs_new[:, :loc_y] .+= y
    locs_new[:, :loc_z] .+= z
    locs_cart2sph!(locs_new)
    locs_cart2pol!(locs_new)

    return locs_new

end

"""
    locs_origin!(locs; <keyword arguments>)

Move locs origin ([0, 0, 0]) along the axes.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `x::Real=0`: move origin along X axis
- `y::Real=0`: move origin along Y axis
- `z::Real=0`: move origin along Z axis

# Returns

Nothing
"""
function locs_origin!(obj::NeuroAnalyzer.NEURO; x::Real=0, y::Real=0, z::Real=0)::Nothing

    obj.locs = locs_origin(obj.locs, x=x, y=y, z=z)

    return nothing

end

"""
    locs_origin!(locs; <keyword arguments>)

Move locs origin ([0, 0, 0]) along the axes.

# Arguments

- `locs::DataFrame`
- `x::Real=0`: move origin along X axis
- `y::Real=0`: move origin along Y axis
- `z::Real=0`: move origin along Z axis

# Returns

Nothing
"""
function locs_origin!(locs::DataFrame; x::Real=0, y::Real=0, z::Real=0)::Nothing

    locs[!, :] = locs_origin(locs, x=x, y=y, z=z)[!, :]

    return nothing

end
