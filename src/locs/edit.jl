"""
    edit_locs(obj; <keyword arguments>)

Edit electrode.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{String, Int64}`: channel number or name
- `x::Union{Real, Nothing}=nothing`: Cartesian X spherical coordinate
- `y::Union{Real, Nothing}=nothing`: Cartesian Y spherical coordinate
- `z::Union{Real, Nothing}=nothing`: Cartesian Z spherical coordinate
- `theta::Union{Real, Nothing}=nothing`: polar planar theta coordinate
- `radius::Union{Real, Nothing}=nothing`: polar planar radius coordinate
- `theta_sph::Union{Real, Nothing}=nothing`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
- `radius_sph::Union{Real, Nothing}=nothing`: spherical radius, the distance from the origin to the point
- `phi_sph::Union{Real, Nothing}=nothing`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees
- `name::String=""`: channel name
- `type::String=""`: channel type

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function edit_locs(obj::NeuroAnalyzer.NEURO; channel::Union{String, Int64}, x::Union{Real, Nothing}=nothing, y::Union{Real, Nothing}=nothing, z::Union{Real, Nothing}=nothing, theta::Union{Real, Nothing}=nothing, radius::Union{Real, Nothing}=nothing, theta_sph::Union{Real, Nothing}=nothing, radius_sph::Union{Real, Nothing}=nothing, phi_sph::Union{Real, Nothing}=nothing, name::String="", type::String="")

    obj_new = deepcopy(obj)
    channel = _get_ch_idx(labels(obj_new), channel)

    name != "" && rename_channel!(obj_new, channel=channel, name=name)
    type != "" && channel_type!(obj_new, channel=channel, type=type)

    x !== nothing && (obj_new.locs[!, :loc_x][channel] = x)
    y !== nothing && (obj_new.locs[!, :loc_y][channel] = y)
    z !== nothing && (obj_new.locs[!, :loc_z][channel] = z)
    theta !== nothing && (obj_new.locs[!, :loc_theta][channel] = theta)
    radius !== nothing && (obj_new.locs[!, :loc_radius][channel] = radius)
    theta_sph !== nothing && (obj_new.locs[!, :loc_theta_sph][channel] = theta_sph)
    radius_sph !== nothing && (obj_new.locs[!, :loc_radius_sph][channel] = radius_sph)
    phi_sph !== nothing && (obj_new.locs[!, :loc_phi_sph][channel] = phi_sph)

    (x !== nothing || y !== nothing || z !== nothing || theta !== nothing || radius !== nothing || theta_sph !== nothing  || radius_sph !== nothing || phi_sph !== nothing) && (obj_new.header.has_locs == true)

    reset_components!(obj_new)
    push!(obj_new.header.history, "edit_locs(OBJ; channel=$channel, x=$x, y=$y, z=$z, theta=$theta, radius=$radius, theta_sph=$theta_sph, radius_sph=$radius_sph, phi_sph=$phi_sph, name=$name, type=$type)")

    return obj_new
end

"""
    edit_locs!(obj; <keyword arguments>)

Edit electrode.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{String, Int64}`: channel number or name
- `x::Union{Real, Nothing}=nothing`: Cartesian X spherical coordinate
- `y::Union{Real, Nothing}=nothing`: Cartesian Y spherical coordinate
- `z::Union{Real, Nothing}=nothing`: Cartesian Z spherical coordinate
- `theta::Union{Real, Nothing}=nothing`: polar planar theta coordinate
- `radius::Union{Real, Nothing}=nothing`: polar planar radius coordinate
- `theta_sph::Union{Real, Nothing}=nothing`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
- `radius_sph::Union{Real, Nothing}=nothing`: spherical radius, the distance from the origin to the point
- `phi_sph::Union{Real, Nothing}=nothing`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees
- `name::String=""`: channel name
- `type::String=""`: channel type
"""
function edit_locs!(obj::NeuroAnalyzer.NEURO; channel::Union{String, Int64}, x::Union{Real, Nothing}=nothing, y::Union{Real, Nothing}=nothing, z::Union{Real, Nothing}=nothing, theta::Union{Real, Nothing}=nothing, radius::Union{Real, Nothing}=nothing, theta_sph::Union{Real, Nothing}=nothing, radius_sph::Union{Real, Nothing}=nothing, phi_sph::Union{Real, Nothing}=nothing, name::String="", type::String="")

    obj_tmp = edit_locs(obj, channel=channel, x=x, y=y, z=z, theta=theta, radius=radius, theta_sph=theta_sph, radius_sph=radius_sph, phi_sph=phi_sph, name=name, type=type)
    obj.header = obj_tmp.header
    obj.locs = obj_tmp.locs
    reset_components!(obj)

    return nothing
end
