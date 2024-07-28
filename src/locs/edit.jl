export edit_locs
export edit_locs!

"""
    edit_locs(obj; <keyword arguments>)

Edit electrode.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `x::Union{Real, Nothing}=nothing`: Cartesian X spherical coordinate
- `y::Union{Real, Nothing}=nothing`: Cartesian Y spherical coordinate
- `z::Union{Real, Nothing}=nothing`: Cartesian Z spherical coordinate
- `theta::Union{Real, Nothing}=nothing`: polar angle
- `radius::Union{Real, Nothing}=nothing`: polar radius
- `theta_sph::Union{Real, Nothing}=nothing`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
- `radius_sph::Union{Real, Nothing}=nothing`: spherical radius, the distance from the origin to the point
- `phi_sph::Union{Real, Nothing}=nothing`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees
- `name::String=""`: channel name
- `type::String=""`: channel type

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function edit_locs(obj::NeuroAnalyzer.NEURO; ch::String, x::Union{Real, Nothing}=nothing, y::Union{Real, Nothing}=nothing, z::Union{Real, Nothing}=nothing, theta::Union{Real, Nothing}=nothing, radius::Union{Real, Nothing}=nothing, theta_sph::Union{Real, Nothing}=nothing, radius_sph::Union{Real, Nothing}=nothing, phi_sph::Union{Real, Nothing}=nothing, name::String="", type::String="")

    obj_new = deepcopy(obj)
    ch = get_channel(obj_new, ch=ch)
    loc_idx = _find_bylabel(obj.locs, labels(obj)[ch])[1]

    @assert length(loc_idx) > 0 "$(labels(obj)[ch]) not found in obj.locs labels."

    name != "" && rename_channel!(obj_new, ch=labels(obj)[ch], name=name)
    type != "" && channel_type!(obj_new, ch=labels(obj)[ch], type=type)

    x !== nothing && (obj_new.locs[loc_idx, :loc_x] = x)
    y !== nothing && (obj_new.locs[loc_idx, :loc_y] = y)
    z !== nothing && (obj_new.locs[loc_idx, :loc_z] = z)
    theta !== nothing && (obj_new.locs[loc_idx, :loc_theta] = theta)
    radius !== nothing && (obj_new.locs[loc_idx, :loc_radius] = radius)
    theta_sph !== nothing && (obj_new.locs[loc_idx, :loc_theta_sph] = theta_sph)
    radius_sph !== nothing && (obj_new.locs[loc_idx, :loc_radius_sph] = radius_sph)
    phi_sph !== nothing && (obj_new.locs[loc_idx, :loc_phi_sph] = phi_sph)

    reset_components!(obj_new)
    push!(obj_new.history, "edit_locs(OBJ; ch=$ch, x=$x, y=$y, z=$z, theta=$theta, radius=$radius, theta_sph=$theta_sph, radius_sph=$radius_sph, phi_sph=$phi_sph, name=$name, type=$type)")

    return obj_new

end

"""
    edit_locs!(obj; <keyword arguments>)

Edit electrode.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel number or name
- `x::Union{Real, Nothing}=nothing`: Cartesian X spherical coordinate
- `y::Union{Real, Nothing}=nothing`: Cartesian Y spherical coordinate
- `z::Union{Real, Nothing}=nothing`: Cartesian Z spherical coordinate
- `theta::Union{Real, Nothing}=nothing`: polar angle
- `radius::Union{Real, Nothing}=nothing`: polar radius
- `theta_sph::Union{Real, Nothing}=nothing`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
- `radius_sph::Union{Real, Nothing}=nothing`: spherical radius, the distance from the origin to the point
- `phi_sph::Union{Real, Nothing}=nothing`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees
- `name::String=""`: channel name
- `type::String=""`: channel type
"""
function edit_locs!(obj::NeuroAnalyzer.NEURO; ch::String, x::Union{Real, Nothing}=nothing, y::Union{Real, Nothing}=nothing, z::Union{Real, Nothing}=nothing, theta::Union{Real, Nothing}=nothing, radius::Union{Real, Nothing}=nothing, theta_sph::Union{Real, Nothing}=nothing, radius_sph::Union{Real, Nothing}=nothing, phi_sph::Union{Real, Nothing}=nothing, name::String="", type::String="")

    obj_new = edit_locs(obj, ch=ch, x=x, y=y, z=z, theta=theta, radius=radius, theta_sph=theta_sph, radius_sph=radius_sph, phi_sph=phi_sph, name=name, type=type)
    obj.locs = obj_new.locs
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end
