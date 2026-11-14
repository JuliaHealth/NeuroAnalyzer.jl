export edit_locs
export edit_locs!

# the interactive editor is opened only when no edit is requested by keyword arguments
function _el_no_kwargs(args...)::Nothing
    all(a -> isnothing(a) || a == "", args) || throw(
        ArgumentError(
            "ch must be specified to edit coordinates, name or type of a channel.",
        ),
    )
    return nothing
end

"""
    edit_locs(obj; <keyword arguments>)

Edit channel location and properties.

Without `ch` (and without any other keyword argument), opens the interactive editor (GLMakie window) for EEG and MEG objects: select a channel from the list or by clicking its electrode, edit its label, type, unit and coordinates, transform all locations (flip, swap X/Y, center, rotate, scale, normalize, move origin, convert between coordinate systems, generate 10-10 positions), load or save a locations file. Every location change can be undone. Coordinates are shown and edited as stored, in their original units.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Nothing}=nothing`: channel name; must resolve to exactly one channel; if `nothing`, open the interactive editor
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

- `NeuroAnalyzer.NEURO`: output NEURO object; for the interactive editor: the edited copy if the changes were applied, otherwise an unchanged copy
"""
function edit_locs(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Nothing} = nothing,
    x::Union{Real, Nothing} = nothing,
    y::Union{Real, Nothing} = nothing,
    z::Union{Real, Nothing} = nothing,
    theta::Union{Real, Nothing} = nothing,
    radius::Union{Real, Nothing} = nothing,
    theta_sph::Union{Real, Nothing} = nothing,
    radius_sph::Union{Real, Nothing} = nothing,
    phi_sph::Union{Real, Nothing} = nothing,
    name::String = "",
    type::String = "",
)::NeuroAnalyzer.NEURO
    if isnothing(ch)
        _el_no_kwargs(x, y, z, theta, radius, theta_sph, radius_sph, phi_sph, name, type)
        obj_new = _edit_locs_gui(obj)
        return isnothing(obj_new) ? deepcopy(obj) : obj_new
    end

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    length(ch) == 1 || throw(ArgumentError("ch must resolve to exactly one channel."))
    ch = ch[1]

    # create new dataset
    obj_tmp = deepcopy(obj)

    loc_idx = _find_bylabel(obj.locs, labels(obj)[ch])
    loc_idx isa Int64 ||
        throw(ArgumentError("$(labels(obj)[ch]) not found in obj.locs labels."))

    # type first: after renaming the old label no longer exists
    type != "" && set_channel_type!(obj_tmp; ch = labels(obj)[ch], type = type)
    name != "" && rename_channel!(obj_tmp; ch = labels(obj)[ch], name = name)

    x !== nothing && (obj_tmp.locs[loc_idx, :loc_x] = x)
    y !== nothing && (obj_tmp.locs[loc_idx, :loc_y] = y)
    z !== nothing && (obj_tmp.locs[loc_idx, :loc_z] = z)
    theta !== nothing && (obj_tmp.locs[loc_idx, :loc_theta] = theta)
    radius !== nothing && (obj_tmp.locs[loc_idx, :loc_radius] = radius)
    theta_sph !== nothing && (obj_tmp.locs[loc_idx, :loc_theta_sph] = theta_sph)
    radius_sph !== nothing && (obj_tmp.locs[loc_idx, :loc_radius_sph] = radius_sph)
    phi_sph !== nothing && (obj_tmp.locs[loc_idx, :loc_phi_sph] = phi_sph)

    push!(
        obj_tmp.history,
        "edit_locs(obj; ch=$ch, x=$x, y=$y, z=$z, theta=$theta, radius=$radius, theta_sph=$theta_sph, radius_sph=$radius_sph, phi_sph=$phi_sph, name=$name, type=$type)",
    )

    return obj_tmp
end

"""
    edit_locs!(obj; <keyword arguments>)

Edit channel location and properties in-place.

Without `ch` (and without any other keyword argument), opens the interactive editor; the changes are written to `obj` when *Apply* is clicked (GLMakie window) for EEG and MEG objects: select a channel from the list or by clicking its electrode, edit its label, type, unit and coordinates, transform all locations (flip, swap X/Y, center, rotate, scale, normalize, move origin, convert between coordinate systems, generate 10-10 positions), load or save a locations file. Every location change can be undone. Coordinates are shown and edited as stored, in their original units.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Nothing}=nothing`: channel name; must resolve to exactly one channel; if `nothing`, open the interactive editor
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

- `Nothing`
"""
function edit_locs!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Nothing} = nothing,
    x::Union{Real, Nothing} = nothing,
    y::Union{Real, Nothing} = nothing,
    z::Union{Real, Nothing} = nothing,
    theta::Union{Real, Nothing} = nothing,
    radius::Union{Real, Nothing} = nothing,
    theta_sph::Union{Real, Nothing} = nothing,
    radius_sph::Union{Real, Nothing} = nothing,
    phi_sph::Union{Real, Nothing} = nothing,
    name::String = "",
    type::String = "",
)::Nothing
    if isnothing(ch)
        _el_no_kwargs(x, y, z, theta, radius, theta_sph, radius_sph, phi_sph, name, type)
        obj_tmp = _edit_locs_gui(obj)
        isnothing(obj_tmp) && return nothing
        obj.header = obj_tmp.header
        obj.data = obj_tmp.data
        obj.locs = obj_tmp.locs
        obj.markers = obj_tmp.markers
        obj.history = obj_tmp.history
        return nothing
    end

    obj_tmp = edit_locs(
        obj;
        ch = ch,
        x = x,
        y = y,
        z = z,
        theta = theta,
        radius = radius,
        theta_sph = theta_sph,
        radius_sph = radius_sph,
        phi_sph = phi_sph,
        name = name,
        type = type,
    )
    # the channel may have been renamed / retyped (header)
    obj.header = obj_tmp.header
    obj.locs = obj_tmp.locs
    obj.history = obj_tmp.history
    obj_tmp = nothing

    return nothing
end
