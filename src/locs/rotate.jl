export locs_rotz
export locs_rotz!
export locs_roty
export locs_roty!
export locs_rotx
export locs_rotx!

"""
    locs_rotz(locs; <keyword arguments>)

Rotate channel locations around the Z axis.

# Arguments

- `locs::DataFrame`
- `a::Real`: angle of rotation (in degrees); positive angle rotates anti-clockwise
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function locs_rotz(locs::DataFrame; a::Real, polar::Bool=true, cart::Bool=true, spherical::Bool=true)

    locs_new = deepcopy(locs)

    if cart
        for idx in eachindex(locs[!, :label])
            locs_new[idx, :loc_x] = locs[idx, :loc_x] * cosd(a) - locs[idx, :loc_y] * sind(a)
            locs_new[idx, :loc_y] = locs[idx, :loc_x] * sind(a) + locs[idx, :loc_y] * cosd(a)
        end
    end

    if spherical
        locs_tmp = deepcopy(locs)
        locs_sph2cart!(locs_tmp)
        for idx in eachindex(locs[!, :label])
            locs_tmp[idx, :loc_x] = locs_tmp[idx, :loc_x] * cosd(a) - locs_tmp[idx, :loc_y] * sind(a)
            locs_tmp[idx, :loc_y] = locs_tmp[idx, :loc_x] * sind(a) + locs_tmp[idx, :loc_y] * cosd(a)
        end
        locs_cart2sph!(locs_tmp)
        locs_new[!, :loc_radius_sph] = locs_tmp[!, :loc_radius_sph]
        locs_new[!, :loc_theta_sph] = locs_tmp[!, :loc_theta_sph]
        locs_new[!, :loc_phi_sph] = locs_tmp[!, :loc_phi_sph]
    end

    polar && locs_new[!, :loc_theta] .+= a % 360

    _locs_round!(locs_new)
    _locs_remove_nans!(locs_new)

    return locs_new

end

"""
    locs_rotz!(locs; <keyword arguments>)

Rotate channel locations in the xy-plane.

# Arguments

- `locs::DataFrame`
- `a::Real`: angle of rotation (in degrees); positive angle rotates anti-clockwise
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates
"""
function locs_rotz!(locs::DataFrame; a::Real, polar::Bool=true, cart::Bool=true, spherical::Bool=true)

    locs[!, :] = locs_rotz(locs, a=a, polar=polar, cart=cart, spherical=spherical)[!, :]

    return nothing

end

"""
    locs_roty(locs; <keyword arguments>)

Rotate channel locations around the Y axis (in the XZ-plane).

# Arguments

- `locs::DataFrame`
- `a::Real`: angle of rotation (in degrees); positive angle rotates clockwise
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function locs_roty(locs::DataFrame; a::Real, polar::Bool=true, cart::Bool=true, spherical::Bool=true)

    locs_new = deepcopy(locs)

    if cart
        for idx in 1:nrow(locs)
            locs_new[idx, :loc_x] = locs[idx, :loc_x] * cosd(a) + locs[idx, :loc_z] * sind(a)
            locs_new[idx, :loc_z] = -locs[idx, :loc_x] * sind(a) + locs[idx, :loc_z] * cosd(a)
        end
    end

    if spherical
        locs_tmp = deepcopy(locs)
        locs_sph2cart!(locs_tmp)
        for idx in 1:nrow(locs)
            locs_tmp[idx, :loc_x] = locs[idx, :loc_x] * cosd(a) + locs[idx, :loc_z] * sind(a)
            locs_tmp[idx, :loc_z] = -locs[idx, :loc_x] * sind(a) + locs[idx, :loc_z] * cosd(a)
        end
        locs_cart2sph!(locs_tmp)
        locs_new[!, :loc_radius_sph] = locs_tmp[!, :loc_radius_sph]
        locs_new[!, :loc_theta_sph] = locs_tmp[!, :loc_theta_sph]
        locs_new[!, :loc_phi_sph] = locs_tmp[!, :loc_phi_sph]
    end

    polar && _warn("For polar coordinates this is a lossy conversion and will be ignored.")

    _locs_round!(locs_new)
    _locs_remove_nans!(locs_new)

    return locs_new

end

"""
    locs_roty!(locs; <keyword arguments>)

Rotate channel locations around the Y axis (in the XZ-plane).

# Arguments

- `locs::DataFrame`
- `a::Real`: angle of rotation (in degrees); positive angle rotates clockwise
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates
"""
function locs_roty!(locs::DataFrame; a::Real, polar::Bool=true, cart::Bool=true, spherical::Bool=true)

    locs[!, :] = locs_roty(locs, a=a, polar=polar, cart=cart, spherical=spherical)[!, :]

    return nothing

end

"""
    locs_rotx(locs; <keyword arguments>)

Rotate channel locations around the X axis (in the YZ-plane).

# Arguments

- `locs::DataFrame`
- `a::Real`: angle of rotation (in degrees); positive angle rotates anti-clockwise
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function locs_rotx(locs::DataFrame; a::Real, polar::Bool=true, cart::Bool=true, spherical::Bool=true)

    locs_new = deepcopy(locs)

    if cart
        for idx in 1:nrow(locs)
            locs_new[idx, :loc_y] = locs[idx, :loc_y] * cosd(a) - locs[idx, :loc_z] * sind(a)
            locs_new[idx, :loc_z] = locs[idx, :loc_y] * sind(a) + locs[idx, :loc_z] * cosd(a)
        end
    end

    if spherical
        locs_tmp = deepcopy(locs)
        locs_sph2cart!(locs_tmp)
        for idx in 1:nrow(locs)
            locs_tmp[idx, :loc_y] = locs[idx, :loc_y] * cosd(a) - locs[idx, :loc_z] * sind(a)
            locs_tmp[idx, :loc_z] = locs[idx, :loc_y] * sind(a) + locs[idx, :loc_z] * cosd(a)
        end
        locs_cart2sph!(locs_tmp)
        locs_new[!, :loc_radius_sph] = locs_tmp[!, :loc_radius_sph]
        locs_new[!, :loc_theta_sph] = locs_tmp[!, :loc_theta_sph]
        locs_new[!, :loc_phi_sph] = locs_tmp[!, :loc_phi_sph]
    end

    polar && _warn("For polar coordinates this is a lossy conversion and will be ignored.")

    _locs_round!(locs_new)
    _locs_remove_nans!(locs_new)

    return locs_new

end

"""
    locs_rotx!(locs; <keyword arguments>)

Rotate channel locations around the X axis (in the YZ-plane).

# Arguments

- `locs::DataFrame`
- `a::Real`: angle of rotation (in degrees); positive angle rotates anti-clockwise
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates
"""
function locs_rotx!(locs::DataFrame; a::Real, polar::Bool=true, cart::Bool=true, spherical::Bool=true)

    locs[!, :] = locs_rotx(locs, a=a, polar=polar, cart=cart, spherical=spherical)[!, :]

    return nothing

end
