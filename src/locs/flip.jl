export locs_flipy
export locs_flipy!
export locs_flipx
export locs_flipx!
export locs_flipz
export locs_flipz!

"""
    locs_flipy(locs; <keyword arguments>)

Flip channel locations along y axis.

# Arguments

- `locs::DataFrame`
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function locs_flipy(locs::DataFrame; polar::Bool=true, cart::Bool=true, spherical::Bool=true)::DataFrame

    locs_new = deepcopy(locs)

    if cart
        locs_new[!, :loc_y] = -locs_new[!, :loc_y]
    end

    if spherical
        locs_tmp = deepcopy(locs)
        locs_sph2cart!(locs_tmp)
        locs_tmp[!, :loc_y] = -locs_tmp[!, :loc_y]
        locs_cart2sph!(locs_tmp)
        locs_new[!, :loc_radius_sph] = locs_tmp[!, :loc_radius_sph]
        locs_new[!, :loc_theta_sph] = locs_tmp[!, :loc_theta_sph]
        locs_new[!, :loc_phi_sph] = locs_tmp[!, :loc_phi_sph]
    end

    if polar
        for idx in eachindex(locs[!, :label])
            t = locs[idx, :loc_theta]
            q = _angle_quadrant(t)
            q == 1 && (t = 360 - t)
            q == 2 && (t = 180 + (180 - t))
            q == 3 && (t = 180 - (t - 180))
            q == 4 && (t = 360 - t)
            t = mod(t, 360)
            locs_new[idx, :loc_theta] = t
        end
    end

    _locs_round!(locs_new)
    _locs_remove_nans!(locs_new)

    return locs_new

end

"""
    locs_flipy!(locs; <keyword arguments>)

Flip channel locations along y axis.

# Arguments

- `locs::DataFrame`
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

Nothing
"""
function locs_flipy!(locs::DataFrame; polar::Bool=true, cart::Bool=true, spherical::Bool=true)::Nothing

    locs[!, :] = locs_flipy(locs, polar=polar, cart=cart, spherical=spherical)[!, :]

    return nothing

end

"""
    locs_flipx(locs; <keyword arguments>)

Flip channel locations along x axis.

# Arguments

- `locs::DataFrame`
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function locs_flipx(locs::DataFrame; polar::Bool=true, cart::Bool=true, spherical::Bool=true)::DataFrame

    locs_new = deepcopy(locs)

    cart && (locs_new[!, :loc_x] = -locs[!, :loc_x])

    if spherical
        locs_tmp = deepcopy(locs)
        locs_sph2cart!(locs_tmp)
        locs_tmp[!, :loc_x] = -locs_tmp[!, :loc_x]
        locs_cart2sph!(locs_tmp)
        locs_new[!, :loc_radius_sph] = locs_tmp[!, :loc_radius_sph]
        locs_new[!, :loc_theta_sph] = locs_tmp[!, :loc_theta_sph]
        locs_new[!, :loc_phi_sph] = locs_tmp[!, :loc_phi_sph]
    end

    if polar
        for idx in 1:nrow(locs)
            t = locs[idx, :loc_theta]
            q = _angle_quadrant(t)
            q == 1 && (t = 90 + (90 - t))
            q == 2 && (t = 90 - (t - 90))
            q == 3 && (t = 270 + (270 - t))
            q == 4 && (t = 270 - (t - 270))
            locs_new[idx, :loc_theta] = t % 360
        end
    end

    _locs_round!(locs_new)
    _locs_remove_nans!(locs_new)

    return locs_new

end

"""
    locs_flipx!(locs; <keyword arguments>)

Flip channel locations along x axis.

# Arguments

- `locs::DataFrame`
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

Nothing
"""
function locs_flipx!(locs::DataFrame; polar::Bool=true, cart::Bool=true, spherical::Bool=true)::Nothing

    locs[!, :] = locs_flipx(locs, polar=polar, cart=cart, spherical=spherical)[!, :]

    return nothing

end

"""
    locs_flipz(locs; <keyword arguments>)

Flip channel locations along z axis.

# Arguments

- `locs::DataFrame`
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function locs_flipz(locs::DataFrame; polar::Bool=true, cart::Bool=true, spherical::Bool=true)::DataFrame

    locs_new = deepcopy(locs)

    if cart
        locs_new[!, :loc_z] = -locs[!, :loc_z]
    end

    if spherical
        locs_tmp = deepcopy(locs)
        locs_sph2cart!(locs_tmp)
        locs_tmp[!, :loc_z] = -locs_tmp[!, :loc_z]
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
    locs_flipz!(locs; <keyword arguments>)

Flip channel locations along z axis.

# Arguments

- `locs::DataFrame`
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

Nothing
"""
function locs_flipz!(locs::DataFrame; polar::Bool=true, cart::Bool=true, spherical::Bool=true)::Nothing

    locs[!, :] = locs_flipz(locs, polar=polar, cart=cart, spherical=spherical)[!, :]

    return nothing

end
