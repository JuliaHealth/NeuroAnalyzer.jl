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

- `locs::DataFrame`: channel location data
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `DataFrame`: modified channel location data
"""
function locs_flipy(
    locs::DataFrame;
    polar::Bool = true,
    cart::Bool = true,
    spherical::Bool = true,
)::DataFrame
    locs_new = deepcopy(locs)

    if cart
        locs_new[!, :loc_y] = -locs_new[!, :loc_y]
    end

    if spherical
        # y → −y: θ → −θ (exact, no conversion)
        locs_new[!, :loc_theta_sph] = _wrap180(-locs[!, :loc_theta_sph])
    end

    if polar
        # mirror across the x axis: θ → −θ
        locs_new[!, :loc_theta] = mod.(-locs[!, :loc_theta], 360)
    end

    _locs_remove_nans!(locs_new)

    return locs_new
end

"""
    locs_flipy!(locs; <keyword arguments>)

Flip channel locations along y axis.

# Arguments

- `locs::DataFrame`: channel location data
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `Nothing`
"""
function locs_flipy!(
    locs::DataFrame; polar::Bool = true, cart::Bool = true, spherical::Bool = true,
)::Nothing
    locs[!, :] = locs_flipy(locs; polar = polar, cart = cart, spherical = spherical)[!, :]

    return nothing
end

"""
    locs_flipx(locs; <keyword arguments>)

Flip channel locations along x axis.

# Arguments

- `locs::DataFrame`: channel location data
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `DataFrame`: modified channel location data
"""
function locs_flipx(
    locs::DataFrame;
    polar::Bool = true,
    cart::Bool = true,
    spherical::Bool = true,
)::DataFrame
    locs_new = deepcopy(locs)

    cart && (locs_new[!, :loc_x] = -locs[!, :loc_x])

    if spherical
        # x → −x: θ → 180° − θ (exact, no conversion)
        locs_new[!, :loc_theta_sph] = _wrap180(180 .- locs[!, :loc_theta_sph])
    end

    if polar
        # mirror across the y axis: θ → 180° − θ
        locs_new[!, :loc_theta] = mod.(180 .- locs[!, :loc_theta], 360)
    end

    _locs_remove_nans!(locs_new)

    return locs_new
end

"""
    locs_flipx!(locs; <keyword arguments>)

Flip channel locations along x axis.

# Arguments

- `locs::DataFrame`: channel location data
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `Nothing`
"""
function locs_flipx!(
    locs::DataFrame;
    polar::Bool = true,
    cart::Bool = true,
    spherical::Bool = true,
)::Nothing
    locs[!, :] = locs_flipx(locs; polar = polar, cart = cart, spherical = spherical)[!, :]

    return nothing
end

"""
    locs_flipz(locs; <keyword arguments>)

Flip channel locations along z axis.

# Arguments

- `locs::DataFrame`: channel location data
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `DataFrame`: modified channel location data
"""
function locs_flipz(
    locs::DataFrame;
    polar::Bool = true,
    cart::Bool = true,
    spherical::Bool = true,
)::DataFrame
    locs_new = deepcopy(locs)

    if cart
        locs_new[!, :loc_z] = -locs[!, :loc_z]
    end

    if spherical
        # z → −z: elevation φ → −φ (exact, no conversion)
        locs_new[!, :loc_phi_sph] = -locs[!, :loc_phi_sph]
    end

    polar &&
        _warn("For polar coordinates this is a lossy conversion and will be ignored.")

    _locs_remove_nans!(locs_new)

    return locs_new
end

"""
    locs_flipz!(locs; <keyword arguments>)

Flip channel locations along z axis.

# Arguments

- `locs::DataFrame`: channel location data
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `Nothing`
"""
function locs_flipz!(
    locs::DataFrame;
    polar::Bool = true,
    cart::Bool = true,
    spherical::Bool = true,
)::Nothing
    locs[!, :] = locs_flipz(locs; polar = polar, cart = cart, spherical = spherical)[!, :]

    return nothing
end
