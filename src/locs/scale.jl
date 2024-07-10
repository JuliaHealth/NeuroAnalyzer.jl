export locs_scale
export locs_scale!
export locs_normalize
export locs_normalize!

"""
    locs_scale(locs; <keyword arguments>)

Scale channel locations.

# Arguments

- `locs::DataFrame`
- `r::Real`: scaling factor
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function locs_scale(locs::DataFrame; r::Real, polar::Bool=true, cart::Bool=true, spherical::Bool=true)

    locs_new = deepcopy(locs)

    polar && (locs_new[!, :loc_radius] .*= r)

    spherical && (locs_new[!, :loc_radius_sph] .*= r)

    if cart
        locs_tmp = deepcopy(locs)
        locs_cart2sph!(locs_tmp)
        locs_tmp[!, :loc_radius_sph] .*= r
        locs_sph2cart!(locs_tmp)
        locs_new[!, :loc_x], locs_new[!, :loc_y], locs_new[!, :loc_z] = locs_tmp[!, :loc_x], locs_tmp[!, :loc_y], locs_tmp[!, :loc_z]
    end

    _locs_round!(locs_new)
    _locs_remove_nans!(locs_new)

    return locs_new

end

"""
    locs_scale!(locs, r, polar, cart, spherical)

Scale channel locations.

# Arguments

- `locs::DataFrame`
- `r::Real`: scaling factor
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates
"""
function locs_scale!(locs::DataFrame; r::Real, polar::Bool=true, cart::Bool=true, spherical::Bool=true)

    locs[!, :] = locs_scale(locs, r=r, polar=polar, cart=cart, spherical=spherical)[!, :]

    return nothing

end

"""
    locs_normalize(locs; <keyword arguments>)

Normalize channel locations to fit the unit sphere.

# Arguments

- `locs::DataFrame`
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function locs_normalize(locs::DataFrame; polar::Bool=true, cart::Bool=true, spherical::Bool=true)

    locs_new = deepcopy(locs)

    if polar
        r1 = locs_new[!, :loc_radius]
        r2 = normalize_minmax(locs_new[!, :loc_radius])
        r = maximum(r2) / maximum(r1)
        locs_new[!, :loc_radius] .*= r
    end

    if spherical
        r1 = locs_new[!, :loc_radius_sph]
        r2 = normalize_minmax(r1)
        r = maximum(r2) / maximum(r1)
        locs_new[!, :loc_radius_sph] .*= r
    end

    if cart
        locs_tmp = deepcopy(locs)
        locs_cart2sph!(locs_tmp)
        r1 = locs_tmp[!, :loc_radius_sph]
        r2 = normalize_minmax(r1)
        r = maximum(r2) / maximum(r1)
        locs_tmp[!, :loc_radius_sph] .*= r
        locs_sph2cart!(locs_tmp)
        locs_new[!, :loc_x], locs_new[!, :loc_y], locs_new[!, :loc_z] = locs_tmp[!, :loc_x], locs_tmp[!, :loc_y], locs_tmp[!, :loc_z]
    end

    _locs_round!(locs_new)
    _locs_remove_nans!(locs_new)

    return locs_new

end

"""
    locs_normalize!(locs; <keyword arguments>)

Normalize channel locations to fit the unit sphere.

# Arguments

- `locs::DataFrame`
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates
"""
function locs_normalize!(locs::DataFrame; polar::Bool=true, cart::Bool=true, spherical::Bool=true)

    locs[!, :] = locs_normalize(locs, polar=polar, cart=cart, spherical=spherical)[!, :]

    return nothing

end
