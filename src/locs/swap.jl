export locs_swapxy
export locs_swapxy!

"""
    locs_swapxy(locs; <keyword arguments>)

Swap channel locations x and y axes: `(x, y) → (−y, x)`, i.e. a 90° rotation about the z axis (so that nose direction +X becomes +Y).

# Arguments

- `locs::DataFrame`
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `DataFrame`: modified channel location data
"""
function locs_swapxy(
    locs::DataFrame;
    polar::Bool = true,
    cart::Bool = true,
    spherical::Bool = true,
)::DataFrame
    locs_new = deepcopy(locs)
    # locs_new = locs_rotz(locs, a=90)

    if cart
        locs_new[!, :loc_x], locs_new[!, :loc_y] = locs[!, :loc_y], locs[!, :loc_x]
        locs_new[!, :loc_x] = -locs_new[!, :loc_x]
    end

    if spherical
        # a 90° rotation about z only changes the horizontal angle (exact, no conversion)
        locs_new[!, :loc_theta_sph] = _wrap180(locs[!, :loc_theta_sph] .+ 90)
    end

    polar && locs_rotz!(locs_new; a = 90, polar = true, cart = false, spherical = false)

    _locs_round!(locs_new)
    _locs_remove_nans!(locs_new)

    return locs_new
end

"""
    locs_swapxy!(locs; <keyword arguments>)

Swap channel locations x and y axes.

# Arguments

- `locs::DataFrame`
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `Nothing`
"""
function locs_swapxy!(
    locs::DataFrame;
    polar::Bool = true,
    cart::Bool = true,
    spherical::Bool = true,
)::Nothing
    locs[!, :] = locs_swapxy(locs; polar = polar, cart = cart, spherical = spherical)[!, :]

    return nothing
end
