export locs_swapxy
export locs_swapxy!

"""
    locs_swapxy(locs; <keyword arguments>)

Swap channel locations x and y axes.

# Arguments

- `locs::DataFrame`
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function locs_swapxy(locs::DataFrame; polar::Bool=true, cart::Bool=true, spherical::Bool=true)::DataFrame

    locs_new = deepcopy(locs)
    # locs_new = locs_rotz(locs, a=90)

    if cart
        locs_new[!, :loc_x], locs_new[!, :loc_y] = locs[!, :loc_y], locs[!, :loc_x]
        locs_new[!, :loc_x] = -locs_new[!, :loc_x]
    end

    if spherical
        locs_tmp = deepcopy(locs)
        locs_sph2cart!(locs_tmp)
        locs_tmp[!, :loc_x], locs_tmp[!, :loc_y] = locs[!, :loc_y], locs[!, :loc_x]
        locs_tmp[!, :loc_x] = -locs_tmp[!, :loc_x]
        locs_cart2sph!(locs_tmp)
        locs_new[!, :loc_radius_sph] = locs_tmp[!, :loc_radius_sph]
        locs_new[!, :loc_theta_sph] = locs_tmp[!, :loc_theta_sph]
        locs_new[!, :loc_phi_sph] = locs_tmp[!, :loc_phi_sph]
    end

    polar && locs_rotz!(locs_new, a=90, polar=true, cart=false, spherical=false)

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

Nothing
"""
function locs_swapxy!(locs::DataFrame; polar::Bool=true, cart::Bool=true, spherical::Bool=true)::Nothing

    locs[!, :] = locs_swapxy(locs, polar=polar, cart=cart, spherical=spherical)[!, :]

    return nothing

end
