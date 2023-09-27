export locs_swapxy
export locs_swapxy!

"""
    locs_swapxy(locs; planar, cart, spherical)

Swap channel locations x and y axes.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function locs_swapxy(locs::DataFrame; planar::Bool=true, cart::Bool=true, spherical::Bool=true)

    locs_new = deepcopy(locs)
    # locs_new = locs_rotz(locs, a=90)

    if cart
        locs_new[!, :loc_x], locs_new[!, :loc_y] = locs[!, :loc_y], locs[!, :loc_x]
    end

    if spherical
        locs_tmp = deepcopy(locs)
        locs_sph2cart!(locs_tmp)
        locs_tmp[!, :loc_x], locs_tmp[!, :loc_y] = locs[!, :loc_y], locs[!, :loc_x]
        locs_cart2sph!(locs_tmp)
        locs_new[!, :loc_radius_sph] = locs_tmp[!, :loc_radius_sph]
        locs_new[!, :loc_theta_sph] = locs_tmp[!, :loc_theta_sph]
        locs_new[!, :loc_phi_sph] = locs_tmp[!, :loc_phi_sph]
    end

    if planar
        locs_tmp = deepcopy(locs)
        locs_pol2cart!(locs_tmp)
        locs_tmp[!, :loc_x], locs_tmp[!, :loc_y] = locs[!, :loc_y], locs[!, :loc_x]
        locs_cart2pol!(locs_tmp)
        locs_new[!, :loc_radius] = locs_tmp[!, :loc_radius]
        locs_new[!, :loc_theta] = locs_tmp[!, :loc_theta]
    end

#=    for idx in eachindex(locs[!, :labels])
        if planar
            t = deg2rad(locs_new[idx, :loc_theta])
            t += pi / 2
            locs_new[idx, :loc_theta] = rad2deg(t)
        end
        if spherical
            locs_new[idx, :loc_x], locs_new[idx, :loc_y] = locs_new[idx, :loc_y], locs_new[idx, :loc_x]
            t = deg2rad(locs_new[idx, :loc_theta_sph])
            t += pi / 2
            locs_new[idx, :loc_theta_sph] = rad2deg(t)
        end
    end=#

    return locs_new
end

"""
    locs_swapxy!(locs; planar, cart, spherical)

Swap channel locations x and y axes.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates
"""
function locs_swapxy!(locs::DataFrame; planar::Bool=true, cart::Bool=true, spherical::Bool=true)
    
    locs[!, :] = locs_swapxy(locs, planar=planar, cart=cart, spherical=spherical)[!, :]

    return nothing

end
