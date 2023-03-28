export locs_swapxy
export locs_swapxy!

"""
    locs_swapxy(locs; planar, spherical)

Swap channel locations x and y axes.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function locs_swapxy(locs::DataFrame; planar::Bool=true, spherical::Bool=true)

    locs_new = deepcopy(locs)

    for idx in eachindex(locs[!, :labels])
        if planar == true
            t = deg2rad(locs_new[idx, :loc_theta])
            t += pi / 2
            locs_new[idx, :loc_theta] = rad2deg(t)
        end
        if spherical == true
            locs_new[idx, :loc_x], locs_new[idx, :loc_y] = locs_new[idx, :loc_y], locs_new[idx, :loc_x]
            t = deg2rad(locs_new[idx, :loc_theta_sph])
            t += pi / 2
            locs_new[idx, :loc_theta_sph] = rad2deg(t)
        end
    end

    return locs_new
end

"""
    locs_swapxy!(locs; planar, spherical)

Swap channel locations x and y axes.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates
"""
function locs_swapxy!(locs::DataFrame; planar::Bool=true, spherical::Bool=true)
    locs[!, :] = locs_swapxy(locs, planar=planar, spherical=spherical)[!, :]
    return nothing
end
