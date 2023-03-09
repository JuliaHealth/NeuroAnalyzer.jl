export locs_flipy
export locs_flipy!
export locs_flipx
export locs_flipx!
export locs_flipz
export locs_flipz!

"""
    locs_flipy(locs; planar, spherical)

Flip channel locations along y axis.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function locs_flipy(locs::DataFrame; planar::Bool=true, spherical::Bool=true)

    locs_new = deepcopy(locs)

    for idx in eachindex(locs[!, :labels])
        if planar == true
            t = locs_new[idx, :loc_theta]
            q = _angle_quadrant(t)
            q == 1 && (t = 360 - t)
            q == 2 && (t = 180 + (180 - t))
            q == 3 && (t = 180 - (t - 180))
            q == 4 && (t = 360 - t)
            t = mod(t, 360)
            locs_new[idx, :loc_theta] = t
        end
        if spherical == true
            locs_new[idx, :loc_y] = -locs_new[idx, :loc_y]
            t = locs_new[idx, :loc_theta_sph]
            q = _angle_quadrant(t)
            q == 1 && (t = 360 - t)
            q == 2 && (t = 180 + (180 - t))
            q == 3 && (t = 180 - (t - 180))
            q == 4 && (t = 360 - t)
            t = mod(t, 360)
            locs_new[idx, :loc_theta_sph] = t
        end            
    end
    # locs_cart2sph!(locs_new)

    return locs_new
end

"""
    locs_flipy!(locs; planar, spherical)

Flip channel locations along y axis.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates
"""
function locs_flipy!(locs::DataFrame; planar::Bool=true, spherical::Bool=true)
    locs[!, :] = locs_flipy(locs, planar=planar, spherical=spherical)[!, :]
    return nothing
end

"""
    locs_flipx(locs; planar, spherical)

Flip channel locations along x axis.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function locs_flipx(locs::DataFrame; planar::Bool=true, spherical::Bool=true)

    locs_new = deepcopy(locs)

    for idx in eachindex(locs[!, :labels])
        if planar == true
            t = locs_new[!, :loc_theta][idx]
            q = _angle_quadrant(t)
            q == 1 && (t = 90 + (90 - t))
            q == 2 && (t = 90 - (t - 90))
            q == 3 && (t = 270 + (270 - t))
            q == 4 && (t = 270 - (t - 270))
            t = mod(t, 360)
            locs_new[!, :loc_theta][idx] = t
        end
        if spherical == true
            locs_new[idx, :loc_x] = -locs_new[idx, :loc_x]
            t = locs_new[idx, :loc_theta_sph]
            q = _angle_quadrant(t)
            q == 1 && (t = 90 + (90 - t))
            q == 2 && (t = 90 - (t - 90))
            q == 3 && (t = 270 + (270 - t))
            q == 4 && (t = 270 - (t - 270))
            t = mod(t, 360)
            locs_new[idx, :loc_theta_sph] = t
        end
    end

    return locs_new
end

"""
    locs_flipx!(locs; planar, spherical)

Flip channel locations along x axis.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates
"""
function locs_flipx!(locs::DataFrame; planar::Bool=true, spherical::Bool=true)
    locs[!, :] = locs_flipx(locs, planar=planar, spherical=spherical)[!, :]
    return nothing
end

"""
    locs_flipz(locs)

Flip channel locations along z axis.

# Arguments

- `locs::DataFrame`

# Returns

- `locs_new::DataFrame`
"""
function locs_flipz(locs::DataFrame)

    locs_new = deepcopy(locs)

    for idx in eachindex(locs[!, :labels])
        locs_new[!, :loc_z][idx] = -locs_new[!, :loc_z][idx]
        locs_new[!, :loc_phi_sph][idx] = -locs_new[!, :loc_phi_sph][idx]
    end
    locs_cart2sph!(locs_new)

    return locs_new
end

"""
    locs_flipz!(locs)

Flip channel locations along z axis.

# Arguments

- `locs::DataFrame`
"""
function locs_flipz!(locs::DataFrame)
    locs[!, :] = locs_flipz(locs)[!, :]
    return nothing
end

