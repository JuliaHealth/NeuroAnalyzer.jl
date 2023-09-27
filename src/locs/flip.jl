export locs_flipy
export locs_flipy!
export locs_flipx
export locs_flipx!
export locs_flipz
export locs_flipz!

"""
    locs_flipy(locs; planar, cart, spherical)

Flip channel locations along y axis.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function locs_flipy(locs::DataFrame; planar::Bool=true, cart::Bool=true, spherical::Bool=true)

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

    if planar
        locs_tmp = deepcopy(locs)
        locs_pol2cart!(locs_tmp)
        locs_tmp[!, :loc_y] = -locs_tmp[!, :loc_y]
        locs_cart2pol!(locs_tmp)
        locs_new[!, :loc_radius] = locs_tmp[!, :loc_radius]
        locs_new[!, :loc_theta] = locs_tmp[!, :loc_theta]
    end

#=    for idx in eachindex(locs[!, :labels])
        if planar == true
            t = locs[idx, :loc_theta]
            q = _angle_quadrant(t)
            q == 1 && (t = 360 - t)
            q == 2 && (t = 180 + (180 - t))
            q == 3 && (t = 180 - (t - 180))
            q == 4 && (t = 360 - t)
            t = mod(t, 360)
            locs_new[idx, :loc_theta] = t
            locs_new[idx, :loc_y] = -locs[idx, :loc_y]
        end

        if spherical == true
            t = locs_new[idx, :loc_theta_sph]
            q = _angle_quadrant(t)
            q == 1 && (t = 360 - t)
            q == 2 && (t = 180 + (180 - t))
            q == 3 && (t = 180 - (t - 180))
            q == 4 && (t = 360 - t)
            t = mod(t, 360)
            locs_new[idx, :loc_theta_sph] = t
            locs_new[idx, :loc_y] = -locs[idx, :loc_y]
        end            
    end
=#

    return locs_new

end

"""
    locs_flipy!(locs; planar, cart, spherical)

Flip channel locations along y axis.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates
"""
function locs_flipy!(locs::DataFrame; planar::Bool=true, cart::Bool=true, spherical::Bool=true)

    locs[!, :] = locs_flipy(locs, planar=planar, cart=cart, spherical=spherical)[!, :]

    return nothing
    
end

"""
    locs_flipx(locs; planar, cart, spherical)

Flip channel locations along x axis.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function locs_flipx(locs::DataFrame; planar::Bool=true, cart::Bool=true, spherical::Bool=true)

    locs_new = deepcopy(locs)

    if cart
        locs_new[!, :loc_x] = -locs[!, :loc_x]
    end

    if spherical
        locs_tmp = deepcopy(locs)
        locs_sph2cart!(locs_tmp)
        locs_tmp[!, :loc_x] = -locs_tmp[!, :loc_x]
        locs_cart2sph!(locs_tmp)
        locs_new[!, :loc_radius_sph] = locs_tmp[!, :loc_radius_sph]
        locs_new[!, :loc_theta_sph] = locs_tmp[!, :loc_theta_sph]
        locs_new[!, :loc_phi_sph] = locs_tmp[!, :loc_phi_sph]
    end

    if planar
        locs_tmp = deepcopy(locs)
        locs_pol2cart!(locs_tmp)
        locs_tmp[!, :loc_x] = -locs_tmp[!, :loc_x]
        locs_cart2pol!(locs_tmp)
        locs_new[!, :loc_radius] = locs_tmp[!, :loc_radius]
        locs_new[!, :loc_theta] = locs_tmp[!, :loc_theta]
    end


#=    for idx in eachindex(locs[!, :labels])
        if planar == true
            t = locs[idx, :loc_theta]
            q = _angle_quadrant(t)
            q == 1 && (t = 90 + (90 - t))
            q == 2 && (t = 90 - (t - 90))
            q == 3 && (t = 270 + (270 - t))
            q == 4 && (t = 270 - (t - 270))
            t = mod(t, 360)
            locs_new[idx, :loc_theta] = t
        end

        if spherical == true
            locs_new[idx, :loc_x] = -locs[idx, :loc_x]
            t = locs[idx, :loc_theta_sph]
            q = _angle_quadrant(t)
            q == 1 && (t = 90 + (90 - t))
            q == 2 && (t = 90 - (t - 90))
            q == 3 && (t = 270 + (270 - t))
            q == 4 && (t = 270 - (t - 270))
            t = mod(t, 360)
            locs_new[idx, :loc_theta_sph] = t
        end
    end
=#
    return locs_new

end

"""
    locs_flipx!(locs; planar, cart, spherical)

Flip channel locations along x axis.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates
"""
function locs_flipx!(locs::DataFrame; planar::Bool=true, cart::Bool=true, spherical::Bool=true)

    locs[!, :] = locs_flipx(locs, planar=planar, cart=cart, spherical=spherical)[!, :]

    return nothing

end

"""
    locs_flipz(locs; planar, cart, spherical)

Flip channel locations along z axis.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function locs_flipz(locs::DataFrame)

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

    if planar
        locs_tmp = deepcopy(locs)
        locs_pol2cart!(locs_tmp)
        locs_tmp[!, :loc_z] = -locs_tmp[!, :loc_z]
        locs_cart2pol!(locs_tmp)
        locs_new[!, :loc_radius] = locs_tmp[!, :loc_radius]
        locs_new[!, :loc_theta] = locs_tmp[!, :loc_theta]
    end

#=    locs_new[!, :loc_z] = -locs[!, :loc_z]
    locs_new[!, :loc_phi_sph] = -locs[!, :loc_phi_sph]
=#
    return locs_new

end

"""
    locs_flipz!(locs; planar, cart, spherical)

Flip channel locations along z axis.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates
"""
function locs_flipz!(locs::DataFrame)

    locs[!, :] = locs_flipz(locs)[!, :]

    return nothing

end

