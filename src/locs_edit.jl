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
    locs_flipz!(eeg)

Flip channel locations along z axis.

# Arguments

- `locs::DataFrame`
"""
function locs_flipz!(locs::DataFrame)
    locs[!, :] = locs_flipz(locs)[!, :]
    return nothing
end

"""
    locs_scale(locs; r, planar, spherical)

Scale channel locations.

# Arguments

- `locs::DataFrame`
- `r::Real`: scaling factor
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function locs_scale(locs::DataFrame; r::Real, planar::Bool=true, spherical::Bool=true)
    locs_new = deepcopy(locs)
    planar == true && (locs_new[!, :loc_radius] .*= r)
    spherical == true && (locs_new[!, :loc_radius_sph] .*= r)
    locs_sph2cart!(locs_new)
    return locs_new
end

"""
    locs_scale!(eeg)

Scale channel locations.

# Arguments

- `locs::DataFrame`
- `r::Real`: scaling factor
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates
"""
function locs_scale!(locs::DataFrame; r::Real, planar::Bool=true, spherical::Bool=true)
    locs[!, :] = locs_scale(locs, r=r, planar=planar, spherical=spherical)[!, :]
    return nothing
end

"""
    locs_rotx(locs; a, planar, spherical)

Rotate channel locations in the x-plane.

# Arguments

- `locs::DataFrame`
- `a::Real`: angle of rotation (in degrees)
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function locs_rotx(locs::DataFrame; a::Real, planar::Bool=true, spherical::Bool=true)
    locs_new = deepcopy(locs)
    planar == true && (locs_new[!, :loc_theta] .+= a)
    spherical == true && (locs_new[!, :loc_theta_sph] .+= a)
    for idx in 1:nrow(locs)
        locs_new[idx, :loc_x] = locs_new[idx, :loc_x] * cosd(a) + locs_new[idx, :loc_y] * sind(a)
        locs_new[idx, :loc_y] = locs_new[idx, :loc_y] * cosd(a) - locs_new[idx, :loc_x] * sind(a)
    end
    return locs_new
end

"""
    locs_rotx!(eeg)

Rotate channel locations in the x-plane.

# Arguments

- `locs::DataFrame`
- `a::Int64`: scaling factor
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates
"""
function locs_rotx!(locs::DataFrame; a::Real, planar::Bool=true, spherical::Bool=true)
    locs[!, :] = locs_rotx(locs, a=a, planar=planar, spherical=spherical)[!, :]
    return nothing
end

"""
    locs_swapxy(locs; planar, spherical)

Swap channel locations x and y axes.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `eeg::NeuroAnalyzer.EEG`
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

"""
    locs_sph2cart(locs)

Convert spherical locations to Cartesian.

# Arguments

- `locs::DataFrame`

# Returns

- `locs_new::DataFrame`
"""
function locs_sph2cart(locs::DataFrame)

    locs_new = deepcopy(locs)

    for idx in eachindex(locs[!, :labels])
        r = locs_new[idx, :loc_radius_sph]
        t = locs_new[idx, :loc_theta_sph]
        p = locs_new[idx, :loc_phi_sph]
        x, y, z = sph2cart(r, t, p)
        locs_new[idx, :loc_x] = x
        locs_new[idx, :loc_y] = y
        locs_new[idx, :loc_z] = z
    end

    return locs_new
end

"""
    locs_sph2cart!(locs)

Convert spherical locations to Cartesian.

# Arguments

- `locs::DataFrame`
"""
function locs_sph2cart!(locs::DataFrame)
    locs[!, :] = locs_sph2cart(locs)[!, :]
    return nothing
end

"""
    locs_cart2sph(locs)

Convert Cartesian locations to spherical.

# Arguments

- `locs::DataFrame`

# Returns

- `locs_new::DataFrame`
"""
function locs_cart2sph(locs::DataFrame)

    locs_new = deepcopy(locs)

    for idx in eachindex(locs[!, :labels])
        x = locs_new[!, :loc_x][idx]
        y = locs_new[!, :loc_y][idx]
        z = locs_new[!, :loc_z][idx]
        r, t, p = cart2sph(x, y, z)
        locs_new[!, :loc_radius_sph][idx] = r
        locs_new[!, :loc_theta_sph][idx] = t
        locs_new[!, :loc_phi_sph][idx] = p
    end

    return locs_new
end

"""
    locs_cart2sph!(locs)

Convert Cartesian locations to spherical.

# Arguments

- `locs::DataFrame`
"""
function locs_cart2sph!(locs::DataFrame)
    locs[!, :] = locs_cart2sph(locs)[!, :]
    return nothing
end

"""
    locs_cart2pol(locs)

Convert Cartesian locations to polar.

# Arguments

- `locs::DataFrame`

# Returns

- `locs_new::DataFrame`
"""
function locs_cart2pol(locs::DataFrame)

    locs_new = deepcopy(locs)

    for idx in eachindex(locs[!, :labels])
        x = locs_new[!, :loc_x][idx]
        y = locs_new[!, :loc_y][idx]
        r, t = cart2pol(x, y)
        locs_new[!, :loc_radius][idx] = r
        locs_new[!, :loc_theta][idx] = t
    end

    return locs_new
end

"""
    locs_cart2pol!(locs)

Convert Cartesian locations to polar.

# Arguments

- `locs::DataFrame`
"""
function locs_cart2pol!(locs::DataFrame)
    locs[!, :] = locs_cart2pol(locs)[!, :]
    return nothing
end

"""
    locs_sph2pol(locs)

Convert spherical locations to polar.

# Arguments

- `locs::DataFrame`

# Returns

- `locs_new::DataFrame`
"""
function locs_sph2pol(locs::DataFrame)

    locs_new = deepcopy(locs)

    for idx in eachindex(locs[!, :labels])
        r_sph = locs_new[!, :loc_radius_sph][idx]
        t_sph = locs_new[!, :loc_theta_sph][idx]
        p_sph = locs_new[!, :loc_phi_sph][idx]
        r_pol, t_pol = sph2pol(r_sph, t_sph, p_sph)
        locs_new[!, :loc_radius][idx] = r_pol
        locs_new[!, :loc_theta][idx] = t_pol
    end

    return locs_new
end

"""
    locs_sph2pol!(locs)

Convert Cartesian locations to polar.

# Arguments

- `locs::DataFrame`
"""
function locs_sph2pol!(locs::DataFrame)
    locs[!, :] = locs_sph2pol(locs)[!, :]
    return nothing
end

"""
    locs_maximize(locs; planar, spherical)

Maximize channel locations to the unit sphere.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=false`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function locs_maximize(locs::DataFrame; planar::Bool=true, spherical::Bool=false)

    locs_new = deepcopy(locs)

    if planar == true
        r1 = locs_new[!, :loc_radius]
        r2 = s_normalize_minmax(locs_new[!, :loc_radius])
        r = maximum(r2) / maximum(r1)
        locs_new[!, :loc_radius] .*= r
    end

    if spherical == true
        r1 = locs_new[!, :loc_radius_sph]
        r2 = s_normalize_minmax(r1)
        r = maximum(r2) / maximum(r1)
        locs_new[!, :loc_radius_sph] .*= r
    end
    locs_sph2cart!(locs_new)

    return locs_new
end

"""
    locs_maximize!(locs; planar, spherical)

Maximize channel locations to the unit sphere.

# Arguments

- `locs::DataFrame`
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=false`: modify spherical coordinates
"""
function locs_maximize!(locs::DataFrame; planar::Bool=true, spherical::Bool=false)
    locs[!, :] = locs_maximize(locs, planar=planar, spherical=spherical)[!, :]
    return nothing
end
