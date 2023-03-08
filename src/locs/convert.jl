export locs_sph2cart
export locs_sph2cart!
export locs_cart2sph
export locs_cart2sph!
export locs_cart2pol
export locs_cart2pol!
export locs_sph2pol
export locs_sph2pol!

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
