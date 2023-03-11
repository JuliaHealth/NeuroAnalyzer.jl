export cart2pol
export pol2cart
export sph2cart
export cart2sph
export sph2pol
export locs_sph2cart
export locs_sph2cart!
export locs_cart2sph
export locs_cart2sph!
export locs_cart2pol
export locs_cart2pol!
export locs_sph2pol
export locs_sph2pol!

"""
    cart2pol(x, y)

Convert Cartesian coordinates to polar.

# Arguments

- `x::Real`
- `y::Real`

# Returns

- `radius::Float64`
- `theta::Float64`
"""
function cart2pol(x::Real, y::Real)

    radius = round(hypot(x, y), digits=3)
    theta = round(atand(y, x), digits=3)

    # q = _angle_quadrant(theta)
    # q == 2 && (theta += 180)
    # q == 3 && (theta += 180)
    # q == 4 && (theta += 360)

    # make theta positive - should not be required
    # theta < 0 && (theta = 360 - abs(theta))

    return radius, theta
end

"""
    pol2cart(radius, theta)

Convert polar coordinates to Cartesian.

# Arguments

- `radius::Real`: polar radius, the distance from the origin to the point, in degrees
- `theta::Real`: polar angle

# Returns

- `x::Float64`
- `y::Float64`
"""
function pol2cart(radius::Real, theta::Real)
    
    # theta = mod(theta, 360)
    x = round(radius * cosd(theta), digits=3)
    y = round(radius * sind(theta), digits=3)

    return x, y
end

"""
    sph2cart(radius, theta, phi)

Convert spherical coordinates to Cartesian.

# Arguments

- `radius::Real`: spherical radius, the distance from the origin to the point
- `theta::Real`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
- `phi::Real`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees

# Returns

- `x::Float64`
- `y::Float64`
- `z::Float64`
"""
function sph2cart(radius::Real, theta::Real, phi::Real)
    
    x = round(radius * sind(90 - phi) * cosd(theta), digits=3)
    y = round(radius * sind(90 - phi) * sind(theta), digits=3)
    z = round(radius * cosd(90 - phi), digits=3)

    return x, y, z
end

"""
    cart2sph(x, y, z)

Convert spherical coordinates to Cartesian.

# Arguments

- `x::Real`
- `y::Real`
- `z::Real`

# Returns

- `radius::Float64`: spherical radius, the distance from the origin to the point
- `theta::Float64`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
- `phi::Float64`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees
"""
function cart2sph(x::Real, y::Real, z::Real)
    
    # radius = sqrt(x^2 + y^2 + z^2)    
    radius = round(hypot(x, y, z), digits=3)
    # theta = tan^-1(x, y)
    theta = round(atand(y, x), digits=3)
    # phi = cos^-1(z / sqrt(x^2 + y^2 + z^2))
    phi = 90 - round(acosd(z / radius), digits=3)

    # make theta positive - should not be required
    # theta < 0 && (theta = 360 - abs(theta))
    
    return radius, theta, phi
end

"""
    sph2pol(radius, theta, phi)

Convert spherical coordinates to polar.

# Arguments

- `radius::Real`: spherical radius, the distance from the origin to the point
- `theta::Real`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
- `phi::Real`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees

# Returns

- `radius::Float64`
- `theta::Float64`
"""
function sph2pol(radius::Real, theta::Real, phi::Real)
    radius = round(radius * abs(cosd(phi)), digits=3)
    return radius, theta
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
