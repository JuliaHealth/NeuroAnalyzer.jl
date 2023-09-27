export cart2pol
export cart2sph
export pol2cart
export pol2sph
export sph2cart
export sph2pol
export locs_cart2pol
export locs_cart2pol!
export locs_cart2sph
export locs_cart2sph!
export locs_pol2cart
export locs_pol2cart!
export locs_pol2sph
export locs_pol2sph!
export locs_sph2cart
export locs_sph2cart!
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

    radius = hypot(x, y)
    theta = atand(y, x)

    return radius, theta

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
- `theta::Float64`: spherical horizontal angle, the angle in the xy plane with respect to the x axis, in degrees
- `phi::Float64`: spherical azimuth angle, the angle with respect to the z axis (elevation), in degrees
"""
function cart2sph(x::Real, y::Real, z::Real)
    
    # radius = sqrt(x^2 + y^2 + z^2)    
    radius = hypot(x, y, z)
    # theta = tan^-1(x, y)
    theta = atand(y, x)
    # phi = cos^-1(z / sqrt(x^2 + y^2 + z^2))
    # phi = round(acosd(z / radius), digits=2)
    phi = 90 - acosd(z / radius)
    
    return radius, theta, phi

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
    
    x = radius * cosd(theta)
    y = radius * sind(theta)

    return x, y

end

"""
    pol2sph(radius, theta)

Convert polar coordinates to spherical.

# Arguments

- `radius::Real`: polar radius, the distance from the origin to the point, in degrees
- `theta::Real`: polar angle

# Returns

- `radius::Float64`: spherical radius, the distance from the origin to the point
- `theta::Float64`: spherical horizontal angle, the angle in the xy plane with respect to the x axis, in degrees
- `phi::Float64`: spherical azimuth angle, the angle with respect to the z axis (elevation), in degrees
"""
function pol2sph(radius::Real, theta::Real)
    
    return radius, theta, 0

end

"""
    sph2cart(radius, theta, phi)

Convert spherical coordinates to Cartesian.

# Arguments

- `radius::Real`: spherical radius, the distance from the origin to the point
- `theta::Real`: spherical horizontal angle, the angle in the xy plane with respect to the x axis, in degrees
- `phi::Real`: spherical azimuth angle, the angle with respect to the z axis (elevation), in degrees

# Returns

- `x::Float64`
- `y::Float64`
- `z::Float64`
"""
function sph2cart(radius::Real, theta::Real, phi::Real)
    
    x = radius * sind(90 - phi) * cosd(theta)
    y = radius * sind(90 - phi) * sind(theta)
    z = radius * cosd(90 - phi)

    return x, y, z

end

"""
    sph2pol(radius, theta, phi)

Convert spherical coordinates to polar.

# Arguments

- `radius::Real`: spherical radius, the distance from the origin to the point
- `theta::Real`: spherical horizontal angle, the angle in the xy plane with respect to the x axis, in degrees
- `phi::Real`: spherical azimuth angle, the angle with respect to the z axis (elevation), in degrees

# Returns

- `radius::Real`: planar radius, the distance from the origin to the point
- `theta::Real`: planar horizontal angle, the angle in the xy plane with respect to the x axis, in degrees
"""
function sph2pol(radius::Real, theta::Real, phi::Real)

    radius = radius * abs(cosd(phi))

    return radius, theta

end

"""
    locs_pol2cart(locs)

Convert polar coordinates to Cartesian.

# Arguments

- `locs::DataFrame`

# Returns

- `locs_new::DataFrame`
"""
function locs_pol2cart(locs::DataFrame)

    locs_new = deepcopy(locs)

    for idx in 1:nrow(locs)
        r = locs_new[idx, :loc_radius]
        t = locs_new[idx, :loc_theta]
        x, y = pol2cart(r, t)
        locs_new[idx, :loc_x] = x
        locs_new[idx, :loc_y] = y
    end

    _locs_norm!(locs_new)
    _locs_round!(locs_new)

    return locs_new

end

"""
    locs_pol2cart!(locs)

Convert polar coordinates to Cartesian.

# Arguments

- `locs::DataFrame`
"""
function locs_pol2cart!(locs::DataFrame)

    locs[!, :] = locs_pol2cart(locs)[!, :]

    return nothing

end

"""
    locs_pol2sph(locs)

Convert polar coordinates to spherical.

# Arguments

- `locs::DataFrame`

# Returns

- `locs_new::DataFrame`
"""
function locs_pol2sph(locs::DataFrame)

    locs_new = deepcopy(locs)

    for idx in 1:nrow(locs)
        r = locs_new[idx, :loc_radius]
        t = locs_new[idx, :loc_theta]
        r, t, p = pol2sph(r, t)
        locs_new[idx, :loc_radius_sph] = r
        locs_new[idx, :loc_theta_sph] = t
        locs_new[idx, :loc_phi_sph] = p
    end

    _locs_round!(locs_new)

    return locs_new

end

"""
    locs_pol2sph!(locs)

Convert polar coordinates to spherical.

# Arguments

- `locs::DataFrame`
"""
function locs_pol2sph!(locs::DataFrame)

    locs[!, :] = locs_pol2sph(locs)[!, :]

    return nothing

end

"""
    locs_sph2cart(locs)

Convert spherical coordinates to Cartesian.

# Arguments

- `locs::DataFrame`

# Returns

- `locs_new::DataFrame`
"""
function locs_sph2cart(locs::DataFrame)

    locs_new = deepcopy(locs)

    for idx in 1:nrow(locs)
        r = locs_new[idx, :loc_radius_sph]
        t = locs_new[idx, :loc_theta_sph]
        p = locs_new[idx, :loc_phi_sph]
        x, y, z = sph2cart(r, t, p)
        locs_new[idx, :loc_x] = x
        locs_new[idx, :loc_y] = y
        locs_new[idx, :loc_z] = z
    end

    # _locs_norm!(locs_new)
    _locs_round!(locs_new)

    return locs_new

end

"""
    locs_sph2cart!(locs)

Convert spherical coordinates to Cartesian.

# Arguments

- `locs::DataFrame`
"""
function locs_sph2cart!(locs::DataFrame)

    locs[!, :] = locs_sph2cart(locs)[!, :]

    return nothing

end

"""
    locs_sph2pol(locs)

Convert spherical coordinates to polar.

# Arguments

- `locs::DataFrame`

# Returns

- `locs_new::DataFrame`
"""
function locs_sph2pol(locs::DataFrame)

    locs_new = deepcopy(locs)

    for idx in 1:nrow(locs)
        r_sph = locs_new[idx, :loc_radius_sph]
        t_sph = locs_new[idx, :loc_theta_sph]
        p_sph = locs_new[idx, :loc_phi_sph]
        r_pol, t_pol = sph2pol(r_sph, t_sph, p_sph)
        locs_new[idx, :loc_radius] = r_pol
        locs_new[idx, :loc_theta] = t_pol
    end

    _locs_round!(locs_new)

    return locs_new

end

"""
    locs_sph2pol!(locs)

Convert Cartesian coordinates to polar.

# Arguments

- `locs::DataFrame`
"""
function locs_sph2pol!(locs::DataFrame)

    locs[!, :] = locs_sph2pol(locs)[!, :]

    return nothing

end

"""
    locs_cart2sph(locs)

Convert Cartesian coordinates to spherical.

# Arguments

- `locs::DataFrame`

# Returns

- `locs_new::DataFrame`
"""
function locs_cart2sph(locs::DataFrame)

    locs_new = deepcopy(locs)

    for idx in 1:nrow(locs)
        x = locs_new[idx, :loc_x]
        y = locs_new[idx, :loc_y]
        z = locs_new[idx, :loc_z]
        r, t, p = cart2sph(x, y, z)
        locs_new[idx, :loc_radius_sph] = r
        locs_new[idx, :loc_theta_sph] = t
        locs_new[idx, :loc_phi_sph] = p
    end

    _locs_round!(locs_new)

    return locs_new
    
end

"""
    locs_cart2sph!(locs)

Convert Cartesian coordinates to spherical.

# Arguments

- `locs::DataFrame`
"""
function locs_cart2sph!(locs::DataFrame)

    locs[!, :] = locs_cart2sph(locs)[!, :]

    return nothing

end

"""
    locs_cart2pol(locs)

Convert Cartesian coordinates to polar.

# Arguments

- `locs::DataFrame`

# Returns

- `locs_new::DataFrame`
"""
function locs_cart2pol(locs::DataFrame)

    locs_new = deepcopy(locs)

    for idx in 1:nrow(locs)
        x = locs_new[idx, :loc_x]
        y = locs_new[idx, :loc_y]
        r, t = cart2pol(x, y)
        locs_new[idx, :loc_radius] = r
        locs_new[idx, :loc_theta] = t
    end

    _locs_round!(locs_new)

    return locs_new

end

"""
    locs_cart2pol!(locs)

Convert Cartesian coordinates to polar.

# Arguments

- `locs::DataFrame`
"""
function locs_cart2pol!(locs::DataFrame)

    locs[!, :] = locs_cart2pol(locs)[!, :]

    return nothing

end
