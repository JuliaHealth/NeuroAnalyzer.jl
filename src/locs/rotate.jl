export locs_rotz
export locs_rotz!
export locs_roty
export locs_roty!
export locs_rotx
export locs_rotx!

"""
    locs_rotz(locs; a, planar, spherical)

Rotate channel locations around the z axis.

# Arguments

- `locs::DataFrame`
- `a::Real`: angle of rotation (in degrees)
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function locs_rotz(locs::DataFrame; a::Real, planar::Bool=true, spherical::Bool=true)

    locs_new = deepcopy(locs)

    for idx in 1:nrow(locs)
        locs_new[idx, :loc_x] = locs[idx, :loc_x] * cosd(a) - locs[idx, :loc_y] * sind(a)
        locs_new[idx, :loc_y] = locs[idx, :loc_x] * sind(a) + locs[idx, :loc_y] * cosd(a)
    end

    planar == true && locs_cart2pol!(locs_new)
    spherical == true && locs_cart2sph!(locs_new)

    return locs_new

end

"""
    locs_rotz!(locs)

Rotate channel locations in the xy-plane.

# Arguments

- `locs::DataFrame`
- `a::Int64`: scaling factor
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates
"""
function locs_rotz!(locs::DataFrame; a::Real, planar::Bool=true, spherical::Bool=true)

    locs[!, :] = locs_rotz(locs, a=a, planar=planar, spherical=spherical)[!, :]

    return nothing
    
end

"""
    locs_roty(locs; a, planar, spherical)

Rotate channel locations in the xz-plane.

# Arguments

- `locs::DataFrame`
- `a::Real`: angle of rotation (in degrees)
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function locs_roty(locs::DataFrame; a::Real, planar::Bool=true, spherical::Bool=true)

    locs_new = deepcopy(locs)

    for idx in 1:nrow(locs)
        locs_new[idx, :loc_x] = locs[idx, :loc_x] * cosd(a) + locs[idx, :loc_z] * sind(a)
        locs_new[idx, :loc_z] = -locs[idx, :loc_x] * sind(a) + locs[idx, :loc_z] * cosd(a)
    end

    planar == true && locs_cart2pol!(locs_new)
    spherical == true && locs_cart2sph!(locs_new)

    return locs_new

end

"""
    locs_roty!(locs)

Rotate channel locations in the xz-plane.

# Arguments

- `locs::DataFrame`
- `a::Int64`: scaling factor
- `planar::Bool=true`: modify planar coordinates
- `spherical::Bool=true`: modify spherical coordinates
"""
function locs_roty!(locs::DataFrame; a::Real, planar::Bool=true, spherical::Bool=true)

    locs[!, :] = locs_roty(locs, a=a, planar=planar, spherical=spherical)[!, :]

    return nothing
    
end

"""
    locs_rotx(locs; a, planar, spherical)

Rotate channel locations in the yz-plane.

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

    for idx in 1:nrow(locs)
        locs_new[idx, :loc_y] = locs[idx, :loc_y] * cosd(a) - locs[idx, :loc_z] * sind(a)
        locs_new[idx, :loc_z] = locs[idx, :loc_y] * sind(a) + locs[idx, :loc_z] * cosd(a)
    end

    planar == true && locs_cart2pol!(locs_new)
    spherical == true && locs_cart2sph!(locs_new)

    return locs_new

end

"""
    locs_rotx!(locs)

Rotate channel locations in the yz-plane.

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
