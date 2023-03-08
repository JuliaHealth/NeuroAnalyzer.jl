export locs_rotx
export locs_rotx!

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
