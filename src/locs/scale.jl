export locs_scale
export locs_scale!
export locs_maximize
export locs_maximize!

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
