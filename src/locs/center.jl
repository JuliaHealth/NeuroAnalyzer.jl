export locs_center
export locs_center!

"""
    locs_center(locs; <keyword arguments>)

Center channel locations at `(0, 0)`

# Arguments

- `locs::DataFrame`: channel location data
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `DataFrame`: modified channel location data
"""
function locs_center(
    locs::DataFrame;
    polar::Bool = true,
    cart::Bool = true,
    spherical::Bool = true,
)::DataFrame
    locs_new = deepcopy(locs)
    lab = lowercase.(locs[!, :label])
    # row of the first channel of `names` present in locs (in order of preference)
    pick(names) = (n = findfirst(in(lab), names); isnothing(n) ? nothing : findfirst(==(names[n]), lab))

    # x offset: a midline channel (Cz preferred, then Fz, Pz)
    cl = pick(["cz", "fz", "pz"])
    isnothing(cl) && throw(ArgumentError("Central line channels (Cz, Fz, Pz) could not be found."))
    x_offset = locs[cl, :loc_x]

    # y offset: a channel of the central coronal line (Cz preferred, then C1, C2)
    cl = pick(["cz", "c1", "c2"])
    y_offset = if isnothing(cl)
        _warn("Cz/C1/C2 electrodes could not be found, cannot center along Y axis.")
        0.0
    else
        locs[cl, :loc_y]
    end

    if cart
        locs_new[!, :loc_x] = locs[!, :loc_x] .- x_offset
        locs_new[!, :loc_y] = locs[!, :loc_y] .- y_offset
    end

    if spherical
        # shift the Cartesian coordinates derived from the spherical ones (unrounded)
        r, t, p = _sph_transform(locs, (x, y, z) -> (x .- x_offset, y .- y_offset, z))
        locs_new[!, :loc_radius_sph], locs_new[!, :loc_theta_sph], locs_new[!, :loc_phi_sph] = r, t, p
    end

    if polar
        # polar coordinates of the shifted Cartesian coordinates
        locs_tmp = deepcopy(locs)
        locs_tmp[!, :loc_x] .-= x_offset
        locs_tmp[!, :loc_y] .-= y_offset
        locs_cart2pol!(locs_tmp)
        locs_new[!, :loc_radius] = locs_tmp[!, :loc_radius]
        locs_new[!, :loc_theta] = locs_tmp[!, :loc_theta]
    end

    _locs_round!(locs_new)
    _locs_remove_nans!(locs_new)

    return locs_new
end

"""
    locs_center!(locs; <keyword arguments>)

Center channel locations at `(0, 0)`.

# Arguments

- `locs::DataFrame`: channel location data
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `Nothing`
"""
function locs_center!(
    locs::DataFrame;
    polar::Bool = true,
    cart::Bool = true,
    spherical::Bool = true,
)::Nothing
    locs[!, :] = locs_center(
        locs;
        polar = polar,
        cart = cart,
        spherical = spherical,
    )[!, :]

    return nothing
end
