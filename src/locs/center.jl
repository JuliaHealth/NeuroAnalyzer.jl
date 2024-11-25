export locs_center
export locs_center!
export locs_origin
export locs_origin!

"""
    locs_center(locs; <keyword arguments>)

Center locs at (0, 0).

# Arguments

- `locs::DataFrame`
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

- `locs_new::DataFrame`
"""
function locs_center(locs::DataFrame; polar::Bool=true, cart::Bool=true, spherical::Bool=true)::DataFrame

    locs_new = deepcopy(locs)
    # locs_new = locs_rotz(locs, a=90)
    # search for central line channels (Fz, Cz, Pz)
    cl = nothing
    x_offset = 0
    findfirst(lowercase.(locs[!, :label]) .== "fz") !== nothing && (cl = findfirst(lowercase.(locs[!, :label]) .== "fz"))
    findfirst(lowercase.(locs[!, :label]) .== "cz") !== nothing && (cl = findfirst(lowercase.(locs[!, :label]) .== "cz"))
    findfirst(lowercase.(locs[!, :label]) .== "pz") !== nothing && (cl = findfirst(lowercase.(locs[!, :label]) .== "pz"))
    @assert cl !== nothing "Central line channels could not be find."
    x_offset = locs[cl, :loc_x]

    cl = nothing
    y_offset = 0
    findfirst(lowercase.(locs[!, :label]) .== "cz") !== nothing && (cl = findfirst(lowercase.(locs[!, :label]) .== "cz"))
    findfirst(lowercase.(locs[!, :label]) .== "c1") !== nothing && (cl = findfirst(lowercase.(locs[!, :label]) .== "c1"))
    findfirst(lowercase.(locs[!, :label]) .== "c2") !== nothing && (cl = findfirst(lowercase.(locs[!, :label]) .== "c2"))
    cl === nothing && _warn("Cz/C1/C2 electrodes could not be find, cannot center along Y axis")
    y_offset = locs[cl, :loc_y]

    if cart
        locs_new[!, :loc_x] = locs[!, :loc_x] .- x_offset
        locs_new[!, :loc_y] = locs[!, :loc_y] .- y_offset
    end

    if spherical
        locs_tmp = deepcopy(locs)
        locs_sph2cart!(locs_tmp)
        locs_new[!, :loc_x] = locs[!, :loc_x] .- x_offset
        locs_new[!, :loc_y] = locs[!, :loc_y] .- y_offset
        locs_cart2sph!(locs_tmp)
        locs_new[!, :loc_radius_sph] = locs_tmp[!, :loc_radius_sph]
        locs_new[!, :loc_theta_sph] = locs_tmp[!, :loc_theta_sph]
        locs_new[!, :loc_phi_sph] = locs_tmp[!, :loc_phi_sph]
    end

    polar && locs_rotz!(locs_new, a=90, polar=true, cart=false, spherical=false)

    _locs_round!(locs_new)
    _locs_remove_nans!(locs_new)

    return locs_new

end

"""
    locs_center!(locs; <keyword arguments>)

Center locs at X=0.

# Arguments

- `locs::DataFrame`
- `polar::Bool=true`: modify polar coordinates
- `cart::Bool=true`: modify Cartesian coordinates
- `spherical::Bool=true`: modify spherical coordinates

# Returns

Nothing
"""
function locs_center!(locs::DataFrame; polar::Bool=true, cart::Bool=true, spherical::Bool=true)::Nothing

    locs[!, :] = locs_center(locs, polar=polar, cart=cart, spherical=spherical)[!, :]

    return nothing

end

"""
    locs_origin(locs; <keyword arguments>)

Move origin ([0, 0, 0]) along the axes.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `x::Real=0`: move origin along X axis
- `y::Real=0`: move origin along Y axis
- `z::Real=0`: move origin along Z axis

# Returns

- `locs_new::DataFrame`
"""
function locs_origin(locs::DataFrame; x::Real=0, y::Real=0, z::Real=0)::DataFrame

    locs_new = deepcopy(locs)
    locs_new[:, :loc_x] .+= x
    locs_new[:, :loc_y] .+= y
    locs_new[:, :loc_z] .+= z
    locs_cart2sph!(locs_new)
    locs_cart2pol!(locs_new)

    return locs_new

end

"""
    locs_origin!(locs; <keyword arguments>)

Move origin ([0, 0, 0]) along the axes.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `x::Real=0`: move origin along X axis
- `y::Real=0`: move origin along Y axis
- `z::Real=0`: move origin along Z axis

# Returns

Nothing
"""
function locs_origin!(obj::NeuroAnalyzer.NEURO; x::Real=0, y::Real=0, z::Real=0)::Nothing

    obj.locs = locs_origin(obj.locs, x=x, y=y, z=z)

    return nothing

end

"""
    locs_origin!(locs; <keyword arguments>)

Move origin ([0, 0, 0]) along the axes.

# Arguments

- `locs::DataFrame`
- `x::Real=0`: move origin along X axis
- `y::Real=0`: move origin along Y axis
- `z::Real=0`: move origin along Z axis

# Returns

Nothing
"""
function locs_origin!(locs::DataFrame; x::Real=0, y::Real=0, z::Real=0)::Nothing

    locs[!, :] = locs_origin(locs, x=x, y=y, z=z)[!, :]

    return nothing

end
