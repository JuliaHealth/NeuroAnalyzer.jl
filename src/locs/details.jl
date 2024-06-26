export locs_details

"""
    locs_details(obj; ch, output)

Return locations of OBJ ch electrode.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, String}`: channel number or name
- `out::Bool=true`: print output if true

# Returns

Named tuple containing:
- `ch::Int64`: channel number
- `label::String`: location label
- `theta::Float64`: polar angle
- `radius::Float64`: polar radius
- `x::Float64`: Cartesian X spherical coordinate
- `y::Float64`: Cartesian Y spherical coordinate
- `z::Float64`: Cartesian Z spherical coordinate
- `theta_sph::Float64`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
- `radius_sph::Float64`: spherical radius, the distance from the origin to the point
- `phi_sph::Float64`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees
"""
function locs_details(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, String}, out::Bool=true)

    @assert _has_locs(obj) "Electrode locations not available, use load_locs() or add_locs() first."

    ch = _get_ch_idx(labels(obj), ch)

    x = obj.locs[ch, :loc_x]
    y = obj.locs[ch, :loc_y]
    z = obj.locs[ch, :loc_z]
    theta_pl = obj.locs[ch, :loc_theta]
    radius_pl = obj.locs[ch, :loc_radius]
    theta_sph = obj.locs[ch, :loc_theta_sph]
    radius_sph = obj.locs[ch, :loc_radius_sph]
    phi_sph = obj.locs[ch, :loc_phi_sph]

    # convert InlineString to String
    l = ""
    for idx in eachindex(obj.locs[ch, :labels])
        l *= obj.locs[ch, :labels][idx]
    end

    if out
        println("Channel: $ch")
        println("  Label: $l")
        println("  theta: $theta_pl (polar)")
        println(" radius: $radius_pl (polar)")
        println("      X: $x (spherical)")
        println("      Y: $y (spherical)")
        println("      Z: $z (spherical)")
        println(" radius: $radius_sph (spherical)")
        println("  theta: $theta_sph (spherical)")
        println("    phi: $phi_sph (spherical)")
    end

    return (ch=ch, label=l, theta_pl=theta_pl, radius_pl=radius_pl, x=x, y=y, z=z, theta_sph=theta_sph, radius_sph=radius_sph, phi_sph=phi_sph)

end
