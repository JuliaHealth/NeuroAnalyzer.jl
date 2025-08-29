export locs_details

"""
    locs_details(obj; <keyword arguments>)

Return channel location details.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `out::Bool=true`: if true, print details

# Returns

Named tuple containing:
- `label::String`: location label
- `theta_pl::Float64`: polar angle
- `radius_pl::Float64`: polar radius
- `x::Float64`: Cartesian X spherical coordinate
- `y::Float64`: Cartesian Y spherical coordinate
- `z::Float64`: Cartesian Z spherical coordinate
- `theta_sph::Float64`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
- `radius_sph::Float64`: spherical radius, the distance from the origin to the point
- `phi_sph::Float64`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees
"""
function locs_details(obj::NeuroAnalyzer.NEURO; ch::String, out::Bool=true)::@NamedTuple{label::String, theta_pl::Float64, radius_pl::Float64, x::Float64, y::Float64, z::Float64, theta_sph::Float64, radius_sph::Float64, phi_sph::Float64}

    ch = intersect(obj.locs[!, :label], [ch])
    locs = Base.filter(:label => in(ch), obj.locs)
    @assert DataFrames.nrow(locs) == 1 "Channel has no location details."

    l = obj.locs[1, :label]
    x = obj.locs[1, :loc_x]
    y = obj.locs[1, :loc_y]
    z = obj.locs[1, :loc_z]
    theta_pl = obj.locs[1, :loc_theta]
    radius_pl = obj.locs[1, :loc_radius]
    theta_sph = obj.locs[1, :loc_theta_sph]
    radius_sph = obj.locs[1, :loc_radius_sph]
    phi_sph = obj.locs[1, :loc_phi_sph]

    if out
        println("  Label: $l")
        println("  Theta: $theta_pl (polar)")
        println(" Radius: $radius_pl (polar)")
        println("      X: $x (spherical)")
        println("      Y: $y (spherical)")
        println("      Z: $z (spherical)")
        println(" Radius: $radius_sph (spherical)")
        println("  Theta: $theta_sph (spherical)")
        println("    Phi: $phi_sph (spherical)")
    end

    return (label=l, theta_pl=theta_pl, radius_pl=radius_pl, x=x, y=y, z=z, theta_sph=theta_sph, radius_sph=radius_sph, phi_sph=phi_sph)

end
