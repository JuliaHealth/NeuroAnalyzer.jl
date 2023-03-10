export locs_details

"""
    locs_details(obj; channel, output)

Return locations of OBJ channel electrode.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, String}`: channel number or name
- `output::Bool=true`: print output if true

# Returns

Named tuple containing:
- `theta::Float64`: polar planar theta coordinate
- `radius::Float64`: polar planar radius coordinate
- `x::Float64`: Cartesian X spherical coordinate
- `y::Float64`: Cartesian Y spherical coordinate
- `z::Float64`: Cartesian Z spherical coordinate
- `theta_sph::Float64`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
- `radius_sph::Float64`: spherical radius, the distance from the origin to the point
- `phi_sph::Float64`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees
"""
function locs_details(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, String}, output::Bool=true)

    obj.header.has_locs == false && throw(ArgumentError("Electrode locations not available, use load_locs() or add_locs() first."))

    channel = _get_ch_idx(labels(obj), channel)

    x = obj.locs[!, :loc_x][channel]
    y = obj.locs[!, :loc_y][channel]
    z = obj.locs[!, :loc_z][channel]
    theta = obj.locs[!, :loc_theta][channel]
    radius = obj.locs[!, :loc_radius][channel]
    theta_sph = obj.locs[!, :loc_theta_sph][channel]
    radius_sph = obj.locs[!, :loc_radius_sph][channel]
    phi_sph = obj.locs[!, :loc_phi_sph][channel]

    if output
        println("Channel: $channel")
        println("  Label: $(labels(obj)[channel])")
        println("  theta: $theta (planar)")
        println(" radius: $radius (planar)")
        println("      X: $x (spherical)")
        println("      Y: $y (spherical)")
        println("      Z: $z (spherical)")
        println(" radius: $radius_sph (spherical)")
        println("  theta: $theta_sph (spherical)")
        println("    phi: $phi_sph (spherical)")
    end
    
    return (theta=theta, radius=radius, x=x, y=y, z=z, theta_sph=theta_sph, radius_sph=radius_sph, phi_sph=phi_sph)
end
