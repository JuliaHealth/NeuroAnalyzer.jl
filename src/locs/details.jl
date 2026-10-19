export locs_details

"""
    locs_details(obj; <keyword arguments>)

Return channel location details.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `out::Bool=true`: if true, print details

# Returns

Named tuple:

- `label::String`: location label
- `theta_pl::Float64`: polar angle
- `radius_pl::Float64`: polar radius
- `x::Float64`: Cartesian x spherical coordinate
- `y::Float64`: Cartesian y spherical coordinate
- `z::Float64`: Cartesian z spherical coordinate
- `theta_sph::Float64`: spherical horizontal angle, the angle in the xy plane with respect to the x-axis, in degrees
- `radius_sph::Float64`: spherical radius, the distance from the origin to the point
- `phi_sph::Float64`: spherical azimuth angle, the angle with respect to the z-axis (elevation), in degrees
"""
function locs_details(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    out::Bool = true,
)::@NamedTuple{
    label::String,
    theta_pl::Float64,
    radius_pl::Float64,
    x::Float64,
    y::Float64,
    z::Float64,
    theta_sph::Float64,
    radius_sph::Float64,
    phi_sph::Float64,
}
    # validate
    length(get_channel(obj; ch = ch)) == 1 ||
        throw(ArgumentError("ch must resolve to exactly one channel."))

    # location row of the channel (labels are matched case-insensitively)
    loc_idx = _find_bylabel(obj.locs, ch)
    loc_idx isa Int64 || throw(ArgumentError("Channel has no location details."))
    loc = obj.locs[loc_idx, :]

    label = loc[:label]
    x = loc[:loc_x]
    y = loc[:loc_y]
    z = loc[:loc_z]
    theta_pl = loc[:loc_theta]
    radius_pl = loc[:loc_radius]
    theta_sph = loc[:loc_theta_sph]
    radius_sph = loc[:loc_radius_sph]
    phi_sph = loc[:loc_phi_sph]

    if out
        # values are rounded for display only; the returned values are exact
        r2(v) = round(v; digits = 2)
        println("  Label: $label")
        println("  Theta: $(r2(theta_pl)) (polar)")
        println(" Radius: $(r2(radius_pl)) (polar)")
        println("      X: $(r2(x)) (Cartesian)")
        println("      Y: $(r2(y)) (Cartesian)")
        println("      Z: $(r2(z)) (Cartesian)")
        println(" Radius: $(r2(radius_sph)) (spherical)")
        println("  Theta: $(r2(theta_sph)) (spherical)")
        println("    Phi: $(r2(phi_sph)) (spherical)")
    end

    return (; label, theta_pl, radius_pl, x, y, z, theta_sph, radius_sph, phi_sph)
end
