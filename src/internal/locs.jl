function _locnorm(x::Union{AbstractVector, Real}, y::Union{AbstractVector, Real})
    xy = s_normalize_minmax(hcat(x, y))
    x = xy[:, 1]
    y = xy[:, 2]
    return x, y
end

function _locnorm(x::Union{AbstractVector, Real}, y::Union{AbstractVector, Real}, z::Union{AbstractVector, Real})
    xyz = s_normalize_minmax(hcat(x, y, z))
    x = xyz[:, 1]
    y = xyz[:, 2]
    z = xyz[:, 3]
    return x, y, z
end

function _round_locs(locs::DataFrame)
    locs[!, :loc_x] = round.(locs[!, :loc_x], digits=3)
    locs[!, :loc_y] = round.(locs[!, :loc_y], digits=3)
    locs[!, :loc_z] = round.(locs[!, :loc_z], digits=3)
    locs[!, :loc_radius] = round.(locs[!, :loc_radius], digits=3)
    locs[!, :loc_theta] = round.(locs[!, :loc_theta], digits=3)
    locs[!, :loc_radius_sph] = round.(locs[!, :loc_radius_sph], digits=3)
    locs[!, :loc_theta_sph] = round.(locs[!, :loc_theta_sph], digits=3)
    locs[!, :loc_phi_sph] = round.(locs[!, :loc_phi_sph], digits=3)
    return locs
end
