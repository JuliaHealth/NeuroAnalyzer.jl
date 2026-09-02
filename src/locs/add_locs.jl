export add_locs
export add_locs!

"""
    add_locs(obj; <keyword arguments>)

Add electrode positions from `locs`.

Electrode locations:

- `labels`          channel label
- `loc_theta`       polar angle
- `loc_radius`      polar radius
- `loc_x`           Cartesian x
- `loc_y`           Cartesian y
- `loc_z`           Cartesian z
- `loc_radius_sph`  spherical radius
- `loc_theta_sph`   spherical horizontal angle
- `loc_phi_sph`     spherical azimuth angle

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `locs::DataFrame`: channel location data

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function add_locs(obj::NeuroAnalyzer.NEURO; locs::DataFrame)::NeuroAnalyzer.NEURO
    no_match = setdiff(labels(obj), locs[!, :label])
    length(no_match) > 0 && _warn(
        "Location$(_pl(no_match)): $(uppercase.(no_match)) could not be found in the LOCS object.",
    )
    locs = Base.filter(:label => in(labels(obj)), locs)

    # create new dataset
    obj_tmp = deepcopy(obj)

    for idx in 1:DataFrames.nrow(locs)
        lidx = findfirst(isequal(locs[idx, :label]), obj_tmp.locs[!, :label])
        isa(lidx, Int64) && (obj_tmp.locs[lidx, :] = locs[idx, :])
    end

    # keep order consistent with labels
    locs_idx = indexin(obj_tmp.locs[:, :label], labels(obj_tmp))
    obj_tmp.locs = obj_tmp.locs[sortperm(locs_idx), :]

    push!(obj_tmp.history, "add_locs(obj, locs)")

    return obj_tmp
end

"""
    add_locs!(obj; <keyword arguments>)

Load electrode positions from `locs` and return `NeuroAnalyzer.NEURO` object attached with channel locations data.

Electrode locations:

- `labels`: channel label
- `loc_theta`: polar angle
- `loc_radius`: polar radius
- `loc_x`: Cartesian X
- `loc_y`: Cartesian Y
- `loc_z`: Cartesian Z
- `loc_radius_sph`: spherical radius
- `loc_theta_sph`: spherical horizontal angle
- `loc_phi_sph`: spherical azimuth angle

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `locs::DataFrame`: channel location data

# Returns

- `Nothing`
"""
function add_locs!(obj::NeuroAnalyzer.NEURO; locs::DataFrame)::Nothing
    obj_tmp = add_locs(obj; locs = locs)
    obj.history = obj_tmp.history
    obj.locs = obj_tmp.locs
    obj_tmp = nothing

    return nothing
end
