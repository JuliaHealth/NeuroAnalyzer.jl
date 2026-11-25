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

Locations are stored as given, in their original units.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `locs::DataFrame`: channel location data

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function add_locs(obj::NeuroAnalyzer.NEURO; locs::DataFrame)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    add_locs!(obj_new; locs = locs)

    return obj_new
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
    no_match = setdiff(labels(obj), locs[!, :label])
    length(no_match) > 0 && _warn(
        "Location$(_pl(no_match)): $(uppercase.(no_match)) could not be found in the LOCS object.",
    )
    locs = Base.filter(:label => in(labels(obj)), locs)

    # update channels that already have a location, add the others
    for idx = 1:DataFrames.nrow(locs)
        lidx = findfirst(isequal(locs[idx, :label]), obj.locs[!, :label])
        if isa(lidx, Int64)
            obj.locs[lidx, :] = locs[idx, :]
        else
            push!(obj.locs, locs[idx, :]; cols = :intersect)
        end
    end

    # keep order consistent with labels
    locs_idx = indexin(obj.locs[:, :label], labels(obj))
    obj.locs = obj.locs[sortperm(locs_idx), :]

    _history!(obj, "add_locs", (; locs))

    return nothing
end
