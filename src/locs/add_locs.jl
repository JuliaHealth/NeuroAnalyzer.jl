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

- `obj::NeuroAnalyzer.NEURO`
- `locs::DataFrame`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function add_locs(obj::NeuroAnalyzer.NEURO; locs::DataFrame)::NeuroAnalyzer.NEURO

    no_match = setdiff(labels(obj), locs[!, :label])
    length(no_match) > 0 && _warn("Location$(_pl(no_match)): $(uppercase.(no_match)) could not be found in the LOCS object.")
    locs = Base.filter(:label => in(labels(obj)), locs)
    # create new dataset
    obj_new = deepcopy(obj)
    for idx in 1:nrow(locs)
        lidx = findfirst(isequal(locs[idx, :label]), obj_new.locs[!, :label])
        isa(lidx, Int64) && (obj_new.locs[lidx, :] = locs[idx, :])
    end

    # add entry to :history field
    push!(obj_new.history, "add_locs(OBJ, locs)")

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

- `obj::NeuroAnalyzer.NEURO`
- `locs::DataFrame`

# Returns

Nothing
"""
function add_locs!(obj::NeuroAnalyzer.NEURO; locs::DataFrame)::Nothing

    obj_new = add_locs(obj, locs=locs)
    obj.history = obj_new.history
    obj.locs = obj_new.locs

    return nothing

end
