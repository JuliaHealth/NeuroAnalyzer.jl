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
function add_locs(obj::NeuroAnalyzer.NEURO; locs::DataFrame)

    f_labels = lowercase.(locs[!, :label])

    e_labels = lowercase.(obj.header.recording[:label])
    no_match = setdiff(e_labels, f_labels)

    length(no_match) > 0 && _warn("Location$(_pl(no_match)): $(uppercase.(no_match)) could not be found in the LOCS object.")

    labels_idx = zeros(Int64, length(e_labels))
    for idx1 in eachindex(e_labels)
        for idx2 in eachindex(f_labels)
            e_labels[idx1] == lowercase.(f_labels)[idx2] && (labels_idx[idx1] = idx2)
        end
    end
    for idx in length(labels_idx):-1:1
        labels_idx[idx] == 0 && deleteat!(labels_idx, idx)
    end

    # create new dataset
    obj_new = deepcopy(obj)
    for idx in labels_idx
        l_idx = findfirst(e_labels .== lowercase.(f_labels)[idx])
        obj_new.locs[l_idx, :loc_radius] = locs[idx, :loc_radius]
        obj_new.locs[l_idx, :loc_theta] = locs[idx, :loc_theta]
        obj_new.locs[l_idx, :loc_x] = locs[idx, :loc_x]
        obj_new.locs[l_idx, :loc_y] = locs[idx, :loc_y]
        obj_new.locs[l_idx, :loc_z] = locs[idx, :loc_z]
        obj_new.locs[l_idx, :loc_radius_sph] = locs[idx, :loc_radius_sph]
        obj_new.locs[l_idx, :loc_theta_sph] = locs[idx, :loc_theta_sph]
        obj_new.locs[l_idx, :loc_phi_sph] = locs[idx, :loc_phi_sph]
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
"""
function add_locs!(obj::NeuroAnalyzer.NEURO; locs::DataFrame)

    obj_new = add_locs(obj, locs=locs)
    obj.history = obj_new.history
    obj.locs = obj_new.locs

    return nothing

end
