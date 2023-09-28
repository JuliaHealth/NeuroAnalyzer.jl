export add_locs
export add_locs!

"""
    add_locs(obj; locs)

Add electrode positions from `locs`. 

Electrode locations:
- `channel`         channel number
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

    f_labels = lowercase.(locs[!, :labels])

    e_labels = lowercase.(obj.header.recording[:labels])
    no_match = setdiff(e_labels, f_labels)

    length(no_match) > 0 && _warn("Labels: $(uppercase.(no_match)) not found in the LOCS object")
    
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
    obj_new.locs = locs

    # add entry to :history field
    push!(obj_new.history, "add_locs(OBJ, locs)")

    return obj_new
    
end

"""
    add_locs!(obj; locs)

Load electrode positions from `locs` and return `NeuroAnalyzer.NEURO` object with metadata: `:channel_locations`, `:loc_theta`, `:loc_radius`, `:loc_x`, `:loc_x`, `:loc_y`, `:loc_radius_sph`, `:loc_theta_sph`, `:loc_phi_sph`. 

Electrode locations:
- `channel`         channel number
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
"""
function add_locs!(obj::NeuroAnalyzer.NEURO; locs::DataFrame)

    obj_new = add_locs(obj, locs=locs)
    obj.header = obj_new.header
    obj.history = obj_new.history
    obj.locs = obj_new.locs

    return nothing

end
