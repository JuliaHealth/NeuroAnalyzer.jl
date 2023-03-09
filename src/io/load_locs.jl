export load_locs
export load_locs!

"""
    load_locs(obj; file_name)

Load electrode positions from `file_name` and return `NeuroAnalyzer.NEURO` object with `locs` data frame. 

Accepted formats:
- CED
- LOCS
- ELC
- TSV
- SFP
- CSD
- GEO
- MAT

Electrode locations:
- `loc_theta`       planar polar angle
- `loc_radius`      planar polar radius
- `loc_x`           spherical Cartesian x
- `loc_y`           spherical Cartesian y
- `loc_z`           spherical Cartesian z
- `loc_radius_sph`  spherical radius
- `loc_theta_sph`   spherical horizontal angle
- `loc_phi_sph`     spherical azimuth angle

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `file_name::String`
- `maximize::Bool=true`: maximize locations after importing

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function load_locs(obj::NeuroAnalyzer.NEURO; file_name::String, maximize::Bool=true)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    length(obj.header[:labels]) > 0 || throw(ArgumentError("OBJ does not contain labels, use add_labels() first."))

    if splitext(file_name)[2] == ".ced"
        locs = locs_import_ced(file_name)
    elseif splitext(file_name)[2] == ".elc"
        locs = locs_import_elc(file_name)
    elseif splitext(file_name)[2] == ".locs"
        locs = locs_import_locs(file_name)
    elseif splitext(file_name)[2] == ".tsv"
        locs = locs_import_tsv(file_name)
    elseif splitext(file_name)[2] == ".sfp"
        locs = locs_import_sfp(file_name)
    elseif splitext(file_name)[2] == ".csd"
        locs = locs_import_csd(file_name)
    elseif splitext(file_name)[2] == ".geo"
        locs = locs_import_geo(file_name)
    elseif splitext(file_name)[2] == ".mat"
        locs = locs_import_mat(file_name)
    else
        throw(ArgumentError("Unknown file format."))
    end

    f_labels = locs[!, :labels]

    loc_theta = float.(locs[!, :loc_theta])
    loc_radius = float.(locs[!, :loc_radius])

    loc_radius_sph = float.(locs[!, :loc_radius_sph])
    loc_theta_sph = float.(locs[!, :loc_theta_sph])
    loc_phi_sph = float.(locs[!, :loc_phi_sph])

    loc_x = float.(locs[!, :loc_x])
    loc_y = float.(locs[!, :loc_y])
    loc_z = float.(locs[!, :loc_z])

    e_labels = lowercase.(obj.header[:labels])
    no_match = setdiff(e_labels, lowercase.(f_labels))
    length(no_match) > 0 && _info("Labels: $(uppercase.(no_match)) not found in $file_name.")

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
    obj_new.header.locations = true
    obj_new.locs = DataFrame(:channel => 1:length(f_labels[labels_idx]),
                                 :labels => f_labels[labels_idx],
                                 :loc_theta => loc_theta[labels_idx],
                                 :loc_radius => loc_radius[labels_idx],
                                 :loc_x => loc_x[labels_idx],
                                 :loc_y => loc_y[labels_idx],
                                 :loc_z => loc_z[labels_idx],
                                 :loc_radius_sph => loc_radius_sph[labels_idx],
                                 :loc_theta_sph => loc_theta_sph[labels_idx],
                                 :loc_phi_sph => loc_phi_sph[labels_idx])

    maximize == true && locs_maximize!(obj_new.locs)

    # add entry to :history field
    push!(obj_new.header.history, "load_locs(OBJ, file_name=$file_name)")

    return obj_new
end

"""
    load_locs!(obj; file_name)

Load electrode positions from `file_name` and return `NeuroAnalyzer.NEURO` object with `locs` data frame. 

Accepted formats:
- CED
- LOCS
- ELC
- TSV
- SFP
- CSD
- GEO
- MAT

Electrode locations:
- `loc_theta`       planar polar angle
- `loc_radius`      planar polar radius
- `loc_x`           spherical Cartesian x
- `loc_y`           spherical Cartesian y
- `loc_z`           spherical Cartesian z
- `loc_radius_sph`  spherical radius
- `loc_theta_sph`   spherical horizontal angle
- `loc_phi_sph`     spherical azimuth angle

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `file_name::String`
- `maximize::Bool=true`: maximize locations after importing
"""
function load_locs!(obj::NeuroAnalyzer.NEURO; file_name::String, maximize::Bool=true)

    obj_tmp = load_locs(obj, file_name=file_name, maximize=maximize)
    obj.locs = obj_tmp.locs
    obj.header.locations = true

    nothing
 end
