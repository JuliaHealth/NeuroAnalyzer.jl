export load_locs
export load_locs!

"""
    load_locs(obj; <keyword arguments>)

Load channel locations from `file_name` and return `NeuroAnalyzer.NEURO` object with `locs` data frame.

Accepted formats:
- CED
- LOCS
- ELC
- TSV
- SFP
- CSD
- GEO
- MAT
- TXT
- DAT
- ASC

Channel locations:

- `loc_theta`: polar angle
- `loc_radius`: polar radius
- `loc_x`: spherical Cartesian x
- `loc_y`: spherical Cartesian y
- `loc_z`: spherical Cartesian z
- `loc_radius_sph`: spherical radius
- `loc_theta_sph`: spherical horizontal angle
- `loc_phi_sph`: spherical azimuth angle

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `file_name::String`: name of the file to load

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function load_locs(obj::NeuroAnalyzer.NEURO; file_name::String)

    @assert isfile(file_name) "File $file_name cannot be loaded."
    @assert length(obj.header.recording[:label]) > 0 "OBJ does not contain labels, use add_label() first."

    _info("Send standard locations for your channels to adam.wysokinski@neuroanalyzer.org")
    _info("Nose direction is set at '+Y'")

    if splitext(file_name)[2] == ".ced"
        locs = import_locs_ced(file_name)
    elseif splitext(file_name)[2] == ".elc"
        locs = import_locs_elc(file_name)
    elseif splitext(file_name)[2] == ".locs"
        locs = import_locs_locs(file_name)
    elseif splitext(file_name)[2] == ".tsv"
        locs = import_locs_tsv(file_name)
    elseif splitext(file_name)[2] == ".sfp"
        locs = import_locs_sfp(file_name)
    elseif splitext(file_name)[2] == ".csd"
        locs = import_locs_csd(file_name)
    elseif splitext(file_name)[2] == ".geo"
        locs = import_locs_geo(file_name)
    elseif splitext(file_name)[2] == ".mat"
        locs = import_locs_mat(file_name)
    elseif splitext(file_name)[2] == ".txt"
        locs = import_locs_txt(file_name)
    elseif splitext(file_name)[2] == ".dat"
        locs = import_locs_dat(file_name)
    elseif splitext(file_name)[2] == ".asc"
        locs = import_locs_asc(file_name)
    else
        @error "Unknown file format."
    end

    # add locations of reference channels
    ref_labels = get_channel(obj, type="ref")
    if length(ref_labels) > 0
        for idx in eachindex(ref_labels)
            (occursin("a", lowercase(ref_labels[idx])) && occursin("1", ref_labels[idx])) && push!(locs, [ref_labels[idx], 1.0, 192.0, -0.92, -0.23, -0.55, 1.10, -165.96, -30.11])
            (occursin("a", lowercase(ref_labels[idx])) && occursin("2", ref_labels[idx])) && push!(locs, [ref_labels[idx], 1.0, -12.0, 0.92, -0.23, -0.55, 1.10, -14.04, -30.11])
            (occursin("m", lowercase(ref_labels[idx])) && occursin("1", ref_labels[idx])) && push!(locs, [ref_labels[idx], 0.95, -173.93, -0.94, -0.1, -0.3, 0.99, -173.93, -17.61])
            (occursin("m", lowercase(ref_labels[idx])) && occursin("2", ref_labels[idx])) && push!(locs, [ref_labels[idx], 0.95, -6.07, 0.94, -0.1, -0.3, 0.99, -6.07, -17.61])
        end
    end

    # add locations of EMG channels
    emg_labels = get_channel(obj, type="emg")
    if length(emg_labels) > 0
        for idx in eachindex(emg_labels)
            occursin("1", emg_labels[idx]) && push!(locs, [emg_labels[idx], 0.99, 135.0, -0.70, 0.70, -1.10, 1.48, 135.00, -48.01])
            occursin("2", emg_labels[idx]) && push!(locs, [emg_labels[idx], 0.99, 45.0, 0.70, 0.70, -1.10, 1.48, 45.00, -48.01])
            # if no numbers, assume that EMG channel is on the right side
            (!occursin("1", emg_labels[idx]) && !occursin("2", emg_labels[idx])) && push!(locs, [emg_labels[idx], 0.99, 45.0, 0.70, 0.70, -1.10, 1.48, 45.00, -48.01])
        end
    end

    # add locations of EOG channels
    eog_labels = get_channel(obj, type="eog")
    if length(eog_labels) > 0
        for idx in eachindex(eog_labels)
            if occursin("1", eog_labels[idx])
                occursin("v", lowercase(eog_labels[idx])) && push!(locs, [eog_labels[idx], 1.01, 149.62, -0.87, 0.51, -0.37, 1.07, 149.62, -20.15])
                occursin("h", lowercase(eog_labels[idx])) && push!(locs, [eog_labels[idx], 1.0, 129.73, -0.64, 0.77, -0.04, 1.00, 129.73, -2.29])
            end
            if occursin("2", eog_labels[idx])
                occursin("v", lowercase(eog_labels[idx])) && push!(locs, [eog_labels[idx], 1.01, 30.38, 0.87, 0.51, -0.37, 1.07, 30.38, -20.15])
                occursin("h", lowercase(eog_labels[idx])) && push!(locs, [eog_labels[idx], 1.0, 50.0, 0.64, 0.77, -0.04, 1.00, 50.27, -2.29])
            end
            # if no numbers, assume that EOG channels are on the right side
            if !occursin("1", eog_labels[idx]) && !occursin("2", eog_labels[idx])
                occursin("v", lowercase(eog_labels[idx])) && push!(locs, [eog_labels[idx], 1.01, 30.38, 0.87, 0.51, -0.37, 1.07, 30.38, -20.15])
                occursin("h", lowercase(eog_labels[idx])) && push!(locs, [eog_labels[idx], 1.0, 50.0, 0.64, 0.77, -0.04, 1.00, 50.27, -2.29])
            end
            # if no V/H indicators, assume that EEG channels are vertical
            if !occursin("v", lowercase(eog_labels[idx])) && !occursin("h", lowercase(eog_labels[idx]))
                occursin("1", eog_labels[idx]) && push!(locs, [eog_labels[idx], 1.01, 149.62, -0.87, 0.51, -0.37, 1.07, 149.62, -20.15])
                occursin("2", eog_labels[idx]) && push!(locs, [eog_labels[idx], 1.0, 35.0, 0.87, 0.51, -0.37, 1.03, 33.52, -21.08])
            end
        end
    end

    no_match = setdiff(labels(obj), locs[!, :label])
    length(no_match) > 0 && _warn("Location$(_pl(no_match)): $(uppercase.(no_match)) could not be found in $file_name")

    # create new dataset
    obj_new = deepcopy(obj)
    obj_new.locs = Base.filter(:label => in(labels(obj)), locs)

    _locs_round!(obj_new.locs)
    _locs_remove_nans!(obj_new.locs)

    push!(obj_new.history, "load_locs(OBJ, file_name=$file_name)")

    return obj_new

end

"""
    load_locs!(obj; <keyword arguments>)

Load channel locations from `file_name` and return `NeuroAnalyzer.NEURO` object with `locs` data frame.

Accepted formats:
- CED
- LOCS
- ELC
- TSV
- SFP
- CSD
- GEO
- MAT
- TXT
- DAT
- ASC

Channel locations:

- `loc_theta`: polar angle
- `loc_radius`: polar radius
- `loc_x`: spherical Cartesian x
- `loc_y`: spherical Cartesian y
- `loc_z`: spherical Cartesian z
- `loc_radius_sph`: spherical radius
- `loc_theta_sph`: spherical horizontal angle
- `loc_phi_sph`: spherical azimuth angle

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `file_name::String`
"""
function load_locs!(obj::NeuroAnalyzer.NEURO; file_name::String)

    obj_tmp = load_locs(obj, file_name=file_name)
    obj.locs = obj_tmp.locs

    return nothing

 end
