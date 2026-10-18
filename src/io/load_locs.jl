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
  - CSV

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

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `file_name::String`: name of the file to load

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function load_locs(obj::NeuroAnalyzer.NEURO; file_name::String)::NeuroAnalyzer.NEURO
    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))
    length(obj.header.recording[:label]) > 0 ||
        throw(ArgumentError("OBJ does not contain labels, use add_label() first."))

    _info(
        "Send standard locations for your channels to adam.wysokinski@neuroanalyzer.org",
    )
    _info("Nose direction is set at '+Y'")

    importers = Dict(
        ".ced" => import_locs_ced, ".elc" => import_locs_elc, ".locs" => import_locs_locs,
        ".tsv" => import_locs_tsv, ".sfp" => import_locs_sfp, ".csd" => import_locs_csd,
        ".geo" => import_locs_geo, ".mat" => import_locs_mat, ".txt" => import_locs_txt,
        ".dat" => import_locs_dat, ".asc" => import_locs_asc, ".csv" => import_locs_csv,
    )
    ext = lowercase(splitext(file_name)[2])
    haskey(importers, ext) || throw(
        ArgumentError(
            "Unknown locations file format \"$ext\"; must be one of: $(Base.join(sort(collect(keys(importers))), ", ")).",
        ),
    )
    locs = importers[ext](file_name)

    # default positions of reference / EMG / EOG channels are added only for
    # channels that the locations file does not already contain
    has_loc(l) = lowercase(l) in lowercase.(locs[!, :label])
    add_loc!(row) = has_loc(row[1]) || push!(locs, row)

    # add locations of reference channels
    ref_labels = get_channel(obj; type = "ref")
    if length(ref_labels) > 0
        for idx in eachindex(ref_labels)
            (
                occursin("a", lowercase(ref_labels[idx])) &&
                occursin("1", ref_labels[idx])
            ) && add_loc!(
                [ref_labels[idx], 1.0, 192.0, -0.92, -0.23, -0.55, 1.1, -165.96, -30.11],
            )
            (
                occursin("a", lowercase(ref_labels[idx])) &&
                occursin("2", ref_labels[idx])
            ) && add_loc!(
                [ref_labels[idx], 1.0, -12.0, 0.92, -0.23, -0.55, 1.1, -14.04, -30.11],
            )
            (
                occursin("m", lowercase(ref_labels[idx])) &&
                occursin("1", ref_labels[idx])
            ) && add_loc!(
                [
                    ref_labels[idx],
                    0.95,
                    -173.93,
                    -0.94,
                    -0.1,
                    -0.3,
                    0.99,
                    -173.93,
                    -17.61,
                ],
            )
            (
                occursin("m", lowercase(ref_labels[idx])) &&
                occursin("2", ref_labels[idx])
            ) && add_loc!(
                [ref_labels[idx], 0.95, -6.07, 0.94, -0.1, -0.3, 0.99, -6.07, -17.61],
            )
        end
    end

    # add locations of EMG channels
    emg_labels = get_channel(obj; type = "emg")
    if length(emg_labels) > 0
        for idx in eachindex(emg_labels)
            occursin("1", emg_labels[idx]) && add_loc!(
                [emg_labels[idx], 0.99, 135.0, -0.7, 0.7, -1.1, 1.48, 135.0, -48.01],
            )
            occursin("2", emg_labels[idx]) && add_loc!(
                [emg_labels[idx], 0.99, 45.0, 0.7, 0.7, -1.1, 1.48, 45.0, -48.01],
            )
            # if no numbers, assume that EMG channel is on the right side
            (!occursin("1", emg_labels[idx]) && !occursin("2", emg_labels[idx])) && add_loc!(
                [emg_labels[idx], 0.99, 45.0, 0.7, 0.7, -1.1, 1.48, 45.0, -48.01],
            )
        end
    end

    # add locations of EOG channels
    eog_labels = get_channel(obj; type = "eog")
    if length(eog_labels) > 0
        for idx in eachindex(eog_labels)
            if occursin("1", eog_labels[idx])
                occursin("v", lowercase(eog_labels[idx])) && add_loc!(
                    [
                        eog_labels[idx],
                        1.01,
                        149.0,
                        -0.87,
                        0.51,
                        -0.37,
                        1.07,
                        149.0,
                        -20.15,
                    ],
                )
                occursin("h", lowercase(eog_labels[idx])) && add_loc!(
                    [
                        eog_labels[idx],
                        1.01,
                        129.73,
                        -0.64,
                        0.77,
                        -0.04,
                        1.0,
                        129.73,
                        -2.29,
                    ],
                )
            end
            if occursin("2", eog_labels[idx])
                occursin("v", lowercase(eog_labels[idx])) && add_loc!(
                    [eog_labels[idx], 1.01, 31.0, 0.87, 0.51, -0.37, 1.07, 31.0, -20.15],
                )
                occursin("h", lowercase(eog_labels[idx])) && add_loc!(
                    [eog_labels[idx], 1.01, 50.0, 0.64, 0.77, -0.04, 1.0, 50.27, -2.29],
                )
            end
            # if no numbers, assume that EOG channels are on the right side
            if !occursin("1", eog_labels[idx]) && !occursin("2", eog_labels[idx])
                occursin("v", lowercase(eog_labels[idx])) && add_loc!(
                    [eog_labels[idx], 1.01, 31.0, 0.87, 0.51, -0.37, 1.07, 31.0, -20.15],
                )
                occursin("h", lowercase(eog_labels[idx])) && add_loc!(
                    [eog_labels[idx], 1.01, 50.0, 0.64, 0.77, -0.04, 1.0, 50.27, -2.29],
                )
            end
            # if no V/H indicators, assume that EEG channels are vertical
            if !occursin("v", lowercase(eog_labels[idx])) &&
               !occursin("h", lowercase(eog_labels[idx]))
                occursin("1", eog_labels[idx]) && add_loc!(
                    [
                        eog_labels[idx],
                        1.01,
                        149.0,
                        -0.87,
                        0.51,
                        -0.37,
                        1.07,
                        149.0,
                        -20.15,
                    ],
                )
                occursin("2", eog_labels[idx]) && add_loc!(
                    [eog_labels[idx], 1.01, 31.0, 0.87, 0.51, -0.37, 1.07, 31.0, -20.15],
                )
            end
        end
    end

    # match location labels to channel labels case-insensitively (e.g. FP1 ↔ Fp1)
    # and use the channel label spelling
    ch_by_lc = Dict(lowercase(l) => l for l in labels(obj))
    locs[!, :label] = [get(ch_by_lc, lowercase(l), l) for l in locs[!, :label]]
    # keep the first location of each channel
    locs = unique(locs, :label)

    no_match = setdiff(labels(obj), locs[!, :label])
    length(no_match) > 0 && _warn(
        "Location$(_pl(no_match)): $(uppercase.(no_match)) could not be found in the LOCS file.",
    )

    # create new dataset
    obj_tmp = deepcopy(obj)

    obj_tmp.locs = Base.filter(:label => in(labels(obj)), locs)

    _locs_round!(obj_tmp.locs)
    _locs_remove_nans!(obj_tmp.locs)

    # keep order consistent with labels
    locs_idx = indexin(obj_tmp.locs[:, :label], labels(obj_tmp))
    obj_tmp.locs = obj_tmp.locs[sortperm(locs_idx), :]

    push!(obj_tmp.history, "load_locs(obj, file_name=$file_name)")

    return obj_tmp
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
  - CSV

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

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `file_name::String`: name of the file to load

# Return

- `Nothing`
"""
function load_locs!(obj::NeuroAnalyzer.NEURO; file_name::String)::Nothing
    obj_tmp = load_locs(obj; file_name = file_name)
    obj.locs = obj_tmp.locs
    obj_tmp = nothing

    return nothing
end
