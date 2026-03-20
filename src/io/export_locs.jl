export export_locs

"""
    export_locs(obj; <keyword arguments>)

Export channel location data from a `NeuroAnalyzer.NEURO` object to a file.

The output format is determined automatically from `file_name`'s extension:

| Extension |                   Format                   |
| --------- | ------------------------------------------ |
| `.csv`    | Full locations table (all columns)         |
| `.ced`    | EEGLAB CED — tab-separated with header     |
| `.locs`   | EEGLAB .locs — tab-separated, no header    |
| `.tsv`    | BIDS-style TSV — tab-separated with header |

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `file_name::String`: output file path; extension determines the format
- `overwrite::Bool=false`: allow overwriting an existing file; throws `ArgumentError` otherwise


# Returns

- `Nothing`

# Throws

- `ArgumentError` if `file_name` already exists and `overwrite=false`
- `ArgumentError` if the file extension is not one of `.csv`, `.ced`, `.locs`, `.tsv`
"""
function export_locs(
        obj::NeuroAnalyzer.NEURO; file_name::String, overwrite::Bool = false
    )::Nothing

    # the .csv branch is only available via the NEURO method (full locs table)
    # all other formats are handled by delegating to the DataFrame method
    if splitext(file_name)[2] == ".csv"
        isfile(file_name) && !overwrite &&
            throw(ArgumentError("File $file_name already exists; use overwrite=true to overwrite."))
        CSV.write(file_name, obj.locs)
    else
        # delegate to the DataFrame method, which handles .ced / .locs / .tsv
        export_locs(obj.locs; file_name, overwrite)
    end

    return nothing
end

"""
    export_locs(locs; <keyword arguments>)

Export a channel-locations `DataFrame` to a file.

The output format is determined automatically from `file_name`'s extension:

| Extension |                   Format                   |
| --------- | ------------------------------------------ |
| `.ced`    | EEGLAB CED — tab-separated with header     |
| `.locs`   | EEGLAB .locs — tab-separated, no header    |
| `.tsv`    | BIDS-style TSV — tab-separated with header |

# Arguments

- `locs::DataFrame`: channel locations table (must contain the standard NeuroAnalyzer location columns)
- `file_name::String`: output file path; extension determines the format
- `overwrite::Bool=false`: allow overwriting an existing file; throws `ArgumentError` otherwise

# Returns

- `Nothing`

# Throws

- `ArgumentError` if `file_name` already exists and `overwrite=false`
- `ArgumentError` if the file extension is not one of `.ced`, `.locs`, `.tsv`
"""
function export_locs(
    locs::DataFrame;
    file_name::String,
    overwrite::Bool = false
)::Nothing

    # guard against accidental overwrites before doing any work
    isfile(file_name) && !overwrite &&
        throw(ArgumentError("File $file_name already exists; use overwrite=true to overwrite."))

    # extract the extension once; used for every branch below
    ext = splitext(file_name)[2]

    # extract all location columns up front
    channels = _find_bylabel(locs, locs[!, :label])
    clabels = locs[!, :label]
    theta = locs[!, :loc_theta]
    radius = locs[!, :loc_radius]
    x = locs[!, :loc_x]
    y = locs[!, :loc_y]
    z = locs[!, :loc_z]
    radius_sph = locs[!, :loc_radius_sph]
    theta_sph = locs[!, :loc_theta_sph]
    phi_sph = locs[!, :loc_phi_sph]

    if ext == ".ced"
        # EEGLAB CED format: tab-separated, column header required.
        df = DataFrame(
            Number = channels,
            labels = clabels,
            theta = theta,
            radius = radius,
            X = x,
            Y = y,
            Z = z,
            sph_theta = theta_sph,
            sph_phi = phi_sph,
            sph_radius = radius_sph,
        )
        CSV.write(file_name, df, delim = "\t", header = true)

    elseif ext == ".locs"
        # EEGLAB .locs format: tab-separated, no column header.
        df = DataFrame(
            Number = channels,
            theta = theta,
            radius = radius,
            labels = clabels,
        )
        CSV.write(file_name, df, delim = "\t", header = false)

    elseif ext == ".tsv"
        # BIDS-style TSV: tab-separated, column header required.
        df = DataFrame(
            labels = clabels,
            x = x,
            y = y,
            z = z,
            theta = theta,
            radius = radius,
            radius_sph = radius_sph,
            theta_sph = theta_sph,
            phi_sph = phi_sph,
        )
        CSV.write(file_name, df, delim = "\t", header = true)

    else
        throw(ArgumentError(
            "Unsupported extension \"$ext\". file_name must end in .ced, .locs, or .tsv."))
    end

    return nothing
end
