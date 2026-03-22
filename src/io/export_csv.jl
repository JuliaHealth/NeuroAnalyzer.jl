export export_csv

"""
    export_csv(obj; <keyword arguments>)

Export a `NeuroAnalyzer.NEURO` object to CSV (and optional sidecar files).

The primary signal data is always written to `file_name`. When the corresponding flag is `true`, companion files are created next to it using the same base name:

|       Option      |       Output file       |
| ----------------- | ----------------------- |
| `header=true`     | `<base>_header.txt`     |
| `epoch_time=true` | `<base>_epoch_time.csv` |
| `markers=true`    | `<base>_markers.csv`    |
| `locs=true`       | `<base>_locs.csv`       |
| `history=true`    | `<base>_history.txt`    |

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `file_name::String`: path of the primary output CSV file
- `names::Bool=true`: write channel names as the CSV header row
- `header::Bool=false`: export recording/subject/experiment metadata to a sidecar text file
- `epoch_time::Bool=false`: export per-epoch time points to a sidecar CSV
- `markers::Bool=false`: export event markers to a sidecar CSV (skipped when the markers table is empty)
- `locs::Bool=false`: export channel locations to a sidecar CSV (skipped when the locations table is empty)
- `history::Bool=false`: export processing history to a sidecar text file (skipped when history is empty)
- `overwrite::Bool=false`: allow overwriting existing files; throws `ArgumentError` otherwise

# Returns

- `Nothing`

# Throws

- `ArgumentError` if any target file already exists and `overwrite=false`
"""
function export_csv(
    obj::NeuroAnalyzer.NEURO;
    file_name::String,
    names::Bool = true,
    header::Bool = false,
    epoch_time::Bool = false,
    markers::Bool = false,
    locs::Bool = false,
    history::Bool = false,
    overwrite::Bool = false
)::Nothing

    # internal guard: throw if the file already exists and overwriting is not permitted
    # defined once here to avoid repeating the condition
    check_overwrite(path) =
        isfile(path) && !overwrite &&
            throw(ArgumentError(
                "File $path already exists; use overwrite=true to overwrite."))

    # derive ALL companion paths upfront from the original base name
    base = replace(file_name, r"\.csv$"i => "")
    header_file = base * "_header.txt"
    epoch_time_file = base * "_epoch_time.csv"
    markers_file = base * "_markers.csv"
    locs_file = base * "_locs.csv"
    history_file = base * "_history.txt"

    # ------------------------------------------------------------------ #
    # signal data                                                        #
    # ------------------------------------------------------------------ #
    check_overwrite(file_name)

    # flatten epochs: (channels × samples × epochs) → (channels × total_samples)
    # `reshape` returns a view here (no copy), so this is allocation-free.
    n_ch, n_s, n_ep = size(obj.data)
    flat = reshape(obj.data, n_ch, n_s * n_ep)

    # prepend the global time axis; transpose so rows correspond to time points.
    col_names = vcat("time", labels(obj))
    df = DataFrame(hcat(obj.time_pts, flat'), col_names)
    CSV.write(file_name, df; writeheader = names)

    # ------------------------------------------------------------------ #
    # metadata header                                                    #
    # ------------------------------------------------------------------ #
    if header
        check_overwrite(header_file)
        open(header_file, "w") do f
            # write key-value pairs from all three metadata sections
            for section in (obj.header.subject, obj.header.recording, obj.header.experiment)
                for (key, value) in section
                    println(f, key, ": ", value)
                end
            end
        end
    end

    # ------------------------------------------------------------------ #
    # epoch time points                                                  #
    # ------------------------------------------------------------------ #
    if epoch_time
        check_overwrite(epoch_time_file)
        CSV.write(epoch_time_file, DataFrame(; epoch_time = obj.epoch_time))
    end

    # ------------------------------------------------------------------ #
    # event markers                                                      #
    # ------------------------------------------------------------------ #
    if markers && DataFrames.nrow(obj.markers) > 0
        check_overwrite(markers_file)
        CSV.write(markers_file, obj.markers)
    end

    # ------------------------------------------------------------------ #
    # channel locations                                                  #
    # ------------------------------------------------------------------ #
    if locs && DataFrames.nrow(obj.locs) > 0
        check_overwrite(locs_file)
        CSV.write(locs_file, obj.locs)
    end

    # ------------------------------------------------------------------ #
    # processing history                                                 #
    # ------------------------------------------------------------------ #
    if history && !isempty(obj.history)
        check_overwrite(history_file)
        open(history_file, "w") do f
            for entry in obj.history
                println(f, entry)
            end
        end
    end

    return nothing

end
