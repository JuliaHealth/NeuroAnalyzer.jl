export zipratio

"""
    zipratio(obj)

Calculate the zip ratio for a NEURO object data as a measure of signal complexity.

The zip ratio is defined as the size of the maximally compressed file (`zip -9`) divided by the size of the raw (uncompressed) CSV export. Lower values indicate lower complexity (higher compressibility).

`zipratio()` requires the `zip` command (Linux/macOS) or `zip.exe` (Windows) to be available on `PATH`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Float64`: ratio of compressed to uncompressed data size, in the range (0, 1]
"""
function zipratio(obj::NeuroAnalyzer.NEURO)::Float64

    # determine the platform-appropriate zip executable name
    zip_cmd = Sys.iswindows() ? "zip.exe" : "zip"
    Sys.which(zip_cmd) === nothing &&
        throw(
            ArgumentError(
                "zip command not found: \"$zip_cmd\". " *
                "Install zip and ensure it is on PATH.",
            ),
        )

    # create a temporary CSV file for the exported signal data
    tmp_path, tmp_io = mktemp()
    close(tmp_io)
    csv_path = tmp_path * ".csv"
    zip_path = tmp_path * ".zip"

    # ensure all temp files are removed regardless of success or failure
    try
        # export the signal data to CSV (no headers/names to minimize
        # non-signal content that could skew the compression ratio)
        export_csv(
            obj;
            file_name = csv_path,
            names = false,
            header = false,
            epoch_time = false,
            markers = false,
            locs = false,
            history = false,
            overwrite = true,
        )

        _info("Compressing exported data to estimate signal complexity")

        # compress at maximum level (-9) and measure the result
        # -q suppresses zip's stdout chatter
        run(`$zip_cmd -9 -q $zip_path $csv_path`)
        zip_size_9 = filesize(zip_path)
        rm(zip_path)

        raw_size = filesize(csv_path)

        # guard against a degenerate empty export
        raw_size == 0 &&
            throw(ArgumentError("Exported CSV is empty; cannot compute zip ratio."))

        zip_ratio = zip_size_9 / raw_size

        return zip_ratio

    finally

        # always remove temp files, even if an exception was thrown above
        isfile(csv_path) && rm(csv_path)
        isfile(zip_path) && rm(zip_path)
        # remove the bare mktemp path that was never used directly
        isfile(tmp_path) && rm(tmp_path)
    end
end
