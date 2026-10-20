export import_alice4

"""
    import_alice4(file_name; <keyword arguments>)

Load an EDF file exported from the Alice 4 Polysomnography System and return a `NeuroAnalyzer.NEURO` object.

Alice 4 EDF files are non-conforming in two ways, both handled by `import_edf`, to which this function delegates after checking the header:
- The `data_records` header field is always `-1` (unknown); the true record count is inferred from file size.
- Channels may have different sampling rates; they are upsampled to the highest rate using Fourier resampling.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: infer channel type from channel label
- `clean_labels::Bool=true`: normalize channel labels (see `import_recording()`)

# Returns

- `NeuroAnalyzer.NEURO`
"""
function import_alice4(file_name::String; detect_type::Bool = true, clean_labels::Bool = true)::NeuroAnalyzer.NEURO
    # validate
    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))

    header = open(file_name, "r") do fid
        buf = UInt8[]
        readbytes!(fid, buf, 256)
        return buf
    end
    length(header) == 256 ||
        throw(ArgumentError("$file_name is not a valid EDF file (header too short)."))

    # bytes 1–8: version; must be 0 for EDF
    tryparse(Int, strip(String(Char.(header[1:8])))) == 0 ||
        throw(ArgumentError("$file_name is not a valid EDF file (version ≠ 0)."))
    # bytes 89–168: local recording identification
    occursin("Alice 4", String(Char.(header[89:168]))) ||
        throw(ArgumentError("$file_name is not a valid Alice 4 EDF file."))

    return import_edf(file_name; detect_type = detect_type, clean_labels = clean_labels)
end
