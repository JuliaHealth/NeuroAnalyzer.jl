export import_recording

"""
    import_recording(file_name; <keyword arguments>)

Load a recording file and return a `NeuroAnalyzer.NEURO` object.

This is a meta-function: it detects the file format from the extension and delegates to the appropriate `import_*()` function. The signal data type (EEG, MEG, fNIRS, …) is auto-detected where possible and stored in `obj.header.recording[:data_type]`.

Supported formats and their extensions:

|             Format             |    Extension(s)   |
| ------------------------------ | ----------------- |
| EDF / EDF+                     | `.edf`            |
| BDF / BDF+                     | `.bdf`            |
| GDF                            | `.gdf`            |
| BrainVision                    | `.vhdr`, `.ahdr`  |
| CSV (plain or gzip-compressed) | `.csv`, `.csv.gz` |
| EEGLAB dataset                 | `.set`            |
| MNE NumPy array                | `.npy`            |
| Extensible Data Format         | `.xdf`            |
| Neurodata Without Borders      | `.nwb`            |
| Neuralynx CSC                  | `.ncs`            |
| Neuromag FIFF                  | `.fif`, `.fiff`   |
| FieldTrip MATLAB               | `.mat`            |
| SNIRF                          | `.snirf`          |
| NIRS                           | `.nirs`           |
| DuoMAG TMS MEP                 | `.ascii`, `.m`    |

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: infer channel type from channel label
- `type::Union{Nothing, Symbol}=nothing`: data type for FieldTrip (`.mat`) files, where auto-detection is not possible; one of `:eeg`, `:meg`, `:nirs`, `:events`
- `sampling_rate::Union{Nothing, Int64}=nothing`: sampling rate for `.npy` files, which store only raw signal data with no embedded metadata
- `n::Int64=0`: subject index to extract from multi-subject SNIRF files

# Returns

- `NeuroAnalyzer.NEURO`: for EEG, MEG, and fNIRS data
- `DataFrame`: when `type = :events` (FieldTrip `.mat` files)

# Throws

- `ArgumentError` if `file_name` does not exist
- `ArgumentError` if the file extension is not recognized
"""
function import_recording(
    file_name::String;
    detect_type::Bool = true,
    type::Union{Nothing, Symbol} = nothing,
    sampling_rate::Union{Nothing, Int64} = nothing,
    n::Int64 = 0
)::NeuroAnalyzer.NEURO

    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))

    ext = lowercase(splitext(file_name)[2])

    # special case: double extension .csv.gz
    if ext == ".gz"
        inner_ext = lowercase(splitext(splitext(file_name)[1])[2])
        inner_ext == ".csv" && return import_csv(file_name, detect_type)
        throw(ArgumentError(
            "Unsupported compressed format \"$inner_ext.gz\" in $file_name."))
    end

    # dispatch table - each branch returns immediately on a match
    ext == ".edf" && return import_edf(file_name, detect_type)
    ext == ".bdf" && return import_bdf(file_name, detect_type)
    ext == ".gdf" && return import_gdf(file_name, detect_type)
    # .vhdr and .ahdr are both BrainVision header formats
    ext in (".vhdr", ".ahdr") && return import_bv(file_name, detect_type)
    ext == ".csv" && return import_csv(file_name, detect_type)
    ext == ".set" && return import_set(file_name, detect_type)
    # .npy stores raw signal data only - sampling rate must be supplied
    ext == ".npy" && return import_npy(file_name; sampling_rate)
    ext == ".xdf" && return import_xdf(file_name)
    ext == ".nwb" && return import_nwb(file_name, detect_type)
    ext == ".ncs" && return import_ncs(file_name)
    # .fif and .fiff are both valid FIFF extensions
    ext in (".fif", ".fiff") && return import_fiff(file_name)
    # .mat requires explicit type hint - auto-detection is not possible
    ext == ".mat" && return import_ft(file_name; type)
    ext == ".snirf" && return import_snirf(file_name; n)
    ext == ".nirs" && return import_nirs(file_name)
    # .ascii and .m are both DuoMAG TMS MEP formats
    ext in (".ascii", ".m") && return import_duomag(file_name)

    throw(ArgumentError(
        "Unsupported file format \"$ext\" in $file_name. " *
        "Supported extensions: .edf, .bdf, .gdf, .vhdr, .ahdr, .csv, .csv.gz, " *
        ".set, .npy, .xdf, .nwb, .ncs, .fif, .fiff, .mat, .snirf, .nirs, .ascii, .m"))

end
