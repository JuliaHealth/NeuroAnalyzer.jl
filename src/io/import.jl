export import_recording

"""
    import_recording(file_name; detect_type)

Load recording file and return `NeuroAnalyzer.NEURO` object. Supported formats:
- EDF/EDF+
- BDF/BDF+
- BrainVision
- CSV
- SET (EEGLAB dataset)
- FIFF

This is a meta-function that triggers appropriate `import_*()` function. File format is detected based on file extension (.edf|.bdf|.vhdr|.csv|.csv.gz|.set|.fif|.fiff). Signal data type (e.g. EEG or MEG is auto-detected) and stored in the `obj.header.recording[:data_type]` field.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `::NeuroAnalyzer.NEURO`
"""
function import_recording(file_name::String; detect_type::Bool=true)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    splitext(file_name)[2] == ".edf" && return import_edf(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".bdf" && return import_bdf(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".vhdr" && return import_bv(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".csv" && return import_csv(file_name, detect_type=detect_type)
    (splitext(file_name)[2] == ".gz" && splitext(splitext(file_name)[1])[2] == ".csv") && return import_csv(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".set" && return import_set(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".fif" && return import_fiff(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".fiff" && return import_fiff(file_name, detect_type=detect_type)
end
