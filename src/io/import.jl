export import_recording

"""
    import_recording(file_name; <keyword arguments>)

Load recording file and return `NeuroAnalyzer.NEURO` object. Supported formats:
- EDF/EDF+
- BDF/BDF+
- GDF
- BrainVision
- CSV
- SET (EEGLAB dataset)
- NWB (Neurodata Without Borders)
- FIFF
- SNIRF
- NIRS

This is a meta-function that triggers appropriate `import_*()` function. File format is detected based on file extension (.edf|.bdf|.gdf|.vhdr|.ahdr|.csv|.csv.gz|.set|.nwb|.fif|.fiff|.snirf|.nirs). Signal data type (e.g. EEG or MEG is auto-detected) and stored in the `obj.header.recording[:data_type]` field.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on channel label
- `n::Int64=0`: subject number to extract in case of multi-subject file

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function import_recording(file_name::String; detect_type::Bool=true, n::Int64=0)::NeuroAnalyzer.NEURO

    @assert isfile(file_name) "File $file_name cannot be loaded."

    splitext(file_name)[2] == ".edf" && return import_edf(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".bdf" && return import_bdf(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".gdf" && return import_gdf(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".vhdr" && return import_bv(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".ahdr" && return import_bv(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".csv" && return import_csv(file_name, detect_type=detect_type)
    (splitext(file_name)[2] == ".gz" && splitext(splitext(file_name)[1])[2] == ".csv") && return import_csv(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".set" && return import_set(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".nwb" && return import_nwb(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".fif" && return import_fiff(file_name)
    splitext(file_name)[2] == ".fiff" && return import_fiff(file_name)
    splitext(file_name)[2] == ".snirf" && return import_snirf(file_name, n=n)
    splitext(file_name)[2] == ".nirs" && return import_nirs(file_name)

end
