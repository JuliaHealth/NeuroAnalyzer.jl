export import_recording

"""
    import_recording(file_name; <keyword arguments>)

Load recording file and return `NeuroAnalyzer.NEURO` object. Supported formats:
- EDF/EDF+
- BDF/BDF+
- GDF
- BrainVision
- CSV or CSV.GZ
- VHDR or AHDR (BrainVision)
- SET (EEGLAB dataset)
- NPY (MNE)
- XDF (Extensible Data Format)
- NWB (Neurodata Without Borders)
- NCS (Neuralinx Continuously Sampled Channels (CSC))
- FIFF (Neuromag)
- MAT (FieldTrip)
- SNIRF
- NIRS
- ASCII or M (DuoMAG TMS MEP)

This is a meta-function that triggers appropriate `import_*()` function. File format is detected based on file extension (.edf|.bdf|.gdf|.vhdr|.ahdr|.csv|.csv.gz|.set|.npy|.xdf|.nwb|.ncs|.fif|.fiff|.mat|.snirf|.nirs|.ascii|.m). Signal data type (e.g. EEG or MEG is auto-detected, if possible) and stored in the `obj.header.recording[:data_type]` field.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on channel label
- `type::Union{Nothing, Symbol}=nothing`: type of imported data (required for FieldTrip files)
    - `:eeg`: EEG
    - `:meg`: MEG
    - `:nirs`: fNIRS
    - `:events`: events
- `n::Int64=0`: subject number to extract in case of multi-subject file (required for SNIRF files)
- `sampling_rate::Union{Nothing, Int64}=nothing`: NPY file contains only signal data, therefore its sampling rate must be provided upon importing

# Returns

- `obj::NeuroAnalyzer.NEURO` - for EEG, MEG, fNIRS data
- `markers::DataFrame` - for events
"""
function import_recording(file_name::String; detect_type::Bool=true, type::Union{Nothing, Symbol}=nothing, sampling_rate::Union{Nothing, Int64}=nothing, n::Int64=0)::NeuroAnalyzer.NEURO

    @assert isfile(file_name) "File $file_name cannot be loaded."

    splitext(file_name)[2] == ".edf" && return import_edf(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".bdf" && return import_bdf(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".gdf" && return import_gdf(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".vhdr" && return import_bv(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".ahdr" && return import_bv(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".csv" && return import_csv(file_name, detect_type=detect_type)
    (splitext(file_name)[2] == ".gz" && splitext(splitext(file_name)[1])[2] == ".csv") && return import_csv(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".set" && return import_set(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".npy" && return import_npy(file_name, sampling_rate=sampling_rate)
    splitext(file_name)[2] == ".xdf" && return import_xdf(file_name)
    splitext(file_name)[2] == ".nwb" && return import_nwb(file_name, detect_type=detect_type)
    splitext(file_name)[2] == ".ncs" && return import_ncs(file_name)
    splitext(file_name)[2] == ".fif" && return import_fiff(file_name)
    splitext(file_name)[2] == ".fiff" && return import_fiff(file_name)
    splitext(file_name)[2] == ".mat" && return import_ft(file_name, type=type)
    splitext(file_name)[2] == ".snirf" && return import_snirf(file_name, n=n)
    splitext(file_name)[2] == ".nirs" && return import_nirs(file_name)
    splitext(file_name)[2] == ".ascii" && return import_duomag(file_name)
    splitext(file_name)[2] == ".m" && return import_duomag(file_name)

end
