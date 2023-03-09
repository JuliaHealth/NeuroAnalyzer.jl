export save
export load

"""
    save(eeg; file_name, overwrite)

Save `eeg` to `file_name` file (HDF5-based).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `file_name::String`: file name
- `overwrite::Bool=false`

# Returns

- `success::Bool`
"""
function save(obj::NeuroAnalyzer.NEURO; file_name::String, overwrite::Bool=false)

    (isfile(file_name) && overwrite == false) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))

    obj.header.recording[:file_name] = file_name

    save_object("/tmp/$(basename(file_name))", eeg)
    obj.header.recording[:file_size_mb] = round(filesize("/tmp/$(basename(file_name))") / 1024, digits=2)
    rm("/tmp/$(basename(file_name))")

    save_object(file_name, eeg)
end

"""
    load(file_name)

Load `eeg` from `file_name` file (HDF5-based).

# Arguments

- `file_name::String`: file name

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function load(file_name::String)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    eeg = load_object(file_name)

    return eeg
end

