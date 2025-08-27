export save
export load

"""
    save(obj; <keyword arguments>)

Save `NeuroAnalyzer.NEURO` object to `file_name` file (HDF5-based).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `file_name::String`: name of the file to save to
- `overwrite::Bool=false`

# Return

Nothing
"""
function save(obj::NeuroAnalyzer.NEURO; file_name::String, overwrite::Bool=false)::Nothing

    @assert !(isfile(file_name) && !overwrite) "File $file_name cannot be saved, to overwrite use overwrite=true."
    @assert lowercase(splitext(file_name)[2]) == ".hdf" "Filename extension must be .hdf"

    obj.header.recording[:file_name] = file_name

    JLD2.save_object("/tmp/$(basename(file_name))", obj)
    obj.header.recording[:file_size_mb] = round(filesize("/tmp/$(basename(file_name))") / 1024, digits=2)
    rm("/tmp/$(basename(file_name))")

    JLD2.save_object(file_name, obj)

    return nothing

end

"""
    load(file_name)

Load `NeuroAnalyzer.NEURO` object from `file_name` file (HDF5-based).

# Arguments

- `file_name::String`: name of the file to load

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function load(file_name::String)::NeuroAnalyzer.NEURO

    @assert isfile(file_name) "File $file_name cannot be loaded."

    obj = JLD2.load_object(file_name)

    _info("Loaded: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(obj.time_pts[end]) s)")

    return obj

end
