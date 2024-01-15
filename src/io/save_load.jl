export save
export load
export save_study
export load_study

"""
    save(obj; file_name, overwrite)

Save `NeuroAnalyzer.NEURO` object to `file_name` file (HDF5-based).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `file_name::String`: name of the file to save to
- `overwrite::Bool=false`
"""
function save(obj::NeuroAnalyzer.NEURO; file_name::String, overwrite::Bool=false)

    @assert !(isfile(file_name) && overwrite == false) "File $file_name cannot be saved, to overwrite use overwrite=true."
    @assert splitext(file_name)[2] == ".hdf" "File name extension should be .HDF."

    obj.header.recording[:file_name] = file_name

    JLD2.save_object("/tmp/$(basename(file_name))", obj)
    obj.header.recording[:file_size_mb] = round(filesize("/tmp/$(basename(file_name))") / 1024, digits=2)
    rm("/tmp/$(basename(file_name))")

    JLD2.save_object(file_name, obj)

end

"""
    load(file_name)

Load `NeuroAnalyzer.NEURO` object from `file_name` file (HDF5-based).

# Arguments

- `file_name::String`: name of the file to load

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function load(file_name::String)

    @assert isfile(file_name) "File $file_name cannot be loaded."

    obj = JLD2.load_object(file_name)

    _info("Loaded: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(obj.time_pts[end]) s)")

    return obj

end

"""
    save_study(obj; file_name, overwrite)

Save `NeuroAnalyzer.STUDY` object to `file_name` file (HDF5-based).

# Arguments

- `obj::NeuroAnalyzer.STUDY`
- `file_name::String`: name of the file to save to
- `overwrite::Bool=false`
"""
function save_study(obj::NeuroAnalyzer.STUDY; file_name::String, overwrite::Bool=false)

    @assert !(isfile(file_name) && overwrite == false) "File $file_name cannot be saved, to overwrite use overwrite=true."

    @assert splitext(file_name)[2] == ".hdf" "File name extension should be .HDF."

    JLD2.save_object(file_name, obj)

end

"""
    load_study(file_name)

Load `NeuroAnalyzer.STUDY` object from `file_name` file (HDF5-based).

# Arguments

- `file_name::String`: name of the file to load

# Returns

- `obj::NeuroAnalyzer.STUDY`
"""
function load_study(file_name::String)

    @assert isfile(file_name) "File $file_name cannot be loaded."

    obj = JLD2.load_object(file_name)

    _info("Loaded study: $(obj_n(obj)) objects")

    return obj

end
