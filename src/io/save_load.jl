export save_neuro
export load_neuro

"""
    save_neuro(obj; <keyword arguments>)

Save `NeuroAnalyzer.NEURO` object to `file_name` file (HDF5-based).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `file_name::String`: name of the file to save to
- `overwrite::Bool=false`

# Return

Nothing
"""
function save_neuro(
    obj::NeuroAnalyzer.NEURO; file_name::String, overwrite::Bool = false,
)::Nothing
    # validate
    (isfile(file_name) && !overwrite) && throw(
        ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."),
    )
    lowercase(splitext(file_name)[2]) == ".hdf" ||
        throw(ArgumentError("Filename extension must be .hdf"))

    obj.header.recording[:file_name] = file_name

    # the file size is stored in the object, so it is measured on a temporary copy first
    tmp_file = tempname() * ".hdf"
    try
        JLD2.save_object(tmp_file, obj)
        obj.header.recording[:file_size_mb] = round(filesize(tmp_file) / 1024^2; digits = 2)
    finally
        rm(tmp_file; force = true)
    end

    JLD2.save_object(file_name, obj)

    return nothing
end

"""
    load_neuro(file_name)

Load `NeuroAnalyzer.NEURO` object from `file_name` file (HDF5-based).

# Arguments

- `file_name::String`: name of the file to load

# Returns

- `NeuroAnalyzer.NEURO`
"""
function load_neuro(file_name::String)::NeuroAnalyzer.NEURO
    # validate
    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    obj = JLD2.load_object(file_name)

    _info(
        "Loaded: " *
        uppercase(obj.header.recording[:data_type]) *
        " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits = 2)) s)",
    )

    return obj
end
