export export_markers

"""
    export_markers(obj; <keyword arguments>)

Export `NeuroAnalyzer.NEURO` object markers to CSV.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `file_name::String`
- `overwrite::Bool=false`

# Returns

- `Nothing`
"""
function export_markers(
        obj::NeuroAnalyzer.NEURO; file_name::String, overwrite::Bool = false
    )::Nothing

    !(!(isfile(file_name) && !overwrite)) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))

    # MARKERS
    if DataFrames.nrow(obj.markers) > 0
        !(!(isfile(file_name) && !overwrite)) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))
        CSV.write(file_name, obj.markers)
    end

    return nothing

end
