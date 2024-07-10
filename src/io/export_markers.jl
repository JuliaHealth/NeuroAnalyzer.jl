export export_markers

"""
    export_markers(obj; <keyword arguments>)

Export `NeuroAnalyzer.NEURO` object markers to CSV.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `file_name::String`
- `overwrite::Bool=false`
"""
function export_markers(obj::NeuroAnalyzer.NEURO; file_name::String, overwrite::Bool=false)

    @assert !(isfile(file_name) && !overwrite) "File $file_name cannot be saved, to overwrite use overwrite=true."

    # MARKERS
    if nrow(obj.markers) > 0
        @assert !(isfile(file_name) && !overwrite) "File $file_name cannot be saved, to overwrite use overwrite=true."
        CSV.write(file_name, obj.markers)
    end

end
