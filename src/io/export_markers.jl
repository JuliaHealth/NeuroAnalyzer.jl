export export_markers

"""
    export_markers(obj; <keyword arguments>)

Export the event markers table of a `NeuroAnalyzer.NEURO` object to a CSV file.

If the markers table is empty, a warning is issued and no file is written.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `file_name::String`: path of the output CSV file
- `overwrite::Bool=false`: allow overwriting an existing file; throws `ArgumentError` otherwise

# Returns

- `Nothing`

# Throws

- `ArgumentError` if `file_name` already exists and `overwrite=false`
"""
function export_markers(
    obj::NeuroAnalyzer.NEURO;
    file_name::String,
    overwrite::Bool = false
)::Nothing

    isfile(file_name) && !overwrite &&
        throw(ArgumentError(
            "File $file_name already exists; use overwrite=true to overwrite."))

    if DataFrames.nrow(obj.markers) == 0
        @warn "No markers to export; $file_name was not written."
        return nothing
    end

    CSV.write(file_name, obj.markers)

    return nothing

end
