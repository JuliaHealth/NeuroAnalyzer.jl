export export_csv

"""
    export_csv(obj; file_name, header, components, markers, overwrite)

Export `NeuroAnalyzer.NEURO` object to CSV.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `file_name::String`
- `header::Bool=false`: export header
- `components::Bool=false`: export components
- `markers::Bool=false`: export event markers
- `locs::Bool=false`: export channel locations
- `overwrite::Bool=false`

# Returns

- `success::Bool`
"""
function export_csv(obj::NeuroAnalyzer.NEURO; file_name::String, header::Bool=false, components::Bool=false, markers::Bool=false, locs::Bool=false, overwrite::Bool=false)

    (isfile(file_name) && overwrite == false) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))

    # DATA
    # unsplit epochs
    out = reshape(obj.data,
                  size(obj.data, 1),
                  size(obj.data, 2) * size(obj.data, 3))
    s = out[:, :, 1]'
    s = hcat(obj.time_pts, s)
    l = vcat("time", labels(obj))
    CSV.write(file_name, DataFrame(s, l))

    # HEADER
    if header
        file_name = replace(file_name, ".csv" => "_header.csv")
        (isfile(file_name) && overwrite == false) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))
        f = open(file_name, "w")
        for (key, value) in obj.header.subject
            println(f, key, ": ", value)
        end
        for (key, value) in obj.header.recording
            println(f, key, ": ", value)
        end
        for (key, value) in obj.header.experiment
            println(f, key, ": ", value)
        end
        close(f)
    end

    # COMPONENTS
    if components
        length(obj.header.component_names) == 0 && throw(ArgumentError("OBJ does not contain components."))
        file_name = replace(file_name, ".csv" => "_components.csv")
        (isfile(file_name) && overwrite == false) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))
        f = open(file_name, "w")
        for c_idx in 1:length(obj.header.component_names)
            println(f, "component: $(obj.header.component_names[c_idx])")
            println(f, obj.components[c_idx])
            println(f, "---")
        end
        close(f)
    end

    # MARKERS
    if markers
        file_name = replace(file_name, ".csv" => "_markers.csv")
        (isfile(file_name) && overwrite == false) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))
        CSV.write(file_name, obj.markers)
    end

    # LOCS
    if locs
        file_name = replace(file_name, ".csv" => "_locs.csv")
        (isfile(file_name) && overwrite == false) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))
        CSV.write(file_name, obj.locs)
    end
end
