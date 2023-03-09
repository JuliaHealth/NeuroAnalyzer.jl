export export_csv

"""
    export_csv(eeg; file_name, header, components, markers, overwrite)

Export EEG data as CSV.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `file_name::String`
- `header::Bool=false`: export header
- `components::Bool=false`: export components
- `markers::Bool=false`: export markers
- `locs::Bool=false`: export locations
- `overwrite::Bool=false`

# Returns

- `success::Bool`
"""
function export_csv(obj::NeuroAnalyzer.NEURO; file_name::String, header::Bool=false, components::Bool=false, markers::Bool=false, locs::Bool=false, overwrite::Bool=false)

    (isfile(file_name) && overwrite == false) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))
    obj.header[:components] == [""] && throw(ArgumentError("EEG does not contain components."))

    # DATA
    # unsplit epochs
    s_merged = reshape(obj.data,
                       size(obj.data, 1),
                       size(obj.data, 2) * size(obj.data, 3))
    s = s_merged[:, :, 1]'
    s = hcat(obj.time_pts, s)
    l = vcat("time", labels(eeg))
    CSV.write(file_name, DataFrame(s, l))

    # HEADER
    if header
        file_name = replace(file_name, ".csv" => "_header.csv")
        (isfile(file_name) && overwrite == false) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))
        f = open(file_name, "w")
        for (key, value) in obj.header
            println(f, key, ": ", value)
        end
        close(f)
    end

    # COMPONENTS
    if components
        file_name = replace(file_name, ".csv" => "_components.csv")
        (isfile(file_name) && overwrite == false) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))
        f = open(file_name, "w")
        for idx in eachindex(obj.header[:components])
            println(f, "component: $(obj.header[:components][idx])")
            println(f, obj.components[idx])
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

