export export_csv

"""
    export_csv(obj; <keyword arguments>)

Export `NeuroAnalyzer.NEURO` object to CSV.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `file_name::String`
- `header::Bool=false`: export header
- `epoch_time::Bool=false`: export epoch time points
- `components::Bool=false`: export components
- `markers::Bool=false`: export event markers
- `locs::Bool=false`: export channel locations
- `history::Bool=false`: export history
- `overwrite::Bool=false`

# Returns

Nothing
"""
function export_csv(obj::NeuroAnalyzer.NEURO; file_name::String, header::Bool=false, epoch_time::Bool=false, components::Bool=false, markers::Bool=false, locs::Bool=false, history::Bool=false, overwrite::Bool=false)::Nothing

    @assert !(isfile(file_name) && !overwrite) "File $file_name cannot be saved, to overwrite use overwrite=true."

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
        file_name = replace(file_name, ".csv" => "_header.txt")
        @assert !(isfile(file_name) && !overwrite) "File $file_name cannot be saved, to overwrite use overwrite=true."
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

    # EPOCH TIME POINST
    if epoch_time
        file_name = replace(file_name, ".csv" => "_epoch_time.csv")
        @assert !(isfile(file_name) && !overwrite) "File $file_name cannot be saved, to overwrite use overwrite=true."
        CSV.write(file_name, obj.epoch_time)
    end

    # COMPONENTS
    if components && length(keys(obj.components)) > 0
        file_name = replace(file_name, ".csv" => "_components.txt")
        @assert !(isfile(file_name) && !overwrite) "File $file_name cannot be saved, to overwrite use overwrite=true."
        f = open(file_name, "w")
        for c_idx in eachindex(keys(obj.components))
            println(f, "component: $(keys(obj.components)[c_idx])")
            println(f, obj.components[c_idx])
            println(f, "---")
        end
        close(f)
    end

    # MARKERS
    if markers && nrow(obj.markers) > 0
        file_name = replace(file_name, ".csv" => "_markers.csv")
        @assert !(isfile(file_name) && !overwrite) "File $file_name cannot be saved, to overwrite use overwrite=true."
        CSV.write(file_name, obj.markers)
    end

    # LOCS
    if locs && nrow(obj.locs) > 0
        file_name = replace(file_name, ".csv" => "_locs.csv")
        @assert (isfile(file_name) && !overwrite) "File $file_name cannot be saved, to overwrite use overwrite=true."
        CSV.write(file_name, obj.locs)
    end

    # HISTORY
    if history && length(obj.history) > 0
        file_name = replace(file_name, ".csv" => "_history.txt")
        for c_idx in eachindex(history)
            println(f, obj.history[c_idx])
        end
    end

    return nothing

end
