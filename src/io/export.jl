export export_recording

# all formats export_recording() can dispatch to
const _EXPORT_FORMATS = (:fiff, :bv, :snirf, :csv, :neuro)

# extension → export format
const _EXPORT_EXT_FORMAT = Dict{String, Symbol}(
    ".fif" => :fiff,
    ".fiff" => :fiff,
    ".vhdr" => :bv,
    ".snirf" => :snirf,
    ".csv" => :csv,
    ".hdf" => :neuro,
)

# sidecar files written by export_csv(): keyword => file suffix
const _CSV_SIDECARS = (
    :header => "_header.txt",
    :epoch_time => "_epoch_time.csv",
    :markers => "_markers.csv",
    :locs => "_locs.csv",
    :history => "_history.txt",
)

"""
    export_recording(obj; <keyword arguments>)

Export a `NeuroAnalyzer.NEURO` object to a file.

This is a meta-function: the format is selected by the file name extension (or by `format`) and the export is delegated to the appropriate `export_*()` function. Other keyword arguments are passed to that function.

|              Format              | Extension(s)    |   `format`  |   Function    |      Data types        |
| -------------------------------- | --------------- | ----------- | ------------- | ---------------------- |
| Neuromag FIFF                    | `.fif`, `.fiff` | `:fiff`     | `export_fiff` | EEG, MEG, ECoG, sEEG   |
| BrainVision                      | `.vhdr`         | `:bv`       | `export_bv`   | EEG, ECoG, sEEG, iEEG, ERP |
| SNIRF                            | `.snirf`        | `:snirf`    | `export_snirf`| NIRS                   |
| CSV                              | `.csv`          | `:csv`      | `export_csv`  | all                    |
| NeuroAnalyzer (JLD2/HDF5)        | `.hdf`          | `:neuro`    | `save_neuro`  | all                    |

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `file_name::String`: name of the file to export to
- `format::Union{Nothing, Symbol}=nothing`: force a specific format (see table) instead of detecting it from the extension; the extension must still be the one required by the exporter
- `overwrite::Bool=false`: if `true`, overwrite existing files
- `kwargs...`: other keyword arguments of the exporter (e.g. `split_size` for FIFF, `markers` for CSV)

# Returns

- `Vector{String}`: names of the written files (FIFF split parts, CSV sidecar files)
"""
function export_recording(
    obj::NeuroAnalyzer.NEURO;
    file_name::String,
    format::Union{Nothing, Symbol} = nothing,
    overwrite::Bool = false,
    kwargs...,
)::Vector{String}
    fmt = if isnothing(format)
        ext = lowercase(splitext(file_name)[2])
        f = get(_EXPORT_EXT_FORMAT, ext, nothing)
        isnothing(f) && throw(
            ArgumentError(
                "Cannot determine export format from extension \"$ext\"; supported: " *
                join(sort(collect(keys(_EXPORT_EXT_FORMAT))), ", ") * " (or use format).",
            ),
        )
        f
    else
        format in _EXPORT_FORMATS || throw(
            ArgumentError("format must be one of: " * join(":" .* string.(_EXPORT_FORMATS), ", ") * "."),
        )
        format
    end

    if fmt === :fiff
        return export_fiff(obj; file_name = file_name, overwrite = overwrite, kwargs...)
    elseif fmt === :bv
        return export_bv(obj; file_name = file_name, overwrite = overwrite, kwargs...)
    elseif fmt === :snirf
        export_snirf(obj; file_name = file_name, overwrite = overwrite, kwargs...)
        return [file_name]
    elseif fmt === :csv
        export_csv(obj; file_name = file_name, overwrite = overwrite, kwargs...)
        base = replace(file_name, r"\.csv$"i => "")
        written = [file_name]
        for (k, sfx) in _CSV_SIDECARS
            get(kwargs, k, false) === true && isfile(base * sfx) && push!(written, base * sfx)
        end
        return written
    else
        isempty(kwargs) || throw(
            ArgumentError("save_neuro() does not accept: " * join(string.(keys(kwargs)), ", ") * "."),
        )
        save_neuro(obj; file_name = file_name, overwrite = overwrite)
        return [file_name]
    end
end
