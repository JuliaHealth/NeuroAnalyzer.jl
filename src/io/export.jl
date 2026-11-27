export export_recording

# all formats export_recording() can dispatch to
const _EXPORT_FORMATS = (:fiff, :bv, :edf, :bdf, :nwb, :snirf, :csv, :neuro)

# data types that can be exported to each signal format (CSV and HDF: all data types)
# - ERP and ERF objects are exported as EEG and MEG signals
# - EDF/BDF: also other physiological signals (sensors, EDA, MEP/EMG); MEG: FIFF only
# - TPT (tapping test): CSV and HDF only
const _EXPORT_DATATYPES = Dict{Symbol, Vector{String}}(
    :fiff  => ["eeg", "ecog", "seeg", "ieeg", "csd", "erp", "meg", "erf"],
    :bv    => ["eeg", "ecog", "seeg", "ieeg", "csd", "erp"],
    :edf   => ["eeg", "ecog", "seeg", "ieeg", "csd", "erp", "sensors", "eda", "mep"],
    :bdf   => ["eeg", "ecog", "seeg", "ieeg", "csd", "erp", "sensors", "eda", "mep"],
    :nwb   => ["eeg", "ecog", "seeg", "ieeg", "csd", "erp"],
    :snirf => ["nirs"],
)

# check the data type of OBJ for an export format
function _check_export_datatype(obj::NeuroAnalyzer.NEURO, fmt::Symbol)::Nothing
    dt = datatype(obj)
    dt in _EXPORT_DATATYPES[fmt] || throw(
        ArgumentError(
            "$(uppercase(dt)) objects cannot be exported to $(uppercase(string(fmt))) " *
            "(data types: $(join(uppercase.(_EXPORT_DATATYPES[fmt]), ", "))); " *
            "all data types can be exported to CSV (export_csv()) and HDF (save_neuro()).",
        ),
    )
    return nothing
end

# extension → export format
const _EXPORT_EXT_FORMAT = Dict{String, Symbol}(
    ".fif" => :fiff,
    ".fiff" => :fiff,
    ".vhdr" => :bv,
    ".edf" => :edf,
    ".bdf" => :bdf,
    ".nwb" => :nwb,
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
| Neuromag FIFF                    | `.fif`, `.fiff` | `:fiff`     | `export_fiff` | EEG, ECoG, sEEG, iEEG, CSD, ERP (as EEG), MEG, ERF (as MEG) |
| BrainVision                      | `.vhdr`         | `:bv`       | `export_bv`   | EEG, ECoG, sEEG, iEEG, CSD, ERP |
| EDF+                             | `.edf`          | `:edf`      | `export_edf`  | EEG, ECoG, sEEG, iEEG, CSD, ERP, sensors, EDA, MEP |
| BDF+                             | `.bdf`          | `:bdf`      | `export_bdf`  | EEG, ECoG, sEEG, iEEG, CSD, ERP, sensors, EDA, MEP |
| Neurodata Without Borders        | `.nwb`          | `:nwb`      | `export_nwb`  | EEG, ECoG, sEEG, iEEG, CSD, ERP |
| SNIRF                            | `.snirf`        | `:snirf`    | `export_snirf`| NIRS                   |
| CSV                              | `.csv`          | `:csv`      | `export_csv`  | all                    |
| NeuroAnalyzer (HDF5)             | `.hdf`          | `:neuro`    | `save_neuro`  | all                    |

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
    elseif fmt === :edf
        export_edf(obj; file_name = file_name, overwrite = overwrite, kwargs...)
        return [file_name]
    elseif fmt === :bdf
        export_bdf(obj; file_name = file_name, overwrite = overwrite, kwargs...)
        return [file_name]
    elseif fmt === :nwb
        export_nwb(obj; file_name = file_name, overwrite = overwrite, kwargs...)
        return [file_name]
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
        bad = setdiff(keys(kwargs), (:compress,))
        isempty(bad) || throw(
            ArgumentError("save_neuro() does not accept: " * join(string.(bad), ", ") * "."),
        )
        save_neuro(obj; file_name = file_name, overwrite = overwrite, kwargs...)
        return [file_name]
    end
end
