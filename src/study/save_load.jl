export save_study
export load_study

# NeuroAnalyzer STUDY file layout (HDF5, see save_study())
const _STUDY_FORMAT = "NeuroAnalyzer STUDY"
const _STUDY_FORMAT_VERSION = 1

"""
    save_study(study; <keyword arguments>)

Save a study (design, cached measures and history; not the recordings) to a `.hdf` file (HDF5). The study is stored as plain tables and arrays, readable by any HDF5 software (layout: root attributes `format` ("NeuroAnalyzer STUDY"), `format_version`, `neuroanalyzer_version`; datasets `/name`, `/notes`, `/root`, `/history`; `/design` table and `/measures` dictionary, encoded as in `save_neuro()`).

# Arguments

- `study::NeuroAnalyzer.STUDY`
- `file_name::String`: name of the file to save to
- `overwrite::Bool=false`

# Returns

- `Nothing`
"""
function save_study(
    study::NeuroAnalyzer.STUDY; file_name::String, overwrite::Bool = false,
)::Nothing
    # validate
    (isfile(file_name) && !overwrite) && throw(
        ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."),
    )
    lowercase(splitext(file_name)[2]) == ".hdf" ||
        throw(ArgumentError("Filename extension must be .hdf"))

    # written to a temporary file: an existing file is replaced only if saving succeeds
    tmp_file = file_name * ".tmp"
    try
        HDF5.h5open(tmp_file, "w") do f
            a = HDF5.attributes(f)
            a["format"] = _STUDY_FORMAT
            a["format_version"] = _STUDY_FORMAT_VERSION
            a["neuroanalyzer_version"] = string(pkgversion(NeuroAnalyzer))
            f["name"] = study.name
            f["notes"] = study.notes
            f["root"] = study.root
            _h5_write(f, "design", study.design)
            _h5_write(f, "measures", study.measures)
            f["history"] = Vector{String}(study.history)
        end
        mv(tmp_file, file_name; force = true)
    finally
        rm(tmp_file; force = true)
    end

    return nothing
end

"""
    load_study(file_name; <keyword arguments>)

Load a study saved by `save_study()` (also files saved with JLD2 by NeuroAnalyzer before v0.26.12).

# Arguments

- `file_name::String`: name of the file to load
- `root::Union{Nothing, String}=nothing`: base folder of the study (e.g. if the study was moved); by default the saved one, or the folder of `file_name` if the saved one does not exist

# Returns

- `NeuroAnalyzer.STUDY`
"""
function load_study(
    file_name::String; root::Union{Nothing, String} = nothing,
)::NeuroAnalyzer.STUDY
    # validate
    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    d = _study_read(file_name)
    r = if !isnothing(root)
        isdir(root) || throw(ArgumentError("Folder $root does not exist."))
        abspath(root)
    elseif isdir(d["root"])
        d["root"]
    else
        _warn(
            "Study folder $(d["root"]) does not exist; the folder of $file_name is used.",
        )
        dirname(abspath(file_name))
    end
    study = NeuroAnalyzer.STUDY(
        d["name"],
        # (studies saved before v0.26.12 have no notes)
        d["notes"],
        r,
        d["design"],
        d["measures"],
        d["history"],
    )
    nf = count(!isfile, study_files(study))
    nf > 0 && _warn("$nf study file(s) not found (see study_files()).")
    _info(
        "Loaded study" * (isempty(study.name) ? "" : " $(study.name)") *
        " ($(DataFrames.nrow(study.design)) recordings; measures: " *
        (isempty(study.measures) ? "none" : join(sort(string.(keys(study.measures))), ", ")) *
        ")",
    )

    return study
end

# contents of a study file: name, notes, root, design, measures (Dict{Symbol, Any}), history
function _study_read(file_name::String)::Dict{String, Any}
    HDF5.ishdf5(file_name) ||
        throw(ArgumentError("$file_name is not a NeuroAnalyzer study (not an HDF5 file)."))
    fmt = HDF5.h5open(file_name, "r") do f
        a = HDF5.attributes(f)
        haskey(a, "format") ? string(read(a["format"])) : ""
    end
    fmt == _NEURO_FORMAT && throw(
        ArgumentError(
            "$file_name does not contain a NeuroAnalyzer study (for NEURO objects use load_neuro()).",
        ),
    )
    if fmt == _STUDY_FORMAT
        return HDF5.h5open(file_name, "r") do f
            v = Int64(HDF5.read_attribute(f, "format_version"))
            v > _STUDY_FORMAT_VERSION && throw(
                ArgumentError(
                    "$file_name was saved by a newer version of NeuroAnalyzer (file format version $v); update NeuroAnalyzer to load it.",
                ),
            )
            hist = read(f["history"])
            Dict{String, Any}(
                "name" => read(f["name"]),
                "notes" => read(f["notes"]),
                "root" => read(f["root"]),
                "design" => _h5_read(f["design"]),
                "measures" => Dict{Symbol, Any}(_h5_read(f["measures"])),
                "history" => isempty(hist) ? String[] : Vector{String}(hist),
            )
        end
    end
    fmt == "" || throw(ArgumentError("$file_name is not a NeuroAnalyzer study ($fmt)."))

    # saved with JLD2 (before v0.26.12)
    d = try
        JLD2.load(file_name)
    catch e
        throw(
            ArgumentError(
                "$file_name is not a NeuroAnalyzer study ($(sprint(showerror, e))).",
            ),
        )
    end
    haskey(d, "neuroanalyzer_study") || throw(
        ArgumentError(
            "$file_name does not contain a NeuroAnalyzer study (for NEURO objects use load_neuro()).",
        ),
    )
    d["measures"] = Dict{Symbol, Any}(Symbol(k) => v for (k, v) in d["measures"])
    # (studies saved before v0.26.12 have no notes)
    d["notes"] = get(d, "notes", "")
    return d
end
