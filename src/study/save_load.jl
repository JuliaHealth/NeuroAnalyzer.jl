export save_study
export load_study

"""
    save_study(study; <keyword arguments>)

Save a study (design, cached measures and history; not the recordings) to a `.hdf` file (HDF5-based). The study is stored as plain tables and arrays.

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

    JLD2.jldsave(
        file_name;
        neuroanalyzer_study = string(pkgversion(NeuroAnalyzer)),
        name = study.name,
        notes = study.notes,
        root = study.root,
        design = study.design,
        measures = Dict{String, Any}(string(k) => v for (k, v) in study.measures),
        history = study.history,
    )

    return nothing
end

"""
    load_study(file_name; <keyword arguments>)

Load a study saved by `save_study()`.

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

    d = JLD2.load(file_name)
    haskey(d, "neuroanalyzer_study") || throw(
        ArgumentError(
            "$file_name does not contain a NeuroAnalyzer study (for NEURO objects use load_neuro()).",
        ),
    )
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
        get(d, "notes", ""),
        r,
        d["design"],
        Dict{Symbol, Any}(Symbol(k) => v for (k, v) in d["measures"]),
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
