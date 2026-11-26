export study_apply!

"""
    study_apply!(study, f; <keyword arguments>)

Process the study recordings: each file is loaded, processed by `f` and saved to `out_dir`; the design table then points to the new files. The recordings are processed one at a time (the processing functions use all threads), so memory use does not grow with the number of recordings. The processing is recorded in the history of each object (see `export_history()`) and of the study.

`f` takes a NEURO object and either modifies it in-place or returns a new NEURO object, e.g. `obj -> (filter!(obj; ...); reference_avg!(obj))` or `obj -> epoch(obj; marker_value = "stim", ep_offset = 0.2, ep_len = 1.0)`.

An error in one recording does not stop the others: its file is not changed and the error is reported in the returned table.

# Arguments

- `study::NeuroAnalyzer.STUDY`: modified in-place
- `f::Function`: processing function
- `out_dir::String`: output folder (relative to the study folder or absolute); files keep their names (and subfolders, if they are in the study folder)
- `overwrite::Bool=false`: if `true`, overwrite existing output files
- `all::Bool=false`: if `true`, process also the recordings excluded from the analyses

# Returns

- `DataFrame`: one row per processed recording: `subject`, `session`, `file` (output), `ok` and `error`

# Notes

The cached measures (`study.measures`) of the previous files are kept, but they are out of date: `study_measure!()` recomputes them.
"""
function study_apply!(
    study::NeuroAnalyzer.STUDY,
    f::Function;
    out_dir::String,
    overwrite::Bool = false,
    all::Bool = false,
)::DataFrame
    rows = all ? collect(1:DataFrames.nrow(study.design)) : findall(study.design.include)
    isempty(rows) && throw(ArgumentError("No recordings to process."))
    od = isabspath(out_dir) ? out_dir : joinpath(study.root, out_dir)

    # output files: the same relative path in the output folder
    src = [_study_path(study, study.design.file[r]) for r in rows]
    dst = map(src) do s
        rel = relpath(abspath(s), study.root)
        startswith(rel, "..") && (rel = basename(s))
        joinpath(od, splitext(rel)[1] * ".hdf")
    end
    allunique(dst) || throw(
        ArgumentError(
            "Output file names are not unique (files from different folders with the same name).",
        ),
    )
    any(abspath.(dst) .== abspath.(src)) &&
        throw(ArgumentError("out_dir must differ from the folder of the input files."))
    if !overwrite
        ex = Base.filter(isfile, dst)
        isempty(ex) || throw(
            ArgumentError(
                "Output files exist: $(join(ex, ", ")); to overwrite use overwrite=true.",
            ),
        )
    end

    ok = falses(length(rows))
    err = fill("", length(rows))
    function _one(k::Int64)
        try
            obj = _silent(() -> load_neuro(src[k]))
            r = f(obj)
            r isa NeuroAnalyzer.NEURO && (obj = r)
            mkpath(dirname(dst[k]))
            _silent(() -> save_neuro(obj; file_name = dst[k], overwrite = true))
            ok[k] = true
        catch e
            err[k] = sprint(showerror, e)
        end
        return nothing
    end
    p = Progress(length(rows); dt = 1, barlen = 20, color = :white, enabled = progress_bar)
    # one recording at a time: the processing functions are multi-threaded
    # (`@threads :static` cannot be nested in a threaded loop)
    for k in eachindex(rows)
        _one(k)
        next!(p)
    end
    finish!(p)

    # the design points to the new files
    for (k, r) in enumerate(rows)
        ok[k] && (study.design.file[r] = _study_rel(study, dst[k]))
    end
    nf = count(.!ok)
    nf > 0 && _warn("$nf of $(length(rows)) recordings could not be processed (see the returned table).")
    push!(
        study.history,
        "study_apply!(study, $(startswith(string(nameof(f)), "#") ? "<anonymous function>" : nameof(f)); out_dir = $(repr(out_dir)), all = $all)" *
        (nf > 0 ? " # $nf failed" : ""),
    )

    return DataFrame(
        :subject => study.design.subject[rows],
        :session => study.design.session[rows],
        :file    => study.design.file[rows],
        :ok      => ok,
        :error   => err,
    )
end
