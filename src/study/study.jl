export study_create
export study_check
export study_files
export study_overview
export study_bids

# required columns of the design table
const _STUDY_COLS = [:subject, :session, :group, :file, :include]

"""
    study_create(files; <keyword arguments>)

Create a study (group-level analysis) from recordings saved as NeuroAnalyzer `.hdf` files (see `save_neuro()`). Only the file names are stored: the recordings are loaded one at a time when they are processed (`study_apply!()`) or when measures are computed (`study_measure!()`).

Each recording is described by a row of the design table (`study.design`): subject, session, group, file, inclusion and any covariates. The conditions are not part of the design: they are the epoch conditions of each recording (epoch metadata, see `epoch()`).

# Arguments

- `files::Vector{String}`: file names (relative to `root` or absolute)
- `subject::Vector{String}=String[]`: subject of each file; by default the BIDS `sub-<label>` part of the file name, otherwise the file name without extension
- `session::Vector{String}=String[]`: session of each file; by default the BIDS `ses-<label>` part of the file name, otherwise `"1"`
- `group::Vector{String}=String[]`: group of each file; by default `"all"`
- `root::String=pwd()`: base folder of the study
- `name::String=""`: study name
- `covariates::NamedTuple=NamedTuple()`: additional columns of the design table, e.g. `(; age = [34, 41, 29])`

# Returns

- `NeuroAnalyzer.STUDY`

# Examples

```julia
study = study_create(["sub-01.hdf", "sub-02.hdf", "sub-03.hdf"];
                     group = ["patient", "control", "patient"],
                     root = "data/my_study")
```
"""
function study_create(
    files::Vector{String};
    subject::Vector{String} = String[],
    session::Vector{String} = String[],
    group::Vector{String} = String[],
    root::String = pwd(),
    name::String = "",
    covariates::NamedTuple = NamedTuple(),
)::NeuroAnalyzer.STUDY
    isempty(files) && throw(ArgumentError("files must not be empty."))
    n = length(files)
    for (v, nm) in ((subject, "subject"), (session, "session"), (group, "group"))
        isempty(v) || length(v) == n ||
            throw(ArgumentError("Length of $nm must equal the number of files ($n)."))
    end

    bn = [splitext(basename(f))[1] for f in files]
    bids(f, key) = (m = match(Regex("(?:^|_)$key-([A-Za-z0-9]+)"), f); isnothing(m) ? nothing : String(m[1]))
    isempty(subject) && (subject = [something(bids(b, "sub"), b) for b in bn])
    isempty(session) && (session = [something(bids(b, "ses"), "1") for b in bn])
    isempty(group) && (group = fill("all", n))

    design = DataFrame(
        :subject => subject,
        :session => session,
        :group   => group,
        :file    => copy(files),
        :include => trues(n),
    )
    for (k, v) in pairs(covariates)
        k in propertynames(design) && throw(ArgumentError("Covariate $k is a design column."))
        length(v) == n ||
            throw(ArgumentError("Length of covariate $k must equal the number of files ($n)."))
        design[!, k] = v
    end

    return study_create(design; root = root, name = name)
end

"""
    study_create(design; <keyword arguments>)

Create a study (group-level analysis) from a design table: one row per recording with the columns `subject`, `session`, `group` and `file` (relative to `root` or absolute); the optional column `include` (`Bool`) excludes recordings from the analyses, other columns are covariates.

# Arguments

- `design::DataFrame`: design table
- `root::String=pwd()`: base folder of the study
- `name::String=""`: study name

# Returns

- `NeuroAnalyzer.STUDY`
"""
function study_create(
    design::DataFrame; root::String = pwd(), name::String = "",
)::NeuroAnalyzer.STUDY
    d = copy(design)
    for c in (:subject, :session, :group, :file)
        c in propertynames(d) || throw(ArgumentError("Design table must have a $c column."))
        d[!, c] = string.(d[!, c])
    end
    :include in propertynames(d) || (d[!, :include] = trues(DataFrames.nrow(d)))
    d[!, :include] = Bool.(d[!, :include])
    DataFrames.nrow(d) > 0 || throw(ArgumentError("Design table must not be empty."))
    # required columns first, covariates after them
    DataFrames.select!(d, _STUDY_COLS, DataFrames.Not(_STUDY_COLS))

    isdir(root) || throw(ArgumentError("Folder $root does not exist."))
    study = NeuroAnalyzer.STUDY(name, abspath(root), d, Dict{Symbol, Any}(), String[])

    # validate
    allunique(study_files(study; all = true)) ||
        throw(ArgumentError("Each file must be listed only once."))
    dup = nonunique(d, [:subject, :session])
    any(dup) && throw(
        ArgumentError(
            "Each subject and session must be listed only once (duplicated: $(join(d.subject[dup] .* "/" .* d.session[dup], ", "))).",
        ),
    )
    missing_files = Base.filter(!isfile, study_files(study; all = true))
    isempty(missing_files) ||
        throw(ArgumentError("Files not found: $(join(missing_files, ", "))."))

    push!(study.history, "study_create(design; root = $(repr(study.root)), name = $(repr(name)))")

    return study
end

"""
    study_bids(root; <keyword arguments>)

Create a study from a BIDS dataset: the selected recordings (see `bids_files()`) are imported with their sidecars (`import_bids()`) and saved as NeuroAnalyzer `.hdf` files in `out_dir` (keeping the BIDS folders and file names), and a study of these files is created. Subjects and sessions are the BIDS entities `sub` and `ses`; the group and the covariates are taken from `participants.tsv`.

A study has one recording per subject and session, so the recordings must be selected (e.g. by `task` or `run`) so that each subject and session has one recording.

# Arguments

- `root::String`: dataset folder (with `dataset_description.json`)
- `out_dir::String=joinpath(root, "derivatives", "neuroanalyzer")`: folder of the `.hdf` files (the root of the study)
- `subject::Union{Nothing, String, Vector{String}}=nothing`: subject(s); `nothing`: all
- `session::Union{Nothing, String, Vector{String}}=nothing`: session(s); `nothing`: all
- `task::Union{Nothing, String, Vector{String}}=nothing`: task(s); `nothing`: all
- `run::Union{Nothing, String, Vector{String}}=nothing`: run(s); `nothing`: all
- `modality::Union{Nothing, String, Vector{String}}=nothing`: modality(ies) (`"eeg"`, `"ieeg"`, `"meg"`, `"nirs"`); `nothing`: all
- `group::Union{Nothing, Symbol}=nothing`: column of `participants.tsv` with the group of each subject (e.g. `:group`); `nothing`: all subjects in group `"all"`
- `covariates::Union{Nothing, Vector{Symbol}}=nothing`: columns of `participants.tsv` added to the design table; `nothing`: all columns
- `overwrite::Bool=false`: if `true`, import again recordings whose `.hdf` file exists; otherwise the existing files are used
- `name::String=""`: study name; `""`: the `Name` of `dataset_description.json`
- `kwargs`: keyword arguments of `import_bids()` (e.g. `events`, `clean_labels`)

# Returns

- `NeuroAnalyzer.STUDY`

# Examples

```julia
study = study_bids("ds002718"; task = "FaceRecognition", covariates = [:age, :sex])
```
"""
function study_bids(
    root::String;
    out_dir::String = joinpath(root, "derivatives", "neuroanalyzer"),
    subject::Union{Nothing, String, Vector{String}} = nothing,
    session::Union{Nothing, String, Vector{String}} = nothing,
    task::Union{Nothing, String, Vector{String}} = nothing,
    run::Union{Nothing, String, Vector{String}} = nothing,
    modality::Union{Nothing, String, Vector{String}} = nothing,
    group::Union{Nothing, Symbol} = nothing,
    covariates::Union{Nothing, Vector{Symbol}} = nothing,
    overwrite::Bool = false,
    name::String = "",
    kwargs...,
)::NeuroAnalyzer.STUDY
    bf = bids_files(root; subject, session, task, run, modality)
    DataFrames.nrow(bf) > 0 || throw(ArgumentError("No recordings selected in $root."))
    bf.session = [isempty(s) ? "1" : s for s in bf.session]

    # one recording per subject and session
    dup = nonunique(bf, [:subject, :session])
    if any(dup)
        k = unique(bf[dup, [:subject, :session]])
        r = bf[[(x.subject, x.session) in zip(k.subject, k.session) for x in eachrow(bf)], :]
        throw(
            ArgumentError(
                "Several recordings of one subject and session; select them by task, run or modality: " *
                join(r.file, ", ") * ".",
            ),
        )
    end

    # participants.tsv columns
    std_cols = ["subject", "session", "task", "acq", "run", "modality", "file"]
    pcols = setdiff(names(bf), std_cols)
    if !isnothing(group)
        string(group) in pcols ||
            throw(ArgumentError("participants.tsv has no column $group."))
    end
    cov = isnothing(covariates) ? setdiff(pcols, isnothing(group) ? String[] : [string(group)]) :
          string.(covariates)
    miss = setdiff(cov, pcols)
    isempty(miss) || throw(
        ArgumentError(
            "participants.tsv has no column$(_pl(miss)) $(join(miss, ", ")) (columns: $(join(pcols, ", "))).",
        ),
    )
    clash = intersect(cov, string.(_STUDY_COLS))
    isempty(clash) || throw(
        ArgumentError("Covariates must not be named $(join(clash, ", ")) (see group)."),
    )

    # import and save the recordings
    mkpath(out_dir)
    hdf = String[]
    ok = trues(DataFrames.nrow(bf))
    for (k, x) in enumerate(eachrow(bf))
        f = joinpath(splitext(x.file)[1] * ".hdf")
        out = joinpath(out_dir, f)
        push!(hdf, f)
        isfile(out) && !overwrite && continue
        try
            obj = _silent(() -> import_bids(joinpath(root, x.file); kwargs...))
            mkpath(dirname(out))
            _silent(() -> save_neuro(obj; file_name = out, overwrite = true))
            _info("Imported $(x.file)")
        catch e
            ok[k] = false
            _warn("$(x.file) cannot be imported: $(sprint(showerror, e))")
        end
    end
    any(ok) || throw(ArgumentError("No recording could be imported."))

    design = DataFrame(
        :subject => bf.subject[ok],
        :session => bf.session[ok],
        :group   => isnothing(group) ? fill("all", count(ok)) :
                    [ismissing(g) ? "n/a" : string(g) for g in bf[ok, string(group)]],
        :file    => hdf[ok],
    )
    for c in cov
        design[!, c] = bf[ok, c]
    end
    if isempty(name)
        dd = joinpath(root, "dataset_description.json")
        isfile(dd) && (name = string(get(JSON.parsefile(dd), "Name", "")))
    end
    study = study_create(design; root = out_dir, name = name)
    study.history[end] = "study_bids($(repr(root)); out_dir = $(repr(out_dir))" *
        join([", $k = $(repr(v))" for (k, v) in
              pairs((; subject, session, task, run, modality, group, covariates))
              if !isnothing(v)]) * ")"
    return study
end

"""
    study_files(study; <keyword arguments>)

Return the full paths of the study files.

# Arguments

- `study::NeuroAnalyzer.STUDY`
- `all::Bool=false`: if `true`, include the recordings excluded from the analyses (`include = false`)

# Returns

- `Vector{String}`
"""
function study_files(study::NeuroAnalyzer.STUDY; all::Bool = false)::Vector{String}
    rows = all ? (1:DataFrames.nrow(study.design)) : findall(study.design.include)
    return [_study_path(study, f) for f in study.design.file[rows]]
end

# full path of a study file
_study_path(study::NeuroAnalyzer.STUDY, f::String)::String =
    isabspath(f) ? f : joinpath(study.root, f)

# file name stored in the design: relative to the root if it is inside the root
function _study_rel(study::NeuroAnalyzer.STUDY, f::String)::String
    r = relpath(abspath(f), study.root)
    return startswith(r, "..") ? abspath(f) : r
end

function Base.show(io::IO, ::MIME"text/plain", study::NeuroAnalyzer.STUDY)
    d = study.design
    inc = d[d.include, :]
    println(io, "NeuroAnalyzer STUDY", isempty(study.name) ? "" : ": $(study.name)")
    println(io, "  root: $(study.root)")
    println(
        io,
        "  recordings: $(DataFrames.nrow(d)) ($(DataFrames.nrow(inc)) included), " *
        "subjects: $(length(unique(inc.subject))), sessions: $(length(unique(inc.session)))",
    )
    grp = [g * " (" * string(count(==(g), unique(inc[!, [:subject, :group]]).group)) * ")"
           for g in unique(inc.group)]
    println(io, "  groups: $(join(grp, ", "))")
    cov = setdiff(propertynames(d), _STUDY_COLS)
    isempty(cov) || println(io, "  covariates: $(join(string.(cov), ", "))")
    ms = sort(collect(keys(study.measures)))
    print(io, "  measures: $(isempty(ms) ? "none" : join(string.(ms), ", "))")
    return nothing
end

Base.show(io::IO, study::NeuroAnalyzer.STUDY) =
    print(io, "STUDY($(repr(study.name)), $(DataFrames.nrow(study.design)) recordings)")

"""
    study_check(study; <keyword arguments>)

Check the study recordings: load each file and compare it with the first included recording (data type, sampling rate, channel labels and, for epoched recordings, the epoch time axis; continuous recordings may differ in length). Measures of the recordings (`study_measure!()`) can be computed only when they are compatible.

# Arguments

- `study::NeuroAnalyzer.STUDY`
- `ch::Union{String, Vector{String}, Regex}="all"`: channels to compare
- `all::Bool=false`: if `true`, check also the recordings excluded from the analyses

# Returns

- `DataFrame`: one row per recording: `subject`, `session`, `file`, `data_type`, `sr`, `channels`, `epoch_len`, `epochs`, `conditions`, `ok` and `issue` (empty if compatible)
"""
function study_check(
    study::NeuroAnalyzer.STUDY;
    ch::Union{String, Vector{String}, Regex} = "all",
    all::Bool = false,
)::DataFrame
    rows = all ? collect(1:DataFrames.nrow(study.design)) : findall(study.design.include)
    isempty(rows) && throw(ArgumentError("No recordings included."))
    out = DataFrame(
        :subject    => study.design.subject[rows],
        :session    => study.design.session[rows],
        :file       => study.design.file[rows],
        :data_type  => fill("", length(rows)),
        :sr         => zeros(Int64, length(rows)),
        :channels   => zeros(Int64, length(rows)),
        :epoch_len  => zeros(Int64, length(rows)),
        :epochs     => zeros(Int64, length(rows)),
        :conditions => fill("", length(rows)),
        :ok         => falses(length(rows)),
        :issue      => fill("", length(rows)),
    )
    ref = nothing
    ref_ep = 0
    for (k, r) in enumerate(rows)
        f = _study_path(study, study.design.file[r])
        obj = try
            _silent(() -> load_neuro(f))
        catch e
            out.issue[k] = "cannot be loaded: " * sprint(showerror, e)
            continue
        end
        chi = try
            _get_ch(obj, ch, false)
        catch e
            out.issue[k] = sprint(showerror, e)
            continue
        end
        out.data_type[k] = datatype(obj)
        out.sr[k] = sr(obj)
        out.channels[k] = length(chi)
        out.epoch_len[k] = epoch_len(obj)
        out.epochs[k] = nepochs(obj)
        out.conditions[k] = join(epoch_conditions(obj), ", ")
        sig = _study_sig(obj, chi)
        if isnothing(ref)
            ref = sig
            ref_ep = nepochs(obj)
            out.ok[k] = true
        else
            # continuous recordings may differ in length: their time axes are compared
            # only when both are epoched
            out.issue[k] = _study_compat(ref, sig; time = ref_ep > 1 && nepochs(obj) > 1)
            out.ok[k] = isempty(out.issue[k])
        end
    end

    return out
end

"""
    study_overview(study; <keyword arguments>)

Overview of the study recordings: which subjects have data of which sessions and conditions, the number of epochs of each recording and condition, and the bad channels of each recording. Each recording is loaded (as in `study_check()`), so the overview does not require measures; plot: `plot_study_overview()`.

# Arguments

- `study::NeuroAnalyzer.STUDY`
- `all::Bool=true`: if `true`, include the recordings excluded from the analyses (`include = false`)

# Returns

Named tuple:

- `recordings::DataFrame`: one row per recording: `subject`, `session`, `group`, `include`, `file`, `epochs` (number of epochs, without the rejected ones), `rejected` (number of rejected epochs), `n_bad` (number of bad channels), `bad_channels` (labels joined by ", "), `ok` (the file could be loaded)
- `epochs::DataFrame`: one row per recording and condition: `subject`, `session`, `group`, `include`, `condition`, `n_epochs` (without the rejected epochs; 0: no epochs of the condition), `n_rejected`
- `conditions::Vector{String}`: conditions of all recordings (`""`: epochs without a condition)
- `bad::BitMatrix`: bad channels, shape (recording, channel)
- `labels::Vector{String}`: channel labels of `bad` (of all recordings, in order of appearance)
"""
function study_overview(study::NeuroAnalyzer.STUDY; all::Bool = true)::NamedTuple
    rows = all ? collect(1:DataFrames.nrow(study.design)) : findall(study.design.include)
    isempty(rows) && throw(ArgumentError("No recordings included."))
    d = study.design[rows, :]
    nr = length(rows)
    ep = Vector{DataFrame}(undef, nr)
    bad = Vector{Vector{String}}(undef, nr)
    ok = falses(nr)
    labels_all = String[]
    for k = 1:nr
        obj = try
            _silent(() -> load_neuro(_study_path(study, d.file[k])))
        catch e
            _warn("$(d.file[k]) cannot be loaded: $(sprint(showerror, e))")
            ep[k] = DataFrame(:condition => String[], :reject => Bool[])
            bad[k] = String[]
            continue
        end
        ok[k] = true
        ep[k] = _epochs_table(obj)
        l = labels(obj)
        b = get(obj.header.recording, :bad_channel, falses(length(l)))
        bad[k] = l[findall(b)]
        for x in l
            x in labels_all || push!(labels_all, x)
        end
    end
    any(ok) || throw(ArgumentError("No recording could be loaded."))
    conds = unique(reduce(vcat, [string.(e.condition) for e in ep]))

    recs = DataFrame(
        :subject      => d.subject,
        :session      => d.session,
        :group        => d.group,
        :include      => d.include,
        :file         => d.file,
        :epochs       => [count(.!e.reject) for e in ep],
        :rejected     => [count(e.reject) for e in ep],
        :n_bad        => length.(bad),
        :bad_channels => [join(b, ", ") for b in bad],
        :ok           => ok,
    )
    epochs = DataFrame(
        :subject    => repeat(d.subject; inner = length(conds)),
        :session    => repeat(d.session; inner = length(conds)),
        :group      => repeat(d.group; inner = length(conds)),
        :include    => repeat(d.include; inner = length(conds)),
        :condition  => repeat(conds, nr),
        :n_epochs   => [count((ep[k].condition .== c) .& .!ep[k].reject) for k = 1:nr for c in conds],
        :n_rejected => [count((ep[k].condition .== c) .& ep[k].reject) for k = 1:nr for c in conds],
    )
    B = falses(nr, length(labels_all))
    for k = 1:nr, b in bad[k]
        B[k, findfirst(==(b), labels_all)] = true
    end

    return (recordings = recs, epochs = epochs, conditions = conds, bad = B, labels = labels_all)
end

# properties of a recording that must match across the study
_study_sig(obj::NeuroAnalyzer.NEURO, chi::Vector{Int64}) = (
    data_type = datatype(obj),
    sr = sr(obj),
    labels = labels(obj)[chi],
    epoch_time = obj.epoch_time,
)

# description of the differences of two recordings ("" if compatible); with
# `time = false`, the epoch time axes are not compared (spectra)
function _study_compat(ref, sig; time::Bool = true)::String
    iss = String[]
    sig.data_type == ref.data_type ||
        push!(iss, "data type $(sig.data_type) (expected: $(ref.data_type))")
    sig.sr == ref.sr ||
        push!(iss, "sampling rate $(sig.sr) Hz (expected: $(ref.sr) Hz; see resample())")
    if sig.labels != ref.labels
        if sort(sig.labels) == sort(ref.labels)
            push!(iss, "channel order differs (see reorder_channel())")
        else
            d1 = setdiff(ref.labels, sig.labels)
            d2 = setdiff(sig.labels, ref.labels)
            push!(
                iss,
                "channels differ" * (isempty(d1) ? "" : "; missing: $(join(d1, ", "))") *
                (isempty(d2) ? "" : "; extra: $(join(d2, ", "))") *
                " (see reorder_channel(), keep_channel())",
            )
        end
    end
    if time && !(length(sig.epoch_time) == length(ref.epoch_time) &&
                 isapprox(sig.epoch_time, ref.epoch_time; atol = 0.5 / ref.sr))
        push!(
            iss,
            "epoch time axis $(round(sig.epoch_time[1]; digits = 3))…$(round(sig.epoch_time[end]; digits = 3)) s, $(length(sig.epoch_time)) samples (expected: $(round(ref.epoch_time[1]; digits = 3))…$(round(ref.epoch_time[end]; digits = 3)) s, $(length(ref.epoch_time)) samples; see epoch())",
        )
    end
    return join(iss, "; ")
end

# run f() without the NeuroAnalyzer info messages and warnings (e.g. of load_neuro());
# reference-counted, so that calls from several threads restore `verbose` only at the end
const _SILENT_LOCK = ReentrantLock()
const _SILENT_N = Ref(0)
const _SILENT_V = Ref(true)
function _silent(f::Function)
    lock(_SILENT_LOCK) do
        _SILENT_N[] == 0 && (_SILENT_V[] = NeuroAnalyzer.verbose)
        _SILENT_N[] += 1
        NeuroAnalyzer.verbose = false
    end
    try
        return f()
    finally
        lock(_SILENT_LOCK) do
            _SILENT_N[] -= 1
            _SILENT_N[] == 0 && (NeuroAnalyzer.verbose = _SILENT_V[])
        end
    end
end
