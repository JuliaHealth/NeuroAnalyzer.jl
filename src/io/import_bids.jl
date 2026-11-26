export import_bids
export bids_files

# modalities (data type folders and suffixes) of BIDS recordings imported by import_bids()
const _BIDS_MODALITIES = ("eeg", "ieeg", "meg", "nirs")

# extensions of BIDS data files (companion files, e.g. .eeg, .vmrk, .fdt, are not listed)
const _BIDS_DATA_EXT = (".edf", ".bdf", ".vhdr", ".set", ".fif", ".snirf", ".nwb")

# BIDS channel types (channels.tsv) → NeuroAnalyzer channel types
const _BIDS_CH_TYPE = Dict(
    "EEG" => "eeg", "ECOG" => "ecog", "SEEG" => "seeg", "DBS" => "seeg",
    "EOG" => "eog", "VEOG" => "eog", "HEOG" => "eog",
    "ECG" => "ecg", "EKG" => "ecg", "EMG" => "emg",
    "TRIG" => "mrk", "REF" => "ref", "GSR" => "eda", "MISC" => "other",
    "MEGMAG" => "mag", "MEGGRADAXIAL" => "grad", "MEGGRADPLANAR" => "grad",
    "NIRSCWAMPLITUDE" => "nirs_int", "NIRSCWOPTICALDENSITY" => "nirs_od",
    "NIRSCWHBO" => "nirs_hbo", "NIRSCWHBR" => "nirs_hbr",
)

# coordinate systems with the x axis towards the nose and the y axis to the left (ALS); all
# other systems are RAS (x to the right, y towards the nose), as NeuroAnalyzer locations
const _BIDS_ALS = ("ctf", "eeglab", "eeglab-hj", "4dbti", "kityokogawa")

"""
    _bids_entities(file_name)

BIDS entities of a file name, e.g. `sub-01_ses-1_task-rest_run-2_eeg.vhdr` → `sub` = `01`, `ses` = `1`, `task` = `rest`, `run` = `2`; the suffix (`eeg`) is returned as `suffix`.

# Returns

- `Dict{String, String}`
"""
function _bids_entities(file_name::AbstractString)::Dict{String, String}
    b = basename(file_name)
    # extension: everything after the first dot (e.g. .nii.gz)
    i = findfirst('.', b)
    stem = isnothing(i) ? b : b[1:(i - 1)]
    parts = split(stem, '_')
    e = Dict{String, String}()
    for p in parts[1:(end - 1)]
        kv = split(p, '-'; limit = 2)
        length(kv) == 2 && (e[kv[1]] = kv[2])
    end
    if length(parts) >= 1
        last = parts[end]
        if occursin('-', last)
            kv = split(last, '-'; limit = 2)
            e[kv[1]] = kv[2]
        else
            e["suffix"] = last
        end
    end
    return e
end

# dataset root: the closest folder with dataset_description.json above the file (nothing: none)
function _bids_root(file_name::AbstractString)::Union{Nothing, String}
    d = dirname(abspath(file_name))
    while true
        isfile(joinpath(d, "dataset_description.json")) && return d
        p = dirname(d)
        p == d && return nothing
        d = p
    end
end

"""
    _bids_sidecars(file_name, suffix, ext)

Sidecar files of a BIDS data file, by the inheritance principle: files named `<entities>_<suffix><ext>` in the folder of the data file or in a folder above it (up to the dataset root), whose entities are a subset of the entities of the data file.

# Returns

- `Vector{String}`: from the least to the most specific (the dataset root first, the folder of the data file last; within a folder, fewer entities first)
"""
function _bids_sidecars(file_name::AbstractString, suffix::String, ext::String)::Vector{String}
    ent = _bids_entities(file_name)
    root = _bids_root(file_name)
    dirs = String[]
    d = dirname(abspath(file_name))
    while true
        pushfirst!(dirs, d)
        (isnothing(root) || d == root) && break
        p = dirname(d)
        p == d && break
        d = p
    end
    out = String[]
    for d in dirs
        cand = Tuple{Int64, String}[]
        for f in readdir(d)
            endswith(f, "_" * suffix * ext) || f == suffix * ext || continue
            isfile(joinpath(d, f)) || continue
            fe = _bids_entities(f)
            get(fe, "suffix", "") == suffix || continue
            delete!(fe, "suffix")
            all(k -> get(ent, k, nothing) == fe[k], keys(fe)) || continue
            push!(cand, (length(fe), joinpath(d, f)))
        end
        sort!(cand; by = first)
        append!(out, last.(cand))
    end
    return out
end

# the most specific sidecar (nothing: none)
function _bids_sidecar(file_name::AbstractString, suffix::String, ext::String)
    s = _bids_sidecars(file_name, suffix, ext)
    return isempty(s) ? nothing : s[end]
end

# JSON sidecars merged by the inheritance principle (more specific files override)
function _bids_json(file_name::AbstractString, suffix::String)::Dict{String, Any}
    js = Dict{String, Any}()
    for f in _bids_sidecars(file_name, suffix, ".json")
        merge!(js, JSON.parsefile(f; dicttype = Dict{String, Any}))
    end
    return js
end

# (BIDS values such as "F" or "T", e.g. sex or outcome codes, are not read as Booleans)
_bids_tsv(f::AbstractString)::DataFrame = CSV.read(
    f, DataFrame;
    delim = '\t', stringtype = String, missingstring = "n/a",
    truestrings = ["true", "True", "TRUE"], falsestrings = ["false", "False", "FALSE"],
)

"""
    import_bids(file_name; <keyword arguments>)

Import a recording of a BIDS dataset (Brain Imaging Data Structure: EEG, iEEG, MEG, fNIRS) with its sidecar files. The data file is imported by `import_recording()` (EDF, BDF, BrainVision, EEGLAB, FIFF, SNIRF, NWB), then the sidecars replace or complete what the data file contains:

- `*_channels.tsv`: channel types (`type`), bad channels (`status`), high-pass and low-pass frequencies (`low_cutoff`, `high_cutoff`) and notch filters (`notch`); channels are matched by name
- `*_events.tsv`: markers (`onset`, `duration`, `trial_type`, `value`; other columns are kept as extra marker columns, see `import_markers()`)
- `*_electrodes.tsv` and `*_coordsystem.json` (EEG, iEEG): channel locations; coordinates of the ALS systems (`CTF`, `EEGLAB`, `4DBti`, `KitYokogawa`: x towards the nose) are rotated to the NeuroAnalyzer orientation (y towards the nose); the units of the file are kept
- `*_optodes.tsv` and `*_coordsystem.json` (NIRS): optode positions (`x`, `y`, `z`, or `template_x`, `template_y`, `template_z` if the positions are `n/a`), rotated as electrode positions; optodes are matched by name (`optode_labels` of the header), or by their order within the sources and detectors (`type`) if no name matches; the unit of the positions (`NIRSCoordinateUnits`) becomes the length unit of the recording
- iEEG recordings (`*_ieeg.*`): data type `"ecog"`, `"seeg"` or `"ieeg"` (both), from the channel types
- `*_<modality>.json`: line frequency (`PowerLineFrequency`), reference (`EEGReference`, `iEEGReference`), task (`TaskName`, `TaskDescription`: experiment name and notes)

Sidecars are found by the BIDS inheritance principle: a sidecar may be in the folder of the data file or in a folder above it (up to the dataset root, the folder with `dataset_description.json`), and applies if its entities (e.g. `task-rest`) are a subset of the entities of the data file; the most specific file is used (JSON files are merged). The subject ID (`sub`) becomes the subject ID of the header.

# Arguments

- `file_name::String`: BIDS data file, e.g. `sub-01/eeg/sub-01_task-rest_eeg.vhdr`
- `channels::Bool=true`: apply `*_channels.tsv`
- `events::Bool=true`: replace the markers by `*_events.tsv`
- `electrodes::Bool=true`: replace the channel locations by `*_electrodes.tsv` (EEG, iEEG)
- `optodes::Bool=true`: replace the optode positions by `*_optodes.tsv` (NIRS)
- `detect_type::Bool=true`: infer channel type from channel label (channel types of `*_channels.tsv` override it)
- `clean_labels::Bool=true`: normalize channel labels (see `import_recording()`); the names in `*_channels.tsv` and `*_electrodes.tsv` are normalized in the same way

# Returns

- `NeuroAnalyzer.NEURO`

# See also

`bids_files()`, `study_bids()`
"""
function import_bids(
    file_name::String;
    channels::Bool = true,
    events::Bool = true,
    electrodes::Bool = true,
    optodes::Bool = true,
    detect_type::Bool = true,
    clean_labels::Bool = true,
)::NeuroAnalyzer.NEURO
    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    ent = _bids_entities(file_name)
    mod = get(ent, "suffix", "")
    mod in _BIDS_MODALITIES || throw(
        ArgumentError(
            "$(basename(file_name)) is not a BIDS data file (suffix must be one of: $(join(_BIDS_MODALITIES, ", "))).",
        ),
    )
    obj = import_recording(file_name; detect_type, clean_labels)
    rec = obj.header.recording
    eeg_like = datatype(obj) in ("eeg", "ecog", "seeg", "ieeg")
    lbl(v) = clean_labels ? _clean_labels(string.(v)) : String.(strip.(string.(v)))

    # modality JSON
    js = _bids_json(file_name, mod)
    if haskey(js, "PowerLineFrequency") && js["PowerLineFrequency"] isa Real &&
       haskey(rec, :line_frequency)
        rec[:line_frequency] = js["PowerLineFrequency"]
    end
    for k in ("EEGReference", "iEEGReference")
        if haskey(js, k) && js[k] isa AbstractString && js[k] != "n/a" && haskey(rec, :reference)
            rec[:reference] = js[k]
        end
    end
    haskey(js, "TaskName") && (obj.header.experiment[:name] = string(js["TaskName"]))
    haskey(js, "TaskDescription") &&
        (obj.header.experiment[:notes] = string(js["TaskDescription"]))
    haskey(ent, "sub") && (obj.header.subject[:id] = ent["sub"])

    # channels.tsv
    f = channels ? _bids_sidecar(file_name, "channels", ".tsv") : nothing
    if !isnothing(f)
        ct = _bids_tsv(f)
        "name" in names(ct) || throw(ArgumentError("$f has no name column."))
        _bids_apply_channels!(obj, ct, lbl(ct.name); types = eeg_like)
    end
    # iEEG: the data type from the channel types (file formats do not distinguish iEEG)
    if mod == "ieeg" && eeg_like
        t = rec[:channel_type]
        ne, ns = count(==("ecog"), t), count(==("seeg"), t)
        rec[:data_type] = ns == 0 && ne > 0 ? "ecog" : ne == 0 && ns > 0 ? "seeg" : "ieeg"
    end

    # events.tsv
    f = events ? _bids_sidecar(file_name, "events", ".tsv") : nothing
    if !isnothing(f)
        ev = _bids_tsv(f)
        ("onset" in names(ev) || "sample" in names(ev)) ||
            throw(ArgumentError("$f has no onset or sample column."))
        obj.markers = sort!(_bids_events(ev, sr(obj)), :start)
    end

    # electrodes.tsv
    f = electrodes && eeg_like ? _bids_sidecar(file_name, "electrodes", ".tsv") : nothing
    if !isnothing(f)
        el = _bids_tsv(f)
        if all(c -> c in names(el), ("name", "x", "y", "z"))
            cs = _bids_coordsystem(file_name)
            obj.locs = _bids_locs(el, lbl(el.name), cs, labels(obj))
            if DataFrames.nrow(obj.locs) > 0
                sys, u = _bids_coord(cs, ("EEG", "iEEG"))
                _set_locs_coord!(obj, sys, u)
            end
        else
            _warn("$f has no name, x, y and z columns: channel locations are not imported.")
        end
    end

    # optodes.tsv
    f = optodes && datatype(obj) == "nirs" ? _bids_sidecar(file_name, "optodes", ".tsv") : nothing
    if !isnothing(f)
        ot = _bids_tsv(f)
        has(c) = all(x -> x in names(ot), c)
        if "name" in names(ot) && (has(("x", "y", "z")) || has(("template_x", "template_y", "template_z")))
            _bids_optodes!(obj, ot, _bids_coordsystem(file_name))
        else
            _warn("$f has no name and x, y, z (or template_x, template_y, template_z) columns: optode positions are not imported.")
        end
    end

    return obj
end

# coordsystem.json merged by the inheritance principle
function _bids_coordsystem(file_name::AbstractString)::Dict{String, Any}
    cs = Dict{String, Any}()
    for j in _bids_sidecars(file_name, "coordsystem", ".json")
        merge!(cs, JSON.parsefile(j; dicttype = Dict{String, Any}))
    end
    return cs
end

# coordinate system and unit ("m", "cm", "mm", "": n/a) of coordsystem.json, for the
# modality prefixes `mods` (e.g. "EEG" → EEGCoordinateSystem, EEGCoordinateUnits)
function _bids_coord(cs::Dict, mods::Tuple)::Tuple{String, String}
    sys, u = "", ""
    for m in mods
        v = get(cs, m * "CoordinateSystem", nothing)
        v isa AbstractString && (sys = string(v))
        v = get(cs, m * "CoordinateUnits", nothing)
        v isa AbstractString && (u = _locs_unit_name(v))
    end
    return sys, u
end

# coordinates in the NeuroAnalyzer orientation (x: right, y: nose) from a BIDS system
function _bids_ras(x::Vector{Float64}, y::Vector{Float64}, sys::String)
    if lowercase(sys) in _BIDS_ALS
        # ALS (x: nose, y: left) → RAS (x: right, y: nose)
        return -y, x
    end
    lowercase(sys) == "other" &&
        _info("Coordinate system \"Other\": the axes are taken as RAS (x: right, y: nose).")
    return x, y
end

# optode positions from optodes.tsv (NIRS)
function _bids_optodes!(obj::NeuroAnalyzer.NEURO, ot::DataFrame, cs::Dict)::Nothing
    rec = obj.header.recording
    ol = string.(rec[:optode_labels])
    src = string.(rec[:src_labels])
    det = string.(rec[:det_labels])
    n = DataFrames.nrow(ot)
    nm = string.(ot.name)
    cols = names(ot)
    val(v) = ismissing(v) ? nothing : v isa Real ? Float64(v) : tryparse(Float64, string(v))
    # position of row i: x, y, z, or the template position if x, y, z are n/a
    function pos(i)
        for (a, b, c) in (("x", "y", "z"), ("template_x", "template_y", "template_z"))
            all(in(cols), (a, b, c)) || continue
            p = val.((ot[i, a], ot[i, b], ot[i, c]))
            any(isnothing, p) || return collect(Float64, p)
        end
        return nothing
    end

    # rows of the optodes: by name, otherwise by the order of the sources and detectors
    row = [findfirst(==(lowercase(l)), lowercase.(nm)) for l in ol]
    if all(isnothing, row)
        typ = "type" in cols ? [lowercase(string(coalesce(t, ""))) for t in ot.type] : fill("", n)
        rs = findall(==("source"), typ)
        rd = findall(==("detector"), typ)
        if isempty(rs) && isempty(rd)
            _warn("No optode of optodes.tsv matches the optodes of the recording: optode positions are not imported.")
            return nothing
        end
        _info("Optode names of optodes.tsv do not match the recording: optodes are matched by their order.")
        row = Union{Nothing, Int64}[
            l in src ? get(rs, findfirst(==(l), src), nothing) :
            l in det ? get(rd, findfirst(==(l), det), nothing) : nothing for l in ol
        ]
    else
        um = Base.filter(x -> !(lowercase(x) in lowercase.(ol)), nm)
        isempty(um) || _info("Optode$(_pl(um)) of optodes.tsv not in the recording: $(join(um, ", ")).")
    end
    P = [isnothing(r) ? nothing : pos(r) for r in row]
    ok = .!isnothing.(P)
    any(ok) || (_warn("optodes.tsv has no positions of the optodes of the recording."); return nothing)
    count(ok) < length(ol) &&
        _warn("Positions of $(length(ol) - count(ok)) optode(s) are not in optodes.tsv: their positions are kept.")

    x = [p[1] for p in P[ok]]
    y = [p[2] for p in P[ok]]
    x, y = _bids_ras(x, y, string(get(cs, "NIRSCoordinateSystem", "")))
    locs = DataFrames.nrow(obj.locs) == 0 ? _initialize_locs() : copy(obj.locs)
    for (j, l) in enumerate(ol[ok])
        p = P[ok][j]
        k = findfirst(==(l), locs.label)
        if isnothing(k)
            push!(locs, (label = l, loc_radius = 0.0, loc_theta = 0.0, loc_x = x[j], loc_y = y[j],
                         loc_z = p[3], loc_radius_sph = 0.0, loc_theta_sph = 0.0,
                         loc_phi_sph = 0.0); cols = :subset)
        else
            locs[k, :loc_x], locs[k, :loc_y], locs[k, :loc_z] = x[j], y[j], p[3]
        end
    end
    # optodes in the order of the header
    o = [something(findfirst(==(l), ol), length(ol) + 1) for l in locs.label]
    locs = locs[sortperm(o), :]
    locs_cart2sph!(locs)
    locs_cart2pol!(locs)
    obj.locs = locs
    sys, u = _bids_coord(cs, ("NIRS",))
    lu = string(get(rec, :length_unit, ""))
    count(ok) < length(ol) && !isempty(lu) && !isempty(u) && lu != u &&
        _warn("Units of optodes.tsv ($u) and of the recording ($lu) differ.")
    _set_locs_coord!(obj, sys, isempty(u) ? lu : u)
    return nothing
end

# apply channels.tsv: types (EEG-like data only), status, filters
function _bids_apply_channels!(
    obj::NeuroAnalyzer.NEURO, ct::DataFrame, names_tsv::Vector{String}; types::Bool,
)::Nothing
    rec = obj.header.recording
    l = labels(obj)
    ll = lowercase.(l)
    nyq = sr(obj) / 2
    cols = names(ct)
    num(x) = x isa Real ? Float64(x) : x isa AbstractString ? tryparse(Float64, x) : nothing
    unmatched = String[]
    for (i, nm) in enumerate(names_tsv)
        k = findfirst(==(nm), l)
        isnothing(k) && (k = findfirst(==(lowercase(nm)), ll))
        if isnothing(k)
            push!(unmatched, nm)
            continue
        end
        if types && "type" in cols && !ismissing(ct[i, :type])
            t = get(_BIDS_CH_TYPE, uppercase(string(ct[i, :type])), "other")
            rec[:channel_type][k] = t
        end
        if "status" in cols && !ismissing(ct[i, :status])
            s = lowercase(string(ct[i, :status]))
            s in ("bad", "good") && (rec[:bad_channel][k] = s == "bad")
        end
        # BIDS: low_cutoff is the high-pass and high_cutoff the low-pass frequency
        if "low_cutoff" in cols && !ismissing(ct[i, :low_cutoff]) && haskey(rec, :highpass)
            v = num(ct[i, :low_cutoff])
            !isnothing(v) && isfinite(v) && v >= 0 && (rec[:highpass][k] = min(v, nyq))
        end
        if "high_cutoff" in cols && !ismissing(ct[i, :high_cutoff]) && haskey(rec, :lowpass)
            v = num(ct[i, :high_cutoff])
            !isnothing(v) && isfinite(v) && v > 0 && (rec[:lowpass][k] = min(v, nyq))
        end
        if "notch" in cols && !ismissing(ct[i, :notch]) && haskey(rec, :notch)
            v = num(ct[i, :notch])
            if isnothing(v) && ct[i, :notch] isa AbstractString
                # list of frequencies, e.g. "[50, 100]"
                v = Float64[parse(Float64, m.match) for m in eachmatch(r"[0-9]*\.?[0-9]+", ct[i, :notch])]
            end
            !isnothing(v) && (rec[:notch][k] = Base.filter(x -> 0 < x < nyq, Float64.(vcat(v))))
        end
    end
    types && (rec[:channel_order] = _sort_channels(rec[:channel_type]))
    isempty(unmatched) || _warn(
        "Channel$(_pl(unmatched)) of channels.tsv not in the recording: $(join(unmatched, ", ")).",
    )
    return nothing
end

# channel locations from electrodes.tsv (positions n/a are left out), in the order of `clabels`
function _bids_locs(
    el::DataFrame, names_el::Vector{String}, cs::Dict, clabels::Vector{String},
)::DataFrame
    sys = ""
    for k in ("EEGCoordinateSystem", "iEEGCoordinateSystem")
        haskey(cs, k) && (sys = string(cs[k]))
    end
    ok = [!ismissing(el[i, :x]) && !ismissing(el[i, :y]) && !ismissing(el[i, :z]) for
          i = 1:DataFrames.nrow(el)]
    val(v) = v isa Real ? Float64(v) : parse(Float64, string(v))
    x = [val(v) for v in el.x[ok]]
    y = [val(v) for v in el.y[ok]]
    z = [val(v) for v in el.z[ok]]
    nm = names_el[ok]
    x, y = _bids_ras(x, y, sys)
    locs = DataFrame(
        :label => nm,
        :loc_radius => zeros(length(nm)),
        :loc_theta => zeros(length(nm)),
        :loc_x => x,
        :loc_y => y,
        :loc_z => z,
        :loc_radius_sph => zeros(length(nm)),
        :loc_theta_sph => zeros(length(nm)),
        :loc_phi_sph => zeros(length(nm)),
    )
    # channels of the recording only, in their order
    lc = lowercase.(clabels)
    idx = [findfirst(==(lowercase(n)), lc) for n in nm]
    keep = .!isnothing.(idx)
    locs = locs[keep, :]
    isempty(locs) && return _initialize_locs()
    locs.label = clabels[Int64.(idx[keep])]
    locs = locs[sortperm(Int64.(idx[keep])), :]
    locs_cart2sph!(locs)
    locs_cart2pol!(locs)
    return locs
end

"""
    bids_files(root; <keyword arguments>)

List the recordings of a BIDS dataset (EEG, iEEG, MEG, fNIRS data files that `import_bids()` can import), with their entities and the participant information of `participants.tsv`.

The folders `derivatives`, `sourcedata`, `code` and `stimuli` and hidden folders are not searched.

# Arguments

- `root::String`: dataset folder (with `dataset_description.json`)
- `subject::Union{Nothing, String, Vector{String}}=nothing`: subject(s) (`sub`, without `sub-`)
- `session::Union{Nothing, String, Vector{String}}=nothing`: session(s) (`ses`)
- `task::Union{Nothing, String, Vector{String}}=nothing`: task(s)
- `run::Union{Nothing, String, Vector{String}}=nothing`: run(s)
- `modality::Union{Nothing, String, Vector{String}}=nothing`: modality(ies): `"eeg"`, `"ieeg"`, `"meg"`, `"nirs"`

# Returns

- `DataFrame`: one row per data file: `subject`, `session` (`""`: none), `task`, `acq`, `run`, `modality`, `file` (relative to `root`), followed by the columns of `participants.tsv` (if present; `participant_id` is not repeated)
"""
function bids_files(
    root::String;
    subject::Union{Nothing, String, Vector{String}} = nothing,
    session::Union{Nothing, String, Vector{String}} = nothing,
    task::Union{Nothing, String, Vector{String}} = nothing,
    run::Union{Nothing, String, Vector{String}} = nothing,
    modality::Union{Nothing, String, Vector{String}} = nothing,
)::DataFrame
    isdir(root) || throw(ArgumentError("Folder $root does not exist."))
    isfile(joinpath(root, "dataset_description.json")) ||
        _warn("$root has no dataset_description.json: it may not be a BIDS dataset.")
    skip = ("derivatives", "sourcedata", "code", "stimuli")
    files = String[]
    for (d, dirs, fs) in walkdir(root)
        if d == root
            Base.filter!(x -> !(x in skip), dirs)
        end
        Base.filter!(x -> !startswith(x, "."), dirs)
        for f in fs
            startswith(f, "sub-") || continue
            any(e -> endswith(lowercase(f), e), _BIDS_DATA_EXT) || continue
            get(_bids_entities(f), "suffix", "") in _BIDS_MODALITIES || continue
            push!(files, relpath(joinpath(d, f), root))
        end
    end
    sort!(files)
    ents = _bids_entities.(files)
    df = DataFrame(
        :subject  => [get(e, "sub", "") for e in ents],
        :session  => [get(e, "ses", "") for e in ents],
        :task     => [get(e, "task", "") for e in ents],
        :acq      => [get(e, "acq", "") for e in ents],
        :run      => [get(e, "run", "") for e in ents],
        :modality => [get(e, "suffix", "") for e in ents],
        :file     => files,
    )
    sel(v) = isnothing(v) ? nothing : v isa String ? [v] : v
    for (c, v) in ((:subject, subject), (:session, session), (:task, task), (:run, run),
                   (:modality, modality))
        isnothing(v) && continue
        df = df[in(sel(v)).(df[!, c]), :]
    end

    # participants.tsv
    pf = joinpath(root, "participants.tsv")
    if isfile(pf) && DataFrames.nrow(df) > 0
        p = _bids_tsv(pf)
        if "participant_id" in names(p)
            p.subject = [replace(string(x), r"^sub-" => "") for x in p.participant_id]
            DataFrames.select!(p, DataFrames.Not(:participant_id))
            for c in setdiff(names(p), ["subject"])
                c in names(df) && DataFrames.rename!(p, c => "participant_" * c)
            end
            df = DataFrames.leftjoin(df, p; on = :subject, order = :left)
        end
    end
    return df
end
