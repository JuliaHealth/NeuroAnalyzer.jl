export export_bids

# BIDS version of the exported datasets
const _BIDS_VERSION = "1.10.0"

# BIDS modality of the data types
const _BIDS_MODALITY = Dict(
    "eeg" => "eeg", "erp" => "eeg", "csd" => "eeg",
    "ecog" => "ieeg", "seeg" => "ieeg", "ieeg" => "ieeg",
    "meg" => "meg", "erf" => "meg",
    "nirs" => "nirs",
)

# data formats of the BIDS modalities (the first one is the default)
const _BIDS_FORMATS = Dict(
    "eeg" => [:bv, :edf, :bdf],
    "ieeg" => [:bv, :edf, :bdf, :nwb],
    "meg" => [:fiff],
    "nirs" => [:snirf],
)
const _BIDS_EXT = Dict(:bv => ".vhdr", :edf => ".edf", :bdf => ".bdf", :nwb => ".nwb",
                       :fiff => ".fif", :snirf => ".snirf")

# NeuroAnalyzer channel types → BIDS channel types (channels.tsv)
const _BIDS_EXPORT_TYPE = Dict(
    "eeg" => "EEG", "csd" => "MISC", "ecog" => "ECOG", "seeg" => "SEEG", "ieeg" => "SEEG",
    "eog" => "EOG", "ecg" => "ECG", "emg" => "EMG", "ref" => "REF", "mrk" => "TRIG",
    "eda" => "GSR", "mag" => "MEGMAG", "grad" => "MEGGRADPLANAR", "meg" => "MEGMAG",
    "nirs_int" => "NIRSCWAMPLITUDE", "nirs_od" => "NIRSCWOPTICALDENSITY",
    "nirs_hbo" => "NIRSCWHBO", "nirs_hbr" => "NIRSCWHBR",
)

# BIDS entity label: letters and digits only
function _bids_label(s::AbstractString, what::String)::String
    l = replace(string(s), r"[^A-Za-z0-9]" => "")
    isempty(l) && throw(ArgumentError("$what $(repr(s)) has no letters or digits (BIDS labels are alphanumeric)."))
    return l
end

# BIDS unit (CMIXF, e.g. uV)
_bids_unit(u::AbstractString)::String =
    isempty(u) ? "n/a" : replace(string(u), "μ" => "u", "µ" => "u", "²" => "^2")

# TSV file with n/a for missing and NaN values
function _bids_write_tsv(file::String, df::DataFrame)::Nothing
    cell(v) = ismissing(v) || (v isa AbstractFloat && isnan(v)) ? "n/a" :
              v isa AbstractString ? (isempty(v) ? "n/a" : replace(v, '\t' => ' ', '\n' => ' ')) :
              string(v)
    open(file, "w") do io
        println(io, join(names(df), '\t'))
        for r in eachrow(df)
            println(io, join([cell(r[c]) for c in names(df)], '\t'))
        end
    end
    return nothing
end

_bids_write_json(file::String, d::AbstractDict)::Nothing =
    (write(file, JSON.json(d; pretty = true) * "\n"); nothing)

"""
    export_bids(obj; <keyword arguments>)
    export_bids(study; <keyword arguments>)

Export a recording, or the included recordings of a study, to a BIDS dataset (Brain Imaging Data Structure: EEG, iEEG, MEG, fNIRS). The data file is written in a format of the BIDS modality of the object, with its sidecars:

| Data type                | Modality | Format (`format`)                                        |
| ------------------------ | -------- | -------------------------------------------------------- |
| EEG, ERP, CSD            | `eeg`    | BrainVision (`:bv`, default), EDF (`:edf`), BDF (`:bdf`) |
| ECoG, sEEG, iEEG         | `ieeg`   | BrainVision (`:bv`, default), EDF, BDF, NWB (`:nwb`)     |
| MEG, ERF                 | `meg`    | FIFF (`:fiff`)                                           |
| NIRS                     | `nirs`   | SNIRF (`:snirf`)                                         |

- `sub-<subject>/[ses-<session>/]<modality>/sub-<subject>[_ses-<session>]_task-<task>[_acq-<acq>][_run-<run>]_<modality>.<ext>`
- `*_<modality>.json`: task, sampling rate, line frequency, reference, channel counts, recording type (continuous or epoched) and duration
- `*_channels.tsv`: channel names, BIDS types, units, high-pass and low-pass frequencies (`low_cutoff`, `high_cutoff`), notch, status (bad channels)
- `*_events.tsv`: the markers (`onset`, `duration`; `trial_type`: the value of stimulus markers, the type of other markers; `value`: the event code of imported BIDS events, otherwise the marker value; `sample`, `channel` and the extra marker columns); read back by `import_bids()`
- `*_events.json`: descriptions of the columns of `*_events.tsv`
- `sub-<subject>[_ses-<session>]_electrodes.tsv` and `*_coordsystem.json` (EEG with channel locations, iEEG: always, `n/a` positions without locations; shared by the recordings of the session): the coordinates in the coordinate system and unit of `locs_coord()` (`"Other"` and `n/a` if unknown; systems with x towards the nose, e.g. EEGLAB, CTF, are written in their own orientation); `sub-<subject>[_ses-<session>]_optodes.tsv` (NIRS)
- dataset: `dataset_description.json` and `README` (created if missing), `participants.tsv` (one row per subject, the row of the subject is added or replaced)

Personal data (subject names) are not written; the recording date is written only into the data file (BrainVision: not stored).

A study (`export_bids(study; …)`) is exported recording by recording; the subjects and sessions are those of the design table (the session entity is left out if all recordings have one session), the group and the covariates are written to `participants.tsv`.

# Arguments

- `obj::NeuroAnalyzer.NEURO` or `study::NeuroAnalyzer.STUDY`
- `root::String`: dataset folder (created if needed)
- `subject::String`: subject label (`obj` only; `sub-` is added)
- `session::String=""`: session label (`obj` only; `""`: none)
- `task::String`: task label
- `acq::String=""`, `run::String=""`: acquisition and run labels (`""`: none)
- `format::Union{Nothing, Symbol}=nothing`: data format (see the table); `nothing`: the default of the modality
- `participant::NamedTuple=NamedTuple()`: columns of `participants.tsv` for the subject (`obj` only), e.g. `(; age = 30, sex = "F", group = "control")`
- `name::String=""`: name of the dataset (`dataset_description.json`, if created)
- `authors::Vector{String}=String[]`: authors of the dataset (`dataset_description.json`, if created)
- `license::String=""`: license of the dataset, e.g. `"CC0"` (`dataset_description.json`, if created)
- `overwrite::Bool=false`: if `true`, overwrite existing files of the recording
- `kwargs`: keyword arguments of the data exporter (e.g. `sex`, `age` of `export_nwb()`)

# Returns

- `Vector{String}`: the written files of the recording (data and sidecars; for a study: of all recordings)
"""
function export_bids(
    obj::NeuroAnalyzer.NEURO;
    root::String,
    subject::String,
    session::String = "",
    task::String,
    acq::String = "",
    run::String = "",
    format::Union{Nothing, Symbol} = nothing,
    participant::NamedTuple = NamedTuple(),
    name::String = "",
    authors::Vector{String} = String[],
    license::String = "",
    overwrite::Bool = false,
    kwargs...,
)::Vector{String}
    dt = datatype(obj)
    haskey(_BIDS_MODALITY, dt) ||
        throw(ArgumentError("$(uppercase(dt)) objects cannot be exported to BIDS."))
    mod = _BIDS_MODALITY[dt]
    fmt = isnothing(format) ? _BIDS_FORMATS[mod][1] : format
    fmt in _BIDS_FORMATS[mod] || throw(
        ArgumentError(
            "format must be one of $(join(":" .* string.(_BIDS_FORMATS[mod]), ", ")) for BIDS $mod.",
        ),
    )
    sub = _bids_label(subject, "subject")
    ses = isempty(session) ? "" : _bids_label(session, "session")
    tsk = _bids_label(task, "task")
    acq_l = isempty(acq) ? "" : _bids_label(acq, "acq")
    run_l = isempty(run) ? "" : _bids_label(run, "run")

    # file names
    dir = joinpath(root, "sub-$sub", (isempty(ses) ? String[] : ["ses-$ses"])..., mod)
    base = "sub-$sub" * (isempty(ses) ? "" : "_ses-$ses") * "_task-$tsk" *
           (isempty(acq_l) ? "" : "_acq-$acq_l") * (isempty(run_l) ? "" : "_run-$run_l")
    # electrodes, optodes and the coordinate system apply to all recordings of the session
    sbase = "sub-$sub" * (isempty(ses) ? "" : "_ses-$ses")
    data_file = joinpath(dir, base * "_$mod" * _BIDS_EXT[fmt])
    (isfile(data_file) && !overwrite) && throw(
        ArgumentError("File $data_file cannot be saved, to overwrite use overwrite=true."),
    )
    mkpath(dir)

    # dataset files
    _bids_dataset!(root, name, authors, license)

    # data file (without personal data)
    o = deepcopy(obj)
    for k in (:first_name, :middle_name, :last_name)
        haskey(o.header.subject, k) && (o.header.subject[k] = "")
    end
    o.header.subject[:id] = "sub-$sub"
    files = if fmt === :bv
        export_bv(o; file_name = data_file, overwrite = true, kwargs...)
    elseif fmt === :edf
        export_edf(o; file_name = data_file, overwrite = true, kwargs...)
        [data_file]
    elseif fmt === :bdf
        export_bdf(o; file_name = data_file, overwrite = true, kwargs...)
        [data_file]
    elseif fmt === :nwb
        export_nwb(o; file_name = data_file, overwrite = true, subject_id = "sub-$sub", kwargs...)
        [data_file]
    elseif fmt === :fiff
        export_fiff(o; file_name = data_file, overwrite = true, kwargs...)
    else
        export_snirf(o; file_name = data_file, overwrite = true, kwargs...)
        [data_file]
    end

    # sidecars
    rec = obj.header.recording
    fs = sr(obj)
    ct = string.(rec[:channel_type])
    clabels = labels(obj)
    ep_n = nepochs(obj)
    btype = [get(_BIDS_EXPORT_TYPE, t, "MISC") for t in ct]

    # *_<mod>.json
    js = Dict{String, Any}(
        "TaskName" => task,
        "SamplingFrequency" => fs,
        "PowerLineFrequency" => (haskey(rec, :line_frequency) && rec[:line_frequency] isa Real &&
                                 rec[:line_frequency] > 0) ? rec[:line_frequency] : "n/a",
        "SoftwareFilters" => "n/a",
        "RecordingDuration" => round(signal_len(obj) * ep_n / fs; digits = 6),
        "RecordingType" => ep_n > 1 ? "epoched" : "continuous",
    )
    ep_n > 1 && (js["EpochLength"] = epoch_len(obj) / fs)
    eq = string(get(rec, :recording, ""))
    isempty(eq) || (js["ManufacturersModelName"] = eq)
    ref = string(get(rec, :reference, ""))
    # channel counts of the BIDS types of channels.tsv
    cnt(t...) = count(in(t), btype)
    if mod == "eeg"
        js["EEGReference"] = isempty(ref) ? "n/a" : ref
        js["EEGChannelCount"] = cnt("EEG")
        js["EOGChannelCount"] = cnt("EOG")
        js["ECGChannelCount"] = cnt("ECG")
        js["EMGChannelCount"] = cnt("EMG")
        js["MISCChannelCount"] = cnt("MISC")
        js["TriggerChannelCount"] = cnt("TRIG")
    elseif mod == "ieeg"
        js["iEEGReference"] = isempty(ref) ? "n/a" : ref
        js["ECOGChannelCount"] = cnt("ECOG")
        js["SEEGChannelCount"] = cnt("SEEG")
        js["EEGChannelCount"] = cnt("EEG")
        js["EOGChannelCount"] = cnt("EOG")
        js["ECGChannelCount"] = cnt("ECG")
        js["EMGChannelCount"] = cnt("EMG")
        js["MiscChannelCount"] = cnt("MISC")
        js["TriggerChannelCount"] = cnt("TRIG")
    elseif mod == "meg"
        js["DewarPosition"] = "n/a"
        js["DigitizedLandmarks"] = false
        js["DigitizedHeadPoints"] = false
        js["MEGChannelCount"] = cnt("MEGMAG", "MEGGRADPLANAR")
        js["EEGChannelCount"] = cnt("EEG")
        js["EOGChannelCount"] = cnt("EOG")
        js["ECGChannelCount"] = cnt("ECG")
    else
        js["NIRSChannelCount"] = count(startswith("nirs_"), ct)
        js["NIRSSourceOptodeCount"] = length(rec[:src_labels])
        js["NIRSDetectorOptodeCount"] = length(rec[:det_labels])
    end
    f = joinpath(dir, base * "_$mod.json")
    _bids_write_json(f, js)
    push!(files, f)

    # *_channels.tsv
    nyq = fs / 2
    chs = DataFrame(:name => clabels, :type => btype)
    if mod == "nirs"
        # source, detector and wavelength of the optical channels
        pairs = rec[:optode_pairs]
        wl = rec[:wavelengths]
        wi = rec[:wavelength_index]
        np = size(pairs, 1)
        chs[!, :source] = [i <= np && pairs[i, 1] > 0 ? rec[:src_labels][pairs[i, 1]] : "n/a"
                           for i in eachindex(clabels)]
        chs[!, :detector] = [i <= np && pairs[i, 2] > 0 ? rec[:det_labels][pairs[i, 2]] : "n/a"
                             for i in eachindex(clabels)]
        chs[!, :wavelength_nominal] = [i <= length(wi) && wi[i] > 0 ? Float64(wl[wi[i]]) : NaN
                                       for i in eachindex(clabels)]
    end
    # (NIRS: units after the optode columns)
    chs[!, :units] = _bids_unit.(string.(rec[:unit]))
    if haskey(rec, :highpass)
        chs[!, :low_cutoff] = Float64.(rec[:highpass])
        chs[!, :high_cutoff] = [l < nyq ? Float64(l) : NaN for l in rec[:lowpass]]
        # (notch is not a NIRS column)
        mod == "nirs" ||
            (chs[!, :notch] = [isempty(n) ? "n/a" : join(string.(n), ",") for n in rec[:notch]])
    end
    chs[!, :status] = [b ? "bad" : "good" for b in rec[:bad_channel]]
    f = joinpath(dir, base * "_channels.tsv")
    _bids_write_tsv(f, chs)
    push!(files, f)

    # *_events.tsv
    m = obj.markers
    if DataFrames.nrow(m) > 0 || ep_n > 1
        ev = DataFrame(
            :onset => Float64.(m.start),
            :duration => Float64.(m.length),
            :trial_type => [t == "stim" ? (isempty(v) ? "stim" : string(v)) : t for
                            (t, v) in zip(m.type, m.value)],
            :value => string.(m.value),
            :sample => round.(Int64, m.start .* fs),
        )
        # event codes of imported BIDS events (`bids_value`) are written back as `value`
        if "bids_value" in names(m)
            for (k, v) in enumerate(m.bids_value)
                (ismissing(v) || (v isa AbstractFloat && isnan(v)) || string(v) in ("", "n/a")) ||
                    (ev.value[k] = string(v))
            end
        end
        any(!isempty, m.channel) && (ev[!, :channel] = [join(c, ",") for c in m.channel])
        for c in _marker_extra_cols(m)
            (c in names(ev) || c == "bids_value") && continue
            ev[!, c] = m[!, c]
        end
        # epochs of epoched objects
        if ep_n > 1
            ept = _epochs_table(obj)
            e = DataFrame(
                :onset => [(k - 1) * epoch_len(obj) / fs for k = 1:ep_n],
                :duration => fill(epoch_len(obj) / fs, ep_n),
                :trial_type => fill("epoch", ep_n),
                :value => [isempty(c) ? "epoch" : string(c) for c in ept.condition],
                :sample => [(k - 1) * epoch_len(obj) for k = 1:ep_n],
            )
            ev = vcat(ev, e; cols = :union)
        end
        sort!(ev, :onset)
        f = joinpath(dir, base * "_events.tsv")
        _bids_write_tsv(f, ev)
        push!(files, f)
        # events.json: descriptions of the columns
        desc = Dict(
            "trial_type" => "Condition (value of stimulus markers) or marker type",
            "value" => "Event code or marker value",
            "sample" => "Onset sample (0-based)",
            "channel" => "Channels of the marker (comma-separated)",
        )
        ej = Dict{String, Any}(c => Dict("Description" => get(desc, c, "Marker column $c"))
                               for c in names(ev) if !(c in ("onset", "duration")))
        f = joinpath(dir, base * "_events.json")
        _bids_write_json(f, ej)
        push!(files, f)
    end

    # electrodes / optodes and the coordinate system (iEEG: always, n/a without locations)
    has_locs = DataFrames.nrow(obj.locs) > 0 &&
               any(r -> r.loc_x != 0 || r.loc_y != 0 || r.loc_z != 0, eachrow(obj.locs))
    if has_locs || mod == "ieeg"
        lc = has_locs ? locs_coord(obj) : (system = "", unit = "", guessed = false)
        sys = isempty(lc.system) ? "Other" : lc.system
        if has_locs
            lbl = string.(obj.locs.label)
            x = Float64.(obj.locs.loc_x)
            y = Float64.(obj.locs.loc_y)
            z = Float64.(obj.locs.loc_z)
            if lowercase(sys) in _BIDS_ALS
                # NeuroAnalyzer orientation (RAS) → ALS (x: nose, y: left)
                x, y = y, -x
            end
        else
            lbl = clabels[in.(ct, Ref(("ecog", "seeg", "ieeg")))]
            isempty(lbl) && (lbl = clabels)
            x = fill(NaN, length(lbl))
            y = copy(x)
            z = copy(x)
        end
        pfx = mod == "ieeg" ? "iEEG" : mod == "nirs" ? "NIRS" : "EEG"
        if mod == "nirs"
            src = string.(rec[:src_labels])
            opt = DataFrame(:name => lbl, :type => [l in src ? "source" : "detector" for l in lbl],
                            :x => x, :y => y, :z => z)
            f = joinpath(dir, sbase * "_optodes.tsv")
            _bids_write_tsv(f, opt)
        else
            el = DataFrame(:name => lbl, :x => x, :y => y, :z => z)
            mod == "ieeg" && (el[!, :size] = fill(NaN, DataFrames.nrow(el)))
            f = joinpath(dir, sbase * "_electrodes.tsv")
            _bids_write_tsv(f, el)
        end
        push!(files, f)
        cs = Dict{String, Any}(
            "$(pfx)CoordinateSystem" => sys,
            "$(pfx)CoordinateUnits" => isempty(lc.unit) ? "n/a" : lc.unit,
        )
        sys == "Other" && (cs["$(pfx)CoordinateSystemDescription"] = has_locs ?
            "Coordinates exported by NeuroAnalyzer (x: right, y: anterior, z: superior)" :
            "Electrode positions not available")
        f = joinpath(dir, sbase * "_coordsystem.json")
        _bids_write_json(f, cs)
        push!(files, f)
    end

    # participants.tsv
    _bids_participant!(root, "sub-$sub", participant)

    return files
end

function export_bids(
    study::NeuroAnalyzer.STUDY;
    root::String,
    task::String,
    acq::String = "",
    run::String = "",
    format::Union{Nothing, Symbol} = nothing,
    name::String = "",
    authors::Vector{String} = String[],
    license::String = "",
    overwrite::Bool = false,
    kwargs...,
)::Vector{String}
    d = study.design[study.design.include, :]
    DataFrames.nrow(d) > 0 || throw(ArgumentError("The study has no included recordings."))
    sessions = length(unique(d.session)) > 1
    subj = [_bids_label(s, "subject") for s in d.subject]
    # subjects with different labels must stay different after removing other characters
    for s in unique(subj)
        length(unique(d.subject[subj .== s])) == 1 || throw(
            ArgumentError("Subjects $(join(unique(d.subject[subj .== s]), ", ")) have the same BIDS label $s."),
        )
    end
    cov = setdiff(propertynames(d), [:subject, :session, :group, :file, :include])
    isempty(name) && (name = study.name)
    files = String[]
    for (k, r) in enumerate(eachrow(d))
        obj = _silent(() -> load_neuro(_study_path(study, r.file)))
        part = NamedTuple{Tuple(vcat(:group, cov))}(Tuple(vcat(r.group, [r[c] for c in cov])))
        append!(files, export_bids(obj; root, subject = subj[k],
                                   session = sessions ? string(r.session) : "",
                                   task, acq, run, format, participant = part, name,
                                   authors, license, overwrite, kwargs...))
    end
    return files
end

# dataset_description.json and README of a BIDS dataset (created if missing)
function _bids_dataset!(root::String, name::String, authors::Vector{String} = String[],
                        license::String = "")::Nothing
    mkpath(root)
    f = joinpath(root, "dataset_description.json")
    if !isfile(f)
        d = Dict{String, Any}(
            "Name" => isempty(name) ? "NeuroAnalyzer BIDS export" : name,
            "BIDSVersion" => _BIDS_VERSION,
            "DatasetType" => "raw",
            "GeneratedBy" => [Dict("Name" => "NeuroAnalyzer", "Version" => string(NeuroAnalyzer.VER))],
        )
        isempty(authors) || (d["Authors"] = authors)
        isempty(license) || (d["License"] = license)
        _bids_write_json(f, d)
    end
    f = joinpath(root, "README")
    isfile(f) || isfile(f * ".md") ||
        write(f, (isempty(name) ? "NeuroAnalyzer BIDS export" : name) * "\n\nDataset exported by NeuroAnalyzer.\n")
    return nothing
end

# add or replace the row of a participant in participants.tsv
function _bids_participant!(root::String, pid::String, cols::NamedTuple)::Nothing
    f = joinpath(root, "participants.tsv")
    p = isfile(f) ? _bids_tsv(f) : DataFrame(:participant_id => String[])
    p.participant_id = string.(p.participant_id)
    # (a named tuple keeps the order of the columns)
    row = merge((; participant_id = pid), cols)
    keep = p.participant_id .!= pid
    p = p[keep, :]
    push!(p, row; cols = :union, promote = true)
    _bids_write_tsv(f, p)
    return nothing
end
