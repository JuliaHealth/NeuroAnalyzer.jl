export export_nwb

# NWB version of the schema cached in exported files (nwb_schema/, see
# scripts/nwb_schema_extract.jl)
const _NWB_EXPORT_VERSION = "2.11.0"

# electrode groups of the channel types: name, description, location (the words are those
# used by import_nwb() to infer the channel types)
const _NWB_GROUPS = Dict(
    "eeg" => ("EEG", "EEG scalp electrodes", "scalp"),
    "seeg" => ("sEEG", "sEEG depth electrodes", "intracranial"),
    "ecog" => ("ECoG", "ECoG grid and strip electrodes", "cortical surface"),
    "ieeg" => ("iEEG", "iEEG intracranial electrodes", "intracranial"),
    "eog" => ("EOG", "EOG electrodes", "periocular"),
    "ecg" => ("ECG", "ECG electrodes", "chest"),
    "emg" => ("EMG", "EMG electrodes", "muscle"),
    "ref" => ("reference", "reference electrodes", "scalp"),
)

# factor converting a voltage unit to volts (nothing: not a voltage)
function _nwb_volts(u::AbstractString)::Union{Nothing, Float64}
    v = _normalize_unit(u)
    v == "V" && return 1.0
    v == "mV" && return 1e-3
    v == "μV" && return 1e-6
    v == "nV" && return 1e-9
    return nothing
end

"""
    export_nwb(obj; <keyword arguments>)

Export a recording to a Neurodata Without Borders (NWB 2.x) file (read by pynwb, MNE-Python via pynwb, MatNWB, `import_nwb()`). The DANDI archive requires the subject ID, sex and age (`subject_id`, `sex`, `age`; checked by nwbinspector). The file contains a copy of the NWB schema (`/specifications`, NWB $(_NWB_EXPORT_VERSION)), so it can be read and validated without the schema installed.

- signal: channels with voltage units (EEG, iEEG, EOG, ECG, EMG, reference) as an `ElectricalSeries` (`acquisition/ElectricalSeries`): 32-bit floats in the units of the channels, `channel_conversion` to volts, gzip compression; channels with other units (e.g. markers, CSD) as one `TimeSeries` per unit (`acquisition/auxiliary_<unit>`)
- electrodes table: one electrode group per channel type (EEG, sEEG, ECoG, iEEG, EOG, ECG, EMG, reference), labels, channel locations (x, y, z: the coordinates of `obj.locs`, in the NeuroAnalyzer orientation, i.e. RAS, and units; `NaN` for channels without a location), filters (high-pass, low-pass, notch), bad channels, reference
- markers: a core `EventsTable` (`events/markers`): `timestamp`, `duration`, `value`, `type`, `channel` (labels joined by ","), extra marker columns
- epoched objects: the epochs are written one after another; `intervals/trials` gives the time, condition and rejection of each epoch
- metadata: session start time (recording date and time, UTC; unknown: 1900-01-01), subject ID and species, the header notes, the recording equipment (device)

# Arguments

- `obj::NeuroAnalyzer.NEURO`: EEG, ECoG, sEEG, iEEG, CSD or ERP object
- `file_name::String`: file name (`.nwb`)
- `session_description::String=""`: description of the session; `""`: the first line of the header notes, or the name of the source file
- `subject_id::String=""`: subject ID; `""`: the subject ID of the header (not written if empty)
- `sex::String="U"`: sex of the subject: `"M"`, `"F"`, `"O"` (other) or `"U"` (unknown)
- `age::String=""`: age of the subject as an ISO 8601 duration, e.g. `"P30Y"` (30 years), `"P18Y/"` (18 years or older); `""`: not written
- `species::String="Homo sapiens"`: species of the subject (`""`: not written)
- `experimenter::Vector{String}=String[]`: experimenters (DANDI form: `"Last, First"`)
- `institution::String=""`, `experiment_description::String=""`, `keywords::Vector{String}=String[]`: written if not empty
- `compress::Int64=4`: gzip level of the signal (0–9; 0: not compressed)
- `overwrite::Bool=false`: if `true`, overwrite an existing file

# Returns

- `Nothing`
"""
function export_nwb(
    obj::NeuroAnalyzer.NEURO;
    file_name::String,
    session_description::String = "",
    subject_id::String = "",
    sex::String = "U",
    age::String = "",
    species::String = "Homo sapiens",
    experimenter::Vector{String} = String[],
    institution::String = "",
    experiment_description::String = "",
    keywords::Vector{String} = String[],
    compress::Int64 = 4,
    overwrite::Bool = false,
)::Nothing
    # validate
    _check_export_datatype(obj, :nwb)
    lowercase(splitext(file_name)[2]) == ".nwb" ||
        throw(ArgumentError("File name extension must be .nwb"))
    (isfile(file_name) && !overwrite) && throw(
        ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."),
    )
    0 <= compress <= 9 || throw(ArgumentError("compress must be in [0, 9]."))
    _check_var(sex, ["M", "F", "O", "U"], "sex")
    isempty(age) || occursin(r"^P(\d+Y)?(\d+M)?(\d+W)?(\d+D)?(T(\d+H)?(\d+M)?(\d+S)?)?(/(P.*)?)?$", age) ||
        throw(ArgumentError("age must be an ISO 8601 duration, e.g. \"P30Y\"."))
    rec = obj.header.recording
    ch_n = nchannels(obj)
    clabels = labels(obj)
    ch_type = string.(rec[:channel_type])
    units = string.(rec[:unit])
    fs = sr(obj)
    nyq = fs / 2
    x = reshape(obj.data, ch_n, :)
    any(v -> isfinite(v) && abs(v) > floatmax(Float32), x) &&
        throw(ArgumentError("Values out of the Float32 range cannot be exported."))

    # channels of the ElectricalSeries (voltage units) and of auxiliary TimeSeries
    volts = [_nwb_volts(u) for u in units]
    el = findall(!isnothing, volts)
    aux = findall(isnothing, volts)
    isempty(el) && throw(ArgumentError("OBJ has no channels with voltage units."))

    # session start (UTC), file creation time
    t0 = _fiff_meas_date(string(get(rec, :recording_date, "")), string(get(rec, :recording_time, "")))
    start = isnothing(t0) ? "1900-01-01T00:00:00+00:00" :
            Dates.format(unix2datetime(t0), "yyyy-mm-ddTHH:MM:SS") * "+00:00"
    isnothing(t0) && _info("Recording date and time unknown: the session start is 1900-01-01.")
    created = Dates.format(Dates.now(Dates.UTC), "yyyy-mm-ddTHH:MM:SS") * "+00:00"
    notes = strip(string(obj.header.notes))
    desc = !isempty(session_description) ? session_description :
           !isempty(notes) ? first(split(notes, '\n')) :
           "NeuroAnalyzer export of " * basename(string(get(rec, :file_name, "")))

    isfile(file_name) && rm(file_name)
    HDF5.h5open(file_name, "w") do f
        # root
        _nwb_attrs!(f, "NWBFile"; nwb_version = _NWB_EXPORT_VERSION)
        f["file_create_date"] = [created]
        _h5_vstr!(f, "identifier", _nwb_uuid())
        _h5_vstr!(f, "session_description", desc)
        _h5_vstr!(f, "session_start_time", start)
        _h5_vstr!(f, "timestamps_reference_time", start)
        for g in ("acquisition", "analysis", "processing", "stimulus/presentation",
                  "stimulus/templates", "general/devices", "general/extracellular_ephys")
            HDF5.create_group(f, g)
        end

        # general: notes, subject, device
        gen = f["general"]
        isempty(notes) || _h5_vstr!(gen, "notes", notes)
        isempty(experimenter) || (gen["experimenter"] = experimenter)
        isempty(institution) || _h5_vstr!(gen, "institution", institution)
        isempty(experiment_description) ||
            _h5_vstr!(gen, "experiment_description", experiment_description)
        isempty(keywords) || (gen["keywords"] = keywords)
        subj = HDF5.create_group(gen, "subject")
        _nwb_attrs!(subj, "Subject")
        sid = isempty(subject_id) ? string(get(obj.header.subject, :id, "")) : subject_id
        isempty(sid) || _h5_vstr!(subj, "subject_id", sid)
        _h5_vstr!(subj, "sex", sex)
        isempty(age) || _h5_vstr!(subj, "age", age)
        isempty(species) || _h5_vstr!(subj, "species", species)
        dev = HDF5.create_group(f["general/devices"], "amplifier")
        equipment = string(get(rec, :recording, ""))
        _nwb_attrs!(dev, "Device"; description = isempty(equipment) ? "recording device" : equipment)

        # electrode groups (one per channel type), linked to the device
        gname = String[]
        gloc = String[]
        ee = f["general/extracellular_ephys"]
        for k in el
            n, d, l = get(_NWB_GROUPS, ch_type[k], ("other", "other electrodes", "unknown"))
            push!(gname, n)
            push!(gloc, l)
            haskey(ee, n) && continue
            g = HDF5.create_group(ee, n)
            _nwb_attrs!(g, "ElectrodeGroup"; description = d, location = l)
            HDF5.API.h5l_create_soft("/general/devices/amplifier", g, "device",
                                     HDF5.API.H5P_DEFAULT, HDF5.API.H5P_DEFAULT)
        end

        # electrodes table
        li = Dict(string(l) => i for (i, l) in enumerate(obj.locs.label))
        coord(i, c) = haskey(li, clabels[i]) ? Float64(obj.locs[li[clabels[i]], c]) : NaN
        function filt(i)
            p = String[]
            rec[:highpass][i] > 0 && push!(p, "HP: $(rec[:highpass][i]) Hz")
            rec[:lowpass][i] < nyq && push!(p, "LP: $(rec[:lowpass][i]) Hz")
            for n in rec[:notch][i]
                push!(p, "N: $n Hz")
            end
            return isempty(p) ? "none" : join(p, " ")
        end
        etab = DataFrame(
            :location => gloc,
            :group => ["/general/extracellular_ephys/" * n for n in gname],
            :group_name => gname,
            :label => clabels[el],
            :x => [coord(i, :loc_x) for i in el],
            :y => [coord(i, :loc_y) for i in el],
            :z => [coord(i, :loc_z) for i in el],
            :filtering => [filt(i) for i in el],
            :bad => Vector{Bool}(rec[:bad_channel][el]),
        )
        ref = string(get(rec, :reference, ""))
        isempty(ref) || (etab[!, :reference] = fill(ref, length(el)))
        _nwb_table!(f, ee, "electrodes", etab;
                    neurodata_type = "ElectrodesTable", namespace = "core",
                    description = "metadata about extracellular electrodes",
                    references = ["group"],
                    col_descriptions = Dict(
                        "location" => "Location of the electrode (channel).",
                        "group" => "Reference to the ElectrodeGroup.",
                        "group_name" => "Name of the ElectrodeGroup.",
                        "label" => "Channel label.",
                        "x" => "x coordinate of the channel location (NeuroAnalyzer orientation, RAS: +x right).",
                        "y" => "y coordinate of the channel location (RAS: +y anterior).",
                        "z" => "z coordinate of the channel location (RAS: +z superior).",
                        "filtering" => "Filters of the channel (HP: high-pass, LP: low-pass, N: notch).",
                        "bad" => "Bad channel.",
                        "reference" => "Reference of the channel.",
                    ))

        # signal: ElectricalSeries (volts via channel_conversion)
        es = HDF5.create_group(f["acquisition"], "ElectricalSeries")
        _nwb_attrs!(es, "ElectricalSeries";
                    description = "signal exported by NeuroAnalyzer ($(uppercase(datatype(obj))))",
                    comments = "no comments")
        d = _nwb_dataset!(es, "data", Float32.(x[el, :]); compress = compress)
        a = HDF5.attributes(d)
        _h5_vstr_attr!(d, "unit", "volts")
        a["conversion"] = 1.0
        a["offset"] = 0.0
        a["resolution"] = -1.0
        es["channel_conversion"] = Float64[volts[k] for k in el]
        HDF5.attributes(es["channel_conversion"])["axis"] = 1
        es["starting_time"] = 0.0
        HDF5.attributes(es["starting_time"])["rate"] = Float64(fs)
        _h5_vstr_attr!(es["starting_time"], "unit", "seconds")
        r = _nwb_dataset!(es, "electrodes", collect(0:(length(el) - 1)))
        _nwb_attrs!(r, "DynamicTableRegion"; namespace = "hdmf-common",
                    description = "electrodes of the signal")
        _nwb_ref_attr!(f, r, "table", "/general/extracellular_ephys/electrodes")

        # auxiliary channels: one TimeSeries per unit
        for u in unique(units[aux])
            ks = aux[units[aux] .== u]
            nm = "auxiliary_" * (isempty(u) ? "no_unit" : replace(u, r"[^A-Za-z0-9]" => "_"))
            g = HDF5.create_group(f["acquisition"], nm)
            _nwb_attrs!(g, "TimeSeries";
                        description = "channels " * join(clabels[ks], ", "),
                        comments = "channel types: " * join(ch_type[ks], ", "))
            d = _nwb_dataset!(g, "data", Float32.(x[ks, :]); compress = compress)
            _h5_vstr_attr!(d, "unit", isempty(u) ? "n/a" : u)
            HDF5.attributes(d)["conversion"] = 1.0
            HDF5.attributes(d)["offset"] = 0.0
            HDF5.attributes(d)["resolution"] = -1.0
            g["starting_time"] = 0.0
            HDF5.attributes(g["starting_time"])["rate"] = Float64(fs)
            _h5_vstr_attr!(g["starting_time"], "unit", "seconds")
        end

        # markers: EventsTable
        if DataFrames.nrow(obj.markers) > 0
            m = obj.markers
            ev = DataFrame(
                :timestamp => Float64.(m.start),
                :duration => Float64.(m.length),
                :value => string.(m.value),
                :type => string.(m.type),
                :channel => [join(c, ",") for c in m.channel],
            )
            for c in _marker_extra_cols(m)
                v = m[!, c]
                T = nonmissingtype(eltype(v))
                ev[!, c] = T <: Real && !(T <: Bool) ? Float64[ismissing(z) ? NaN : z for z in v] :
                           T <: Bool ? Bool[coalesce(z, false) for z in v] :
                           String[ismissing(z) ? "" : string(z) for z in v]
            end
            HDF5.create_group(f, "events")
            t = _nwb_table!(f, f["events"], "markers", ev;
                            neurodata_type = "EventsTable", namespace = "core",
                            description = "markers exported by NeuroAnalyzer",
                            col_descriptions = Dict(
                                "value" => "Marker value.", "type" => "Marker type.",
                                "channel" => "Channels of the marker (labels joined by \",\"; empty: all channels)."))
            _h5_vstr_attr!(t, "source_description", "NeuroAnalyzer markers")
            _nwb_attrs!(t["timestamp"], "TimestampVectorData";
                        description = "The time that each event occurred, in seconds, from the session start time.",
                        unit = "seconds")
            _nwb_attrs!(t["duration"], "DurationVectorData";
                        description = "The duration of each event, in seconds. A value of NaN can be used for events without a duration.",
                        unit = "seconds")
        end

        # epochs: trials
        if nepochs(obj) > 1
            ep = DataFrames.nrow(obj.epochs) == nepochs(obj) ? obj.epochs :
                 _epochs_table(obj)
            el_s = epoch_len(obj) / fs
            tr = DataFrame(
                :start_time => [(k - 1) * el_s for k = 1:nepochs(obj)],
                :stop_time => [k * el_s for k = 1:nepochs(obj)],
                :condition => string.(ep.condition),
                :reject => Vector{Bool}(ep.reject),
            )
            HDF5.create_group(f, "intervals")
            _nwb_table!(f, f["intervals"], "trials", tr;
                        neurodata_type = "TimeIntervals", namespace = "core",
                        description = "epochs of the recording",
                        col_descriptions = Dict(
                            "start_time" => "Start time of epoch, in seconds",
                            "stop_time" => "Stop time of epoch, in seconds",
                            "condition" => "Condition of the epoch.",
                            "reject" => "Rejected epoch."))
        end

        # cached schema
        spec = HDF5.create_group(f, "specifications")
        root = joinpath(NeuroAnalyzer.PATH, "nwb_schema")
        for ns in readdir(root)
            isdir(joinpath(root, ns)) || continue
            for v in readdir(joinpath(root, ns))
                g = HDF5.create_group(spec, "$ns/$v")
                for fn in readdir(joinpath(root, ns, v))
                    endswith(fn, ".json") || continue
                    _h5_vstr!(g, fn[1:(end - 5)], read(joinpath(root, ns, v, fn), String))
                end
            end
        end
        _nwb_ref_attr!(f, f, ".specloc", "/specifications")
    end

    return nothing
end
