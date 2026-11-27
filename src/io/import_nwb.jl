export import_nwb
export nwb_series

"""
    import_nwb(file_name; <keyword arguments>)

Load a signal of a Neurodata Without Borders (NWB 2.x) file and return `NeuroAnalyzer.NEURO` object.

The signal is an NWB `ElectricalSeries`, found by its type anywhere in the file: in `acquisition` (raw data) or in `processing` (e.g. LFP in `processing/ecephys`). If the file contains several series (e.g. EEG and EOG, or raw and LFP data), the first `ElectricalSeries` in `acquisition` is imported, unless another series is selected with `series` (see `nwb_series()`); a plain `TimeSeries` (e.g. auxiliary channels) can also be selected.

Signal:

- data: scaled to the unit of the series (`data × channel_conversion × conversion + offset`); volts are converted to μV
- sampling rate: the rate of the series, or the median interval of its timestamps; irregular timestamps (intervals differing by more than 1% from the median) are resampled to a regular time axis at the nearest integer rate (linear interpolation)

Electrodes table:

- channel labels: the `label` column (otherwise `ch_1`, `ch_2`, …)
- channel types: from the electrode group of each channel (its name, description and location: e.g. depth electrode, shank → sEEG; grid, strip → ECoG; intracranial → iEEG; scalp, EEG → EEG; EOG, ECG, EMG), labels of EOG, ECG and EMG channels take precedence, labels of reference channels (e.g. A1, A2) for scalp electrodes only (in depth electrodes they are contacts of shaft A) (`detect_type`); the data type is ECoG, sEEG, iEEG (both or unspecified intracranial channels) or EEG
- channel locations: the `x`, `y`, `z` columns (channels without finite coordinates have no location); the NWB schema defines the axes as +x posterior, +y inferior, +z right, but human recordings usually store RAS coordinates (e.g. MNI): see `locs_axes`; the unit is not defined by NWB (`locs_coord()` guesses it)
- bad channels: a Bool column `bad` (or `good`)
- filters: the `filtering` column (or the filtering of the series; band ranges such as `0.1-100 Hz` or `0.5 to 1000 Hz` are read as high-pass and low-pass frequencies); reference: the `reference` column (if all channels have the same reference)

Metadata: recording date and time (session start time, UTC), subject ID (`general/subject/subject_id` or the file identifier); subject (sex, age, species, description) and session (description, experiment description, experimenter, institution, lab, session ID, notes) in the header notes.

Events (`events = true`), with times relative to the start of the imported series (markers outside the signal are left out):

- `intervals/trials` and other `TimeIntervals` tables: markers of type `stim`, length: `stop_time − start_time`, value: the `condition` (or `trial_type`, `stimulus`, `label`, `type`) column or the table name; other columns are kept as extra marker columns
- `intervals/epochs`: markers of type `epoch`; `intervals/invalid_times`: bad segments (type `BAD`)
- `EventsTable` (NWB core): markers (`timestamp`, `duration`; value: the first text column or the table name; type: the `type` column, otherwise `response` for tables named after responses or button presses, `stim` for the others; channels: the `channel` column, labels joined by ",")
- `AnnotationSeries`: markers of type `annotation`
- ndx-events 0.2 types: `Events` (value: the name), `LabeledEvents` and `TTLs` (value: the label of the code, otherwise the code), `AnnotatedEventsTable` (one marker per event time, value: the label)
- `stimulus/presentation`: stimulus times (series with timestamps) as markers (value: the series name)

BIDS sidecars of NWB recordings (`*_ieeg.json`, `*_events.tsv`, …) are applied by `import_bids()`.

# Arguments

- `file_name::String`: name of the file to load
- `series::Union{Nothing, String}=nothing`: series to import: its path (e.g. `"/processing/ecephys/LFP/v_trace"`) or name (e.g. `"v_trace"`); `nothing`: the only `ElectricalSeries`, or the first one in `acquisition`
- `detect_type::Bool=true`: detect channel types (electrode groups, labels); otherwise all channels are EEG
- `events::Bool=true`: import events as markers
- `locs_axes::Symbol=:ras`: axes of the electrode coordinates: `:ras` (x: right, y: anterior, z: superior, e.g. MNI) or `:pir` (NWB schema: x: posterior, y: inferior, z: right)

# Returns

- `NeuroAnalyzer.NEURO`

# See also

`nwb_series()`, `import_bids()`

# References

1. Rübel O et al. The Neurodata Without Borders ecosystem for neurophysiological data science. eLife 2022;11:e78362.
"""
function import_nwb(
    file_name::String;
    series::Union{Nothing, String} = nothing,
    detect_type::Bool = true,
    events::Bool = true,
    locs_axes::Symbol = :ras,
)::NeuroAnalyzer.NEURO
    # validate
    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    detect_format(file_name) === :nwb ||
        throw(ArgumentError("$file_name is not an NWB file."))
    _check_var(locs_axes, [:ras, :pir], "locs_axes")

    r = HDF5.h5open(file_name, "r") do f
        sel = _nwb_select_series(f, series)
        s = _nwb_series(f, f[sel.path])
        n = ndims(s.data) == 1 ? 1 : size(s.data, 1)
        el = !isnothing(s.electrodes) && !isempty(s.electrodes_table) &&
             length(s.electrodes) == n ?
             _nwb_electrodes(f, s.electrodes_table, s.electrodes; detect_type, locs_axes) :
             nothing
        meta = _nwb_metadata(f)
        mrk = events ? _nwb_markers(f) : nothing
        (series = s, path = sel.path, electrodes = el, meta = meta, markers = mrk)
    end
    s = r.series

    # data: (channels, time points), volts → μV
    x = s.data
    ndims(x) == 1 && (x = reshape(x, 1, :))
    eltype(x) <: Real || throw(ArgumentError("Series $(r.path) does not contain numeric data."))
    x = Float64.(x)
    lowercase(s.unit) in ("volts", "volt", "v") && (x .*= 1e6)

    # sampling rate; irregular timestamps resampled to a regular time axis
    t0 = s.starting_time
    if !isnothing(s.rate)
        fs = s.rate
    else
        ts = s.timestamps
        (isnothing(ts) || length(ts) < 2) &&
            throw(ArgumentError("Series $(r.path) has neither a sampling rate nor timestamps."))
        dt = diff(ts)
        all(>(0), dt) || throw(ArgumentError("Timestamps of series $(r.path) are not increasing."))
        fs = 1 / median(dt)
        if maximum(abs.(dt .- median(dt))) > 0.01 * median(dt)
            # (resampled to the nearest integer rate)
            fs = max(1, round(Int64, fs))
            x = _nwb_resample(x, ts, fs)
            _info("Irregular timestamps: signal resampled to $fs Hz.")
        end
    end
    data = reshape(x, size(x, 1), size(x, 2), 1)
    data, sampling_rate = _integer_rate(data, fs)
    time_pts = _time_axis(size(data, 2) * size(data, 3), sampling_rate)
    epoch_time = _time_axis(size(data, 2), sampling_rate)

    # channels
    ch_n = size(data, 1)
    el = r.electrodes
    clabels = !isnothing(el) && !isnothing(el.labels) ? el.labels : ["ch_$idx" for idx = 1:ch_n]
    ch_type = if !isnothing(el) && detect_type
        el.types
    elseif detect_type
        _set_channel_types(clabels, "eeg")
    else
        repeat(["eeg"], ch_n)
    end
    data_type = _nwb_datatype(ch_type)
    units = [_ch_units(t) for t in ch_type]
    bad = isnothing(el) ? zeros(Bool, ch_n) : el.bad
    prefiltering = _nwb_filtering.(isnothing(el) || isnothing(el.filtering) ? repeat([s.filtering], ch_n) : el.filtering)
    reference = isnothing(el) ? "" : el.reference

    m = r.meta
    sub = _create_subject(;
        id = m.subject_id, first_name = "", middle_name = "", last_name = "",
        handedness = "", head_circumference = -1, weight = -1, height = -1,
    )
    rec = _create_recording_eeg(;
        data_type = data_type,
        file_name = file_name,
        file_size_mb = round(filesize(file_name) / 1024^2; digits = 2),
        file_type = "NWB",
        recording = m.device,
        recording_date = m.date,
        recording_time = m.time,
        channel_type = ch_type,
        channel_order = _sort_channels(ch_type),
        reference = reference,
        clabels = clabels,
        transducers = repeat([""], ch_n),
        units = units,
        prefiltering = prefiltering,
        line_frequency = 50,
        sampling_rate = sampling_rate,
        gain = ones(ch_n),
        bad_channels = bad,
    )
    hdr = _create_header(;
        subject = sub,
        recording = rec,
        notes = _header_notes(m.notes..., "Series" => r.path, "Series description" => s.description),
    )

    # markers: times relative to the start of the series, within the signal
    markers = DataFrame(
        :type => String[], :start => Float64[], :length => Float64[],
        :value => String[], :channel => Vector{String}[],
    )
    if !isnothing(r.markers) && DataFrames.nrow(r.markers) > 0
        mk = r.markers
        mk.start .-= t0
        dur = size(data, 2) / sampling_rate
        inside = (mk.start .>= 0) .& (mk.start .<= dur)
        out = count(.!inside)
        out > 0 && _info("$out event(s) outside the imported signal are left out.")
        markers = sort!(mk[inside, :], :start)
    end

    obj = NeuroAnalyzer.NEURO(hdr, String[], markers, _initialize_locs(), time_pts, epoch_time, data)
    _initialize_locs!(obj)
    if !isnothing(el) && DataFrames.nrow(el.locs) > 0
        obj.locs = el.locs
    end

    _info(
        "Imported: " * uppercase(obj.header.recording[:data_type]) *
        " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits = 2)) s) from $(r.path)",
    )

    return obj
end

"""
    nwb_series(file_name)

List the signals of a Neurodata Without Borders (NWB) file: its `ElectricalSeries` and other `TimeSeries` with numeric data (e.g. auxiliary or behavioral signals), to select one with `import_nwb(; series)`.

# Arguments

- `file_name::String`

# Returns

- `DataFrame`: one row per series (`ElectricalSeries` first, those in `acquisition` before the others): `path`, `name`, `type` (neurodata type), `channels`, `samples`, `rate` (Hz; `NaN`: timestamps), `unit`, `description`
"""
function nwb_series(file_name::String)::DataFrame
    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    detect_format(file_name) === :nwb ||
        throw(ArgumentError("$file_name is not an NWB file."))
    return HDF5.h5open(_nwb_series_list, file_name, "r")
end

# neurodata types of the signals listed by nwb_series()
const _NWB_SIGNAL_TYPES = ["ElectricalSeries", "TimeSeries"]

# signals of an NWB file (see nwb_series())
function _nwb_series_list(f::HDF5.File)::DataFrame
    o = _nwb_objects(f; types = _NWB_SIGNAL_TYPES)
    out = DataFrame(
        :path => String[], :name => String[], :type => String[], :channels => Int64[],
        :samples => Int64[], :rate => Float64[], :unit => String[], :description => String[],
    )
    for row in eachrow(o)
        g = f[row.path]
        try
            haskey(g, "data") || continue
            d = g["data"]
            sz = size(d)
            T = eltype(d)
            unit = string(something(_nwb_attr(d, "unit"), ""))
            close(d)
            (T <: Real && !isempty(sz)) || continue
            rate = NaN
            if haskey(g, "starting_time")
                st = g["starting_time"]
                rate = Float64(something(_nwb_attr(st, "rate"), NaN))
                close(st)
            end
            desc = string(something(_nwb_attr(g, "description"), ""))
            # NWB (time, channels) is read as (channels, time)
            ch, n = length(sz) == 1 ? (1, sz[1]) : (prod(sz[1:(end - 1)]), sz[end])
            push!(out, (row.path, basename(row.path), row.type, ch, n, rate, unit, desc))
        finally
            close(g)
        end
    end
    # ElectricalSeries first, acquisition before other groups (otherwise in file order)
    key = [(t != "ElectricalSeries", !startswith(p, "/acquisition/")) for
           (t, p) in zip(out.type, out.path)]
    return out[sortperm(key; alg = Base.Sort.DEFAULT_STABLE), :]
end

# series to import: given by its path or name, otherwise the first ElectricalSeries
function _nwb_select_series(f::HDF5.File, series::Union{Nothing, String})::DataFrameRow
    l = _nwb_series_list(f)
    DataFrames.nrow(l) > 0 || throw(ArgumentError("The file contains no signals (TimeSeries)."))
    avail() = join(["$(x.path) ($(x.type), $(x.channels) ch)" for x in eachrow(l)], ", ")
    if isnothing(series)
        es = findall(==("ElectricalSeries"), l.type)
        isempty(es) && throw(
            ArgumentError("The file contains no ElectricalSeries; select a series (available: $(avail())).",),
        )
        length(es) > 1 && _info(
            "$(length(es)) ElectricalSeries in the file, importing $(l.path[es[1]]); others: " *
            join(l.path[es[2:end]], ", ") * " (select with series).",
        )
        return l[es[1], :]
    end
    p = startswith(series, "/") ? series : "/" * series
    k = findfirst(==(p), l.path)
    if isnothing(k)
        kk = findall(==(series), l.name)
        length(kk) > 1 && throw(
            ArgumentError("Several series are named $series: $(join(l.path[kk], ", ")); give the path."),
        )
        k = isempty(kk) ? nothing : kk[1]
    end
    isnothing(k) && throw(ArgumentError("Series $series not found (available: $(avail()))."))
    return l[k, :]
end

# date ("yyyy-mm-dd") and UTC time ("HH:MM:SS") of an ISO 8601 date-time
function _nwb_datetime(t::String)::Tuple{String, String}
    m = match(
        r"^(\d{4})-(\d{2})-(\d{2})(?:[T ](\d{2}):(\d{2}):(\d{2})(?:\.\d+)?(Z|[+-]\d{2}:?\d{2})?)?",
        strip(t),
    )
    isnothing(m) && return ("", "")
    y, mo, d = parse.(Int64, m.captures[1:3])
    isnothing(m.captures[4]) &&
        return (string(y, "-", lpad(mo, 2, '0'), "-", lpad(d, 2, '0')), "")
    dt = Dates.DateTime(y, mo, d, parse.(Int64, m.captures[4:6])...)
    tz = m.captures[7]
    if !isnothing(tz) && tz != "Z"
        sgn = tz[1] == '-' ? -1 : 1
        dt -= Dates.Minute(sgn * (60 * parse(Int64, tz[2:3]) + parse(Int64, tz[(end - 1):end])))
    end
    return (Dates.format(dt, "yyyy-mm-dd"), Dates.format(dt, "HH:MM:SS"))
end

# ---------------------------------------------------------------------------
# electrodes
# ---------------------------------------------------------------------------

"""
    _nwb_electrodes(f, table, rows; <keyword arguments>)

Channel metadata of a series from the rows `rows` of its electrodes table (see `import_nwb()`).

# Returns

Named tuple: `labels` (`nothing`: no `label` column), `types`, `bad`, `filtering` (`nothing`: no `filtering` column), `reference` (`""`: none or several), `locs` (channel locations of the channels with finite coordinates)
"""
function _nwb_electrodes(
    f::HDF5.File,
    table::AbstractString,
    rows::Vector{Int64};
    detect_type::Bool,
    locs_axes::Symbol,
)::NamedTuple
    e = _nwb_table(f, f[table])[rows, :]
    n = DataFrames.nrow(e)
    cols = names(e)
    labels = "label" in cols ? [string(strip(string(l))) for l in e.label] : nothing
    lab = isnothing(labels) ? ["ch_$i" for i = 1:n] : labels

    # electrode group of each channel: its name, description and location
    gtext = fill("", n)
    if "group" in cols
        cache = Dict{String, String}()
        for (i, p) in enumerate(e.group)
            p isa AbstractString || continue
            gtext[i] = get!(cache, p) do
                g = f[p]
                t = join([basename(p), string(something(_nwb_attr(g, "description"), "")),
                          string(something(_nwb_attr(g, "location"), ""))], " ")
                close(g)
                t
            end
        end
    elseif "group_name" in cols
        gtext = string.(e.group_name)
    end
    types = detect_type ? [_nwb_ch_type(lab[i], gtext[i]) for i = 1:n] : fill("eeg", n)

    bad = if "bad" in cols && eltype(e.bad) <: Union{Bool, Missing}
        Bool[coalesce(b, false) for b in e.bad]
    elseif "good" in cols && eltype(e.good) <: Union{Bool, Missing}
        Bool[!coalesce(g, true) for g in e.good]
    else
        zeros(Bool, n)
    end
    filtering = "filtering" in cols ? string.(e.filtering) : nothing
    reference = "reference" in cols && length(unique(e.reference)) == 1 ?
                string(e.reference[1]) : ""

    # channel locations: finite x, y, z
    locs = _initialize_locs()
    if all(c -> c in cols, ("x", "y", "z"))
        val(v) = v isa Real ? Float64(v) : NaN
        x, y, z = val.(e.x), val.(e.y), val.(e.z)
        ok = isfinite.(x) .& isfinite.(y) .& isfinite.(z)
        if any(ok)
            if locs_axes === :pir
                # NWB schema (x: posterior, y: inferior, z: right) → RAS
                x, y, z = z, -x, -y
            end
            k = count(ok)
            locs = DataFrame(
                :label => lab[ok], :loc_radius => zeros(k), :loc_theta => zeros(k),
                :loc_x => x[ok], :loc_y => y[ok], :loc_z => z[ok],
                :loc_radius_sph => zeros(k), :loc_theta_sph => zeros(k), :loc_phi_sph => zeros(k),
            )
            locs_cart2sph!(locs)
            locs_cart2pol!(locs)
        end
    end
    return (labels = labels, types = types, bad = bad, filtering = filtering,
            reference = reference, locs = locs)
end

# channel type from the label and the electrode group (name, description, location): labels
# of EOG, ECG, EMG and marker channels take precedence; labels of reference channels (e.g.
# A1, A2) only for non-intracranial groups (A1, A2, … are contacts of shaft A of depth
# electrodes)
function _nwb_ch_type(label::AbstractString, group::AbstractString)::String
    tl = _set_channel_types([string(label)], "")[1]
    tl in ("eog", "ecg", "emg", "mrk") && return tl
    g = lowercase(group)
    occursin(r"\beog\b|ocular", g) && return "eog"
    occursin(r"\b(ecg|ekg)\b", g) && return "ecg"
    occursin(r"\bemg\b", g) && return "emg"
    occursin(r"seeg|stereo|depth|shank|shaft", g) && return "seeg"
    occursin(r"ecog|grid|strip|subdural|electrocortico", g) && return "ecog"
    occursin(r"ieeg|intracranial", g) && return "ieeg"
    (tl == "ref" || occursin(r"\breference\b", g)) && return "ref"
    occursin(r"\beeg\b|scalp|\bcap\b|10-10|10-20|10-5", g) && return "eeg"
    return isempty(tl) ? "eeg" : tl
end

# data type of a recording from its channel types
function _nwb_datatype(types::Vector{String})::String
    ne, ns, ni = count(==("ecog"), types), count(==("seeg"), types), count(==("ieeg"), types)
    ne + ns + ni == 0 && return "eeg"
    ni == 0 && ns == 0 && return "ecog"
    ni == 0 && ne == 0 && return "seeg"
    return "ieeg"
end

# filter description of NWB (free text, e.g. "0.1-100 Hz", "Passband, 0.5 to 1000 Hz",
# "HP: 0.1 Hz, LP: 70 Hz"): band ranges as "HP: … Hz LP: … Hz" (read by _prefilter_freqs())
function _nwb_filtering(s::AbstractString)::String
    m = match(r"([0-9]*\.?[0-9]+)\s*(?:hz)?\s*(?:-|–|to)\s*([0-9]*\.?[0-9]+)\s*hz"i, s)
    isnothing(m) && return string(s)
    return "HP: $(m[1]) Hz LP: $(m[2]) Hz"
end

# signal (channels × samples) at irregular timestamps `ts` resampled to a regular time axis at
# `fs` Hz from ts[1] (linear interpolation)
function _nwb_resample(x::AbstractMatrix, ts::Vector{Float64}, fs::Real)::Matrix{Float64}
    n = floor(Int64, (ts[end] - ts[1]) * fs) + 1
    t = ts[1] .+ (0:(n - 1)) ./ fs
    out = Matrix{Float64}(undef, size(x, 1), n)
    @inbounds for j = 1:n
        k = clamp(searchsortedlast(ts, t[j]), 1, length(ts) - 1)
        w = (t[j] - ts[k]) / (ts[k + 1] - ts[k])
        for i in axes(x, 1)
            out[i, j] = (1 - w) * x[i, k] + w * x[i, k + 1]
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# metadata
# ---------------------------------------------------------------------------

# recording date and time, subject ID, device and header notes (pairs) of an NWB file
function _nwb_metadata(f::HDF5.File)::NamedTuple
    function rd(p)
        haskey(f, p) || return ""
        v = read(f[p])
        return v isa AbstractVector ? join(string.(v), "; ") : string(v)
    end
    t = rd("session_start_time")
    isempty(t) && (t = first(split(rd("file_create_date"), "; ")))
    date, time = isempty(t) ? ("", "") : _nwb_datetime(String(t))
    subj = rd("general/subject/subject_id")
    isempty(subj) && (subj = rd("identifier"))
    dev = _nwb_objects(f; types = ["Device"])
    device = DataFrames.nrow(dev) == 0 ? "" : begin
        g = f[dev.path[1]]
        d = string(something(_nwb_attr(g, "description"), ""))
        close(g)
        isempty(d) ? basename(dev.path[1]) : basename(dev.path[1]) * " (" * d * ")"
    end
    subject = join(Base.filter(!isempty, [
        (v = rd("general/subject/sex"); isempty(v) ? "" : "sex $v"),
        (v = rd("general/subject/age"); isempty(v) ? "" : "age $v"),
        (v = rd("general/subject/species"); isempty(v) ? "" : "species $v"),
        rd("general/subject/description"),
    ]), ", ")
    notes = [
        "" => rd("session_description"),
        "Experiment" => rd("general/experiment_description"),
        "Session" => rd("general/session_id"),
        "Experimenter" => rd("general/experimenter"),
        "Institution" => rd("general/institution"),
        "Lab" => rd("general/lab"),
        "Subject" => subject,
        "Notes" => rd("general/notes"),
    ]
    return (date = date, time = time, subject_id = subj, device = device, notes = notes)
end

# ---------------------------------------------------------------------------
# events
# ---------------------------------------------------------------------------

# value of an extra marker column (ragged values joined)
_nwb_cell(v) = v isa AbstractVector ? join(string.(v), ",") : v isa NamedTuple ? string(v) : v

# markers (fixed columns + extra columns) from vectors
function _nwb_mk(type, start, len, value, extra::DataFrame = DataFrame())::DataFrame
    n = length(start)
    m = DataFrame(
        :type => type isa AbstractVector ? string.(type) : fill(string(type), n),
        :start => Float64.(start),
        :length => Float64[isfinite(l) ? l : 0.0 for l in len],
        :value => string.(value),
        :channel => [String[] for _ = 1:n],
    )
    for c in names(extra)
        c in names(m) && continue
        m[!, c] = _nwb_cell.(extra[!, c])
    end
    return m
end

"""
    _nwb_markers(f)

Markers of the events of an NWB file (see `import_nwb()`), times relative to the session start.

# Returns

- `DataFrame`
"""
function _nwb_markers(f::HDF5.File)::DataFrame
    out = DataFrame(:type => String[], :start => Float64[], :length => Float64[],
                    :value => String[], :channel => Vector{String}[])
    objs = _nwb_objects(f)
    for o in eachrow(objs)
        name = basename(o.path)
        m = try
            if o.type == "TimeIntervals"
                t = _nwb_table(f, f[o.path])
                DataFrames.nrow(t) == 0 && continue
                skip = ["id", "start_time", "stop_time", "timeseries", "timeseries_index"]
                tags = "tags" in names(t) ? [join(string.(x), ",") for x in t.tags] : fill("", DataFrames.nrow(t))
                len = t.stop_time .- t.start_time
                if name == "invalid_times"
                    _nwb_mk("BAD", t.start_time, len, [isempty(g) ? "BAD" : "BAD_" * g for g in tags])
                elseif name == "epochs"
                    _nwb_mk("epoch", t.start_time, len, [isempty(g) ? "epoch" : g for g in tags],
                            t[:, setdiff(names(t), vcat(skip, "tags"))])
                else
                    vc = findfirst(c -> c in names(t), ["condition", "trial_type", "stimulus", "label", "type"])
                    val = isnothing(vc) ? fill(name, DataFrames.nrow(t)) :
                          string.(t[!, ["condition", "trial_type", "stimulus", "label", "type"][vc]])
                    _nwb_mk("stim", t.start_time, len, val, t[:, setdiff(names(t), skip)])
                end
            elseif o.type == "EventsTable"
                t = _nwb_table(f, f[o.path])
                (DataFrames.nrow(t) == 0 || !("timestamp" in names(t))) && continue
                len = "duration" in names(t) ? Float64.(t.duration) : zeros(DataFrames.nrow(t))
                other = setdiff(names(t), ["id", "timestamp", "duration"])
                vc = findfirst(c -> eltype(t[!, c]) <: AbstractString, other)
                val = isnothing(vc) ? fill(name, DataFrames.nrow(t)) : t[!, other[vc]]
                # marker types: the type column (e.g. written by export_nwb()), otherwise
                # from the table name
                typ = "type" in names(t) ? string.(t.type) :
                      occursin(r"response|button|press|key", lowercase(name)) ? "response" : "stim"
                mk = _nwb_mk(typ, t.timestamp, len, val, t[:, setdiff(other, ["type", "channel"])])
                # channels of the markers (labels joined by ",")
                if "channel" in names(t)
                    mk.channel = [isempty(c) ? String[] : String.(split(c, ",")) for c in string.(t.channel)]
                end
                mk
            elseif o.type == "AnnotationSeries"
                g = f[o.path]
                ts = haskey(g, "timestamps") ? Float64.(read(g["timestamps"])) : Float64[]
                d = haskey(g, "data") ? string.(read(g["data"])) : String[]
                close(g)
                _nwb_mk("annotation", ts, zeros(length(ts)), d)
            elseif o.type == "Events"
                ts = Float64.(read(f[o.path * "/timestamps"]))
                _nwb_mk("stim", ts, zeros(length(ts)), fill(name, length(ts)))
            elseif o.type in ("LabeledEvents", "TTLs")
                ts = Float64.(read(f[o.path * "/timestamps"]))
                d = f[o.path * "/data"]
                codes = Int64.(read(d))
                labs = string.(something(_nwb_attr(d, "labels"), String[]))
                close(d)
                lab(c) = 0 <= c < length(labs) && !isempty(labs[c + 1]) ? labs[c + 1] : string(c)
                _nwb_mk("stim", ts, zeros(length(ts)), lab.(codes))
            elseif o.type == "AnnotatedEventsTable"
                t = _nwb_table(f, f[o.path])
                ("event_times" in names(t) && "label" in names(t)) || continue
                other = setdiff(names(t), ["id", "event_times", "label", "event_times_index"])
                rows = [(i, tm) for i = 1:DataFrames.nrow(t) for tm in t.event_times[i]]
                isempty(rows) && continue
                ri = first.(rows)
                _nwb_mk([_marker_type(string(t.label[i]); default = "annotation") for i in ri],
                        last.(rows), zeros(length(rows)), string.(t.label[ri]), t[ri, other])
            elseif startswith(o.path, "/stimulus/presentation/") && haskey(f, o.path * "/timestamps")
                # stimulus times (stimulus signals with a sampling rate: see import_nwb(; series))
                ts = Float64.(read(f[o.path * "/timestamps"]))
                _nwb_mk("stim", ts, zeros(length(ts)), fill(name, length(ts)))
            else
                continue
            end
        catch err
            _warn("Events of $(o.path) cannot be read: $(sprint(showerror, err))")
            continue
        end
        _markers_append!(out, m)
    end
    return out
end
