export import_dtee

# DigiTrack (Elmiko) native exam: an .eef XML header lists the signal (.1) and event (.evx) files

# attribute value of an XML tag string (`nothing` if absent)
function _xml_attr(tag::AbstractString, name::AbstractString)::Union{Nothing, String}
    m = match(Regex("\\b" * name * "=\"([^\"]*)\""), tag)
    return isnothing(m) ? nothing : String(m.captures[1])
end

# resolve the .eef file: an .eef file or a folder containing exactly one .eef file
function _dtee_eef(path::String)::String
    if isdir(path)
        eef = Base.filter(f -> lowercase(splitext(f)[2]) == ".eef", readdir(path))
        isempty(eef) &&
            throw(ArgumentError("Folder $path does not contain a DigiTrack .eef file."))
        length(eef) > 1 && throw(
            ArgumentError(
                "Folder $path contains $(length(eef)) .eef files; pass one of them: $(join(eef, ", ")).",
            ),
        )
        return joinpath(path, eef[1])
    end
    isfile(path) || throw(ArgumentError("File $path cannot be loaded."))
    return path
end

# DigiTrack export manifest (.dtEE) in the exam folder: compressed or encrypted exports are not supported
function _dtee_check_manifest(dir::String)::Nothing
    for f in readdir(dir)
        lowercase(splitext(f)[2]) == ".dtee" || continue
        fp = joinpath(dir, f)
        isfile(fp) || continue
        txt = read(fp, String)
        occursin("ELMIKO_EXPORT_FILE", txt) || continue
        m_c = match(r"<compression\b[^>]*\btype=\"([^\"]*)\"", txt)
        m_e = match(r"<encription\b[^>]*\bvalue=\"([^\"]*)\"", txt)
        (!isnothing(m_c) && m_c.captures[1] != "0") && throw(
            ArgumentError("$fp: compressed DigiTrack exports are not supported."),
        )
        (!isnothing(m_e) && m_e.captures[1] != "0") && throw(
            ArgumentError("$fp: encrypted DigiTrack exports are not supported."),
        )
    end
    return nothing
end

# UTC time (μs since the Unix epoch) of the local date-time `t`, using the UTC offset of an
# event wall-clock time (e.g. "2024-04-09T09:37:56.151+02:00"); `nothing` if no offset found
function _dtee_utc_us(t::DateTime, evx::AbstractString)::Union{Nothing, Int64}
    m = match(r"<info\b[^>]*\btime=\"[^\"]*([+-])(\d{2}):(\d{2})\"", evx)
    isnothing(m) && return nothing
    off = (parse(Int64, m.captures[2]) * 60 + parse(Int64, m.captures[3])) * 60_000_000
    m.captures[1] == "-" && (off = -off)
    return round(Int64, datetime2unix(t) * 1e6) - off
end

"""
    import_dtee(file_name; <keyword arguments>)

Load a DigiTrack (Elmiko) native EEG exam and return a `NeuroAnalyzer.NEURO` object.

A DigiTrack exam consists of an XML header (`.eef`), a binary signal file (`.1`) and an XML event file (`.evx`); the names of the signal and event files are read from the `.eef` file. Either the `.eef` file or a folder containing it (e.g. an exported `.dtEE` folder) may be given.

The signal file stores 16-bit samples of all channels interleaved, followed by a marker channel. Channel labels, types (EEG, ECG, EOG, EMG) and per-channel resolutions are read from the binary header; signals are converted to μV. The marker channel is imported as a `"Marker"` channel (type `"mrk"`) with the raw 16-bit event codes.

Events are imported as markers. DigiTrack does not store the signal during pauses of the recording, so event times are shifted by the duration of the preceding pauses (pause and resume are placed at the same time point). Event values (e.g. elapsed time of hyperventilation or photic stimulation frequency) are appended to the event name (e.g. `"Stimulation_HyperventilationIsGoing 10"`). Events outside the recording (e.g. impedance checks before the start) are skipped.

Compressed or encrypted exports are not supported.

# Arguments

- `file_name::String`: path to the `.eef` file or to the folder containing it
- `detect_type::Bool=true`: infer the type of channels of unknown type from their labels
- `clean_labels::Bool=true`: normalize channel labels (see `import_recording()`)

# Returns

- `NeuroAnalyzer.NEURO`
"""
function import_dtee(
    file_name::String;
    detect_type::Bool = true,
    clean_labels::Bool = true,
)::NeuroAnalyzer.NEURO
    # ------------------------------------------------------------------ #
    # .eef header                                                         #
    # ------------------------------------------------------------------ #
    eef_file = _dtee_eef(file_name)
    dir = dirname(abspath(eef_file))
    eef = read(eef_file, String)
    occursin("<ElmikoExamFile", eef) ||
        throw(ArgumentError("$eef_file is not a DigiTrack .eef file."))
    _dtee_check_manifest(dir)

    # signal (eeg) and event file names
    comp_file(type) = begin
        m = match(
            Regex(
                "<component\\b[^>]*\\btype=\"" * type *
                "\"[^>]*>.*?<file\\b[^>]*>([^<]+)</file>",
                "s",
            ),
            eef,
        )
        isnothing(m) ? nothing : joinpath(dir, strip(m.captures[1]))
    end
    sig_file = comp_file("eeg")
    evx_file = comp_file("events")
    (isnothing(sig_file) || !isfile(sig_file)) && throw(
        ArgumentError(
            "Signal file $(isnothing(sig_file) ? "(not listed in $eef_file)" : sig_file) not found.",
        ),
    )

    m = match(r"<patient\b[^>]*>", eef)
    sex = isnothing(m) ? nothing : _xml_attr(m.match, "sex")
    age = isnothing(m) ? nothing : _xml_attr(m.match, "age")
    m = match(r"<device\b[^>]*>", eef)
    device = isnothing(m) ? "" : something(_xml_attr(m.match, "name"), "")
    m = match(r"<application\b[^>]*type=\"recorder\"[^>]*>", eef)
    recorder = isnothing(m) ? "" :
        strip(something(_xml_attr(m.match, "name"), "") * " " *
              something(_xml_attr(m.match, "version"), ""))

    # ------------------------------------------------------------------ #
    # binary header of the signal file                                    #
    # ------------------------------------------------------------------ #
    rd16(b, o) = ltoh(reinterpret(UInt16, b[(o + 1):(o + 2)])[1])
    hdr = open(io -> read(io, 4800), sig_file)
    (length(hdr) >= 0x46 && rd16(hdr, 0) == 0xffff && rd16(hdr, 16) == 0xffff) ||
        throw(ArgumentError("$sig_file is not a DigiTrack signal file."))
    sampling_rate = Int64(rd16(hdr, 4))
    ch_n = Int64(rd16(hdr, 6))
    hdr_size = Int64(rd16(hdr, 10))
    (sampling_rate > 0 && ch_n > 0 && hdr_size >= 0x46 + 26 * ch_n) ||
        throw(ArgumentError("$sig_file: invalid DigiTrack header."))
    length(hdr) < hdr_size && (hdr = open(io -> read(io, hdr_size), sig_file))

    # start: date (year, day, month) and time (minute, hour, hundredths, second)
    start = try
        DateTime(
            rd16(hdr, 18), hdr[22], hdr[21], hdr[24], hdr[23], hdr[26], 10 * hdr[25],
        )
    catch
        nothing
    end

    # channel records: type, label, resolution (μV / bit)
    clabels = Vector{String}(undef, ch_n)
    ch_type = Vector{String}(undef, ch_n)
    gain = zeros(ch_n)
    types = Dict(0x01 => "eeg", 0x02 => "ecg", 0x03 => "eog", 0x04 => "emg")
    for idx = 1:ch_n
        o = 0x46 + 26 * (idx - 1)
        rec = hdr[(o + 1):(o + 26)]
        lbl = rec[9:16]
        z = findfirst(==(0x00), lbl)
        clabels[idx] = strip(String(isnothing(z) ? lbl : lbl[1:(z - 1)]))
        ch_type[idx] = get(types, rec[1], "other")
        gain[idx] = Float64(ltoh(reinterpret(Float32, rec[19:22])[1]))
    end

    # ------------------------------------------------------------------ #
    # signal data: ch_n × uint16 + marker (16 bit) per sample             #
    # ------------------------------------------------------------------ #
    data_bytes = filesize(sig_file) - hdr_size
    frame = 2 * (ch_n + 1)
    data_bytes % frame == 0 || throw(
        ArgumentError(
            "$sig_file: data size is not a multiple of $(ch_n + 1) channels × 2 bytes (compressed export?).",
        ),
    )
    n_samples = data_bytes ÷ frame
    raw = Vector{UInt16}(undef, (ch_n + 1) * n_samples)
    open(sig_file) do io
        seek(io, hdr_size)
        read!(io, raw)
    end
    raw .= ltoh.(raw)
    raw = reshape(raw, ch_n + 1, n_samples)
    data = zeros(ch_n + 1, n_samples, 1)
    @inbounds for s = 1:n_samples, ch = 1:ch_n
        data[ch, s, 1] = (Float64(raw[ch, s]) - 32768.0) * gain[ch]
    end
    @views data[ch_n + 1, :, 1] .= raw[ch_n + 1, :]

    clabels = clean_labels ? _clean_labels(clabels) : clabels
    if detect_type && any(==("other"), ch_type)
        guess = _set_channel_types(clabels, "other")
        ch_type = [t == "other" ? guess[i] : t for (i, t) in enumerate(ch_type)]
    end
    push!(clabels, "Marker")
    push!(ch_type, "mrk")
    push!(gain, 1.0)
    units, _ = _resolve_units(fill("μV", ch_n + 1), ch_type)
    ch_n += 1

    # ------------------------------------------------------------------ #
    # events (.evx)                                                       #
    # ------------------------------------------------------------------ #
    duration = n_samples / sampling_rate
    ev_t = Int64[]
    ev_name = String[]
    ev_value = String[]
    evx = ""
    if !isnothing(evx_file) && isfile(evx_file)
        evx = read(evx_file, String)
        for mt in eachmatch(r"<type\b([^>]*)>(.*?)</type>"s, evx)
            sid = something(_xml_attr(mt.captures[1], "strId"), _xml_attr(mt.captures[1], "name"), "event")
            for me in eachmatch(r"<event\b([^>]*?)(?:/>|>(.*?)</event>)"s, mt.captures[2])
                t = _xml_attr(me.captures[1], "time")
                isnothing(t) && continue
                v = ""
                if !isnothing(me.captures[2])
                    mi = match(r"<info\b[^>]*\btype=\"1\"[^>]*>", me.captures[2])
                    !isnothing(mi) && (v = something(_xml_attr(mi.match, "value"), ""))
                end
                push!(ev_t, parse(Int64, t))
                push!(ev_name, sid)
                push!(ev_value, v)
            end
        end
    elseif !isnothing(evx_file)
        _warn("Event file $evx_file not found; markers not imported.")
    end

    m_start = Float64[]
    m_name = String[]
    m_value = String[]
    if !isempty(ev_t)
        p = sortperm(ev_t)
        ev_t, ev_name, ev_value = ev_t[p], ev_name[p], ev_value[p]
        # time of the first sample: the exam start event, otherwise the header start time
        i0 = findfirst(==("Technical_ExamStart"), ev_name)
        t0 = if !isnothing(i0)
            ev_t[i0]
        else
            u = isnothing(start) ? nothing : _dtee_utc_us(start, evx)
            if isnothing(u)
                _warn("Recording start not found in the event file; markers are relative to the first event.")
                ev_t[1]
            else
                u
            end
        end
        # pauses are not stored: shift events by the accumulated pause time
        gap = 0
        pause_t = nothing
        for (t, name, v) in zip(ev_t, ev_name, ev_value)
            tt = isnothing(pause_t) ? t : pause_t
            pos = (tt - t0 - gap) / 1e6
            if name == "Technical_ExamPause"
                isnothing(pause_t) && (pause_t = t)
            elseif name == "Technical_ExamResume" && !isnothing(pause_t)
                gap += t - pause_t
                pause_t = nothing
            end
            (-0.5 / sampling_rate <= pos <= duration) || continue
            push!(m_start, max(pos, 0.0))
            push!(m_name, name)
            push!(m_value, isempty(v) ? name : name * " " * v)
        end
        length(m_start) < length(ev_t) &&
            _info("$(length(ev_t) - length(m_start)) events outside the recording skipped")
    end
    # pauses of the exam are discontinuities of the signal
    markers = DataFrame(
        :type => [
            n in ("Technical_ExamPause", "Technical_ExamResume") ? "boundary" : "annotation"
            for n in m_name
        ],
        :start => m_start,
        :length => zeros(length(m_start)),
        :value => m_value,
        :channel => [String[] for _ = 1:length(m_start)],
    )

    # ------------------------------------------------------------------ #
    # assemble NEURO object                                               #
    # ------------------------------------------------------------------ #
    time_pts = _time_axis(n_samples, sampling_rate)
    epoch_time = _time_axis(n_samples, sampling_rate)
    file_size_mb = round(filesize(sig_file) / 1024^2; digits = 2)

    notes = String[]
    !isnothing(sex) && push!(notes, "sex: $sex")
    !isnothing(age) && push!(notes, "age: $age")
    isempty(recorder) || push!(notes, "recorder: $recorder")

    s = _create_subject(;
        id = "", first_name = "", middle_name = "", last_name = "",
        head_circumference = -1, handedness = "", weight = -1, height = -1,
    )
    r = _create_recording_eeg(;
        data_type = "eeg",
        file_name = eef_file,
        file_size_mb = file_size_mb,
        file_type = "DigiTrack",
        recording = device,
        recording_date = isnothing(start) ? "" : Dates.format(start, "yyyy-mm-dd"),
        recording_time = isnothing(start) ? "" :
            Dates.format(start, "HH:MM:SS") * "." *
            lpad(Dates.millisecond(start) ÷ 10, 2, '0'),
        channel_type = ch_type,
        channel_order = _sort_channels(ch_type),
        reference = _detect_montage(clabels, ch_type, "eeg"),
        clabels = clabels,
        transducers = repeat([""], ch_n),
        units = units,
        prefiltering = repeat([""], ch_n),
        line_frequency = 50, # TODO: make this a keyword argument
        sampling_rate = sampling_rate,
        gain = gain,
        bad_channels = zeros(Bool, ch_n),
    )
    hdr_obj = _create_header(; subject = s, recording = r, notes = join(notes, "\n"))

    locs = _initialize_locs()
    obj = NeuroAnalyzer.NEURO(hdr_obj, String[], markers, locs, time_pts, epoch_time, data)
    _initialize_locs!(obj)

    _info(
        "Imported: " *
        uppercase(obj.header.recording[:data_type]) *
        " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj))" *
        "; $(round(obj.time_pts[end], digits = 2)) s)",
    )

    return obj
end
