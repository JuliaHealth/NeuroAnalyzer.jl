export export_edf
export export_bdf

"""
    export_edf(obj; <keyword arguments>)

Export a recording to EDF+ (European Data Format, 16-bit integers; one of the formats of BIDS EEG and iEEG; read by EDFbrowser, EEGLAB, MNE-Python, FieldTrip and most clinical software). See `export_bdf()` for 24-bit BDF+ files.

- signal: each channel is scaled to the 16-bit range between its minimum and maximum (physical minimum and maximum of the header), so the resolution is (maximum − minimum) / 65535 in the unit of the channel; `NaN` values are written as 0
- channels: labels (up to 16 characters), units (e.g. `uV`), transducers and filters (`HP:0.5Hz LP:40Hz N:50Hz`, read back by `import_edf()`)
- markers (EDF+): annotations with the start, length and value of each marker (the value of `BAD` markers starts with `BAD_`, `boundary` markers are written as `boundary`); the marker types, channels and extra columns are not stored by the format (`import_edf()` gives markers the type `annotation`, `BAD` markers the type `BAD`)
- epoched objects: the epochs are written one after another, with an annotation `boundary` at the start of each epoch
- subject and recording: subject ID and name (patient field), recording date and time, equipment; use `anonymize()` to remove personal data before exporting
- the signal is written in data records of `record_duration` seconds; the last record is padded with zeros (the length of the file is a multiple of the record duration)
- channel locations and bad channels are not stored by the format

# Arguments

- `obj::NeuroAnalyzer.NEURO`: EEG, ECoG, sEEG, iEEG, CSD, ERP, sensors, EDA or MEP object (MEG: `export_fiff()`, NIRS: `export_snirf()`)
- `file_name::String`: file name (`.edf`)
- `plus::Bool=true`: write EDF+ (with annotations); `false`: plain EDF (no markers)
- `record_duration::Real=1`: duration of a data record in seconds (`sampling rate × record_duration` must be an integer)
- `overwrite::Bool=false`: if `true`, overwrite an existing file

# Returns

- `Nothing`
"""
function export_edf(
    obj::NeuroAnalyzer.NEURO;
    file_name::String,
    plus::Bool = true,
    record_duration::Real = 1,
    overwrite::Bool = false,
)::Nothing
    _export_edf(obj; file_name, bdf = false, plus, record_duration, overwrite)
    return nothing
end

"""
    export_bdf(obj; <keyword arguments>)

Export a recording to BDF+ (BioSemi Data Format, 24-bit integers; one of the formats of BIDS EEG and iEEG). As `export_edf()`, with a resolution of (maximum − minimum) / 16777215 in the unit of the channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: EEG, ECoG, sEEG, iEEG, CSD, ERP, sensors, EDA or MEP object (MEG: `export_fiff()`, NIRS: `export_snirf()`)
- `file_name::String`: file name (`.bdf`)
- `plus::Bool=true`: write BDF+ (with annotations); `false`: plain BDF (no markers)
- `record_duration::Real=1`: duration of a data record in seconds (`sampling rate × record_duration` must be an integer)
- `overwrite::Bool=false`: if `true`, overwrite an existing file

# Returns

- `Nothing`
"""
function export_bdf(
    obj::NeuroAnalyzer.NEURO;
    file_name::String,
    plus::Bool = true,
    record_duration::Real = 1,
    overwrite::Bool = false,
)::Nothing
    _export_edf(obj; file_name, bdf = true, plus, record_duration, overwrite)
    return nothing
end

# EDF header text: printable ASCII only (other characters as "?"), padded or cut to `n` bytes
function _edf_field(s::AbstractString, n::Int64)::String
    t = String([isascii(c) && isprint(c) ? c : '?' for c in string(s)])
    return rpad(first(t, n), n)
end

# number with at most `n` characters (plain decimal), rounded down (`down = true`: the written
# value is ≤ `x`), up (`false`: ≥ `x`) or to the nearest value (`nothing`)
function _edf_num(x::Real, n::Int64 = 8; down::Union{Nothing, Bool} = nothing)::String
    for d = 6:-1:0
        k = 10^d
        v = (isnothing(down) ? round : down ? floor : ceil)(BigInt, big(Float64(x)) * k)
        a = abs(v)
        s = string(a ÷ k)
        if d > 0
            frac = rstrip(lpad(string(a % k), d, '0'), '0')
            isempty(frac) || (s *= "." * frac)
        end
        v < 0 && (s = "-" * s)
        length(s) <= n && return s
    end
    throw(ArgumentError("Value $x cannot be written in $n characters (EDF/BDF header)."))
end

# patient and recording fields of EDF+ (subfields separated by spaces, spaces in a subfield
# replaced by "_", unknown: "X")
_edf_sub(s::AbstractString)::String = isempty(strip(s)) ? "X" : replace(strip(s), r"\s+" => "_")

function _export_edf(
    obj::NeuroAnalyzer.NEURO;
    file_name::String,
    bdf::Bool,
    plus::Bool,
    record_duration::Real,
    overwrite::Bool,
)::Nothing
    fmt = bdf ? "BDF" : "EDF"
    ext = bdf ? ".bdf" : ".edf"
    # validate
    _check_export_datatype(obj, bdf ? :bdf : :edf)
    lowercase(splitext(file_name)[2]) == ext ||
        throw(ArgumentError("File name extension must be $ext"))
    (isfile(file_name) && !overwrite) && throw(
        ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."),
    )
    record_duration > 0 || throw(ArgumentError("record_duration must be > 0."))
    fs = sr(obj)
    spr = fs * record_duration
    isapprox(spr, round(spr); atol = 1e-9) && spr >= 1 || throw(
        ArgumentError("sampling rate × record_duration must be an integer (sampling rate: $fs Hz)."),
    )
    spr = round(Int64, spr)
    rec = obj.header.recording
    ch_n = nchannels(obj)
    clabels = labels(obj)
    units = string.(rec[:unit])
    bps = bdf ? 3 : 2
    dmin, dmax = bdf ? (-8388608, 8388607) : (-32768, 32767)

    # signal: all epochs one after another, padded to whole records
    x = Float64.(reshape(obj.data, ch_n, :))
    n = size(x, 2)
    n_rec = cld(n, spr)
    if any(isnan, x)
        _warn("NaN values are written as 0.")
        x[isnan.(x)] .= 0
    end
    any(isinf, x) && throw(ArgumentError("Infinite values cannot be exported."))

    # physical range of each channel (written with at most 8 characters, outwards)
    pmin = String[]
    pmax = String[]
    for i = 1:ch_n
        lo, hi = extrema(view(x, i, :))
        lo == hi && (lo -= 1; hi += 1)
        push!(pmin, _edf_num(lo; down = true))
        push!(pmax, _edf_num(hi; down = false))
    end
    plo = parse.(Float64, pmin)
    phi = parse.(Float64, pmax)

    # channel labels (at most 16 characters, unique)
    ll = [first(String([isascii(c) && isprint(c) ? c : '?' for c in l]), 16) for l in clabels]
    any(length.(clabels) .> 16) && _warn("Channel labels longer than 16 characters are cut.")
    allunique(ll) || throw(ArgumentError("Channel labels are not unique after cutting to 16 characters."))

    # filters: high-pass, low-pass and notch frequencies
    nyq = fs / 2
    function prefilt(i)
        p = String[]
        rec[:highpass][i] > 0 && push!(p, "HP:" * _edf_num(rec[:highpass][i], 10) * "Hz")
        rec[:lowpass][i] < nyq && push!(p, "LP:" * _edf_num(rec[:lowpass][i], 10) * "Hz")
        for f in rec[:notch][i]
            push!(p, "N:" * _edf_num(f, 10) * "Hz")
        end
        return join(p, " ")
    end
    transducers = haskey(rec, :transducers) && length(rec[:transducers]) == ch_n ?
                  string.(rec[:transducers]) : fill("", ch_n)

    # ------------------------------------------------------------------ #
    # annotations (EDF+/BDF+): time-stamped annotation lists per record  #
    # ------------------------------------------------------------------ #
    tals = [String[] for _ = 1:n_rec]
    if plus
        # (onsets with up to 6 decimals; "+" sign required)
        tnum(t) = (s = _edf_num(t, 20); startswith(s, "-") ? s : "+" * s)
        ann = DataFrame(:start => Float64[], :length => Float64[], :text => String[])
        if nepochs(obj) > 1
            for e = 1:nepochs(obj)
                push!(ann, ((e - 1) * epoch_len(obj) / fs, 0.0, "boundary"))
            end
        end
        n_ch = 0
        for m in eachrow(obj.markers)
            isempty(m.channel) || (n_ch += 1)
            # (boundary markers: "boundary", read back as boundary markers)
            txt = m.type == "boundary" ? "boundary" : _fiff_an_desc(m.type, string(m.value))
            push!(ann, (m.start, m.length, txt))
        end
        n_ch > 0 && _info("Channels of $n_ch marker(s) are not stored by $fmt+.")
        isempty(_marker_extra_cols(obj.markers)) ||
            _info("Extra marker columns ($(join(_marker_extra_cols(obj.markers), ", "))) are not exported.")
        for a in eachrow(ann)
            r = clamp(floor(Int64, a.start / record_duration) + 1, 1, n_rec)
            txt = replace(a.text, '\x14' => ' ', '\x15' => ' ', '\0' => ' ')
            dur = a.length > 0 ? "\x15" * _edf_num(a.length, 20) : ""
            push!(tals[r], tnum(a.start) * dur * "\x14" * txt * "\x14\0")
        end
    end
    # time-keeping TAL first in each record
    ann_bytes = [codeunits("+" * _edf_num((r - 1) * record_duration, 20) * "\x14\x14\0" *
                           join(tals[r])) for r = 1:n_rec]
    ann_spr = plus ? max(1, cld(maximum(length.(ann_bytes)), bps)) : 0

    # ------------------------------------------------------------------ #
    # header                                                             #
    # ------------------------------------------------------------------ #
    ns = ch_n + (plus ? 1 : 0)
    subj = obj.header.subject
    name = strip(string(get(subj, :last_name, "")) * " " * string(get(subj, :first_name, "")))
    sid = string(get(subj, :id, ""))
    t0 = _fiff_meas_date(string(get(rec, :recording_date, "")), string(get(rec, :recording_time, "")))
    dt = isnothing(t0) ? nothing : unix2datetime(t0)
    months = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]
    equipment = string(get(rec, :recording, ""))
    patient, recording = if plus
        (
            join([_edf_sub(sid), "X", "X", _edf_sub(name)], " "),
            join(["Startdate",
                  isnothing(dt) ? "X" :
                  lpad(day(dt), 2, '0') * "-" * months[month(dt)] * "-" * string(year(dt)),
                  "X", "X", _edf_sub(equipment)], " "),
        )
    else
        (strip(sid * " " * name), equipment)
    end
    date_s = isnothing(dt) ? "01.01.85" :
             lpad(day(dt), 2, '0') * "." * lpad(month(dt), 2, '0') * "." * lpad(year(dt) % 100, 2, '0')
    time_s = isnothing(dt) ? "00.00.00" :
             lpad(hour(dt), 2, '0') * "." * lpad(minute(dt), 2, '0') * "." * lpad(second(dt), 2, '0')
    isnothing(dt) && _info("Recording date and time unknown: 01.01.85 00.00.00 is written.")

    ann_label = bdf ? "BDF Annotations" : "EDF Annotations"
    h = IOBuffer()
    if bdf
        write(h, UInt8(0xff))
        write(h, _edf_field("BIOSEMI", 7))
    else
        write(h, _edf_field("0", 8))
    end
    write(h, _edf_field(patient, 80))
    write(h, _edf_field(recording, 80))
    write(h, date_s, time_s)
    write(h, _edf_field(string(256 * (ns + 1)), 8))
    write(h, _edf_field(plus ? fmt * "+C" : (bdf ? "24BIT" : ""), 44))
    write(h, _edf_field(string(n_rec), 8))
    write(h, _edf_field(_edf_num(record_duration), 8))
    write(h, _edf_field(string(ns), 4))
    sigf(f) = (for i = 1:ch_n; write(h, f(i)); end; plus && write(h, f(0)))
    sigf(i -> _edf_field(i == 0 ? ann_label : ll[i], 16))
    sigf(i -> _edf_field(i == 0 ? "" : transducers[i], 80))
    sigf(i -> _edf_field(i == 0 ? "" : replace(units[i], "μ" => "u", "µ" => "u", "²" => "2", "°" => "deg"), 8))
    sigf(i -> _edf_field(i == 0 ? "-1" : pmin[i], 8))
    sigf(i -> _edf_field(i == 0 ? "1" : pmax[i], 8))
    sigf(i -> _edf_field(string(dmin), 8))
    sigf(i -> _edf_field(string(dmax), 8))
    sigf(i -> _edf_field(i == 0 ? "" : prefilt(i), 80))
    sigf(i -> _edf_field(string(i == 0 ? ann_spr : spr), 8))
    sigf(i -> _edf_field("", 32))
    hdr = take!(h)
    length(hdr) == 256 * (ns + 1) || error("Internal error: EDF header size.")

    # ------------------------------------------------------------------ #
    # data records                                                       #
    # ------------------------------------------------------------------ #
    # digital values: (x − pmin) / (pmax − pmin) × (dmax − dmin) + dmin
    g = [(dmax - dmin) / (phi[i] - plo[i]) for i = 1:ch_n]
    dig = Matrix{Int32}(undef, ch_n, n_rec * spr)
    @inbounds for i = 1:ch_n, j = 1:(n_rec * spr)
        v = j <= n ? x[i, j] : 0.0
        dig[i, j] = clamp(round(Int32, (v - plo[i]) * g[i] + dmin), Int32(dmin), Int32(dmax))
    end
    open(file_name, "w") do io
        write(io, hdr)
        buf = Vector{UInt8}(undef, bps * spr)
        for r = 1:n_rec
            cols = ((r - 1) * spr + 1):(r * spr)
            for i = 1:ch_n
                for (k, j) in enumerate(cols)
                    v = dig[i, j]
                    o = bps * (k - 1)
                    buf[o + 1] = UInt8(v & 0xff)
                    buf[o + 2] = UInt8((v >> 8) & 0xff)
                    bdf && (buf[o + 3] = UInt8((v >> 16) & 0xff))
                end
                write(io, buf)
            end
            if plus
                a = zeros(UInt8, bps * ann_spr)
                a[1:length(ann_bytes[r])] .= ann_bytes[r]
                write(io, a)
            end
        end
    end

    return nothing
end
