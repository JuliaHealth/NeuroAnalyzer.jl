export export_fiff

# NeuroAnalyzer channel type → FIFF channel kind (FIFFV_*_CH); other types → misc (502)
const _FIFF_EXPORT_KIND = Dict{String, Int64}(
    "eeg" => 2, "ref" => 2, "csd" => 2,
    "meg" => 1, "mag" => 1, "grad" => 1, "meg_ref" => 301,
    "stim" => 3, "mrk" => 3,
    "eog" => 202, "emg" => 302, "ecg" => 402, "resp" => 602, "bio" => 102,
    "seeg" => 802, "ieeg" => 802, "dbs" => 803, "ecog" => 902,
)

# NeuroAnalyzer unit → (FIFF unit code, factor to SI)
const _FIFF_EXPORT_UNIT = Dict{String, Tuple{Int64, Float64}}(
    "fT" => (112, 1e-15), "pT" => (112, 1e-12), "nT" => (112, 1e-9), "μT" => (112, 1e-6),
    "T" => (112, 1.0),
    "fT/cm" => (201, 1e-13), "fT/m" => (201, 1e-15), "T/m" => (201, 1.0),
    "nV" => (107, 1e-9), "μV" => (107, 1e-6), "mV" => (107, 1e-3), "V" => (107, 1.0),
    "μV/m²" => (117, 1e-6), "V/m²" => (117, 1.0),
    "" => (-1, 1.0),
)

# electrode-like FIFF kinds: EEG, sEEG, DBS, ECoG
const _FIFF_ELECTRODE_KINDS = (2, 802, 803, 902)

# NeuroAnalyzer MEG coil class → FIFF coil type, used only when no FIFF coil id is stored
const _FIFF_DEFAULT_COIL = Dict{String, Int64}(
    "mag" => 3022,    # vv_mag_t1
    "pgrad" => 3012,  # vv_planar_t1
    "agrad" => 2001,  # axial_grad_5cm
)

# measurement date (UTC, Unix seconds) from NeuroAnalyzer date and time strings; nothing if not parseable
function _fiff_meas_date(d::AbstractString, t::AbstractString)::Union{Nothing, Int64}
    d = strip(d)
    t = strip(t)
    isempty(d) && return nothing
    m = match(r"^(\d{1,4})[-./](\d{1,2})[-./](\d{1,4})$", d)
    isnothing(m) && return nothing
    a, b, c = parse.(Int64, m.captures)
    # yyyy-mm-dd or dd-mm-yyyy / dd.mm.yy (EDF: yy 85-99 → 19yy, otherwise 20yy)
    y, mo, dd = length(m.captures[1]) == 4 ? (a, b, c) : (c, b, a)
    y < 100 && (y += y >= 85 ? 1900 : 2000)
    hh, mm, ss = 0, 0, 0
    mt = match(r"^(\d{1,2})[:.](\d{1,2})[:.](\d{1,2})", t)
    !isnothing(mt) && ((hh, mm, ss) = parse.(Int64, mt.captures))
    dt = try
        DateTime(y, mo, dd, hh, mm, ss)
    catch
        return nothing
    end
    return round(Int64, datetime2unix(dt))
end

# low-pass and high-pass frequencies from a prefiltering string "LP: x Hz; HP: y Hz"
function _fiff_prefiltering(s::AbstractString)::Tuple{Union{Nothing, Float64}, Union{Nothing, Float64}}
    num(r) = (m = match(r, s); isnothing(m) ? nothing : tryparse(Float64, m.captures[1]))
    lp = num(r"LP:?\s*([0-9.]+)"i)
    hp = num(r"HP:?\s*([0-9.]+)"i)
    return (isnothing(lp) || lp <= 0 ? nothing : lp), (isnothing(hp) ? nothing : hp)
end

# factor converting NeuroAnalyzer locs coordinates to meters
function _fiff_locs_scale(xyz::AbstractMatrix{<:Real}, unit::Symbol)::Float64
    unit === :m && return 1.0
    unit === :cm && return 1e-2
    unit === :mm && return 1e-3
    r = median([sqrt(sum(abs2, row)) for row in eachrow(xyz)])
    if r < 0.5
        return 1.0
    elseif r <= 1.5
        _info("Channel locations look normalized (unit sphere); scaled to a 95 mm head radius")
        return 0.095 / r
    elseif r < 30
        _info("Channel locations look like centimeters; converted to meters")
        return 1e-2
    else
        _info("Channel locations look like millimeters; converted to meters")
        return 1e-3
    end
end

"""
    export_fiff(obj; <keyword arguments>)

Export an EEG, MEG, ECoG or sEEG object to a FIFF raw file (Elekta/MEGIN Neuromag Functional Image File Format), readable by MNE-Python (`mne.io.read_raw_fif()`), Brainstorm, FieldTrip and `import_fiff()`.

Written:

- signals as 32-bit floats in SI units (V, T, T/m), sample-major data buffers of 1 s
- channel information: FIFF channel kind (from channel type), coil type, unit and the 12-element `loc` vector (sensor position and orientation); channel names longer than 15 characters are stored as MNE-Python does (short unique names in `ch_info` + full names in per-channel blocks)
- MEG geometry stored by `import_fiff()` (`:coil_loc`, `:coil_id`, `:dev_head_t`): MEG sensors in the device frame, EEG electrodes in the head frame; without it, positions come from `OBJ.locs` (head frame, meters) and MEG sensors get default coil types and orientations
- EEG/sEEG/ECoG electrode positions also as digitization points; cHPI coil frequencies and digitized positions
- SSP projectors (inactive), bad channels, sampling rate, line frequency, low-pass/high-pass (from the prefiltering string), measurement date, subject (ID, names, handedness, weight, height), recording description and experiment name
- markers as MNE annotations (onset relative to the first sample, duration, description = marker value, or the type if the value is empty; `BAD` markers are written as `BAD_<value>` and `boundary` markers as `BAD boundary`, which MNE-Python treats as bad segments); marker channels as the channels of the annotations (`ch_names` in MNE-Python); marker types, values and extra marker columns as the annotation extras (`extras` in MNE-Python), read back by `import_fiff()`
- MEG `first_time` as `first_samp`; IAS (MaxShield) objects as MaxShield raw data

Recordings larger than `split_size` are split as MNE-Python does: `name_raw.fif`, `name_raw-1.fif`, `name_raw-2.fif`, … (each part links to the next one; `import_fiff()` and MNE-Python read them as one recording).

The object must be continuous (one epoch).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `file_name::String`: file name (`.fif` or `.fiff`; MNE-Python expects raw files to end with `raw.fif`, `_meg.fif` or `_eeg.fif`)
- `locs_unit::Symbol=:auto`: unit of `OBJ.locs` coordinates, used when no FIFF geometry is stored:
    - `:auto`: guessed from the median distance from the origin (< 0.5 → m, 0.5-1.5 → unit sphere scaled to a 95 mm head radius, < 30 → cm, otherwise mm)
    - `:m`, `:cm`, `:mm`
- `split_size::Real=2.0`: maximum size of one file in GB (at most 2 GB)
- `overwrite::Bool=false`: if `true`, overwrite existing files

# Returns

- `Vector{String}`: names of the written files (one per split part)

# References

1. Elekta Neuromag: Functional Image File Format Description. FIFF v1.3, 2011.
2. MNE-Python FIFF writer (mne/_fiff/write.py, mne/io/base.py), used as the reference implementation.
"""
function export_fiff(
    obj::NeuroAnalyzer.NEURO;
    file_name::String,
    locs_unit::Symbol = :auto,
    split_size::Real = 2.0,
    overwrite::Bool = false,
)::Vector{String}
    # ------------------------------------------------------------------ #
    # validate                                                           #
    # ------------------------------------------------------------------ #
    _check_datatype(obj, ["eeg", "meg", "ecog", "seeg"])
    nepochs(obj) == 1 || throw(ArgumentError("OBJ must be continuous (one epoch)."))
    lowercase(splitext(file_name)[2]) in (".fif", ".fiff") ||
        throw(ArgumentError("File name extension must be .fif or .fiff"))
    locs_unit in (:auto, :m, :cm, :mm) ||
        throw(ArgumentError("locs_unit must be :auto, :m, :cm or :mm."))
    0 < split_size <= 2.0 || throw(ArgumentError("split_size must be in (0, 2] GB."))

    rec = obj.header.recording
    ch_n = nchannels(obj)
    n_samples = size(obj.data, 2)
    sr = rec[:sampling_rate]
    clabels = String.(labels(obj))
    ch_type = String.(rec[:channel_type])
    units = String.(rec[:unit])

    # ------------------------------------------------------------------ #
    # channel kinds, units, coil types                                   #
    # ------------------------------------------------------------------ #
    kinds = [get(_FIFF_EXPORT_KIND, t, 502) for t in ch_type]
    others = unique(ch_type[kinds .== 502])
    isempty(others) || _info("Channels of type $(join(others, ", ")) are exported as FIFF misc channels")

    unit_code = Vector{Int64}(undef, ch_n)
    scale = Vector{Float64}(undef, ch_n)
    for idx = 1:ch_n
        u = get(_FIFF_EXPORT_UNIT, units[idx], nothing)
        if isnothing(u)
            _warn("Channel $(clabels[idx]): unit $(units[idx]) has no FIFF equivalent; exported without unit and scaling.")
            u = (-1, 1.0)
        end
        unit_code[idx], scale[idx] = u
    end

    has_geom = haskey(rec, :coil_loc) && size(rec[:coil_loc]) == (ch_n, 12)
    coil_id_ok = haskey(rec, :coil_id) && length(rec[:coil_id]) == ch_n
    coil_cls = haskey(rec, :coil_type) && length(rec[:coil_type]) == ch_n ? String.(rec[:coil_type]) : fill("", ch_n)

    coil = zeros(Int64, ch_n)
    meg_default = false
    for idx = 1:ch_n
        k = kinds[idx]
        if k in (1, 301)
            if coil_id_ok && rec[:coil_id][idx] != 0
                coil[idx] = rec[:coil_id][idx]
            else
                cls = ch_type[idx] == "mag" ? "mag" : coil_cls[idx] in keys(_FIFF_DEFAULT_COIL) ? coil_cls[idx] : "pgrad"
                coil[idx] = _FIFF_DEFAULT_COIL[cls]
                meg_default = true
            end
        elseif coil_id_ok
            # FIFF coil type stored by import_fiff()
            coil[idx] = rec[:coil_id][idx]
        elseif k in _FIFF_ELECTRODE_KINDS
            coil[idx] = ch_type[idx] == "csd" ? 6 : 1
        end
    end
    meg_default && _warn("No FIFF coil types stored for some MEG channels; Vectorview coil types are assumed.")

    # ------------------------------------------------------------------ #
    # sensor geometry (loc: r0, ex, ey, ez)                              #
    # ------------------------------------------------------------------ #
    loc = zeros(Float64, ch_n, 12)
    if has_geom
        loc .= rec[:coil_loc]
    elseif DataFrames.nrow(obj.locs) > 0
        l = obj.locs
        xyz = Matrix{Float64}(l[:, [:loc_x, :loc_y, :loc_z]])
        nz = [any(!iszero, r) for r in eachrow(xyz)]
        if any(nz)
            f = _fiff_locs_scale(xyz[nz, :], locs_unit)
            for idx = 1:ch_n
                j = findfirst(isequal(clabels[idx]), l[!, :label])
                (isnothing(j) || !nz[j]) && continue
                loc[idx, 1:3] .= xyz[j, :] .* f
                if kinds[idx] in (1, 301)
                    # unknown orientation: radial coil normal, tangential ex/ey
                    ez = normalize(loc[idx, 1:3])
                    ex = normalize(cross(abs(ez[3]) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0], ez))
                    loc[idx, 4:12] .= vcat(ex, cross(ez, ex), ez)
                end
            end
            any(in((1, 301)), kinds) && _warn(
                "No FIFF sensor geometry stored; MEG sensor positions taken from OBJ.locs (head frame) with radial orientations - not suitable for source modeling.",
            )
        end
    end

    # ------------------------------------------------------------------ #
    # names                                                              #
    # ------------------------------------------------------------------ #
    short = _fiff_short_names(clabels)
    long_names = short != clabels
    name_map = Dict(zip(clabels, short))

    # ------------------------------------------------------------------ #
    # metadata                                                           #
    # ------------------------------------------------------------------ #
    meas_date = _fiff_meas_date(string(get(rec, :recording_date, "")), string(get(rec, :recording_time, "")))
    lp, hp = isempty(rec[:prefiltering]) ? (nothing, nothing) : _fiff_prefiltering(string(rec[:prefiltering][1]))
    !isnothing(lp) && lp > sr / 2 && (lp = nothing)
    line_freq = get(rec, :line_frequency, 0)
    first_samp = round(Int64, get(rec, :first_time, 0.0) * sr)
    maxshield = get(rec, :maxshield, false) === true
    dev_head_t = get(rec, :dev_head_t, Matrix{Float64}(undef, 0, 0))
    bads = clabels[get(rec, :bad_channel, falses(ch_n))]
    hpi = get(rec, :hpi, Dict{Symbol, Any}())

    subj = obj.header.subject
    hand = Dict("right" => 1, "left" => 2, "ambidextrous" => 3)

    # SSP projectors: rows of ssp_data, columns = channels marked in ssp_channels
    ssp = Tuple{String, Vector{String}, Vector{Float64}}[]
    if haskey(rec, :ssp_data) && !isempty(rec[:ssp_data]) && length(get(rec, :ssp_channels, Bool[])) == ch_n
        sch = clabels[rec[:ssp_channels]]
        if size(rec[:ssp_data], 2) == length(sch)
            for row in axes(rec[:ssp_data], 1)
                push!(ssp, (string(rec[:ssp_labels][row]), sch, rec[:ssp_data][row, :]))
            end
        else
            _warn("SSP projectors do not match the SSP channels; not exported.")
        end
    end

    # markers → annotations (onsets relative to the first sample, as MNE-Python reads them
    # when the annotations block stores no meas_date)
    mrk = obj.markers
    # description = value (or type, if the value is empty); bad segments are written as
    # "BAD_…" and discontinuities as "BAD boundary", which MNE-Python treats as bad segments
    an_desc = [_fiff_an_desc(t, v) for (t, v) in zip(string.(mrk[!, :type]), string.(mrk[!, :value]))]
    # marker channels (MNE: JSON list of channel name lists, empty: all channels); written only
    # if a marker refers to some channels, as MNE-Python does
    an_ch = ""
    if any(!isempty, mrk[!, :channel])
        lbl_set = Set(clabels)
        mch = [String[c for c in ch if c in lbl_set] for ch in mrk[!, :channel]]
        any(length.(mch) .!= length.(mrk[!, :channel])) &&
            _warn("Marker channels that are not in the object are not exported.")
        an_ch = JSON.json(mch)
    end
    # marker type and value as the annotation extras of MNE-Python (FIFF free list), so that
    # import_fiff() restores the markers exactly (the description may have a BAD_ prefix)
    # (with the extra marker columns; names reserved by MNE-Python are skipped)
    ex_cols = Base.filter(
        c -> !(c in ("onset", "duration", "description", "ch_names", "hed_string")),
        _marker_extra_cols(mrk),
    )
    length(ex_cols) < length(_marker_extra_cols(mrk)) &&
        _warn("Extra marker columns with names reserved by MNE-Python are not exported.")
    ex_val(x) = ismissing(x) || isnothing(x) ? nothing : (x isa Union{AbstractString, Real} ? x : string(x))
    an_extra = JSON.json([
        merge(
            Dict{String, Any}("type" => string(mrk[i, :type]), "value" => string(mrk[i, :value])),
            Dict{String, Any}(c => ex_val(mrk[i, c]) for c in ex_cols),
        ) for i = 1:DataFrames.nrow(mrk)
    ])

    # ------------------------------------------------------------------ #
    # output files                                                       #
    # ------------------------------------------------------------------ #
    buf_len = sr
    bytes_per_buf = 16 + 4 * ch_n * buf_len
    limit = floor(Int64, split_size * 1024^3) - 4096
    est_info = 2048 + ch_n * (112 + (long_names ? 300 : 0)) + 64 * DataFrames.nrow(mrk) + ncodeunits(an_ch) + ncodeunits(an_extra) +
               sum(x -> 4 * length(x[3]) + 32 * length(x[2]), ssp; init = 0) + 1024
    max_bufs = max(1, (limit - est_info) ÷ bytes_per_buf)
    n_bufs = cld(n_samples, buf_len)
    n_parts = cld(n_bufs, max_bufs)
    base, ext = splitext(file_name)
    part_names = [i == 1 ? file_name : base * "-$(i - 1)" * ext for i = 1:n_parts]
    for f in part_names
        isfile(f) && !overwrite &&
            throw(ArgumentError("File $f cannot be saved, to overwrite use overwrite=true."))
    end
    n_parts > 1 && _info("Data are split into $n_parts files")

    machid = (rand(Int32), rand(Int32))
    now_s = round(Int64, datetime2unix(Dates.now(Dates.UTC)))
    t_id = (Int32(now_s), Int32(0))

    tag(s) = _find_fiff_tag(s)
    blk(s) = _find_fiff_block(s)

    for part = 1:n_parts
        b1 = (part - 1) * max_bufs + 1
        b2 = min(part * max_bufs, n_bufs)
        s1 = (b1 - 1) * buf_len + 1
        s2 = min(b2 * buf_len, n_samples)

        open(part_names[part], "w") do io
            _fiff_start_file(io, machid, t_id)
            _fiff_start_block(io, blk("meas"))
            _fiff_write_id(io, tag("block_id"), machid, t_id)

            # -------------------------- meas_info --------------------------
            _fiff_start_block(io, blk("meas_info"))

            # cHPI coils
            if !isempty(hpi) && haskey(hpi, :freqs)
                _fiff_start_block(io, blk("hpi_meas"))
                _fiff_write_int(io, tag("hpi_ncoil"), length(hpi[:freqs]))
                for (c, fr) in enumerate(hpi[:freqs])
                    _fiff_start_block(io, blk("hpi_coil"))
                    _fiff_write_int(io, tag("hpi_coil_no"), hpi[:coil_no][c])
                    _fiff_write_float(io, tag("hpi_coil_freq"), fr)
                    _fiff_end_block(io, blk("hpi_coil"))
                end
                _fiff_end_block(io, blk("hpi_meas"))
            end

            # digitization (head frame): HPI coils, EEG electrodes
            dig_hpi = !isempty(hpi) && haskey(hpi, :dig) ? hpi[:dig] : zeros(0, 3)
            el = [idx for idx = 1:ch_n if kinds[idx] in _FIFF_ELECTRODE_KINDS && any(!iszero, loc[idx, 1:3])]
            if size(dig_hpi, 1) > 0 || !isempty(el)
                _fiff_start_block(io, blk("isotrak"))
                for c in axes(dig_hpi, 1)
                    _fiff_write_dig_point(io, 2, hpi[:coil_no][c], dig_hpi[c, :])
                end
                for (n, idx) in enumerate(el)
                    _fiff_write_dig_point(io, 3, n, loc[idx, 1:3])
                end
                _fiff_end_block(io, blk("isotrak"))
            end

            # device → head transform (FIFFV_COORD_DEVICE = 1, FIFFV_COORD_HEAD = 4)
            size(dev_head_t) == (4, 4) && _fiff_write_coord_trans(io, 1, 4, dev_head_t)

            r = string(get(rec, :recording, ""))
            isempty(r) || _fiff_write_string(io, tag("description"), r)
            pn = string(get(obj.header.experiment, :name, ""))
            isempty(pn) || _fiff_write_string(io, tag("proj_name"), pn)
            isnothing(meas_date) || _fiff_write_int(io, tag("meas_date"), [meas_date, 0])
            _fiff_write_int(io, tag("nchan"), ch_n)
            _fiff_write_float(io, tag("sfreq"), sr)
            isnothing(lp) || _fiff_write_float(io, tag("lowpass"), lp)
            isnothing(hp) || _fiff_write_float(io, tag("highpass"), hp)
            line_freq isa Real && line_freq > 0 && _fiff_write_float(io, tag("line_freq"), line_freq)
            _fiff_write_int(io, tag("data_pack"), _FIFFT_FLOAT)
            # re-referenced EEG: stops MNE-Python from adding an average reference projector
            isempty(string(get(rec, :reference, ""))) || _fiff_write_int(io, tag("mne_custom_ref"), 1)

            # subject
            sid = string(get(subj, :id, ""))
            fields = [
                ("subj_his_id", string(get(subj, :id, ""))),
                ("subj_last_name", string(get(subj, :last_name, ""))),
                ("subj_first_name", string(get(subj, :first_name, ""))),
                ("subj_middle_name", string(get(subj, :middle_name, ""))),
            ]
            h = get(hand, lowercase(string(get(subj, :handedness, ""))), nothing)
            w = get(subj, :weight, -1)
            ht = get(subj, :height, -1)
            if any(!isempty(last(x)) for x in fields) || !isnothing(h) || (w isa Real && w > 0) || (ht isa Real && ht > 0)
                _fiff_start_block(io, blk("subject"))
                id_int = tryparse(Int32, sid)
                isnothing(id_int) || _fiff_write_int(io, tag("subj_id"), id_int)
                for (t, v) in fields
                    isempty(v) || _fiff_write_string(io, tag(t), v)
                end
                isnothing(h) || _fiff_write_int(io, tag("subj_hand"), h)
                w isa Real && w > 0 && _fiff_write_float(io, tag("subj_weight"), w)
                ht isa Real && ht > 0 && _fiff_write_float(io, tag("subj_height"), ht)
                _fiff_end_block(io, blk("subject"))
            end

            # SSP projectors (FIFFB_PROJ = 313, FIFFB_PROJ_ITEM = 314)
            if !isempty(ssp)
                _fiff_start_block(io, blk("xfit_proj"))
                for (desc, names, v) in ssp
                    _fiff_start_block(io, blk("xfit_proj_item"))
                    _fiff_write_name_list(io, tag("proj_item_ch_name_list"), [name_map[n] for n in names])
                    _fiff_write_string(io, tag("name"), desc)
                    _fiff_write_int(io, tag("nchan"), length(names))
                    # FIFFV_PROJ_ITEM_FIELD
                    _fiff_write_int(io, tag("proj_item_kind"), 1)
                    _fiff_write_float(io, tag("proj_item_time"), 0.0)
                    _fiff_write_int(io, tag("proj_item_nvec"), 1)
                    _fiff_write_int(io, tag("mne_proj_item_active"), 0)
                    _fiff_write_float_matrix(io, tag("proj_item_vectors"), reshape(v, 1, :))
                    _fiff_end_block(io, blk("xfit_proj_item"))
                end
                _fiff_end_block(io, blk("xfit_proj"))
            end

            # channels
            for idx = 1:ch_n
                _fiff_write_ch_info(
                    io;
                    scan_no = idx, log_no = idx, kind = kinds[idx], range = 1.0, cal = 1.0,
                    coil_type = coil[idx], loc = loc[idx, :], unit = unit_code[idx],
                    unit_mul = 0, ch_name = short[idx],
                )
            end
            if long_names
                for idx = 1:ch_n
                    _fiff_start_block(io, blk("ch_info"))
                    _fiff_write_int(io, tag("ch_scan_no"), idx)
                    _fiff_write_int(io, tag("ch_logical_no"), idx)
                    _fiff_write_int(io, tag("ch_kind"), kinds[idx])
                    _fiff_write_float(io, tag("ch_range"), 1.0)
                    _fiff_write_float(io, tag("ch_cal"), 1.0)
                    _fiff_write_int(io, tag("ch_coil_type"), coil[idx])
                    _fiff_write_float(io, tag("ch_pos"), loc[idx, :])
                    _fiff_write_int(io, tag("ch_unit"), unit_code[idx])
                    _fiff_write_int(io, tag("ch_unit_mul"), 0)
                    _fiff_write_string(io, tag("ch_dacq_name"), clabels[idx])
                    # FIFFV_COORD_DEVICE (MEG) / FIFFV_COORD_HEAD
                    _fiff_write_int(io, tag("ch_coord_frame"), kinds[idx] in (1, 301) ? 1 : 4)
                    _fiff_end_block(io, blk("ch_info"))
                end
            end

            # bad channels
            if !isempty(bads)
                _fiff_start_block(io, blk("bad_channels"))
                _fiff_write_name_list(io, tag("mne_ch_name_list"), [name_map[b] for b in bads])
                _fiff_end_block(io, blk("bad_channels"))
            end

            _fiff_end_block(io, blk("meas_info"))

            # -------------------------- annotations -------------------------
            if DataFrames.nrow(mrk) > 0
                _fiff_start_block(io, blk("mne_annotations"))
                _fiff_write_float(io, tag("mne_baseline_min"), mrk[!, :start])
                _fiff_write_float(io, tag("mne_baseline_max"), mrk[!, :start] .+ mrk[!, :length])
                _fiff_write_name_list(io, tag("description"), an_desc)
                isempty(an_ch) || _fiff_write_string(io, tag("mne_epochs_drop_log"), an_ch)
                _fiff_write_string(io, tag("free_list"), an_extra)
                _fiff_end_block(io, blk("mne_annotations"))
            end

            # -------------------------- raw data ----------------------------
            raw_blk = maxshield ? blk("maxshield_raw_data") : blk("raw_data")
            _fiff_start_block(io, raw_blk)
            fs = first_samp + s1 - 1
            fs != 0 && _fiff_write_int(io, tag("first_sample"), fs)
            if part > 1
                _fiff_start_block(io, blk("ref"))
                # FIFFV_ROLE_PREV_FILE
                _fiff_write_int(io, tag("ref_role"), 1)
                _fiff_write_string(io, tag("ref_file_name"), basename(part_names[part - 1]))
                _fiff_write_int(io, tag("ref_file_num"), part - 2)
                _fiff_end_block(io, blk("ref"))
            end
            k_buf = tag("data_buffer")
            for s = s1:buf_len:s2
                e = min(s + buf_len - 1, s2)
                # sample-major: all channels of sample 1, then sample 2, ...
                b = Float32.(view(obj.data, :, s:e, 1) .* scale)
                _fiff_write_tag(io, k_buf, _FIFFT_FLOAT, _fiff_be_bytes(vec(b)))
            end
            if part < n_parts
                _fiff_start_block(io, blk("ref"))
                # FIFFV_ROLE_NEXT_FILE
                _fiff_write_int(io, tag("ref_role"), 2)
                _fiff_write_string(io, tag("ref_file_name"), basename(part_names[part + 1]))
                _fiff_write_int(io, tag("ref_file_num"), part)
                _fiff_end_block(io, blk("ref"))
            end
            _fiff_end_block(io, raw_blk)

            _fiff_end_block(io, blk("meas"))
            _fiff_end_file(io)
        end
    end

    return part_names
end
