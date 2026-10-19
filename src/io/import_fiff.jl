export load_fiff
export import_fiff
export import_locs_fiff

# ---------------------------------------------------------------------------
# internal helpers
# ---------------------------------------------------------------------------

# all decoded values of a tag, optionally restricted to a block type and/or block instance
# fiff_object entries: (index, tag, block, data, instance)
function _fiff_getall(
    fo::Vector{Any},
    tag::String;
    block::Union{Nothing, String} = nothing,
    inst::Union{Nothing, Int64} = nothing,
)::Vector{Any}
    return Any[
        x[4] for x in fo if x[2] == tag &&
        (isnothing(block) || x[3] == block) &&
        (isnothing(inst) || x[5] == inst)
    ]
end

function _fiff_get(
    fo::Vector{Any},
    tag::String;
    block::Union{Nothing, String} = nothing,
    inst::Union{Nothing, Int64} = nothing,
)
    v = _fiff_getall(fo, tag; block = block, inst = inst)
    return isempty(v) ? nothing : v[1]
end

# instance numbers of all occurrences of a block type, in file order
_fiff_instances(fo::Vector{Any}, block::String)::Vector{Int64} =
    unique(Int64[x[5] for x in fo if x[3] == block])

# colon-separated FIFF name list; MNE escapes literal colons as {COLON}
_fiff_name_list(s)::Vector{String} =
    (isnothing(s) || isempty(s)) ? String[] :
    [replace(String(v), "{COLON}" => ":") for v in split(s, ':')]

_fiff_vec(x) = isnothing(x) ? Any[] : (x isa AbstractVector ? x : [x])

_fiff_nan2nothing(x) = (x isa Real && isnan(x)) ? nothing : x

# channel labels as used in NEURO objects: ch_name (or kind + logical number when empty)
# normalized with _clean_meg_labels()
function _fiff_labels(mi::Dict{Symbol, Any})::Vector{String}
    labels = [
        isempty(n) ? uppercase(_find_fiff_chtype(ch.kind)) * " " * string(ch.log_no) : n for
        (n, ch) in zip(mi[:ch_names], mi[:ch_info])
    ]
    return _clean_meg_labels(labels)
end

# least-squares sphere fit to points (n × 3); returns center and radius
function _fiff_fit_sphere(p::AbstractMatrix{<:Real})::Tuple{Vector{Float64}, Float64}
    size(p, 1) >= 4 || throw(ArgumentError("At least 4 points are required to fit a sphere."))
    a = hcat(2 .* p, ones(size(p, 1)))
    b = vec(sum(p .^ 2; dims = 2))
    sol = a \ b
    c = Vector{Float64}(sol[1:3])
    return c, sqrt(sol[4] + sum(abs2, c))
end

# head-shape digitization points (head frame, m):
# EEG (kind 3) and extra (kind 4) points, nose/face region (z < 0, y > 0) excluded as in MNE-Python
_fiff_hsp(mi::Dict{Symbol, Any})::Vector{Vector{Float64}} = [
    [d[3], d[4], d[5]] for
    d in mi[:dig] if d isa Tuple && d[1] in (3, 4) && !(d[5] < 0 && d[4] > 0)
]

# head origin (center of a sphere fitted to head-shape points, head frame, m);
# empty if fewer than 10 points are available
function _fiff_head_origin(mi::Dict{Symbol, Any})::Vector{Float64}
    hsp = _fiff_hsp(mi)
    length(hsp) >= 10 || return Float64[]
    return first(_fiff_fit_sphere(permutedims(reduce(hcat, hsp))))
end

# cHPI coil frequencies (Hz) and digitized coil positions (head frame, m): coils sorted by
# number; digitized points sorted by identifier and matched by the HPI digitization order
# (as in MNE-Python), otherwise by coil number
function _fiff_hpi(mi::Dict{Symbol, Any})::Dict{Symbol, Any}
    # (coil number, frequency)
    coils = mi[:hpi_coils]
    dig = [(d[2], [d[3], d[4], d[5]]) for d in mi[:dig] if d isa Tuple && d[1] == 2]
    (isempty(coils) || isempty(dig)) && return Dict{Symbol, Any}()
    sort!(coils; by = first)
    sort!(dig; by = first)
    order = get(mi, :hpi_dig_order, Int64[])
    pos = Vector{Vector{Float64}}()
    if length(order) == length(coils) && all(in(1:length(dig)), order)
        pos = [dig[o][2] for o in order]
    else
        for (no, _) in coils
            i = findfirst(d -> d[1] == no, dig)
            isnothing(i) && break
            push!(pos, dig[i][2])
        end
    end
    if length(pos) != length(coils)
        # coil numbers do not match digitization identifiers: fall back to order
        length(dig) == length(coils) || begin
            _warn("Number of digitized HPI coils ($(length(dig))) differs from cHPI coils ($(length(coils))); cHPI information not stored.")
            return Dict{Symbol, Any}()
        end
        _warn("HPI coil numbers do not match digitization identifiers; coils matched by order.")
        pos = [d[2] for d in dig]
    end
    return Dict{Symbol, Any}(
        :coil_no => [first(c) for c in coils],
        :freqs => Float64[last(c) for c in coils],
        :dig => permutedims(reduce(hcat, pos)),
    )
end

# device → head rigid transform (rotation, translation) or nothing
function _fiff_dev_head_t(mi::Dict{Symbol, Any})
    candidates = Any[mi[:dev_head_t]]
    hr = mi[:hpi][:hpi_result]
    for (k, v) in hr
        startswith(string(k), "coord_trans") && push!(candidates, v)
    end
    for t in candidates
        t isa Tuple && length(t) == 6 || continue
        # FIFFV_COORD_DEVICE = 1, FIFFV_COORD_HEAD = 4
        t[1] == 1 && t[2] == 4 && return (t[3], t[4])
        t[1] == 4 && t[2] == 1 && return (t[5], t[6])
    end
    return nothing
end

"""
    _fiff_locs(mi, labels)

Build a NeuroAnalyzer `locs` DataFrame from FIFF channel information.

- MEG sensor positions (`ch_info.loc[1:3]`, device frame) are transformed to the head frame using the device → head transform; without it (e.g. empty-room data) they stay in the device frame.
- EEG/sEEG/ECoG/DBS electrode positions are stored in the head frame.
- Reference MEG sensors and channels without a position (stim, misc, ...) are skipped.
- With at least 10 head-shape digitisation points (extra and EEG points, nose region excluded, as in MNE-Python), coordinates are centred on the fitted head sphere and divided by its radius: the head surface (EEG electrodes) is at radius ≈ 1 and MEG sensors are at their true relative distance. Otherwise coordinates are centred on a sphere fitted to the channel positions and scaled so that the outermost MEG sensor has radius 1.2 and the outermost electrode radius 1.0.

Both FIFF head and device frames are RAS (x → right, y → anterior, z → superior).
"""
function _fiff_locs(mi::Dict{Symbol, Any}, labels::Vector{String})::DataFrame
    chs = mi[:ch_info]
    t = _fiff_dev_head_t(mi)
    meg_kinds = ("meg",)
    el_kinds = ("eeg", "seeg", "ecog", "dbs")
    has_meg = any(ch -> _find_fiff_chtype(ch.kind) in meg_kinds, chs)
    has_meg && isnothing(t) && _info(
        "No device → head transform found; MEG sensor locations are given in the device frame.",
    )

    idx = Int64[]
    pos = Vector{Vector{Float64}}()
    is_meg = Bool[]
    for (i, ch) in enumerate(chs)
        r0 = ch.loc[1:3]
        (all(isfinite, r0) && !all(iszero, r0)) || continue
        k = _find_fiff_chtype(ch.kind)
        if k in meg_kinds
            !isnothing(t) && (r0 = t[1] * r0 .+ t[2])
        elseif !(k in el_kinds)
            continue
        end
        push!(idx, i)
        push!(pos, r0)
        push!(is_meg, k in meg_kinds)
    end

    l = DataFrame(
        :label => String[],
        :loc_radius => Float64[],
        :loc_theta => Float64[],
        :loc_x => Float64[],
        :loc_y => Float64[],
        :loc_z => Float64[],
        :loc_radius_sph => Float64[],
        :loc_theta_sph => Float64[],
        :loc_phi_sph => Float64[],
    )
    if isempty(idx)
        _info("FIFF data contain no channel positions.")
        return l
    end
    # n × 3
    p = permutedims(reduce(hcat, pos))

    # head-shape points: kind 3 = EEG, 4 = extra; exclude nose/face (z < 0, y > 0)
    hsp = _fiff_hsp(mi)
    if length(hsp) >= 10 && (!has_meg || !isnothing(t))
        # head sphere: head surface at radius 1, sensors at their true relative distance
        c, r_head = _fiff_fit_sphere(permutedims(reduce(hcat, hsp)))
        p = (p .- c') ./ r_head
    else
        # no usable head shape: center on the channels
        # outermost MEG sensor at radius 1.2 and outermost electrode at radius 1.0
        c = size(p, 1) >= 4 ? first(_fiff_fit_sphere(p)) : zeros(3)
        p = p .- c'
        r = vec(sqrt.(sum(p .^ 2; dims = 2)))
        for (mask, r_out) in ((is_meg, 1.2), (.!is_meg, 1.0))
            any(mask) || continue
            rmax = maximum(r[mask])
            rmax > 0 && (p[mask, :] .*= r_out / rmax)
        end
        _info("Fewer than 10 head-shape points: channel locations scaled to fixed radii (MEG 1.2, EEG 1.0).")
    end

    for (j, i) in enumerate(idx)
        push!(l, (labels[i], 0.0, 0.0, p[j, 1], p[j, 2], p[j, 3], 0.0, 0.0, 0.0))
    end
    locs_cart2sph!(l)
    locs_cart2pol!(l)

    return l
end

"""
    _fiff_raw_data(buf, fb, insts, raw_inst, nchan)

Assemble the continuous data matrix (channels × samples) from the data buffers of one raw-data block. Handles `first_sample`, `data_skip` (in buffers) and `data_skip_samp` (in samples). Skips before the first buffer advance `first_samp` (as MNE-Python does); later skips are filled with zeros.

Data are returned uncalibrated (raw ADC units).
"""
function _fiff_raw_data(
    buf::Vector{Vector{UInt8}},
    fb::Matrix{Int64},
    insts::Vector{Int64},
    raw_inst::Int64,
    nchan::Int64,
)::Tuple{Matrix{Float64}, Int64, String}
    k_first = _find_fiff_tag("first_sample")
    k_skip = _find_fiff_tag("data_skip")
    k_skip_samp = _find_fiff_tag("data_skip_samp")
    k_buf = _find_fiff_tag("data_buffer")

    first_samp = 0
    pending_skip = 0
    seen_buffer = false
    orig_format = ""
    # (tag index or 0 = zeros, n samples)
    segs = Tuple{Int64, Int64}[]

    for i in findall(isequal(raw_inst), insts)
        kind = fb[i, 2]
        if kind == k_first && !seen_buffer
            first_samp = _i32i64(buf[i][1:4])
        elseif kind == k_skip
            pending_skip += _i32i64(buf[i][1:4])
        elseif kind == k_skip_samp
            n = _i32i64(buf[i][1:4])
            n > 0 && push!(segs, (0, n))
        elseif kind == k_buf
            (fb[i, 3] & 0xFFFF0000) != 0 &&
                throw(ArgumentError("Matrix-coded FIFF data buffers are not supported."))
            t = fb[i, 3] & 0x0000FFFF
            nb = t in (2, 16) ? 2 : t in (3, 4) ? 4 : t == 5 ? 8 : 0
            nb == 0 && throw(
                ArgumentError(
                    "FIFF data buffers of type $(_find_fiff_dt(t)) are not supported; please send this file to adam.wysokinski@neuroanalyzer.org",
                ),
            )
            isempty(orig_format) && (orig_format = _find_fiff_dt(t))
            sz = length(buf[i])
            sz % (nb * nchan) == 0 || throw(
                ArgumentError(
                    "FIFF data buffer size ($sz bytes) is not a multiple of $nchan channels × $nb bytes.",
                ),
            )
            ns = sz ÷ (nb * nchan)
            if pending_skip > 0
                if seen_buffer
                    push!(segs, (0, pending_skip * ns))
                else
                    first_samp += pending_skip * ns
                end
                pending_skip = 0
            end
            push!(segs, (i, ns))
            seen_buffer = true
        end
    end

    seen_buffer || throw(ArgumentError("FIFF raw data block contains no data buffers."))

    data = zeros(Float64, nchan, sum(last, segs))
    col = 0
    for (i, ns) in segs
        if i > 0
            t = fb[i, 3] & 0x0000FFFF
            T = t in (2, 16) ? Int16 : t == 3 ? Int32 : t == 4 ? Float32 : Float64
            # each buffer is sample-major: all channels of sample 1, then sample 2, ...
            @views data[:, (col + 1):(col + ns)] .= reshape(_be(T, buf[i]), nchan, ns)
        end
        col += ns
    end

    return data, first_samp, orig_format
end

# path of the next part of a split FIFF file (or nothing)
function _fiff_next_file(fo::Vector{Any}, file_name::String)::Union{Nothing, String}
    for inst in _fiff_instances(fo, "ref")
        role = _fiff_get(fo, "ref_role"; inst = inst)
        # FIFFV_ROLE_NEXT_FILE
        role == 2 || continue
        name = _fiff_get(fo, "ref_file_name"; inst = inst)
        num = _fiff_get(fo, "ref_file_num"; inst = inst)
        if !isnothing(name) && !isempty(name)
            return joinpath(dirname(file_name), basename(name))
        elseif !isnothing(num)
            # MNE default naming: xxx_raw.fif → xxx_raw-1.fif → xxx_raw-2.fif ...
            m = match(r"^(.*?)(-\d+)?(\.fif(?:f)?)$"i, file_name)
            isnothing(m) && return nothing
            return m.captures[1] * "-$num" * m.captures[3]
        end
    end
    return nothing
end

"""
    load_fiff(file_name)

Load an Elekta/MEGIN Neuromag FIFF (Functional Image File Format) file and return the parsed FIFF structure.

Tags are decoded according to their data type, so unknown tag kinds do not abort reading. Raw data are read from `raw_data`, `continuous_data` or MaxShield (`maxshield_raw_data`) blocks; skips and `first_sample` are handled. 

Gzip-compressed files (`.fif.gz`) are not supported - decompress them first.

# Arguments

- `file_name::String`: path to the `.fif` or `.fiff` file
- `read_data::Bool=true`: read raw data buffers; if `false`, only measurement info is read (`:raw_data` is `nothing`) and files without raw data (e.g. evoked) are accepted

# Returns

- `Dict{Symbol, Any}`: FIFF structure (`:meas_info`, `:raw_data`)
- `Vector{Any}`: FIFF object (one `(index, tag, block, data, block_instance)` tuple per tag; data buffers are not decoded here)
- `Matrix{Int64}`: FIFF blocks (tag position, kind, data type, size, block level, block type)

# References

1. Elekta Neuromag: Functional Image File Format Description. FIFF v1.3, 2011.
2. MNE-Python FIFF reader (mne/_fiff), used as the reference implementation.
"""
function load_fiff(
    file_name::String;
    read_data::Bool = true,
)::Tuple{Dict{Symbol, Any}, Vector{Any}, Matrix{Int64}}
    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    endswith(lowercase(file_name), ".gz") && throw(
        ArgumentError("Gzip-compressed FIFF files are not supported; decompress $file_name first."),
    )

    buf, fiff_blocks, insts = open(file_name, "r") do fid
        # the first tag must be file_id (id_t, 20 bytes)
        tag_kind, tag_type, tag_size, _, _ = try
            _read_fiff_tag(fid)
        catch
            throw(ArgumentError("$file_name: first FIFF tag cannot be read."))
        end
        tag_kind == _find_fiff_tag("file_id") ||
            throw(ArgumentError("$file_name is not a FIFF file."))
        (tag_type == 31 && tag_size == 20) ||
            throw(ArgumentError("$file_name is not a valid FIFF file (malformed file_id)."))
        return _create_fiff_block(fid)
    end

    # ------------------------------------------------------------------ #
    # decode all tags (by data type)                                     #
    # ------------------------------------------------------------------ #
    k_buf = _find_fiff_tag("data_buffer")
    fiff_object = Any[]
    for idx in eachindex(buf)
        kind = fiff_blocks[idx, 2]
        dtype = fiff_blocks[idx, 3]
        d = kind == k_buf ? nothing : _fiff_decode(dtype, buf[idx])
        # code → name conversions for enumerated tags
        if d isa Integer
            if kind == 210
                d = _find_fiff_aspect(d)
            elseif kind == 264
                d = _find_fiff_sss_job(d)
            elseif kind == 252
                d = _find_fiff_chtype(d)
            elseif kind == 256
                d = _find_fiff_unit(d)
            elseif kind == 3411
                d = _find_fiff_proj_item(d)
            elseif kind == 280
                d = _find_fiff_gantry_type(d)
            end
        end
        push!(
            fiff_object,
            (idx, _find_fiff_tag(kind), _find_fiff_block(fiff_blocks[idx, 6]), d, insts[idx]),
        )
    end

    # ------------------------------------------------------------------ #
    # measurement info                                                   #
    # ------------------------------------------------------------------ #
    # main meas_info block = first one that defines nchan (processing
    # history records may contain other blocks of the same type)
    mi_inst = nothing
    for inst in _fiff_instances(fiff_object, "meas_info")
        if !isnothing(_fiff_get(fiff_object, "nchan"; inst = inst))
            mi_inst = inst
            break
        end
    end
    isnothing(mi_inst) &&
        throw(ArgumentError("$file_name: no measurement info (meas_info block with nchan) found."))
    mi(tag) = _fiff_get(fiff_object, tag; inst = mi_inst)

    nchan = mi("nchan")
    ch_info = _fiff_getall(fiff_object, "ch_info"; inst = mi_inst)
    length(ch_info) == nchan ||
        throw(ArgumentError("$file_name: nchan = $nchan but $(length(ch_info)) ch_info records found."))

    # channel names: ch_info struct holds max. 15 characters
    # MNE stores longer names in per-channel ch_info blocks (FIFF_NAME tag)
    ch_names = String[ch.ch_name for ch in ch_info]
    ext_names = _fiff_getall(fiff_object, "name"; block = "ch_info")
    if length(ext_names) == nchan && all(x -> x isa String && !isempty(x), ext_names)
        ch_names = String.(ext_names)
    end

    meas_date = mi("meas_date")
    if !isnothing(meas_date)
        md = _fiff_vec(meas_date)
        secs = Int64(md[1])
        usecs = length(md) > 1 ? Int64(md[2]) : 0
        # (0, 2^31 - 1) is MNE's "no date" sentinel
        meas_date =
            (secs == 0 && usecs == typemax(Int32)) ? nothing :
            unix2datetime(secs + usecs / 1e6)
    end

    meas_info = Dict{Symbol, Any}(
        :nchan => nchan,
        :sfreq => mi("sfreq"),
        :lowpass => _fiff_nan2nothing(mi("lowpass")),
        :highpass => _fiff_nan2nothing(mi("highpass")),
        :line_freq => _fiff_nan2nothing(mi("line_freq")),
        :gantry_angle => mi("gantry_angle"),
        :data_pack => mi("data_pack"),
        :meas_date => meas_date,
        :experimenter => mi("experimenter"),
        :description => mi("description"),
        :ch_info => ch_info,
        :ch_names => ch_names,
        :bad_chs => _fiff_get(fiff_object, "bad_chs"; block = "meas_info"),
    )
    isnothing(meas_info[:sfreq]) && throw(ArgumentError("$file_name: sampling frequency not found."))

    # bad channels (MNE convention: bad_channels block, colon-separated names)
    bads = String[]
    for s in _fiff_getall(fiff_object, "mne_ch_name_list"; block = "bad_channels")
        append!(bads, _fiff_name_list(s))
    end
    meas_info[:bads] = unique(bads)

    meas_info[:subject_info] = _pack_fiff_blocks(
        fiff_object,
        "subject",
        [
            "subj_id", "subj_his_id", "subj_last_name", "subj_first_name",
            "subj_middle_name", "subj_birth_day", "subj_sex", "subj_hand",
            "subj_weight", "subj_height", "subj_comment",
        ],
    )
    meas_info[:project_info] = Dict{Symbol, Any}(
        Symbol(f) => _fiff_get(fiff_object, f; block = "project") for
        f in ["proj_id", "proj_name", "proj_aim", "proj_persons", "proj_comment"]
    )

    # digitization / HPI
    meas_info[:dig] = _fiff_getall(fiff_object, "dig_point"; block = "isotrak")
    meas_info[:hpi] = Dict(
        :hpi_result => _pack_fiff_blocks(
            fiff_object,
            "hpi_result",
            [
                "dig_point", "hpi_digitization_order", "hpi_coils_used",
                "hpi_coil_moments", "hpi_fit_goodness", "hpi_fit_good_limit",
                "hpi_fit_dist_limit", "hpi_fit_accept", "coord_trans",
            ],
        ),
        :hpi_coil => _pack_fiff_blocks(
            fiff_object,
            "hpi_coil",
            ["hpi_coil_no", "epoch", "hpi_slopes", "hpi_corr_coeff", "hpi_coil_freq", "event_bits"],
        ),
        :hpi_subsystem =>
            _pack_fiff_blocks(fiff_object, "hpi_subsystem", ["hpi_ncoil", "event_channel"]),
    )
    meas_info[:dev_head_t] = _fiff_get(fiff_object, "coord_trans"; inst = mi_inst)

    # cHPI coils: (coil number, frequency) from hpi_coil blocks that define a frequency
    hpi_coils = Tuple{Int64, Float64}[]
    for inst in _fiff_instances(fiff_object, "hpi_coil")
        f = _fiff_get(fiff_object, "hpi_coil_freq"; inst = inst)
        f isa Real || continue
        no = _fiff_get(fiff_object, "hpi_coil_no"; inst = inst)
        push!(hpi_coils, (no isa Integer ? Int64(no) : length(hpi_coils) + 1, Float64(f)))
    end
    # the same coils may be described in several hpi_meas blocks; keep the first
    meas_info[:hpi_coils] = unique(first, hpi_coils)
    # digitization order of the coils (last initial HPI fit, as in MNE-Python): coil i (in
    # coil number order) is the order[i]-th digitized HPI point
    meas_info[:hpi_dig_order] = Int64[]
    for inst in _fiff_instances(fiff_object, "hpi_result")
        o = _fiff_get(fiff_object, "hpi_digitization_order"; inst = inst)
        isnothing(o) || (meas_info[:hpi_dig_order] = Int64.(vec(collect(o))))
    end

    # SSP projectors - one entry per projection item block
    ssp = Dict{Symbol, Any}[]
    for inst in _fiff_instances(fiff_object, "xfit_proj_item")
        g = tag -> _fiff_get(fiff_object, tag; inst = inst)
        desc = g("description")
        isnothing(desc) && (desc = g("name"))
        vecs = g("proj_item_vectors")
        vecs isa AbstractVector && (vecs = reshape(vecs, 1, :))
        push!(
            ssp,
            Dict{Symbol, Any}(
                :description => isnothing(desc) ? "" : desc,
                :kind => g("proj_item_kind"),
                :nvec => g("proj_item_nvec"),
                :active => g("mne_proj_item_active"),
                :ch_names => _fiff_name_list(g("proj_item_ch_name_list")),
                :vectors => vecs,
            ),
        )
    end
    meas_info[:ssp] = ssp

    # DACQ events table
    events = Dict{Symbol, Any}[]
    for inst in _fiff_instances(fiff_object, "events")
        push!(
            events,
            Dict{Symbol, Any}(
                :channels => _fiff_get(fiff_object, "event_channels"; inst = inst),
                :list => _fiff_get(fiff_object, "event_list"; inst = inst),
            ),
        )
    end
    meas_info[:events] = events

    # MNE annotations
    annotations = Dict{Symbol, Any}[]
    for inst in _fiff_instances(fiff_object, "mne_annotations")
        g = tag -> _fiff_get(fiff_object, tag; inst = inst)
        onset = Float64.(_fiff_vec(g("mne_baseline_min")))
        offset = Float64.(_fiff_vec(g("mne_baseline_max")))
        desc = _fiff_name_list(g("description"))
        (length(onset) == length(offset) == length(desc)) || begin
            _warn("Inconsistent MNE annotations block; ignored.")
            continue
        end
        push!(
            annotations,
            Dict{Symbol, Any}(
                :onset => onset,
                :duration => offset .- onset,
                :description => desc,
                :has_orig_time => !isnothing(g("meas_date")),
            ),
        )
    end
    meas_info[:annotations] = annotations

    dacq_pars = _pack_fiff_blocks(fiff_object, "dacq_pars", ["dacq_pars", "dacq_stim"])
    if dacq_pars[:dacq_pars] isa AbstractString
        lines = rstrip.(split(dacq_pars[:dacq_pars], '\n'))
        dacq_pars[:dacq_pars] = split.(Base.filter(!isempty, lines), ' ')
    end
    meas_info[:dacq_pars] = dacq_pars

    # ------------------------------------------------------------------ #
    # raw data                                                           #
    # ------------------------------------------------------------------ #
    maxshield = false
    raw_inst = nothing
    for b in ("raw_data", "continuous_data", "maxshield_raw_data")
        i = _fiff_instances(fiff_object, b)
        if !isempty(i)
            raw_inst = i[1]
            maxshield = b == "maxshield_raw_data"
            length(i) > 1 && _warn("$file_name contains $(length(i)) $b blocks; only the first one is read.")
            break
        end
    end
    read_data && isnothing(raw_inst) && throw(
        ArgumentError(
            "$file_name contains no continuous raw data (evoked / epochs FIFF files are not supported by this importer).",
        ),
    )
    # IAS (MaxShield) raw data are flagged here; import_fiff() decides whether to accept them
    meas_info[:maxshield] = maxshield

    raw_data = if read_data
        data, first_samp, orig_format =
            _fiff_raw_data(buf, fiff_blocks, insts, raw_inst, nchan)
        Dict{Symbol, Any}(
            :data => data,
            :first_samp => first_samp,
            :first_time => first_samp / meas_info[:sfreq],
            :orig_format => orig_format,
            :next_file => _fiff_next_file(fiff_object, file_name),
        )
    else
        nothing
    end

    fiff = Dict{Symbol, Any}(:meas_info => meas_info, :raw_data => raw_data)

    return fiff, fiff_object, fiff_blocks
end

"""
    import_fiff(file_name; <keyword arguments>)

Load an Elekta/MEGIN Neuromag FIFF raw file (MEG and/or EEG) and return a `NeuroAnalyzer.NEURO` object.

Supports Vectorview/TRIUX raw and MaxFilter-processed files, MNE-Python-saved raw files (including long channel names, bad channels and annotations) and split files (`*_raw.fif`, `*_raw-1.fif`, ...).

Signals are calibrated as `raw × range × cal` (as in MNE-Python) and converted to fT, fT/cm and μV.

Files recorded with internal active shielding (IAS, "MaxShield") are refused by default, as in MNE-Python: the compensation coils distort the measured field and the data must be processed with SSS (e.g. MaxFilter) before analysis. 

When accepted with `allow_maxshield = true`, the recording header contains `:maxshield => true`.

# Arguments

- `file_name::String`: path to the `.fif` or `.fiff` file (first part of a split recording)
- `allow_maxshield::Bool=false`: import IAS (MaxShield) raw data instead of throwing an error
- `channel_locs::Symbol=:fiff`: source of channel locations in `OBJ.locs`:
    - `:fiff`: sensor/electrode positions stored in the file (see `import_locs_fiff()`); falls back to the bundled 306-channel layout if the file has no positions
    - `:layout`: bundled flattened 306-channel Vectorview/TRIUX layout (`meg_306flattened.csv`)
    - `:none`: no locations

# Returns

- `NeuroAnalyzer.NEURO`
"""
function import_fiff(
    file_name::String;
    allow_maxshield::Bool = false,
    channel_locs::Symbol = :fiff,
)::NeuroAnalyzer.NEURO
    # validate
    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    detect_format(file_name) === :fiff ||
        throw(ArgumentError("$file_name is not FIFF file."))
    channel_locs in (:fiff, :layout, :none) ||
        throw(ArgumentError("channel_locs must be :fiff, :layout or :none."))

    fiff, _, _ = load_fiff(file_name)
    mi = fiff[:meas_info]
    maxshield = mi[:maxshield]
    data = fiff[:raw_data][:data]
    first_samp = fiff[:raw_data][:first_samp]

    # ------------------------------------------------------------------ #
    # split files                                                        #
    # ------------------------------------------------------------------ #
    next_file = fiff[:raw_data][:next_file]
    expected = first_samp + size(data, 2)
    visited = Set([abspath(file_name)])
    while !isnothing(next_file)
        if abspath(next_file) in visited
            _warn("Circular split-file reference to $next_file; stopping.")
            break
        end
        if !isfile(next_file)
            _warn("Next part of split FIFF file not found: $next_file; data are incomplete.")
            break
        end
        push!(visited, abspath(next_file))
        part, _, _ = load_fiff(next_file)
        (
            part[:meas_info][:nchan] == mi[:nchan] &&
            part[:meas_info][:sfreq] == mi[:sfreq] &&
            part[:meas_info][:ch_names] == mi[:ch_names]
        ) || throw(ArgumentError("$next_file does not match channels / sampling rate of $file_name."))
        part[:raw_data][:first_samp] != expected && _warn(
            "$next_file starts at sample $(part[:raw_data][:first_samp]), expected $expected; parts are concatenated without gap filling.",
        )
        maxshield |= part[:meas_info][:maxshield]
        data = hcat(data, part[:raw_data][:data])
        expected = part[:raw_data][:first_samp] + size(part[:raw_data][:data], 2)
        next_file = part[:raw_data][:next_file]
    end

    # ------------------------------------------------------------------ #
    # IAS (MaxShield)                                                    #
    # ------------------------------------------------------------------ #
    if maxshield
        allow_maxshield || throw(
            ArgumentError(
                "$file_name contains raw MaxShield data (internal active shielding was on). These data must be processed with SSS (e.g. MaxFilter) before analysis. Use allow_maxshield = true to import them anyway.",
            ),
        )
        _warn(
            "$file_name contains raw MaxShield data (internal active shielding was on); without SSS processing these data are not suitable for analysis.",
        )
    end

    # ------------------------------------------------------------------ #
    # sampling rate                                                      #
    # ------------------------------------------------------------------ #
    sfreq = mi[:sfreq]
    sampling_rate = round(Int64, sfreq)
    isapprox(sfreq, sampling_rate; atol = 1e-6) || _warn(
        "Sampling rate $sfreq Hz is not an integer; stored as $sampling_rate Hz - time axis and marker positions will be slightly off.",
    )

    # ------------------------------------------------------------------ #
    # channels                                                           #
    # ------------------------------------------------------------------ #
    chs = mi[:ch_info]
    ch_n = mi[:nchan]
    data = reshape(data, ch_n, :, 1)

    units = Vector{String}(undef, ch_n)
    ch_type = Vector{String}(undef, ch_n)
    fiff_kind = Vector{String}(undef, ch_n)
    coil_type = Vector{String}(undef, ch_n)
    clabels = _fiff_labels(mi)

    any(ch -> ch.unit_mul != 0, chs) && _info(
        "Some channels define a non-zero unit multiplier; as in MNE-Python it is not applied to the data.",
    )

    for (idx, ch) in enumerate(chs)
        fiff_kind[idx] = _find_fiff_chtype(ch.kind)
        ch_type[idx] = fiff_kind[idx]
        ch_type[idx] == "stim" && (ch_type[idx] = "mrk")
        ch_type[idx] in (
            "ias", "syst", "exci", "quat0", "quat1", "quat2", "quat3", "quat4",
            "quat5", "quat6", "hpi_goodness", "hpi_error", "hpi_movement",
            "dipole_wave", "goodness_fit",
        ) && (ch_type[idx] = "other")
        startswith(ch_type[idx], "unknown_") && (ch_type[idx] = "other")
        coil_type[idx] = _find_fiff_coiltype(ch.coil_type)
        units[idx] = _find_fiff_unit(ch.unit)
        # FIFF "none" / "unitless" → no unit
        units[idx] in ("none", "unitless") && (units[idx] = "")
        # calibration: physical value = raw × range × cal
        @views data[idx, :, 1] .*= (ch.range * ch.cal)
    end

    # convert to standard units
    for idx = 1:ch_n
        if units[idx] == "T"
            @views data[idx, :, 1] .*= 1.0e15
            units[idx] = "fT"
        elseif units[idx] == "T/m"
            @views data[idx, :, 1] .*= (1.0e15 / 100)
            units[idx] = "fT/cm"
        elseif units[idx] == "V" && ch_type[idx] in ("eeg", "emg", "eog", "ref")
            @views data[idx, :, 1] .*= 1.0e6
            units[idx] = "μV"
        end
    end

    # classify MEG sensors (reference sensors are not counted)
    planar = (
        "nm_122", "vv_planar_w", "vv_planar_t1", "vv_planar_t2", "vv_planar_t3",
        "vv_planar_t4",
    )
    axial = (
        "axial_grad_5cm", "ctf_grad", "magnes_grad", "kit_grad", "baby_grad",
        "artemis123_grad", "kriss_grad", "compumedics_adult_grad",
        "compumedics_pediatric_grad",
    )
    mags = (
        "point_magnetometer", "nm_24", "vv_mag_w", "vv_mag_t1", "vv_mag_t2",
        "vv_mag_t3", "vv_mag_t4", "magnes_mag", "baby_mag", "quspin_zfopm_mag",
        "quspin_zfopm_mag2", "fieldline_opm_mag_gen1", "kernel_opm_mag_gen1",
    )
    magnetometers = Int64[]
    gradiometers = Int64[]
    eeg_chs = Int64[]
    for idx = 1:ch_n
        ct = coil_type[idx]
        if fiff_kind[idx] == "meg"
            if ct in planar
                coil_type[idx] = "pgrad"
                ch_type[idx] = "grad"
                push!(gradiometers, idx)
            elseif ct in axial
                coil_type[idx] = "agrad"
                ch_type[idx] = "grad"
                push!(gradiometers, idx)
            elseif ct in mags
                coil_type[idx] = "mag"
                ch_type[idx] = "mag"
                push!(magnetometers, idx)
            else
                _warn("Channel $(clabels[idx]): unknown MEG coil type $ct.")
            end
        elseif fiff_kind[idx] == "eeg"
            push!(eeg_chs, idx)
        end
    end

    # bad channels
    bad_channels = zeros(Bool, ch_n)
    for b in mi[:bads]
        i = findfirst(isequal(b), mi[:ch_names])
        isnothing(i) ? _warn("Bad channel $b not found in channel list.") : (bad_channels[i] = true)
    end
    !isnothing(mi[:bad_chs]) && _warn(
        "Legacy bad_chs tag found and ignored (as in MNE-Python); please send this file to adam.wysokinski@neuroanalyzer.org",
    )

    # ------------------------------------------------------------------ #
    # markers                                                            #
    # ------------------------------------------------------------------ #
    m_id = String[]
    m_start = Float64[]
    m_len = Float64[]
    m_val = String[]
    m_ch = Int64[]

    # DACQ event table: (sample, before, after) triplets
    for ev in mi[:events]
        isnothing(ev[:list]) && continue
        lst = _fiff_vec(ev[:list])
        length(lst) % 3 == 0 || (_warn("Malformed FIFF event list; ignored."); continue)
        evts = permutedims(reshape(Int64.(lst), 3, :))
        # sample numbers are normally absolute (include first_samp)
        offset = all(>=(first_samp), evts[:, 1]) ? first_samp : 0
        for row in axes(evts, 1)
            push!(m_id, "mrk")
            push!(m_start, (evts[row, 1] - offset) / sfreq)
            push!(m_len, 0.0)
            push!(m_val, string(evts[row, 3]))
            push!(m_ch, 0)
        end
    end

    # MNE annotations: onsets are relative to meas_date when it is stored, i.e. they include first_samp / sfreq
    for an in mi[:annotations]
        t0 = an[:has_orig_time] ? first_samp / sfreq : 0.0
        for k in eachindex(an[:onset])
            push!(m_id, "annotation")
            push!(m_start, an[:onset][k] - t0)
            push!(m_len, an[:duration][k])
            push!(m_val, an[:description][k])
            push!(m_ch, 0)
        end
    end

    markers = DataFrame(
        :id => m_id,
        :start => m_start,
        :length => m_len,
        :value => m_val,
        :channel => m_ch,
    )
    sort!(markers, :start)

    # ------------------------------------------------------------------ #
    # SSP projectors                                                     #
    # ------------------------------------------------------------------ #
    ssp_labels = String[]
    ssp_channels = zeros(Bool, ch_n)
    # full-length (ch_n) vectors
    ssp_rows = Vector{Vector{Float64}}()

    for (item_idx, item) in enumerate(mi[:ssp])
        vecs = item[:vectors]
        isnothing(vecs) && continue
        names = _clean_meg_labels(item[:ch_names])
        if size(vecs, 2) != length(names)
            _warn("SSP item $item_idx: $(size(vecs, 2)) vector elements but $(length(names)) channel names; ignored.")
            continue
        end
        ch_idx = [findfirst(isequal(l), clabels) for l in names]
        any(isnothing, ch_idx) && _warn(
            "SSP item $item_idx: $(count(isnothing, ch_idx)) channel(s) not found in the recording; their weights are dropped.",
        )
        for row in axes(vecs, 1)
            full = zeros(Float64, ch_n)
            for (j, ci) in enumerate(ch_idx)
                isnothing(ci) && continue
                full[ci] = real(vecs[row, j])
                ssp_channels[ci] = true
            end
            push!(ssp_rows, full)
            label = isempty(item[:description]) ? "PCA-v$(lpad(length(ssp_rows), 2, '0'))" : item[:description]
            size(vecs, 1) > 1 && (label *= "-$row")
            push!(ssp_labels, label)
        end
    end
    ssp_data = if isempty(ssp_rows)
        Matrix{Float64}(undef, 0, 0)
    else
        # rows = projection vectors, columns = SSP channels in recording order
        permutedims(reduce(hcat, ssp_rows))[:, ssp_channels]
    end

    # ------------------------------------------------------------------ #
    # sensor geometry                                                    #
    # ------------------------------------------------------------------ #
    # ch_info.loc: r0, ex, ey, ez (12 values, metres / unit vectors);
    # MEG in the device frame, EEG in the head frame
    coil_loc = Matrix{Float64}(permutedims(reduce(hcat, [ch.loc for ch in chs])))
    t = _fiff_dev_head_t(mi)
    dev_head_t = if isnothing(t)
        Matrix{Float64}(undef, 0, 0)
    else
        [t[1] t[2]; 0.0 0.0 0.0 1.0]
    end
    coil_id = Int64[ch.coil_type for ch in chs]
    head_origin = _fiff_head_origin(mi)
    hpi_info = _fiff_hpi(mi)

    # ------------------------------------------------------------------ #
    # data type                                                          #
    # ------------------------------------------------------------------ #
    data_type = if !isempty(magnetometers) || !isempty(gradiometers)
        "meg"
    elseif !isempty(eeg_chs)
        "eeg"
    else
        _warn("Cannot determine data type from channel list; defaulting to \"meg\".")
        "meg"
    end

    # ------------------------------------------------------------------ #
    # time axes                                                          #
    # ------------------------------------------------------------------ #
    n_samples = size(data, 2) * size(data, 3)
    time_pts = _time_axis(n_samples, sampling_rate)
    epoch_time = _time_axis(size(data, 2), sampling_rate)

    # ------------------------------------------------------------------ #
    # subject / recording metadata                                       #
    # ------------------------------------------------------------------ #
    get_field(d, k) = (v = get(d, k, nothing); isnothing(v) ? "" : string(v))
    get_num(d, k) = (v = get(d, k, nothing); v isa Real ? v : -1)
    si = mi[:subject_info]

    hand = get(si, :subj_hand, nothing)
    handedness =
        hand == 1 ? "right" : hand == 2 ? "left" : hand == 3 ? "ambidextrous" :
        isnothing(hand) ? "" : string(hand)

    subj_id = get_field(si, :subj_id)
    isempty(subj_id) && (subj_id = get_field(si, :subj_his_id))

    date = mi[:meas_date]
    rec_d, rec_t = if isnothing(date)
        "", ""
    else
        (
            string(Dates.day(date)) * "-" * string(Dates.month(date)) * "-" *
            string(Dates.year(date)),
            string(Dates.hour(date)) * ":" * string(Dates.minute(date)) * ":" *
            string(Dates.second(date)),
        )
    end

    lp = isnothing(mi[:lowpass]) ? 0 : round(mi[:lowpass]; digits = 1)
    hp = isnothing(mi[:highpass]) ? 0 : round(mi[:highpass]; digits = 1)
    lf = isnothing(mi[:line_freq]) ? 0 : round(mi[:line_freq]; digits = 1)

    # ------------------------------------------------------------------ #
    # assemble NEURO object                                              #
    # ------------------------------------------------------------------ #
    file_size_mb = round(filesize(file_name) / 1024^2; digits = 2)

    s = _create_subject(;
        id = subj_id,
        first_name = get_field(si, :subj_first_name),
        middle_name = get_field(si, :subj_middle_name),
        last_name = get_field(si, :subj_last_name),
        head_circumference = -1,
        handedness = handedness,
        weight = get_num(si, :subj_weight),
        height = get_num(si, :subj_height),
    )

    recording = isnothing(mi[:description]) ? "" : string(mi[:description])
    prefiltering = repeat(["LP: $lp Hz; HP: $hp Hz"], ch_n)

    r = if data_type == "meg"
        _create_recording_meg(;
            data_type = data_type,
            file_name = file_name,
            file_size_mb = file_size_mb,
            file_type = "FIFF",
            recording = recording,
            recording_date = rec_d,
            recording_time = rec_t,
            recording_notes = "",
            channel_type = ch_type,
            channel_order = _sort_channels(ch_type),
            reference = "",
            clabels = clabels,
            units = units,
            prefiltering = prefiltering,
            line_frequency = lf,
            sampling_rate = sampling_rate,
            magnetometers = magnetometers,
            gradiometers = gradiometers,
            coil_type = coil_type,
            bad_channels = bad_channels,
            ssp_labels = ssp_labels,
            ssp_channels = ssp_channels,
            ssp_data = ssp_data,
            maxshield = maxshield,
            coil_loc = coil_loc,
            dev_head_t = dev_head_t,
            coil_id = coil_id,
            head_origin = head_origin,
            hpi = hpi_info,
            first_time = first_samp / sampling_rate,
        )
    else
        isempty(ssp_rows) ||
            _info("SSP projectors (e.g. EEG average reference) are not stored for EEG-only data.")
        _create_recording_eeg(;
            data_type = data_type,
            file_name = file_name,
            file_size_mb = file_size_mb,
            file_type = "FIFF",
            recording = recording,
            recording_date = rec_d,
            recording_time = rec_t,
            recording_notes = "",
            channel_type = ch_type,
            channel_order = _sort_channels(ch_type),
            reference = "",
            clabels = clabels,
            transducers = repeat([""], ch_n),
            units = units,
            prefiltering = prefiltering,
            line_frequency = lf,
            sampling_rate = sampling_rate,
            gain = ones(ch_n),
            bad_channels = bad_channels,
        )
    end

    e = _create_experiment(;
        name = get_field(mi[:project_info], :proj_name),
        notes = "",
        design = "",
    )
    hdr = _create_header(; subject = s, recording = r, experiment = e)

    locs = _initialize_locs()
    obj = NeuroAnalyzer.NEURO(hdr, String[], markers, locs, time_pts, epoch_time, data)
    _initialize_locs!(obj)
    is_306 = length(magnetometers) == 102 && length(gradiometers) == 204
    l = channel_locs === :fiff ? _fiff_locs(mi, clabels) : nothing
    if !isnothing(l) && nrow(l) > 0
        add_locs!(obj; locs = l)
    elseif channel_locs !== :none && is_306
        # the bundled layout describes the 306-channel Vectorview/TRIUX helmet only
        channel_locs === :fiff && _info("Using the bundled 306-channel Vectorview/TRIUX layout.")
        add_locs!(obj; locs = import_locs_csv(joinpath(NeuroAnalyzer.res_path, "meg_306flattened.csv")))
    elseif channel_locs !== :none
        _info("No channel locations available; OBJ.locs left empty.")
    end

    _info(
        "Imported: " *
        uppercase(obj.header.recording[:data_type]) *
        " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj))" *
        "; $(round(obj.time_pts[end], digits = 2)) s)",
    )

    return obj
end

"""
    import_locs_fiff(file_name)

Load channel locations from a FIFF file (raw, evoked or any file with measurement info).

MEG sensor positions are transformed from the device to the head frame (when the device → head transform is present), EEG electrode positions are taken from the head frame.

Coordinates are centered on a sphere fitted to the head-shape digitization (or to the channel positions) and scaled so that the outermost channel has radius 1.2. Reference MEG sensors and channels without positions are omitted.

Sensor orientations (coil normals) are not stored.

# Arguments

- `file_name::String`: path to the `.fif` or `.fiff` file

# Returns

- `DataFrame`: channel locations (`label`, `loc_radius`, `loc_theta`, `loc_x`, `loc_y`, `loc_z`, `loc_radius_sph`, `loc_theta_sph`, `loc_phi_sph`); labels match those created by `import_fiff()`
"""
function import_locs_fiff(file_name::String)::DataFrame
    detect_format(file_name) === :fiff ||
        throw(ArgumentError("$file_name is not FIFF file."))
    fiff, _, _ = load_fiff(file_name; read_data = false)
    mi = fiff[:meas_info]
    return _fiff_locs(mi, _fiff_labels(mi))
end
