export load_fiff
export import_fiff

"""
    load_fiff(file_name)

Load an Elekta-Neuromag FIFF (Functional Image File Format) file and return the parsed FIFF structure.

# Arguments

- `file_name::String`: path to the `.fif` or `.fiff` file

# Returns

- `Dict{Symbol, Dict{Any, Any}}`: FIFF structure
- `Vector{Any}`: FIFF object
- `Matrix{Int64}`: FIFF blocks

# Throws

- `ArgumentError` if the file cannot be opened or is not a valid FIFF file

# References

1. Elekta Neuromag: Functional Image File Format Description. FIFF v1.3, 2011.
"""
function load_fiff(
    file_name::String
)::Tuple{
    Dict{Symbol, Dict{Any, Any}},
    Vector{Any},
    Matrix{Int64}
}

    fiff_object, fiff_blocks, buf = open(file_name, "r") do fid

        # verify first tag is file_id (required by the FIFF spec).
        tag_kind, tag_type, tag_size, data, tag_next =
            try
                _read_fiff_tag(fid)
            catch
                throw(ArgumentError("$file_name: first FIFF tag cannot be read."))
            end
        tag_kind == _find_fiff_tag("file_id") ||
            throw(ArgumentError("$file_name is not a FIFF file."))
        tag_size == 20 ||
            throw(ArgumentError("$file_name is not a valid FIFF file (unexpected file_id size)."))

        buf, fiff_blocks = _create_fiff_block(fid)
        buf, fiff_blocks # returned from do block
    end
    # file closed here in all cases

    fiff_object = Any[]
    @inbounds for block_idx in eachindex(buf)
        tag_type = fiff_blocks[block_idx, 2]
        tag_dt   = fiff_blocks[block_idx, 3]
        buf_tmp  = @views buf[block_idx]
        d = nothing

        if tag_type in [107, 108]
            # void

        elseif tag_type in [210]
            d = _find_fiff_aspect(_i32i64(buf_tmp))

        elseif tag_type in [264]
            d = _find_fiff_sss_job(_i32i64(buf_tmp))

        elseif tag_type in [270, 271, 273, 274, 275, 302, 800, 3415]
            d = @views _fiff_matrix(tag_dt, buf_tmp)

        elseif tag_type in [115]
            d = @views _i32i64(buf_tmp)
            d = d == 1 ? "prev_file" : "next_file"

        elseif tag_type in [234]
            # dig_string
            kind  = @views _i32i64(buf_tmp[1:4])
            ident = @views _i32i64(buf_tmp[5:8])
            np    = @views _i32f64(buf_tmp[9:12])
            rr    = Float64[]
            [push!(rr, @views _f32f64(buf_tmp[idx:(idx+3)])) for idx in 13:4:length(buf_tmp)]
            d = (kind, ident, np, rr)

        elseif tag_type in [255]
            # ch_pos: 12 × 4-byte floats starting at byte 5
            coil_type = @views _i32i64(buf_tmp[1:4])
            r0_1 = @views _f32f64(buf_tmp[5:8])
            r0_2 = @views _f32f64(buf_tmp[9:12])
            r0_3 = @views _f32f64(buf_tmp[13:16])
            ex_1 = @views _f32f64(buf_tmp[17:20])
            ex_2 = @views _f32f64(buf_tmp[21:24])
            ex_3 = @views _f32f64(buf_tmp[25:28])
            ey_1 = @views _f32f64(buf_tmp[29:32])
            ey_2 = @views _f32f64(buf_tmp[33:36])
            ey_3 = @views _f32f64(buf_tmp[37:40])
            ez_1 = @views _f32f64(buf_tmp[41:44])
            ez_2 = @views _f32f64(buf_tmp[45:48])
            ez_3 = @views _f32f64(buf_tmp[49:52])
            d = (coil_type, r0_1, r0_2, r0_3,
                 ex_1, ex_2, ex_3, ey_1, ey_2, ey_3, ez_1, ez_2, ez_3)

        elseif tag_type in [404, 405, 406]
            d = @views _i32i64(buf_tmp)
            tag_type == 404 && (d = unix2datetime(d))

        elseif tag_type in [242]
            d = @views _ui32i32(buf_tmp)

        elseif tag_type in [222]
            # coord_trans: rotation (3×3 Float32), translation (3 Float32),
            # inverse rotation (3×3 Float32), inverse translation (3 Float32).
            # FIX: original used _i16f64 (reads 2-byte Int16) for all fields
            # here — these are 4-byte Float32; must use _f32f64.
            from = @views _i32i64(buf_tmp[1:4])
            to   = @views _i32i64(buf_tmp[5:8])
            rot  = zeros(3, 3)
            for row in 1:3
                d_row = Float64[]
                [push!(d_row, @views _f32f64(buf_tmp[(8 + col):(8 + col + 3)])) for col in 1:4:9]
                rot[row, :] = d_row
            end
            move = Float64[]
            [push!(move, @views _f32f64(buf_tmp[(44 + i):(44 + i + 3)])) for i in 1:4:9]
            invrot = zeros(3, 3)
            for row in 1:3
                d_row = Float64[]
                [push!(d_row, @views _f32f64(buf_tmp[(56 + col):(56 + col + 3)])) for col in 1:4:9]
                invrot[row, :] = d_row
            end
            invmove = Float64[]
            [push!(invmove, @views _f32f64(buf_tmp[(92 + i):(92 + i + 3)])) for i in 1:4:9]
            d = (from, to, rot, move, invrot, invmove)

        elseif tag_type in [300]
            # data_buffer
            df = @views _find_fiff_dt(tag_dt)
            d  = Float64[]
            if df == "dau_pack16" || df == "int16"
                [push!(d, @views _i16f64(buf_tmp[idx:(idx+1)])) for idx in 1:2:length(buf_tmp)]
            elseif df == "int32"
                [push!(d, @views _i32f64(buf_tmp[idx:(idx+3)])) for idx in 1:4:length(buf_tmp)]
            elseif df == "float"
                [push!(d, @views _f32f64(buf_tmp[idx:(idx+3)])) for idx in 1:4:length(buf_tmp)]
            else
                _warn("Data type $df not implemented; please send this file to adam.wysokinski@neuroanalyzer.org")
            end

        elseif tag_type in [213]
            # dig_point
            kind  = @views _i32i64(buf_tmp[1:4])
            ident = @views _i32i64(buf_tmp[5:8])
            r_1   = @views _f32f64(buf_tmp[9:12])
            r_2   = @views _f32f64(buf_tmp[13:16])
            r_3   = @views _f32f64(buf_tmp[17:20])
            d = (kind, ident, r_1, r_2, r_3)

        elseif tag_type in [252]
            d = @views _find_fiff_chtype(_i32i64(buf_tmp))
        elseif tag_type in [256]
            d = _find_fiff_unit(_i32i64(buf_tmp))
        elseif tag_type in [3411]
            d = _find_fiff_proj_item(_i32i64(buf_tmp))
        elseif tag_type in [3416]
            d = _find_fiff_proj_by(_i32i64(buf_tmp))
        elseif tag_type in [280]
            d = _find_fiff_gantry_type(_i32i64(buf_tmp))

        elseif tag_type in [203]
            # ch_info (75-byte record)
            scan_no   = @views _i32i64(buf_tmp[1:4])
            log_no    = @views _i32i64(buf_tmp[5:8])
            kind      = @views _find_fiff_chtype(_i32i64(buf_tmp[9:12]))
            range_val = @views _f32f64(buf_tmp[13:16])
            cal       = @views _f32f64(buf_tmp[17:20])
            coil_type = @views _i32i64(buf_tmp[21:24])
            r0_1 = @views _f32f64(buf_tmp[25:28]); r0_2 = @views _f32f64(buf_tmp[29:32])
            r0_3 = @views _f32f64(buf_tmp[33:36])
            ex_1 = @views _f32f64(buf_tmp[37:40]); ex_2 = @views _f32f64(buf_tmp[41:44])
            ex_3 = @views _f32f64(buf_tmp[45:48])
            ey_1 = @views _f32f64(buf_tmp[49:52]); ey_2 = @views _f32f64(buf_tmp[53:56])
            ey_3 = @views _f32f64(buf_tmp[57:60])
            ez_1 = @views _f32f64(buf_tmp[61:64]); ez_2 = @views _f32f64(buf_tmp[65:68])
            ez_3 = @views _f32f64(buf_tmp[69:72])
            unit     = @views _find_fiff_unit(_i32i64(buf_tmp[73:76]))
            unit_mul = @views _find_fiff_mul(_i32i64(buf_tmp[77:80]))
            d = (scan_no, log_no, kind, range_val, cal, coil_type,
                 r0_1, r0_2, r0_3, ex_1, ex_2, ex_3,
                 ey_1, ey_2, ey_3, ez_1, ez_2, ez_3, unit, unit_mul)

        elseif tag_type in [214]
            # ch_pos_vec (obsolete)

        elseif tag_type in [
                111, 112, 113, 114, 118, 150, 151, 205, 206, 212, 227, 233,
                237, 258, 281, 401, 402, 403, 409, 410, 501, 502, 503, 504,
                602, 2020, 3102, 3300, 3406, 3407, 3417, 3501]
            d = @views _v2s(string.(Char.(buf_tmp)))

        elseif tag_type in [282]
            d = _i8i8(buf_tmp[1:4])

        elseif tag_type in [
                101, 102, 104, 105, 106, 117, 200, 202, 204, 207, 208, 209,
                211, 216, 217, 221, 228, 230, 231, 245, 250, 251, 257, 263,
                266, 267, 268, 277, 278, 301, 303, 400, 500, 701, 702, 703,
                2004, 2010, 2012, 2014, 2023, 2032, 3013, 3104, 3414]
            d = @views _i32i64(buf_tmp[1:4])
            tag_type == 204 && (d = unix2datetime(d))

        elseif tag_type in [220, 232, 225, 246, 247, 269, 304, 305, 600, 601, 603, 3413]
            d = Int64[]
            [push!(d, @views _i32i64(buf_tmp[idx:(idx+3)])) for idx in 1:4:length(buf_tmp)]

        elseif tag_type in [276]
            d = @views _f32f64(buf_tmp[1:4])

        elseif tag_type in [
                201, 218, 219, 223, 229, 235, 236, 240, 241, 243, 244,
                253, 254, 272, 279, 407, 408, 2005, 2009, 2011, 2013,
                2015, 2016, 2018, 2019, 2024, 2040, 3109, 3113, 3405, 3412]
            d = @views _f32f64(buf_tmp[1:4])

        elseif tag_type in [215, 224, 226, 265]
            d = Float64[]
            [push!(d, @views _f32f64(buf_tmp[idx:(idx+3)])) for idx in 1:4:length(buf_tmp)]

        elseif tag_type in [100, 103, 109, 110, 116, 120]
            # id_t
            fiff_v_major     = @views _i16i64(buf_tmp[1:2])
            fiff_v_minor     = @views _i16i64(buf_tmp[3:4])
            mach_id1         = @views _i32i64(buf_tmp[5:8])
            mach_id2         = @views _i32i64(buf_tmp[9:12])
            time_sec         = @views _i32i64(buf_tmp[13:16])
            id_creation_date = unix2datetime(time_sec)
            time_usec        = @views _i32i64(buf_tmp[17:20])
            d = (fiff_v_major, fiff_v_minor, mach_id1, mach_id2, id_creation_date, time_usec)

        else
            _warn("Tag $tag_type not implemented; please send this file to adam.wysokinski@neuroanalyzer.org")
        end

        push!(fiff_object, (
            block_idx,
            _find_fiff_tag(fiff_blocks[block_idx, 2]),
            _find_fiff_block(fiff_blocks[block_idx, end]),
            d,
        ))

    end

    # ------------------------------------------------------------------ #
    # process blocks into structured Dicts                               #
    # ------------------------------------------------------------------ #
    bidx, btypes = _get_blocks(fiff_blocks)

    project_info = Dict()
    for f in ["proj_id", "proj_name"]
        tmp = fiff_object[[fiff_object[idx][2] for idx in eachindex(fiff_object)] .== f]
        push!(project_info, Symbol(f) => isempty(tmp) ? nothing : tmp[1][4])
    end

    ref = Dict()
    for f in ["ref_role", "ref_file_id", "ref_file_num", "ref_file_name", "ref_block_id"]
        tmp = fiff_object[[fiff_object[idx][2] for idx in eachindex(fiff_object)] .== f]
        push!(ref, Symbol(f) => isempty(tmp) ? nothing : tmp[1][4])
    end
    !isnothing(ref[:ref_role]) && _warn(
        "Data split into multiple files is not implemented yet; please send this file to adam.wysokinski@neuroanalyzer.org")

    meas_info = Dict()
    for f in ["sfreq", "lowpass", "highpass", "data_pack", "line_freq", "gantry_angle", "bad_chs"]
        tmp = fiff_object[[fiff_object[idx][2] for idx in eachindex(fiff_object)] .== f]
        push!(meas_info, Symbol(f) => isempty(tmp) ? nothing : tmp[1][4])
    end

    ch_info = Dict()
    block_idx = [fiff_object[idx][2] for idx in eachindex(fiff_object)] .== "ch_info"
    tmp = fiff_object[block_idx][[fiff_object[block_idx][idx][2] for idx in eachindex(fiff_object[block_idx])] .== "ch_info"]
    if length(tmp) == 1
        push!(ch_info, Symbol("ch_info") => tmp[1][4])
    elseif length(tmp) > 1
        [push!(ch_info, Symbol("ch_info_$n") => tmp[n][4]) for n in eachindex(tmp)]
    else
        push!(ch_info, Symbol("ch_info") => nothing)
    end
    push!(meas_info, :ch_info => ch_info)

    subject_info = _pack_fiff_blocks(fiff_object, "subject", [
        "subj_id", "subj_his_id", "subj_last_name", "subj_first_name",
        "subj_middle_name", "subj_birth_day", "subj_sex", "subj_hand",
        "subj_weight", "subj_height", "comment"])

    hpi_result = _pack_fiff_blocks(fiff_object, "hpi_result", [
        "dig_point", "hpi_digitization_order", "hpi_coils_used",
        "hpi_coil_moments", "hpi_fit_goodness", "hpi_fit_good_limit",
        "hpi_fit_dist_limit", "hpi_fit_accept", "coord_trans"])
    hpi_coil   = _pack_fiff_blocks(fiff_object, "hpi_coil",
        ["hpi_coil_no", "epoch", "hpi_slopes", "hpi_corr_coeff", "hpi_coil_freq"])
    isotrak    = _pack_fiff_blocks(fiff_object, "isotrak", ["dig_point"])
    hpi = Dict(:hpi_result => hpi_result, :hpi_coil => hpi_coil, :isotrak => isotrak)

    hpi_coil2      = _pack_fiff_blocks(fiff_object, "hpi_coil", ["event_bits"])
    hpi_subsystem  = _pack_fiff_blocks(fiff_object, "hpi_subsystem", ["hpi_ncoil", "event_channel"])
    hpi_subsystem  = Dict(:hpi_subsystem => hpi_subsystem, :hpi_coil => hpi_coil2)

    xfit_proj_item = _pack_fiff_blocks(fiff_object, "xfit_proj_item", [
        "description", "proj_item_kind", "nchan", "proj_item_ch_name_list",
        "proj_item_time", "proj_item_nvec", "proj_item_vectors"])
    nchan     = _pack_fiff_blocks(fiff_object, "xfit_proj_item", ["nchan"])
    xfit_proj = Dict(:nchan => nchan, :xfit_proj_item => xfit_proj_item)

    events   = _pack_fiff_blocks(fiff_object, "events", ["event_channels", "event_list"])
    dacq_pars = _pack_fiff_blocks(fiff_object, "dacq_pars", ["dacq_pars", "dacq_stim"])
    dacq_pars[:dacq_pars] = split(dacq_pars[:dacq_pars], "\n")
    dacq_pars[:dacq_pars][end] == "" &&
        deleteat!(dacq_pars[:dacq_pars], lastindex(dacq_pars[:dacq_pars]))
    dacq_pars[:dacq_pars] = split.(rstrip.(dacq_pars[:dacq_pars]), ' ')

    merge!(meas_info, _pack_fiff_blocks(fiff_object, "meas_info",
        ["experimenter", "description", "meas_date"]))

    push!(meas_info, :subject_info => subject_info)
    push!(meas_info, :project_info => project_info)
    push!(meas_info, :hpi          => hpi)
    push!(meas_info, :ssp          => xfit_proj)
    push!(meas_info, :events       => events)
    push!(meas_info, :dacq_pars    => dacq_pars)
    push!(meas_info, :nchan        => length(ch_info))
    meas_info[:sfreq] = round(Int64, meas_info[:sfreq])

    raw_data_tmp1 = sort(collect(_pack_fiff_blocks(fiff_object, "raw_data", ["data_buffer"])))
    buffer_length = length(raw_data_tmp1[1][2])
    raw_data_flat = Float64[]
    for item in raw_data_tmp1
        append!(raw_data_flat, item[2])
    end

    raw_data = Dict(:raw_data => raw_data_flat)
    merge!(raw_data, _pack_fiff_blocks(fiff_object, "raw_data",
        ["first_samp", "data_skip", "data_skip_samp"]))

    isnothing(raw_data[:data_skip]) && (raw_data[:data_skip] = 0)
    isnothing(raw_data[:data_skip_samp]) && (raw_data[:data_skip_samp] = 0)

    if isnothing(raw_data[:first_samp])
        raw_data[:first_samp] = 0
    else
        raw_data[:first_samp] /= meas_info[:sfreq]
    end

    raw_data[:data_skip] > 0 && _warn("data_skip not implemented; please send this file to adam.wysokinski@neuroanalyzer.org")
    raw_data[:data_skip_samp] > 0 && _warn("data_skip_samp not implemented; please send this file to adam.wysokinski@neuroanalyzer.org")

    fiff = Dict(:meas_info => meas_info, :raw_data => raw_data)

    return fiff, fiff_object, fiff_blocks

end

"""
    import_fiff(file_name)

Load an Elekta-Neuromag FIFF file (MEG or EEG) and return a `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: path to the `.fif` or `.fiff` file

# Returns

- `NeuroAnalyzer.NEURO`

# Throws

- `ArgumentError` if the file does not exist
"""
function import_fiff(file_name::String)::NeuroAnalyzer.NEURO

    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))

    fiff, _, _ = load_fiff(file_name)

    sampling_rate = fiff[:meas_info][:sfreq]
    ch_n = fiff[:meas_info][:nchan]
    data = @views reshape(fiff[:raw_data][:raw_data], ch_n, :, 1)

    units = repeat([""], ch_n)
    ch_type = repeat([""], ch_n)
    coil_type = repeat([""], ch_n)
    clabels = repeat([""], ch_n)

    @inbounds for (k, v) in fiff[:meas_info][:ch_info]
        ch = v[1]
        ch_type[ch] = v[3]
        ch_type[ch] == "stim" && (ch_type[ch] = "mrk")
        ch_type[ch] in ("ias", "syst") && (ch_type[ch] = "other")
        range_val = v[4]
        cal = v[5]
        coil_type[ch] = _find_fiff_coiltype(v[6])
        units[ch] = v[19]
        unit_mul = v[20] != "none" ? 10^v[20] : 1
        clabels[ch] = uppercase(v[3]) * " " * string(v[2])
        @views data[ch, :, 1] .*= (range_val * cal * unit_mul)
    end

    clabels = _clean_meg_labels(clabels)

    # convert to standard units
    @inbounds for ch_idx in 1:ch_n
        if units[ch_idx] == "T"
            @views data[ch_idx, :, 1] .*= 1e15;  units[ch_idx] = "fT"
        elseif units[ch_idx] == "T/m"
            @views data[ch_idx, :, 1] .*= (1e15 / 100);  units[ch_idx] = "fT/cm"
        elseif units[ch_idx] == "T/cm"
            @views data[ch_idx, :, 1] .*= 1e15;  units[ch_idx] = "fT/cm"
        elseif units[ch_idx] == "V" && ch_type[ch_idx] in ("eeg", "emg", "eog", "ref")
            @views data[ch_idx, :, 1] .*= 1e6;   units[ch_idx] = "μV"
        end
    end

    # classify coil types
    magnetometers = Int64[];  gradiometers = Int64[];  eeg_chs = Int64[]
    @inbounds for ch_idx in 1:ch_n
        ct = coil_type[ch_idx]
        if ct in ("vv_planar_w", "vv_planar_t1", "vv_planar_t2", "vv_planar_t3")
            coil_type[ch_idx] = "pgrad";  ch_type[ch_idx] = "grad";  push!(gradiometers, ch_idx)
        elseif ct in ("magnes_grad",)
            coil_type[ch_idx] = "grad";   ch_type[ch_idx] = "grad";  push!(gradiometers, ch_idx)
        elseif ct in ("axial_grad_5cm", "ctf_grad")
            coil_type[ch_idx] = "agrad";  ch_type[ch_idx] = "grad";  push!(gradiometers, ch_idx)
        elseif ct in ("point_magnetometer", "vv_mag_w", "vv_mag_t1",
                      "vv_mag_t2", "vv_mag_t3", "magnes_mag")
            coil_type[ch_idx] = "mag";    ch_type[ch_idx] = "mag";   push!(magnetometers, ch_idx)
        elseif ct == "eeg"
            push!(eeg_chs, ch_idx)
        end
    end

    bad_channels = zeros(Bool, ch_n)
    !isnothing(fiff[:meas_info][:bad_chs]) && _warn(
        "bad_channels tag not implemented; please send this file to adam.wysokinski@neuroanalyzer.org")

    # HPI digitization points
    id = Int64[]
    p = Int64[]
    x = Float64[]
    y = Float64[]
    z = Float64[]
    for (k, v) in fiff[:meas_info][:hpi][:isotrak]
        push!(id, parse(Int64, match(r"(.+)_(\d+)", string(k)).captures[2]))
        push!(p, v[2]);  push!(x, v[3]);  push!(y, v[4]);  push!(z, v[5])
    end
    hpi = sort(DataFrame(:id => id, :p => p, :x => x, :y => y, :z => z), :id)

    # ------------------------------------------------------------------ #
    # build markers from event table                                     #
    # ------------------------------------------------------------------ #
    events_ch  = fiff[:meas_info][:events][:event_channels]
    event_list = fiff[:meas_info][:events][:event_list]

    markers = if isnothing(event_list) || isempty(event_list)
        DataFrame(:id => String[], :start => Float64[],
                  :length => Float64[], :value => String[], :channel => Int64[])
    else
        evts = reshape(event_list, 3, :)' # columns: sample, before, after
        n_ev = size(evts, 1)
        DataFrame(
            :id => fill("mrk", n_ev),
            :start => Float64.(evts[:, 1]) ./ sampling_rate,
            :length => fill(0.0, n_ev),
            :value => string.(evts[:, 3]), # "after" encodes the event code
            :channel => fill(0, n_ev),
        )
    end

    # SSP projectors
    ssp_labels = String[]
    ssp_channels = zeros(Bool, ch_n)
    ssp_data = Matrix{Float64}(undef, 0, 0)

    if :ssp in keys(fiff[:meas_info])
        ssp = fiff[:meas_info][:ssp]
        if !isnothing(ssp[:nchan][:nchan])
            xfit_proj_item = ssp[:xfit_proj_item]
            xfit_proj_vecs = Vector{Vector{Float64}}()
            idx = 1
            while Symbol("proj_item_vectors_$idx") in keys(xfit_proj_item)
                push!(xfit_proj_vecs, xfit_proj_item[Symbol("proj_item_vectors_$idx")][:])
                idx += 1
            end
            raw_ssp_labels = string.(split(xfit_proj_item[:proj_item_ch_name_list_1], ':'))
            raw_ssp_labels = _clean_meg_labels(raw_ssp_labels)
            ssp_data       = zeros(length(xfit_proj_vecs), length(xfit_proj_vecs[1]))
            for (i, v) in enumerate(xfit_proj_vecs)
                ssp_data[i, :] = v
            end
            ssp_sorting_idx = [findfirst(isequal(l), clabels) for l in raw_ssp_labels]
            ssp_channels[filter(!isnothing, ssp_sorting_idx)] .= true
            ssp_labels = ["PCA-v$(lpad(i, 2, '0'))" for i in eachindex(xfit_proj_vecs)]
        end
    end

    # ------------------------------------------------------------------ #
    # data type detection                                                  #
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
    n_samples  = size(data, 2) * size(data, 3)
    time_pts   = round.(range(0; step = 1/sampling_rate, length = n_samples);  digits = 4)
    epoch_time = round.(range(0; step = 1/sampling_rate, length = size(data,2)); digits = 4)

    # subject / recording metadata
    get_field(d, k) = begin v = d[k]; isnothing(v) ? "" : string(v) end
    get_num(d, k)   = begin v = d[k]; isnothing(v) ? -1  : v         end
    si = fiff[:meas_info][:subject_info]

    date = fiff[:meas_info][:meas_date]
    rec_d, rec_t = if isnothing(date)
        "", ""
    else
        (string(Dates.day(date))  * "-" * string(Dates.month(date))  * "-" * string(Dates.year(date)),
         string(Dates.hour(date)) * ":" * string(Dates.minute(date)) * ":" * string(Dates.second(date)))
    end

    lp = isnothing(fiff[:meas_info][:lowpass])   ? 0 : round(fiff[:meas_info][:lowpass];   digits=1)
    hp = isnothing(fiff[:meas_info][:highpass])  ? 0 : round(fiff[:meas_info][:highpass];  digits=1)
    lf = isnothing(fiff[:meas_info][:line_freq]) ? 0 : round(fiff[:meas_info][:line_freq]; digits=1)

    # ------------------------------------------------------------------ #
    # assemble NEURO object                                              #
    # ------------------------------------------------------------------ #
    file_size_mb = round(filesize(file_name) / 1024^2; digits = 2)

    s = _create_subject(
        id = get_field(si, :subj_id),
        first_name = get_field(si, :subj_first_name),
        middle_name = get_field(si, :subj_middle_name),
        last_name = get_field(si, :subj_last_name),
        head_circumference = -1,
        handedness = get_field(si, :subj_hand),
        weight = get_num(si,   :subj_weight),
        height = get_num(si,   :subj_height))
    r = _create_recording_meg(
        data_type = data_type,
        file_name = file_name,
        file_size_mb = file_size_mb,
        file_type = "FIFF",
        recording = get_field(fiff[:meas_info], :description),
        recording_date = rec_d,
        recording_time = rec_t,
        recording_notes = "",
        channel_type = ch_type,
        channel_order = _sort_channels(ch_type),
        reference = "",
        clabels = clabels,
        units = units,
        prefiltering = repeat(["LP: $lp Hz; HP: $hp Hz"], ch_n),
        line_frequency = lf,
        sampling_rate = sampling_rate,
        magnetometers = magnetometers,
        gradiometers = gradiometers,
        coil_type = coil_type,
        bad_channels = bad_channels,
        ssp_labels = ssp_labels,
        ssp_channels = ssp_channels,
        ssp_data = ssp_data)
    e   = _create_experiment(name = get_field(fiff[:meas_info][:project_info], :proj_name),
                             notes = "",
                             design = "")
    hdr = _create_header(subject = s, recording = r, experiment = e)

    locs = _initialize_locs()
    obj  = NeuroAnalyzer.NEURO(hdr, String[], markers, locs, time_pts, epoch_time, data)
    _initialize_locs!(obj)
    l = import_locs_csv(joinpath(NeuroAnalyzer.res_path, "meg_306flattened.csv"))
    add_locs!(obj; locs = l)

    _info("Imported: " *
        uppercase(obj.header.recording[:data_type]) *
        " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj))" *
        "; $(round(obj.time_pts[end]; digits=2)) s)")

    return obj

end