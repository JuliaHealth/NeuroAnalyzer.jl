export load_fiff
export import_fiff

"""
    load_fiff(file_name)

Load Elekta-Neuromag 306 FIFF (Functional Image File Format) file (MEG, EEG) and return FIFF object.

# Arguments

- `file_name::String`: name of the file to load

# Returns

- `fiff::Dict{Symbol, Dict{Any, Any}}`
- `fiff_object::Vector{Any}`
- `fiff_blocks::Matrix{Int64}`

# Source

1. Elekta Neuromag: Functional Image File Format Description. FIFF version 1.3. March 2011
"""
function load_fiff(file_name::String)::Tuple{Dict{Symbol, Dict{Any, Any}}, Vector{Any}, Matrix{Int64}}

    fid = nothing
    try
        fid = open(file_name, "r")
    catch
        error("File $file_name cannot be loaded.")
    end

    tag_kind = nothing
    tag_type = nothing
    tag_size = nothing
    data = nothing
    tag_next = nothing

    # check file_id tag
    try
        tag_kind, tag_type, tag_size, data, tag_next = _read_fiff_tag(fid)
    catch
        error("File $file_name first tag cannot be read.")
    end
    @assert tag_kind == _find_fiff_tag("file_id") "File $file_name is not a FIFF file."
    @assert tag_size == 20 "File $file_name is not a FIFF file."

    buf, fiff_blocks = _create_fiff_block(fid)

    close(fid)

    fiff_object = Any[]
    @inbounds for block_idx in eachindex(buf)
        tag_type = fiff_blocks[block_idx, 2]
        tag_dt = fiff_blocks[block_idx, 3]
        buf_tmp = @views buf[block_idx]
        d = nothing
        if tag_type in [107, 108]
            # void
        elseif tag_type in [210]
            # aspect
            d = _find_fiff_aspect(_i32i64(buf_tmp))
        elseif tag_type in [264]
            # sss_job
            d = _find_fiff_sss_job(_i32i64(buf_tmp))
        elseif tag_type in [270, 271, 273, 274, 275, 302, 800, 3415]
            # sss_cal_chans
            # sss_cal_cors
            # sss_base_in
            # sss_base_out
            # sss_base_virt
            # epoch
            # decoupler_matrix
            # proj_item_vectors
            d = @views NeuroAnalyzer._fiff_matrix(tag_dt, buf_tmp)
        elseif tag_type in [115]
            # role
            d = @views _i32i64(buf_tmp)
            d = d == 1 ? "prev_file" : "next_file"
        elseif tag_type in [234]
            # dig_string
            kind = @views _i32i64(buf_tmp[1:4])
            ident = @views _i32i64(buf_tmp[5:8])
            np = @views _i32f64(buf_tmp[9:12])
            rr = Float64[]
            [push!(rr, @views _f32f64(buf_tmp[idx:(idx + 3)])) for idx in 13:4:length(buf_tmp)]
            d = (kind, ident, np, rr)
        elseif tag_type in [255]
            # ch_pos
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
            ez_3 = @views _f32f64(buf_tmp[59:52])
            d = (coil_type, r0_1, r0_2, r0_3, ex_1, ex_2, ex_3, ey_1, ey_2, ey_3, ez_1, ez_2, ez_3)
        elseif tag_type in [404, 405, 406]
            # dob, sex, handedness
            d = @views _i32i64(buf_tmp)
            tag_type == 404 && (d = unix2datetime(d))
        elseif tag_type in [242]
            # hpi_mask
            d = @views _ui32i32(buf_tmp)
        elseif tag_type in [222]
            # coord_trans
            from = @views _i32i64(buf_tmp[1:4])
            to = @views _i32i64(buf_tmp[5:8])
            rot = zeros(3, 3)
            for idx1 in 1:3
                d = Float64[]
                [push!(d, @views _i16f64(buf_tmp[(8 + idx2):(8 + idx2 + 3)])) for idx2 in 1:4:9]
                rot[idx1, :] = d
            end
            move = Float64[]
            [push!(move, @views _i16f64(buf_tmp[(44 + idx):(44 + idx + 3)])) for idx in 1:4:9]
            invrot = zeros(3, 3)
            for idx1 in 1:3
                d = Float64[]
                [push!(d, @views _i16f64(buf_tmp[(56 + idx2):(56 + idx2 + 3)])) for idx2 in 1:4:9]
                invrot[idx1, :] = d
            end
            invmove = Float64[]
            [push!(invmove, @views _i16f64(buf_tmp[(92 + idx):(92 + idx + 3)])) for idx in 1:4:9]
            d = (from, to, rot, move, invrot, invmove)
        elseif tag_type in [300]
            # data_buffer
            df = @views _find_fiff_dt(tag_dt)
            d = Float64[]
            if df == "dau_pack16" || df == "int16"
                [push!(d, @views _i16f64(buf_tmp[idx:(idx + 1)])) for idx in 1:2:length(buf_tmp)]
            elseif df == "int32"
                [push!(d, @views _i32f64(buf_tmp[idx:(idx + 3)])) for idx in 1:4:length(buf_tmp)]
            elseif df == "float"
                [push!(d, @views _f32f64(buf_tmp[idx:(idx + 3)])) for idx in 1:4:length(buf_tmp)]
            else
                _warn("Data type $df is not implemented yet; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org")
            end
        elseif tag_type in [213]
            # dig_point
            kind = @views _i32i64(buf_tmp[1:4])
            ident = @views _i32i64(buf_tmp[5:8])
            r_1 = @views _f32f64(buf_tmp[9:12])
            r_2 = @views _f32f64(buf_tmp[13:16])
            r_3 = @views _f32f64(buf_tmp[17:20])
            d = (kind, ident, r_1, r_2, r_3)
        elseif tag_type in [252]
            # ch_type
            d =  @views _find_fiff_chtype(_i32i64(buf_tmp))
        elseif tag_type in [256]
            # ch_unit
            d = _find_fiff_unit(_i32i64(buf_tmp))
        elseif tag_type in [3411]
            # proj_item
            d = _find_fiff_proj_item(_i32i64(buf_tmp))
        elseif tag_type in [3416]
            # proj_item
            d = _find_fiff_proj_by(_i32i64(buf_tmp))
        elseif tag_type in [280]
            # gantry_type
            d = _find_fiff_gantry_type(_i32i64(buf_tmp))
        elseif tag_type in [203]
            # ch_info
            scan_no = @views _i32i64(buf_tmp[1:4])
            log_no = @views _i32i64(buf_tmp[5:8])
            kind = @views _find_fiff_chtype(_i32i64(buf_tmp[9:12]))
            range = @views _f32f64(buf_tmp[13:16])
            cal = @views _f32f64(buf_tmp[17:20])
            coil_type = @views _i32i64(buf_tmp[21:24])
            r0_1 = @views _f32f64(buf_tmp[25:28])
            r0_2 = @views _f32f64(buf_tmp[29:32])
            r0_3 = @views _f32f64(buf_tmp[33:36])
            ex_1 = @views _f32f64(buf_tmp[37:40])
            ex_2 = @views _f32f64(buf_tmp[41:44])
            ex_3 = @views _f32f64(buf_tmp[45:48])
            ey_1 = @views _f32f64(buf_tmp[49:52])
            ey_2 = @views _f32f64(buf_tmp[53:56])
            ey_3 = @views _f32f64(buf_tmp[57:60])
            ez_1 = @views _f32f64(buf_tmp[61:64])
            ez_2 = @views _f32f64(buf_tmp[65:68])
            ez_3 = @views _f32f64(buf_tmp[69:72])
            unit = @views _find_fiff_unit(_i32i64(buf_tmp[73:76]))
            unit_mul = @views _find_fiff_mul(_i32i64(buf_tmp[77:80]))
            d = (scan_no, log_no, kind, range, cal, coil_type, r0_1, r0_2, r0_3, ex_1, ex_2, ex_3, ey_1, ey_2, ey_3, ez_1, ez_2, ez_3, unit, unit_mul)
        elseif tag_type in [214]
            # ch_pos_vec
            # obsolete
        elseif tag_type in [111, 112, 113, 114, 118, 150, 151, 205, 206, 212, 227, 233, 237, 258, 281, 401, 402, 403, 409, 410, 501, 502, 503, 504, 602, 2020, 3102, 3300, 3406, 3407, 3417, 3501]
            # string
            d = @views _v2s(string.(Char.(buf_tmp)))
        elseif tag_type in [282]
            # int8
            d = _i8i8(buf_tmp[1:4])
        elseif tag_type in [101, 102, 104, 105, 106, 117, 200, 202, 204, 207, 208, 209, 211, 216, 217, 221, 228, 230, 231, 245, 250, 251, 257, 263, 266, 267, 268, 277, 278, 301, 303, 400, 500, 701, 702, 703, 2004, 2010, 2012, 2014, 2023, 2032, 3013, 3104, 3414]
            # int32
            d = @views _i32i64(buf_tmp[1:4])
            if tag_type == 204
                d = unix2datetime(d)
            end
        elseif tag_type in [220, 232, 225, 246, 247, 269, 304, 305, 600, 601, 603, 3413]
            # int32*
            d = Int64[]
            [push!(d, @views _i32i64(buf_tmp[idx:(idx + 3)])) for idx in 1:4:length(buf_tmp)]
        elseif tag_type in [276]
            # double
            d = @views _f32f64(buf_tmp[1:4])
        elseif tag_type in [201, 218, 219, 223, 229, 235, 236, 240, 241, 243, 244, 253, 254, 272, 279, 407, 408, 2005, 2009, 2011, 2013, 2015, 2016, 2018, 2019, 2024, 2040, 3109, 3113, 3405, 3412]
            # float
            d = @views _f32f64(buf_tmp[1:4])
        elseif tag_type in [215, 224, 226, 265]
            # float*
            d = Float64[]
            [push!(d, @views _f32f64(buf_tmp[idx:(idx + 3)])) for idx in 1:4:length(buf_tmp)]
        elseif tag_type in [100, 103, 109, 110, 116, 120]
            # id_t
            fiff_v_major = @views _i16i64(buf_tmp[1:2])
            fiff_v_minor = @views _i16i64(buf_tmp[3:4])
            mach_id1 = @views _i32i64(buf_tmp[5:8])
            mach_id2 = @views _i32i64(buf_tmp[9:12])
            time_sec = @views _i32i64(buf_tmp[13:16])
            id_creation_date = unix2datetime(time_sec)
            time_usec = @views _i32i64(buf_tmp[17:20])
            d = (fiff_v_major, fiff_v_minor, mach_id1, mach_id2, id_creation_date, time_usec)
        else
            _warn("$tag_type is not implemented yet; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org")
        end
        push!(fiff_object, (block_idx, _find_fiff_tag(fiff_blocks[block_idx, 2]), _find_fiff_block(fiff_blocks[block_idx, end]), d))
    end

    # process blocks
    bidx, btypes = _get_blocks(fiff_blocks)
    project_info = Dict()
    fields = ["proj_id", "proj_name"]
    @inbounds for f in fields
        tmp = fiff_object[[fiff_object[idx][2] for idx in eachindex(fiff_object)] .== f]
        if length(tmp) != 0
            push!(project_info, Symbol(f)=>tmp[1][4])
        else
            push!(project_info, Symbol(f)=>nothing)
        end
    end

    ref = Dict()
    fields = ["ref_role", "ref_file_id", "ref_file_num", "ref_file_name", "ref_block_id"]
    @inbounds for f in fields
        tmp = fiff_object[[fiff_object[idx][2] for idx in eachindex(fiff_object)] .== f]
        if length(tmp) != 0
            push!(ref, Symbol(f)=>tmp[1][4])
        else
            push!(ref, Symbol(f)=>nothing)
        end
    end
    !isnothing(ref[:ref_role]) && _warn("Data split into multiple files is not implemented yet; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org")

    meas_info = Dict()
    fields = ["sfreq", "lowpass", "highpass", "data_pack", "line_freq", "gantry_angle", "bad_chs"]
    @inbounds for f in fields
        tmp = fiff_object[[fiff_object[idx][2] for idx in eachindex(fiff_object)] .== f]
        if length(tmp) != 0
            push!(meas_info, Symbol(f)=>tmp[1][4])
        else
            push!(meas_info, Symbol(f)=>nothing)
        end
    end

    ch_info = Dict()
    block_idx = [fiff_object[idx][2] for idx in eachindex(fiff_object)] .== "ch_info"
    tmp = fiff_object[block_idx][[fiff_object[block_idx][idx][2] for idx in eachindex(fiff_object[block_idx])] .== "ch_info"]
    if length(tmp) == 1
        push!(ch_info, Symbol("ch_info")=>tmp[1][4])
    elseif length(tmp) > 1
        [push!(ch_info, Symbol("ch_info" * "_$n")=>tmp[n][4]) for n in eachindex(tmp)]
    else
        push!(ch_info, Symbol("ch_info")=>nothing)
    end
    push!(meas_info, :ch_info=>ch_info)

    subject_info = Dict()
    fields = ["subj_id", "subj_his_id", "subj_last_name", "subj_first_name", "subj_middle_name", "subj_birth_day", "subj_sex", "subj_hand", "subj_weight", "subj_height", "comment"]
    subject_info = _pack_fiff_blocks(fiff_object, "subject", fields)

    fields = ["dig_point", "hpi_digitization_order", "hpi_coils_used", "hpi_coil_moments", "hpi_fit_goodness", "hpi_fit_good_limit", "hpi_fit_dist_limit", "hpi_fit_accept", "coord_trans"]
    hpi_result = _pack_fiff_blocks(fiff_object, "hpi_result", fields)
    fields = ["hpi_coil_no", "epoch", "hpi_slopes", "hpi_corr_coeff", "hpi_coil_freq"]
    hpi_coil = _pack_fiff_blocks(fiff_object, "hpi_coil", fields)
    isotrak = _pack_fiff_blocks(fiff_object, "isotrak", ["dig_point"])
    hpi = Dict(:hpi_result=>hpi_result, :hpi_coil=>hpi_coil, :isotrak=>isotrak)

    hpi_coil = _pack_fiff_blocks(fiff_object, "hpi_coil", ["event_bits"])
    hpi_subsystem = _pack_fiff_blocks(fiff_object, "hpi_subsystem", ["hpi_ncoil", "event_channel"])
    hpi_subsystem = Dict(:hpi_subsystem=>hpi_subsystem, :hpi_coil=>hpi_coil)

    fields = ["description", "proj_item_kind", "nchan", "proj_item_ch_name_list", "proj_item_time", "proj_item_nvec", "proj_item_vectors"]
    xfit_proj_item = _pack_fiff_blocks(fiff_object, "xfit_proj_item", fields)
    nchan = _pack_fiff_blocks(fiff_object, "xfit_proj_item", ["nchan"])
    xfit_proj = Dict(:nchan=>nchan, :xfit_proj_item=>xfit_proj_item)

    fields = ["event_channels", "event_list"]
    events = _pack_fiff_blocks(fiff_object, "events", fields)

    fields = ["dacq_pars", "dacq_stim"]
    dacq_pars = _pack_fiff_blocks(fiff_object, "dacq_pars", fields)
    dacq_pars[:dacq_pars] = split(dacq_pars[:dacq_pars], "\n")
    dacq_pars[:dacq_pars][end] == "" && deleteat!(dacq_pars[:dacq_pars], length(dacq_pars[:dacq_pars]))
    dacq_pars[:dacq_pars] = split.(rstrip.(dacq_pars[:dacq_pars]), ' ')

    fields = ["experimenter", "description", "meas_date"]
    merge!(meas_info, _pack_fiff_blocks(fiff_object, "meas_info", fields))

    push!(meas_info, :subject_info=>subject_info)
    push!(meas_info, :project_info=>project_info)
    push!(meas_info, :hpi=>hpi)
    push!(meas_info, :ssp=>xfit_proj)
    push!(meas_info, :events=>events)
    push!(meas_info, :dacq_pars=>dacq_pars)
    push!(meas_info, :nchan=>length(ch_info))
    meas_info[:sfreq] = round(Int64, meas_info[:sfreq])

    fields = ["data_buffer"]
    raw_data_tmp1 = sort(collect(_pack_fiff_blocks(fiff_object, "raw_data", fields)))
    buffer_length = length(raw_data_tmp1[1][2])
    raw_data_tmp2 = Float64[]
    for idx in eachindex(raw_data_tmp1)
        append!(raw_data_tmp2, raw_data_tmp1[idx][2])
    end
    raw_data = Dict()
    push!(raw_data, :raw_data=>raw_data_tmp2)
    fields = ["first_samp", "data_skip", "data_skip_samp"]
    merge!(raw_data, _pack_fiff_blocks(fiff_object, "raw_data", fields))
    isnothing(raw_data[:data_skip]) && (raw_data[:data_skip] = 0)
    isnothing(raw_data[:data_skip_samp]) && (raw_data[:data_skip_samp] = 0)
    if isnothing(raw_data[:first_samp])
        raw_data[:first_samp] = 0
    else
        raw_data[:first_samp] /= meas_info[:sfreq]
    end
    if raw_data[:data_skip] > 0
        #raw_data[:raw_data] = vcat(raw_data[:raw_data], zeros(buffer_length * raw_data[:data_skip]))
        _warn("data_skip is not implemented yet; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org")
    elseif raw_data[:data_skip_samp] > 0
        #raw_data[:raw_data] = vcat(raw_data[:raw_data], zeros(raw_data[:data_skip_samp]))
        _warn("data_skip_samp is not implemented yet; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org")
    end

    fiff = Dict(:meas_info=>meas_info, :raw_data=>raw_data)

    return fiff, fiff_object, fiff_blocks

end

"""
    import_fiff(file_name)

Load Elekta-Neuromag 306 FIFF (Functional Image File Format) file (MEG, EEG) and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function import_fiff(file_name::String)::NeuroAnalyzer.NEURO

    @assert isfile(file_name) "File $file_name cannot be loaded."

    fiff, _, _ = load_fiff(file_name)

    sampling_rate = fiff[:meas_info][:sfreq]
    ch_n = fiff[:meas_info][:nchan]
    data = @views reshape(fiff[:raw_data][:raw_data], fiff[:meas_info][:nchan], :, 1)

    units = repeat([""], ch_n)
    ch_type = repeat([""], ch_n)
    coil_type = repeat([""], ch_n)
    clabels = repeat([""], ch_n)
    @inbounds for (k, v) in fiff[:meas_info][:ch_info]
        ch = v[1]
        ch_type[ch] = v[3]
        ch_type[ch] == "stim" && (ch_type[ch] = "mrk")
        ch_type[ch] == "ias" && (ch_type[ch] = "other")
        ch_type[ch] == "syst" && (ch_type[ch] = "other")
        range = v[4]
        cal = v[5]
        coil_type[ch] = _find_fiff_coiltype(v[6])
        units[ch] = v[19]
        unit_mul = 10^v[20]
        clabels[ch] = uppercase(v[3]) * " " * string(v[2])
        @views data[ch, :, 1] .*= (range * cal * unit_mul)
    end

    clabels = _clean_meg_labels(clabels)

    # convert data to standard units
    @inbounds for ch_idx in 1:ch_n
        if units[ch_idx] == "T"
            @views data[ch_idx, :, 1] .*= 10^15
            units[ch_idx] = "fT"
        elseif units[ch_idx] == "T/m"
            @views data[ch_idx, :, 1] .*= (10^15 / 100)
            units[ch_idx] = "fT/cm"
        elseif units[ch_idx] == "T/cm"
            @views data[ch_idx, :, 1] .*= 10^15
            units[ch_idx] = "fT/cm"
        elseif units[ch_idx] == "V" && ch_type[ch_idx] in ["eeg", "emg", "eog", "ref"]
            @views data[ch_idx, :, 1] .*= 10^6
            units[ch_idx] = "μV"
        end
    end

    # coil types
    magnetometers = Int64[]
    gradiometers = Int64[]
    eeg = Int64[]
    @inbounds for ch_idx in 1:ch_n
        if coil_type[ch_idx] in ["vv_planar_w", "vv_planar_t1", "vv_planar_t2", "vv_planar_t3"]
            coil_type[ch_idx] = "pgrad"
            ch_type[ch_idx] = "grad"
            push!(gradiometers, ch_idx)
        elseif coil_type[ch_idx] in ["magnes_grad"]
            coil_type[ch_idx] = "grad"
            ch_type[ch_idx] = "grad"
            push!(gradiometers, ch_idx)
        elseif coil_type[ch_idx] in ["axial_grad_5cm", "ctf_grad"]
            coil_type[ch_idx] = "agrad"
            ch_type[ch_idx] = "grad"
            push!(gradiometers, ch_idx)
        elseif coil_type[ch_idx] in ["point_magnetometer", "vv_mag_w", "vv_mag_t1", "vv_mag_t2", "vv_mag_t3", "magnes_mag"]
            coil_type[ch_idx] = "mag"
            push!(magnetometers, ch_idx)
            ch_type[ch_idx] = "mag"
        elseif coil_type[ch_idx] in ["eeg"]
            push!(eeg, ch_idx)
        end
    end

    bad_channels = zeros(Bool, ch_n, 1)
    !isnothing(fiff[:meas_info][:bad_chs]) && _warn("bad_channels tag is not implemented yet; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org")

    # HPI
    # fiff[:meas_info][:hpi]
    id = Int64[]
    p = Int64[]
    x = Float64[]
    y = Float64[]
    z = Float64[]
    for (k, v) in fiff[:meas_info][:hpi][:isotrak]
        push!(id, parse(Int64, match(r"(.+)_(\d+)", string(k)).captures[2]))
        push!(p, v[2])
        push!(x, v[3])
        push!(y, v[4])
        push!(z, v[5])
    end
    hpi = DataFrame(:id=>id, :p=>p, :x=>x, :y=>y, :z=>z)
    hpi = sort(hpi, :id)

    # markers
    # fiff[:meas_info][:events]
    events_ch = fiff[:meas_info][:events][:event_channels]
    # number of sample, before, after
    events = reshape(fiff[:meas_info][:events][:event_list], 3, :)'

    markers = DataFrame(:id=>String[],
                        :start=>Float64[],
                        :length=>Float64[],
                        :value=>String[],
                        :channel=>Int64[])

    # MaxShield

    # MRI

    # SSP
    ssp_labels = String[]
    ssp_channels = Bool[]
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
            # we assume that the content of all ch_name_list is the same; this might not be true for some recordings
            ssp_labels = string.(split(xfit_proj_item[:proj_item_ch_name_list_1], ':'))
            ssp_labels = _clean_meg_labels(ssp_labels)
            ssp_channels = zeros(Bool, ch_n)
            ssp_sorting_idx = Int64[]
            ssp_data = zeros(length(xfit_proj_vecs), length(xfit_proj_vecs[1]))
            @inbounds for idx in eachindex(xfit_proj_vecs)
                ssp_data[idx, :] = xfit_proj_vecs[idx]
            end
            for idx in eachindex(ssp_labels)
                push!(ssp_sorting_idx, findfirst(isequal(ssp_labels[idx]), clabels))
            end
            ssp_channels[ssp_sorting_idx] .= true
            ssp_labels = Vector{String}()
            for idx in eachindex(xfit_proj_vecs)
                push!(ssp_labels, "PCA-v$(lpad(idx, 2, '0'))")
            end
        end
    end

    # create signal details
    time_pts = round.(collect(0:1/sampling_rate:size(data, 2) * size(data, 3) / sampling_rate)[1:end-1], digits=3)
    epoch_time = round.((collect(0:1/sampling_rate:size(data, 2) / sampling_rate))[1:end-1], digits=3)

    file_size_mb = round(filesize(file_name) / 1024^2, digits=2)

    data_type = ""
    if occursin("meg", lowercase(fiff[:meas_info][:experimenter]))
        data_type = "meg"
    elseif occursin("eeg", lowercase(fiff[:meas_info][:experimenter]))
        data_type = "eeg"
    else
        _warn("Unknown data type: $(lowercase(fiff[:meas_info][:experimenter])).")
    end

    id = fiff[:meas_info][:subject_info][:subj_id]
    id = isnothing(id) ? "" : string(id)
    first_name = fiff[:meas_info][:subject_info][:subj_first_name]
    isnothing(first_name) && (first_name = "")
    middle_name = fiff[:meas_info][:subject_info][:subj_middle_name]
    isnothing(middle_name) && (middle_name = "")
    last_name = fiff[:meas_info][:subject_info][:subj_last_name]
    isnothing(last_name) && (last_name = "")
    handedness = fiff[:meas_info][:subject_info][:subj_hand]
    handedness = isnothing(handedness) ? "" : string(handedness)
    weight = fiff[:meas_info][:subject_info][:subj_weight]
    isnothing(weight) && (weight = -1)
    height = fiff[:meas_info][:subject_info][:subj_height]
    isnothing(height) && (height = -1)
    recording = fiff[:meas_info][:description] 
    isnothing(recording) && (recording = "")
    project_info = fiff[:meas_info][:project_info][:proj_name]
    isnothing(project_info) && (project_info = "")

    date = fiff[:meas_info][:meas_date]
    if isnothing(date)
        rec_d = ""
        rec_t = ""
    else
        rec_d = string(Dates.day(date)) * "-" * string(Dates.month(date)) * "-" * string(Dates.year(date))
        rec_t = string(Dates.hour(date)) * ":" * string(Dates.minute(date)) * ":" * string(Dates.second(date))
    end

    lp = fiff[:meas_info][:lowpass]
    if isnothing(lp)
        lp = 0
    else
        lp = round(lp, digits=1)
    end
    hp = fiff[:meas_info][:highpass]
    if isnothing(hp)
        hp = 0
    else
        hp = round(hp, digits=1)
    end
    lf = fiff[:meas_info][:line_freq]
    if isnothing(lf)
        lf = 0
    else
        lf = round(lf, digits=1)
    end

    s = _create_subject(id=id,
                        first_name=first_name,
                        middle_name=middle_name,
                        last_name=last_name,
                        head_circumference=-1,
                        handedness=handedness,
                        weight=weight,
                        height=height)
    r = _create_recording_meg(data_type=data_type,
                              file_name=file_name,
                              file_size_mb=file_size_mb,
                              file_type="FIFF",
                              recording=recording,
                              recording_date=rec_d,
                              recording_time=rec_t,
                              recording_notes="",
                              channel_type=ch_type,
                              channel_order=_sort_channels(ch_type),
                              reference="",
                              clabels=clabels,
                              units=units,
                              prefiltering=repeat(["LP: $lp Hz; HP: $hp Hz"], ch_n),
                              line_frequency=lf,
                              sampling_rate=sampling_rate,
                              magnetometers=magnetometers,
                              gradiometers=gradiometers,
                              coil_type=coil_type,
                              bad_channels=bad_channels,
                              ssp_labels=ssp_labels,
                              ssp_channels=ssp_channels,
                              ssp_data=ssp_data)
    e = _create_experiment(name=project_info, notes="", design="")
    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = String[]

    locs = _initialize_locs()
    obj = NeuroAnalyzer.NEURO(hdr, time_pts, epoch_time, data, components, markers, locs, history)
    _initialize_locs!(obj)
    l = import_locs_csv(joinpath(NeuroAnalyzer.res_path, "meg_306flattened.csv"))
    add_locs!(obj, locs=l)

    _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits=3)) s)")

    return obj

end
