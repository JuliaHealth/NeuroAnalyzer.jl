export import_fiff

"""
    import_fiff(file_name; <keyword arguments>)

Load FIFF (Functional Image File Format) file and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `obj::NeuroAnalyzer.NEURO`

# Source

Elekta Neuromag: Functional Image File Format Description. FIFF version 1.3. March 2011
"""
function import_fiff(file_name::String; detect_type::Bool=true)

    @assert isfile(file_name) "File $file_name cannot be loaded."

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
        d = nothing
        if tag_type in [107, 108]
            # void
        elseif tag_type in [210]
            # aspect
            d = _find_fiff_aspect(Int64(ntoh.(reinterpret(Int32, buf[block_idx])[1])))
        elseif tag_type in [264]
            # sss_job
            d = _find_fiff_sss_job(Int64(ntoh.(reinterpret(Int32, buf[block_idx])[1])))
        elseif tag_type in [270, 271, 273, 274, 275, 302, 800, 3415]
            # sss_cal_chans
            # sss_cal_cors
            # sss_base_in
            # sss_base_out
            # sss_base_virt
            # epoch
            # decoupler_matrix
            # proj_item_vectors
            d = NeuroAnalyzer._fiff_matrix(fiff_blocks[block_idx, 3], buf[block_idx])
        elseif tag_type in [115]
            # role
            d = Int64(ntoh.(reinterpret(Int32, buf[block_idx])[1]))
            d = d == 1 ? "prev_file" : "next_file"
        elseif tag_type in [234]
            # dig_string
            kind = Int64(ntoh.(reinterpret(Int32, buf[block_idx][1:4])[1]))
            ident = Int64(ntoh.(reinterpret(Int32, buf[block_idx][5:8])[1]))
            np = Float64(ntoh.(reinterpret(Int32, buf[block_idx][9:12])[1]))
            rr = Float64[]
            [push!(rr, Float64(ntoh.(reinterpret(Float32, buf[block_idx][idx:(idx + 3)])[1]))) for idx in 13:4:length(buf[block_idx])]
            d = (kind, ident, np, rr)
        elseif tag_type in [255]
            # ch_pos
            coil_type = Int64(ntoh.(reinterpret(Int32, buf[block_idx][1:4])[1]))
            r0_1 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][5:8])[1]))
            r0_2 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][9:12])[1]))
            r0_3 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][13:16])[1]))
            ex_1 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][17:20])[1]))
            ex_2 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][21:24])[1]))
            ex_3 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][25:28])[1]))
            ey_1 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][29:32])[1]))
            ey_2 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][33:36])[1]))
            ey_3 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][37:40])[1]))
            ez_1 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][41:44])[1]))
            ez_2 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][45:48])[1]))
            ez_3 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][59:52])[1]))
            d = (coil_type, r0_1, r0_2, r0_3, ex_1, ex_2, ex_3, ey_1, ey_2, ey_3, ez_1, ez_2, ez_3)
        elseif tag_type in [404, 405, 406]
            # dob, sex, handedness
            d = Int64(ntoh.(reinterpret(Int32, buf[block_idx])[1]))
            tag_type == 404 && (d = unix2datetime(d))
        elseif tag_type in [242]
            # hpi_mask
            d = ntoh.(reinterpret(UInt32, buf[block_idx])[1])
        elseif tag_type in [222]
            # coord_trans
            from = Int64(ntoh.(reinterpret(Int32, buf[block_idx][1:4])[1]))
            to = Int64(ntoh.(reinterpret(Int32, buf[block_idx][5:8])[1]))
            rot = zeros(3, 3)
            for idx1 in 1:3
                d = Float64[]
                [push!(d, Float64(ntoh.(reinterpret(Int16, buf[block_idx][(8 + idx2):(8 + idx2 + 3)])[1]))) for idx2 in 1:4:9]
                rot[idx1, :] = d
            end
            move = Float64[]
            [push!(move, Float64(ntoh.(reinterpret(Int16, buf[block_idx][(44 + idx):(44 + idx + 3)])[1]))) for idx in 1:4:9]
            invrot = zeros(3, 3)
            for idx1 in 1:3
                d = Float64[]
                [push!(d, Float64(ntoh.(reinterpret(Int16, buf[block_idx][(56 + idx2):(56 + idx2 + 3)])[1]))) for idx2 in 1:4:9]
                invrot[idx1, :] = d
            end
            invmove = Float64[]
            [push!(invmove, Float64(ntoh.(reinterpret(Int16, buf[block_idx][(92 + idx):(92 + idx + 3)])[1]))) for idx in 1:4:9]
            d = (from, to, rot, move, invrot, invmove)
        elseif tag_type in [300]
            # data_buffer
            df = _find_fiff_dt(fiff_blocks[block_idx, 3])
            d = Float64[]
            if df == "dau_pack16" || df == "int16"
                [push!(d, Float64(ntoh.(reinterpret(Int16, buf[block_idx][idx:(idx + 1)])[1]))) for idx in 1:2:length(buf[block_idx])]
            elseif df == "int32"
                [push!(d, Float64(ntoh.(reinterpret(Int32, buf[block_idx][idx:(idx + 3)])[1]))) for idx in 1:4:length(buf[block_idx])]
            elseif df == "float"
                [push!(d, Float64(ntoh.(reinterpret(Float32, buf[block_idx][idx:(idx + 3)])[1]))) for idx in 1:4:length(buf[block_idx])]
            else
                _warn("Data type $df is not implemented yet; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org")
            end
        elseif tag_type in [213]
            # dig_point
            kind = Int64(ntoh.(reinterpret(Int32, buf[block_idx][1:4])[1]))
            ident = Int64(ntoh.(reinterpret(Int32, buf[block_idx][5:8])[1]))
            r_1 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][9:12])[1]))
            r_2 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][13:16])[1]))
            r_3 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][17:20])[1]))
            d = (kind, ident, r_1, r_2, r_3)
        elseif tag_type in [252]
            # ch_type
            d =  _find_fiff_chtype(Int64(ntoh.(reinterpret(Int32, buf[block_idx])[1])))
        elseif tag_type in [256]
            # ch_unit
            d = _find_fiff_unit(Int64(ntoh.(reinterpret(Int32, buf[block_idx])[1])))
        elseif tag_type in [3411]
            # proj_item
            d = _find_fiff_proj_item(Int64(ntoh.(reinterpret(Int32, buf[block_idx])[1])))
        elseif tag_type in [3416]
            # proj_item
            d = _find_fiff_proj_by(Int64(ntoh.(reinterpret(Int32, buf[block_idx])[1])))
        elseif tag_type in [280]
            # gantry_type
            d = _find_fiff_gantry_type(Int64(ntoh.(reinterpret(Int32, buf[block_idx])[1])))
        elseif tag_type in [203]
            # ch_info
            scan_no = Int64(ntoh.(reinterpret(Int32, buf[block_idx][1:4])[1]))
            log_no = Int64(ntoh.(reinterpret(Int32, buf[block_idx][5:8])[1]))
            kind = _find_fiff_chtype(Int64(ntoh.(reinterpret(Int32, buf[block_idx][9:12])[1])))
            range = Float64(ntoh.(reinterpret(Float32, buf[block_idx][13:16])[1]))
            cal = Float64(ntoh.(reinterpret(Float32, buf[block_idx][17:20])[1]))
            coil_type = Int64(ntoh.(reinterpret(Int32, buf[block_idx][21:24])[1]))
            r0_1 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][25:28])[1]))
            r0_2 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][29:32])[1]))
            r0_3 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][33:36])[1]))
            ex_1 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][37:40])[1]))
            ex_2 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][41:44])[1]))
            ex_3 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][45:48])[1]))
            ey_1 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][49:52])[1]))
            ey_2 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][53:56])[1]))
            ey_3 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][57:60])[1]))
            ez_1 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][61:64])[1]))
            ez_2 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][65:68])[1]))
            ez_3 = Float64(ntoh.(reinterpret(Float32, buf[block_idx][69:72])[1]))
            unit = _find_fiff_unit(Int64(ntoh.(reinterpret(Int32, buf[block_idx][73:76])[1])))
            unit_mul = _find_fiff_mul(Int64(ntoh.(reinterpret(Int32, buf[block_idx][77:80])[1])))
            d = (scan_no, log_no, kind, range, cal, coil_type, r0_1, r0_2, r0_3, ex_1, ex_2, ex_3, ey_1, ey_2, ey_3, ez_1, ez_2, ez_3, unit, unit_mul)
        elseif tag_type in [214]
            # ch_pos_vec
            # obsolete
        elseif tag_type in [111, 112, 113, 114, 118, 150, 151, 205, 206, 212, 227, 233, 237, 258, 281, 401, 402, 403, 409, 410, 501, 502, 503, 504, 602, 2020, 3102, 3300, 3406, 3407, 3417, 3501]
            # string
            d = _v2s(string.(Char.(buf[block_idx])))
        elseif tag_type in [282]
            # int8
            d = ntoh.(reinterpret(Int8, buf[block_idx][1:4])[1])
        elseif tag_type in [101, 102, 104, 105, 106, 117, 200, 202, 204, 207, 208, 209, 211, 216, 217, 221, 228, 230, 231, 245, 250, 251, 257, 263, 266, 267, 268, 277, 278, 301, 303, 400, 500, 701, 702, 703, 2004, 2010, 2012, 2014, 2023, 2032, 3013, 3104, 3414]
            # int32
            d = Int64(ntoh.(reinterpret(Int32, buf[block_idx][1:4])[1]))
            if tag_type == 204
                d = unix2datetime(d)
            end
        elseif tag_type in [220, 232, 225, 246, 247, 269, 304, 305, 600, 601, 603, 3413]
            # int32*
            d = Int64[]
            [push!(d, Int64(ntoh.(reinterpret(Int32, buf[block_idx][idx:(idx + 3)])[1]))) for idx in 1:4:length(buf[block_idx])]
        elseif tag_type in [276]
            # double
            d = Float64(ntoh.(reinterpret(Float64, buf[block_idx][1:4])[1]))
        elseif tag_type in [201, 218, 219, 223, 229, 235, 236, 240, 241, 243, 244, 253, 254, 272, 279, 407, 408, 2005, 2009, 2011, 2013, 2015, 2016, 2018, 2019, 2024, 2040, 3109, 3113, 3405, 3412]
            # float
            d = Float64(ntoh.(reinterpret(Float32, buf[block_idx][1:4])[1]))
        elseif tag_type in [215, 224, 226, 265]
            # float*
            d = Float64[]
            [push!(d, Float64(ntoh.(reinterpret(Float32, buf[block_idx][idx:(idx + 3)])[1]))) for idx in 1:4:length(buf[block_idx])]
        elseif tag_type in [100, 103, 109, 110, 116, 120]
            # id_t
            fiff_v_major = Int64(ntoh.(reinterpret(Int16, buf[block_idx][1:2])[1]))
            fiff_v_minor = Int64(ntoh.(reinterpret(Int16, buf[block_idx][3:4])[1]))
            mach_id1 = Int64(ntoh.(reinterpret(Int32, buf[block_idx][5:8])[1]))
            mach_id2 = Int64(ntoh.(reinterpret(Int32, buf[block_idx][9:12])[1]))
            time_sec = Int64(ntoh.(reinterpret(Int32, buf[block_idx][13:16])[1]))
            id_creation_date = unix2datetime(time_sec)
            time_usec = Int64(ntoh.(reinterpret(Int32, buf[block_idx][17:20])[1]))
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
    for f in fields
        tmp = fiff_object[[fiff_object[idx][2] for idx in eachindex(fiff_object)] .== f]
        if length(tmp) != 0
            push!(project_info, Symbol(f)=>tmp[1][4])
        else
            push!(project_info, Symbol(f)=>nothing)
        end
    end

    ref = Dict()
    fields = ["ref_role", "ref_file_id", "ref_file_num", "ref_file_name", "ref_block_id"]
    for f in fields
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
    for f in fields
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

    raw_data[:raw_data] = reshape(raw_data[:raw_data], meas_info[:nchan], :, 1)
    
    ch_n = meas_info[:nchan]

    units = repeat([""], ch_n)
    ch_type = repeat([""], ch_n)
    coil_type = repeat([""], ch_n)
    clabels = repeat([""], ch_n)
    for (k, v) in meas_info[:ch_info]
        ch = v[1]
        ch_type[ch] = v[3]
        range = v[4]
        cal = v[5]
        coil_type[ch] = _find_fiff_coiltype(v[6])
        units[ch] = v[19]
        unit_mul = 10^v[20]
        clabels[ch] = uppercase(v[3]) * " " * string(v[2])
        raw_data[:raw_data][ch, :] .*= (range * cal * unit_mul)
    end

    # convert data to standard units
    for ch_idx in 1:ch_n
        if units[ch_idx] == "T"
            raw_data[:raw_data][ch_idx, :, :] .*= 10^15
            units[ch_idx] = "fT"
        end
        if units[ch_idx] == "T/m"
            raw_data[:raw_data][ch_idx, :, :] .*= (10^15 / 100)
            units[ch_idx] = "fT/cm"
        end
        if units[ch_idx] == "V" && ch_type[ch_idx] in ["eeg", "emg", "eog", "ref"]
            raw_data[:raw_data][ch_idx, :, :] .*= 10^6
            units[ch_idx] = "μV"
        end 
    end

    # coil types
    magnetometers = Int64[]
    gradiometers = Int64[]
    eeg = Int64[]
    for ch_idx in 1:ch_n
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

    fiff = Dict(:meas_info=>meas_info, :raw_data=>raw_data)

    data = fiff[:raw_data][:raw_data]
    sampling_rate = fiff[:meas_info][:sfreq]
    bad_channels = zeros(Bool, ch_n, 1)
    !isnothing(fiff[:meas_info][:bad_chs]) && _warn("bad_channels tag is not implemented yet; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org")

    # locs
    # SSS, PCA
    # events
    # HPI
    # MaxShield

    # create signal details
    time_pts = round.(collect(0:1/sampling_rate:size(data, 2) * size(data, 3) / sampling_rate)[1:end-1], digits=3)
    epoch_time = round.((collect(0:1/sampling_rate:size(data, 2) / sampling_rate))[1:end-1], digits=3)

    file_size_mb = round(filesize(file_name) / 1024^2, digits=2)

    s = _create_subject(id="",
                        first_name="",
                        middle_name="",
                        last_name="",
                        head_circumference=-1,
                        handedness="",
                        weight=-1,
                        height=-1)
    r = _create_recording_meg(data_type=lowercase(fiff[:meas_info][:experimenter]),
                              file_name=file_name,
                              file_size_mb=file_size_mb,
                              file_type="FIFF",
                              recording="",
                              recording_date=isnothing(fiff[:meas_info][:meas_date]) ? "" : string(Dates.day(date)) * "-" * string(Dates.month(date)) * "-" * string(Dates.year(date)),
                              recording_time=isnothing(fiff[:meas_info][:meas_date]) ? "" : string(Dates.hour(date)) * ":" * string(Dates.minute(date)) * ":" * string(Dates.second(date)),
                              recording_notes="",
                              channel_type=ch_type,
                              channel_order=_sort_channels(ch_type),
                              reference="",
                              clabels=clabels,
                              units=units,
                              prefiltering=repeat(["LP: $(fiff[:meas_info][:lowpass]) Hz; HP: $(fiff[:meas_info][:highpass]) Hz"], ch_n),
                              line_frequency=round(Int64, fiff[:meas_info][:line_freq]),
                              sampling_rate=sampling_rate,
                              magnetometers=magnetometers,
                              gradiometers=gradiometers,
                              coil_type=coil_type,
                              bad_channels=bad_channels)
    e = _create_experiment(name="", notes="", design="")

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = String[]

    markers = DataFrame()

    locs = _initialize_locs()

    obj = NeuroAnalyzer.NEURO(hdr, time_pts, epoch_time, data, components, markers, locs, history)

    _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits=2)) s)")

    return obj, fiff, fiff_object

end
