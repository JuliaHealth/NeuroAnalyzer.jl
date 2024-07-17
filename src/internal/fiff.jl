_find_fiff_tag(t::String) = fiff_tags[:id][findfirst(isequal(t), fiff_tags[:tag])]
_find_fiff_tag(n::Int64) = fiff_tags[:tag][findfirst(isequal(n), fiff_tags[:id])]
_find_fiff_block(t::String) = fiff_blocks[:id][findfirst(isequal(t), fiff_blocks[:block])]
_find_fiff_block(n::Int64) = fiff_blocks[:block][findfirst(isequal(n), fiff_blocks[:id])]
_find_fiff_dt(n::Int64) = fiff_data_type[:name][findfirst(isequal(n & 0x00000FFF), fiff_data_type[:id])]
_find_fiff_unit(n::Int64) = fiff_units[:unit][findfirst(isequal(n), fiff_units[:id])]
_find_fiff_mul(n::Int64) = fiff_multipliers[:id][findfirst(isequal(n), fiff_multipliers[:id])]
_find_fiff_chtype(n::Int64) = fiff_channel_type[:channel_type][findfirst(isequal(n), fiff_channel_type[:id])]
_find_fiff_gantry_type(n::Int64) = fiff_gantry_type[:gantry_type][findfirst(isequal(n), fiff_gantry_type[:id])]
_find_fiff_dacq_system(n::Int64) = fiff_dacq_system[:dacq_system][findfirst(isequal(n), fiff_dacq_system[:id])]
_find_fiff_proj_item(n::Int64) = fiff_proj_item[:proj_item][findfirst(isequal(n), fiff_proj_item[:id])]
_find_fiff_proj_by(n::Int64) = fiff_proj_by[:proj_by][findfirst(isequal(n), fiff_proj_by[:id])]
_find_fiff_coiltype(n::Int64) = fiff_coil_type[:coil_type][findfirst(isequal(n), fiff_coil_type[:id])]

# data type
fiff_data_type= Dict(:id=>[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 11, 13, 14, 16, 20, 21, 23, 30, 31, 32, 33, 34, 35],
               :name=>["void", "byte", "int16", "int32", "float", "double", "julian", "uint16", "uint32", "uint64", "string", "ascii", "int64", "dau_pack13", "dau_pack14", "dau_pack16", "complex_float", "complex_double", "old_pack", "ch_info_struct", "id_struct", "dir_entry_struct", "dig_point_struct", "ch_pos_struct", "coord_trans_struct", "old_pack"],
               :data_type=>["void_t", "byte_t", "int16_t", "int32_t", "float_t", "double_t", "julian_t", "uint16_t", "uint32_t", "uint64_t", "byte_t", "byte_t", "int64_t", "dau_pack13_t", "dau_pack14_t", "dau_pack16_t", "complex_float_t", "complex_double_t", "old_pack_tvariable", "ch_info_t", "id_t", "dir_entry_t", "dig_point_t", "ch_pos_t", "coord_trans_t", "old_pack_tvariable"],
               :size=>[1, 1, 2, 4, 4, 8, 8, 2, 4, 8, 1, 1, 8, 2, 2, 2, 8, 16, 80, 20, 16, 20, 52, 80])
    
# units
fiff_units = Dict(:id=>[-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 201, 202],
                  :unit=>["none", "unitless", "m", "kg", "S", "A", "k", "mol", "rad", "sr", "cd", "Hz", "N", "Pa", "J", "W", "C", "V", "F", "Ω", "Mho", "Wb", "T", "H", "C", "lm", "lx", "T/m", "Am"])

# value multipliers
fiff_multipliers = Dict(:id=>[18, 15, 12, 9, 6, 3, 2, 1, 0, -1, -2, -3, -6, -9, -12, -15, -18],
                        :multiplier=>["e","pet","t","gig","meg","k","h","da","none","d","c","m","mu","n","p","f","a"])

# id tags
fiff_tags = Dict(:id=>[100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 120, 150, 151, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 233, 234, 235, 236, 237, 240, 241, 242, 243, 244, 245, 246, 247, 250, 251, 252, 253, 254, 255, 256, 257, 258, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 300, 301, 302, 303, 304, 305, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 500, 501, 502, 503, 504, 600, 601, 602, 603, 701, 702, 703, 800, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2030, 2031, 2032, 2040, 2041, 2100, 2200, 3000, 3001, 3002, 3003, 3101, 3102, 3103, 3104, 3105, 3106, 3107, 3108, 3109, 3110, 3111, 3112, 3113, 3201, 3300, 3401, 3402, 3403, 3404, 3405, 3406, 3407, 3408, 3411, 3412, 3413, 3414, 3415, 3416, 3417, 3501],
                 :tag=>["file_id", "dir_pointer", "dir", "block_id", "block_start", "block_end", "free_list", "free_block", "nop", "parent_file_id", "parent_block_id", "block_name", "block_version", "creator", "modifier", "ref_role", "ref_file_id", "ref_file_num", "ref_file_name", "ref_block_id", "dacq_pars", "dacq_stim", "nchan", "sfreq", "data_pack", "ch_info", "meas_date", "subject", "description", "nave", "first_sample", "last_sample", "aspect_kind", "ref_event", "experimenter", "dig_point", "ch_pos_vec", "hpi_slopes", "hpi_ncoil", "req_event", "req_limit", "lowpass", "bad_chs", "artef_removal", "coord_trans", "highpass", "ch_cals_vec", "hpi_bad_chs", "hpi_corr_coeff", "event_comment", "no_samples", "first_time", "subave_size", "subave_first", "name", "dig_string", "line_freq", "hpi_coil_freq", "signal_channel", "hpi_coil_moments", "hpi_fit_goodness", "hpi_fit_accept", "hpi_fit_good_limit", "hpi_fit_dist_limit", "hpi_coil_no", "hpi_coils_used", "hpi_digitization_order", "ch_scan_no", "ch_logical_no", "ch_kind", "ch_range", "ch_cal", "ch_pos", "ch_unit", "ch_unit_mul", "ch_dacq_name", "sss_frame", "sss_job", "sss_origin", "sss_ord_in", "sss_ord_out", "sss_nmag", "sss_components", "sss_cal_chans", "sss_cal_corrs", "sss_st_corr", "sss_base_in", "sss_base_out", "sss_base_virt", "sss_norm", "sss_iterate", "sss_nfree", "sss_st_length", "gantry_type", "gantry_model", "gantry_angle", "data_buffer", "data_skip", "epoch", "data_skip_samp", "data_buffer2", "time_stamp", "subj_id", "subj_first_name", "subj_middle_name", "subj_last_name", "subj_birth_day", "subj_sex", "subj_hand", "subj_weight", "subj_height", "subj_comment", "subj_his_id", "proj_id", "proj_name", "proj_aim", "proj_persons", "proj_comment", "event_channels", "event_list", "event_channel", "event_bits", "squid_bias", "squid_offset", "squid_gate", "decoupler_matrix", "volume_type", "mri_source_format", "mri_pixel_encoding", "mri_pixel_data_offset", "mri_pixel_scale", "mri_pixel_data", "mri_pixel_overlay_encoding", "mri_pixel_overlay_data", "mri_bounding_box", "mri_width", "mri_width_m", "mri_height", "mri_height_m", "mri_depth", "mri_depth_m", "mri_thickness", "mri_scene_aim", "mri_calibration_scale", "mri_calibration_offset", "mri_orig_source_path", "mri_orig_source_format", "mri_orig_pixel_encoding", "mri_orig_pixel_data_offset", "mri_time", "mri_voxel_data", "mri_voxel_encoding", "voxel_nchannels", "mri_diffusion_weight", "mri_diffusion_param", "mri_mrilab_setup", "mri_seg_region_id", "conductor_model_kind", "sphere_origin", "modelsphere_coord_frame", "sphere_layers", "bem_surf_id", "bem_surf_name", "bem_surf_nnode", "bem_surf_ntri", "bem_surf_nodes", "bem_surf_triangles", "bem_surf_normals", "bem_surf_curvs", "bem_surf_curv_values", "bem_pot_solution", "bem_approx", "bem_coord_frame", "bem_sigma", "source_dipole", "beamformer_instructions", "xfit_lead_products", "xfit_map_products", "xfit_grad_map_products", "xfit_vol_integration", "xfit_integration_radius", "xfit_conductor_model_name", "xfit_conductor_model_trans_name", "xfit_cont_surf_type", "proj_item_kind", "proj_item_time", "proj_item_ign_chs", "proj_item_nvec", "proj_item_vectors", "proj_item_definition", "proj_item_ch_name_list", "xplotter_layout"])

# block types
fiff_blocks = Dict(:id=>[999, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 114, 115, 116, 117, 118, 119, 120, 121, 122, 200, 201, 202, 203, 204, 205, 206, 300, 310, 311, 312, 313, 314, 315, 359, 400, 500, 501, 502, 503, 504, 505, 510, 900, 901],
                  :block=>["root", "meas", "meas_info", "raw_data", "processed_data", "evoked", "aspect", "subject", "isotrak", "hpi_meas", "hpi_result", "hpi_coil", "project", "continuous_data", "void", "events", "index", "dacq_pars", "ref", "maxshield_raw_data", "maxshield_aspect", "hpi_subsystem", "phantom_subsystem", "structural_data", "volume_data", "volume_slice", "scenery", "scene", "mri_seg", "mri_seg_region", "sphere", "bem", "bem_surf", "conductor_model", "xfit_proj", "xfit_proj_item", "xfit_aux", "bad_channels", "vol_info", "data_correction", "channels_decoupler", "sss_info", "sss_cal_adjust", "sss_st_info", "sss_bases", "maxshield", "processing_history", "processing_record"])

# channel type
fiff_channel_type = Dict(:id=>[1, 2, 3, 102, 201, 202, 301, 302, 402, 502, 602, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 900, 910, 920, 1000, 1001],
                         :channel_type=>["meg", "eeg", "stim", "bio", "mcg", "eog", "meg_ref", "emg", "ecg", "misc", "resp", "quat0", "quat1", "quat2", "quat3", "quat4", "quat5", "quat6", "hpi_goodness", "hpi_error", "hpi_movement", "syst", "ias", "exci", "dipole_wave", "goodness_fit"])

# coil types
fiff_coil_type = Dict(:id=>[0, 1, 2, 3, 4, 5, 200, 1000, 2000, 2001, 3011, 3012, 3013, 3014, 3021, 3022, 3023, 3024, 4001, 4002, 5001],
                      :coil_type=>["none", "eeg", "nm_122", "nm_24", "nm_mcg_axial", "eeg_bipolar", "dipole", "mcg_42", "point_magnetometer", "axial_grad_5cm", "vv_planar_w", "vv_planar_t1", "vv_planar_t2", "vv_planar_t3", "vv_mag_w", "vv_mag_t1", "vv_mag_t2", "vv_mag_t3", "magnes_mag", "magnes_grad", "ctf_grad"])

# gantry
fiff_gantry_type = Dict(:id=>[0, 1, 2],
                        :gantry_type=>["fixed", "uni_axial", "free"])

# acquisition system
fiff_dacq_system = Dict(:id=>[0, 1, 2, 3, 4],
                        :dacq_system=>["dau", "vxi", "rpu", "orion", "triux"])

# proj_item
fiff_proj_item = Dict(:id=>[0, 1, 2, 3, 4, 5, 10],
                      :proj_item=>["none", "field", "dip_fix", "dip_rot", "homog_grad", "homog_field", "eeg_avref"])

# proj_bt
fiff_proj_by = Dict(:id=>[0, 1],
                      :proj_by=>["complement", "space"])

function _read_fiff_tag(fid::IOStream)
    tag = reinterpret(Int32, read(fid, 4 * sizeof(Int32)))
    tag .= ntoh.(tag)
    tag_kind = tag[1]
    data_type = tag[2]
    data_size = tag[3]
    tag_next = tag[4]
    data = read(fid, data_size)
    tag_next > 0 && seek(fid, tag_next)
    return tag_kind, data_type, data_size, data, tag_next
end

function _get_fiff_block_type(fid::IOStream, tag::Tuple{Int64, Int64, Int64, Int64, Vector{UInt8}, Int64})
    seek(fid, tag[1] + 16)
    buf = zeros(UInt8, tag[4])
    readbytes!(fid, buf, tag[4])
    return reinterpret(Int32, reverse(buf))
end

function _create_fiff_block(fid::IOStream)

    # read tags
    seek(fid, 0)

    # tags: position in file, tag_id, data_type, data_size, next
    tags = Vector{Tuple{Int64, Int64, Int64, Int64, Vector{UInt8}, Int64}}()
    tag_next = nothing
    while tag_next != -1
        current_position = position(fid)
        tag_kind, tag_type, tag_size, data, tag_next = NeuroAnalyzer._read_fiff_tag(fid)
        push!(tags, (current_position, tag_kind, tag_type, tag_size, data, tag_next))
    end

    _info("Number of tags: $(length(tags))")

    seek(fid, 0)

    # create list of tag IDs
    tag_pos = [tags[idx][1] for idx in eachindex(tags)]
    tag_ids = [tags[idx][2] for idx in eachindex(tags)]
    tag_type = [tags[idx][3] for idx in eachindex(tags)]
    tag_size = [tags[idx][4] for idx in eachindex(tags)]

    # create block structure
    block_level = ones(Int64, length(tags))
    block_type = Vector{Int64}()
    block_type_current = 999
    for tag_idx in eachindex(tags)
        if tags[tag_idx][2] == _find_fiff_tag("block_start")
            block_level[tag_idx:end] .+= 1
            block_type_current = _get_fiff_block_type(fid, tags[tag_idx])[]
            push!(block_type, block_type_current)
        elseif tags[tag_idx][2] == _find_fiff_tag("block_end")
            push!(block_type, block_type_current)
            block_level[tag_idx:end] .-= 1
        else
            push!(block_type, block_type_current)
        end
    end
    d = [tags[idx][5] for idx in eachindex(tags)]
    return d, hcat(tag_pos, tag_ids, tag_type, tag_size, block_level, block_type)
end

function _view_fiff_block(fb::Matrix{Int64})
    indent = ""
    for tag_idx in 1:size(fb, 1)
        fb[tag_idx, 2] == 104 && (indent = repeat("  ", length(indent) + 1))
        println(indent * "$(fb[tag_idx, 2]) [$(fb[tag_idx, 5])]")
        fb[tag_idx, 2] == 105 && (indent = repeat("  ", length(indent) - 1))
    end
end

function _get_fiff_data(d::Vector{Any}, id::Int64)
    data = Any[]
    for idx in eachindex(d)
        if d[idx][2] == NeuroAnalyzer._find_fiff_tag(id) && length(d[idx]) == 4
            push!(data, d[idx][4])
        end
    end
    return data
end

function _get_blocks(b::Matrix{Int64})
    levels = unique(b[:, 5])
    bidx = Vector{Int64}[]
    for idx in eachindex(levels)
        push!(bidx, findall(isequal(levels[idx]), b[:, 5]))
    end
    btmp = Int64[]
    for idx in eachindex(bidx)
        push!(btmp, length(bidx[idx]))
    end
    btypes = Int64[]
    for idx in eachindex(btmp)
        push!(btypes, unique(b[:, 6][bidx[idx][:]])[1])
    end
    return bidx, btypes
end

function _pack_fiff_blocks(fiff_object::Vector{Any}, block::String, fields::Vector{String})
    block_idx = [fiff_object[idx][3] for idx in eachindex(fiff_object)] .== block
    d = Dict()
    for f in fields
        tmp = fiff_object[block_idx][[fiff_object[block_idx][idx][2] for idx in eachindex(fiff_object[block_idx])] .== f]
        if length(tmp) == 1
            push!(d, Symbol(f)=>tmp[1][4])
        elseif length(tmp) > 1
            for n in eachindex(tmp)
                push!(d, Symbol(f * "_$(lpad(n, length(string(length(tmp))), '0'))")=>tmp[n][4])
            end
        else
            push!(d, Symbol(f)=>nothing)
        end
    end
    return d
end