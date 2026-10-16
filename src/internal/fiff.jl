# ---------------------------------------------------------------------------
# FIFF (Elekta/MEGIN Neuromag Functional Image File Format) low-level reader
#
# All FIFF data are big-endian. Tag header = 4 × Int32: kind, type, size, next.
# Tag payloads are decoded by their DATA TYPE (low 16 bits of `type`), not by
# their kind; the high 16 bits encode matrix storage (dense / sparse CCS / RCS).
#
# Reference implementation used for cross-checking: MNE-Python (mne/_fiff/).
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# byte-reinterpretation helpers (scalar)
# ---------------------------------------------------------------------------
_i16i64(x)::Int64 = Int64(ntoh(reinterpret(Int16, x)[1]))
_i32i64(x)::Int64 = Int64(ntoh(reinterpret(Int32, x)[1]))
_f16f64(x)::Float64 = Float64(ntoh(reinterpret(Float16, x)[1]))
_f32f64(x)::Float64 = Float64(ntoh(reinterpret(Float32, x)[1]))
_f64f64(x)::Float64 = Float64(ntoh(reinterpret(Float64, x)[1]))
_i16f64(x)::Float64 = Float64(ntoh(reinterpret(Int16, x)[1]))
_i32f64(x)::Float64 = Float64(ntoh(reinterpret(Int32, x)[1]))
_i8i8(x)::Int8 = ntoh(reinterpret(Int8, x)[1])
_ui32i32(x)::Int32 = ntoh(reinterpret(UInt32, x)[1])
_i32i32(x)::Int32 = ntoh(reinterpret(Int32, x)[1])

# vectorized big-endian → host conversion; `b` must be a contiguous byte vector
# whose length is a multiple of sizeof(T)
_be(::Type{T}, b::AbstractVector{UInt8}) where {T} = ntoh.(reinterpret(T, Vector{UInt8}(b)))

# FIFF strings are not guaranteed to be NUL-terminated; cut at the first NUL
function _fiff_string(b::AbstractVector{UInt8})::String
    i = findfirst(iszero, b)
    return String(Vector{UInt8}(isnothing(i) ? b : b[1:(i - 1)]))
end

# ---------------------------------------------------------------------------
# FIFF lookup helpers
# Integer → name lookups never throw: unknown codes (newer systems, vendor
# extensions) are returned as "unknown_<code>" so that a single unrecognised
# tag cannot abort reading of the whole file.
# ---------------------------------------------------------------------------
function _fiff_lookup(tbl::Dict, key::Symbol, val::Symbol, k)
    idx = findfirst(isequal(k), tbl[key])
    return isnothing(idx) ? nothing : tbl[val][idx]
end
function _fiff_name(tbl::Dict, key::Symbol, val::Symbol, n::Integer)::String
    v = _fiff_lookup(tbl, key, val, Int64(n))
    return isnothing(v) ? "unknown_$n" : v
end
function _fiff_id(tbl::Dict, key::Symbol, val::Symbol, t::String)::Int64
    v = _fiff_lookup(tbl, val, key, t)
    isnothing(v) && throw(ArgumentError("Unknown FIFF identifier: $t."))
    return v
end

_find_fiff_tag(t::String)::Int64 = _fiff_id(fiff_tags, :id, :tag, t)
_find_fiff_tag(n::Integer)::String = _fiff_name(fiff_tags, :id, :tag, n)
_find_fiff_block(t::String)::Int64 = _fiff_id(fiff_blocks, :id, :block, t)
_find_fiff_block(n::Integer)::String = _fiff_name(fiff_blocks, :id, :block, n)
_find_fiff_dt(n::Integer)::String = _fiff_name(fiff_data_type, :id, :name, n & 0x0000FFFF)
_find_fiff_unit(n::Integer)::String = _fiff_name(fiff_units, :id, :unit, n)
_find_fiff_mul(n::Integer)::String = _fiff_name(fiff_multipliers, :id, :multiplier, n)
_find_fiff_chtype(n::Integer)::String =
    _fiff_name(fiff_channel_type, :id, :channel_type, n)
_find_fiff_gantry_type(n::Integer)::String =
    _fiff_name(fiff_gantry_type, :id, :gantry_type, n)
_find_fiff_dacq_system(n::Integer)::String =
    _fiff_name(fiff_dacq_system, :id, :dacq_system, n)
_find_fiff_proj_item(n::Integer)::String = _fiff_name(fiff_proj_item, :id, :proj_item, n)
_find_fiff_proj_by(n::Integer)::String = _fiff_name(fiff_proj_by, :id, :proj_by, n)
_find_fiff_coiltype(n::Integer)::String = _fiff_name(fiff_coil_type, :id, :coil_type, n)
_find_fiff_aspect(n::Integer)::String = _fiff_name(fiff_aspect, :id, :aspect, n)
_find_fiff_sss_job(n::Integer)::String = _fiff_name(fiff_sss_job, :id, :sss_job, n)

# build a two-column lookup Dict from (id => name) pairs - keeps ids and names
# aligned by construction
_fiff_table(pairs, key::Symbol, val::Symbol) =
    Dict(key => [first(p) for p in pairs], val => [last(p) for p in pairs])

# ---------------------------------------------------------------------------
# data types (FIFFT_*)
# ---------------------------------------------------------------------------
fiff_data_type = Dict(
    :id => [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 16, 20, 21, 23, 30, 31, 32, 33, 34, 35, 36, 37],
    :name => [
        "void", "byte", "int16", "int32", "float", "double", "julian",
        "uint16", "uint32", "uint64", "string", "int64",
        "dau_pack13", "dau_pack14", "dau_pack16", "complex_float",
        "complex_double", "old_pack", "ch_info_struct", "id_struct",
        "dir_entry_struct", "dig_point_struct", "ch_pos_struct",
        "coord_trans_struct", "dig_string_struct", "stream_segment_struct",
    ],
    :size => [0, 1, 2, 4, 4, 8, 4, 2, 4, 8, 1, 8, 2, 2, 2, 8, 16, 0, 96, 20, 16, 20, 52, 104, 0, 0],
)

const _FIFFT_MATRIX = 0x40000000
const _FIFFT_SPARSE_CCS = 0x00100000
const _FIFFT_SPARSE_RCS = 0x00200000

# ---------------------------------------------------------------------------
# physical units (FIFF_UNIT_*)
# ---------------------------------------------------------------------------
fiff_units = _fiff_table(
    [
        -1 => "none", 0 => "unitless", 1 => "m", 2 => "kg", 3 => "s", 4 => "A",
        5 => "K", 6 => "mol", 7 => "rad", 8 => "sr", 9 => "cd", 10 => "mol/m³",
        101 => "Hz", 102 => "N", 103 => "Pa", 104 => "J", 105 => "W", 106 => "C",
        107 => "V", 108 => "F", 109 => "Ω", 110 => "S", 111 => "Wb", 112 => "T",
        113 => "H", 114 => "°C", 115 => "lm", 116 => "lx", 117 => "V/m²",
        118 => "s²", 201 => "T/m", 202 => "Am", 203 => "Am/m²", 204 => "Am/m³",
        210 => "px",
    ],
    :id,
    :unit,
)

# value multipliers (FIFF_UNITM_*), stored as power-of-10 exponents
fiff_multipliers = _fiff_table(
    [
        18 => "exa", 15 => "peta", 12 => "tera", 9 => "giga", 6 => "mega",
        3 => "k", 2 => "h", 1 => "da", 0 => "none", -1 => "d", -2 => "c",
        -3 => "m", -6 => "μ", -9 => "n", -12 => "p", -15 => "f", -18 => "a",
    ],
    :id,
    :multiplier,
)

# ---------------------------------------------------------------------------
# tag kinds (FIFF_*) - MNE-Python _fiff/constants.py
# ---------------------------------------------------------------------------
fiff_tags = _fiff_table(
    [
        # file structure
        100 => "file_id", 101 => "dir_pointer", 102 => "dir", 103 => "block_id",
        104 => "block_start", 105 => "block_end", 106 => "free_list",
        107 => "free_block", 108 => "nop", 109 => "parent_file_id",
        110 => "parent_block_id", 111 => "block_name", 112 => "block_version",
        113 => "creator", 114 => "modifier", 115 => "ref_role",
        116 => "ref_file_id", 117 => "ref_file_num", 118 => "ref_file_name",
        120 => "ref_block_id",
        150 => "dacq_pars", 151 => "dacq_stim",
        # device / helium (TRIUX)
        152 => "device_type", 153 => "device_model", 154 => "device_serial",
        155 => "device_site", 156 => "he_level_raw", 157 => "helium_level",
        158 => "orig_file_guid", 159 => "utc_offset",
        # measurement info
        200 => "nchan", 201 => "sfreq", 202 => "data_pack", 203 => "ch_info",
        204 => "meas_date", 205 => "subject", 206 => "description", 207 => "nave",
        208 => "first_sample", 209 => "last_sample", 210 => "aspect_kind",
        211 => "ref_event", 212 => "experimenter", 213 => "dig_point",
        214 => "ch_pos_vec", 215 => "hpi_slopes", 216 => "hpi_ncoil",
        217 => "req_event", 218 => "req_limit", 219 => "lowpass", 220 => "bad_chs",
        221 => "artef_removal", 222 => "coord_trans", 223 => "highpass",
        224 => "ch_cals_vec", 225 => "hpi_bad_chs", 226 => "hpi_corr_coeff",
        227 => "event_comment", 228 => "no_samples", 229 => "first_time",
        230 => "subave_size", 231 => "subave_first", 233 => "name",
        234 => "dig_string", 235 => "line_freq", 236 => "hpi_coil_freq",
        237 => "signal_channel", 240 => "hpi_coil_moments",
        241 => "hpi_fit_goodness", 242 => "hpi_fit_accept",
        243 => "hpi_fit_good_limit", 244 => "hpi_fit_dist_limit",
        245 => "hpi_coil_no", 246 => "hpi_coils_used",
        247 => "hpi_digitization_order",
        # channel info (single-field form)
        250 => "ch_scan_no", 251 => "ch_logical_no", 252 => "ch_kind",
        253 => "ch_range", 254 => "ch_cal", 255 => "ch_pos", 256 => "ch_unit",
        257 => "ch_unit_mul", 258 => "ch_dacq_name",
        # SSS
        263 => "sss_frame", 264 => "sss_job", 265 => "sss_origin",
        266 => "sss_ord_in", 267 => "sss_ord_out", 268 => "sss_nmag",
        269 => "sss_components", 270 => "sss_cal_chans", 271 => "sss_cal_corrs",
        272 => "sss_st_corr", 273 => "sss_base_in", 274 => "sss_base_out",
        275 => "sss_base_virt", 276 => "sss_norm", 277 => "sss_iterate",
        278 => "sss_nfree", 279 => "sss_st_length",
        # gantry
        280 => "gantry_type", 281 => "gantry_model", 282 => "gantry_angle",
        # data
        300 => "data_buffer", 301 => "data_skip", 302 => "epoch",
        303 => "data_skip_samp", 304 => "data_buffer2", 305 => "time_stamp",
        # extended channel info
        350 => "ch_coil_type", 351 => "ch_coord_frame",
        # subject
        400 => "subj_id", 401 => "subj_first_name", 402 => "subj_middle_name",
        403 => "subj_last_name", 404 => "subj_birth_day", 405 => "subj_sex",
        406 => "subj_hand", 407 => "subj_weight", 408 => "subj_height",
        409 => "subj_comment", 410 => "subj_his_id",
        # project
        500 => "proj_id", 501 => "proj_name", 502 => "proj_aim",
        503 => "proj_persons", 504 => "proj_comment",
        # events
        600 => "event_channels", 601 => "event_list", 602 => "event_channel",
        603 => "event_bits",
        # SQUID
        701 => "squid_bias", 702 => "squid_offset", 703 => "squid_gate",
        800 => "decoupler_matrix",
        1101 => "ref_path",
        # MRI
        2001 => "volume_type", 2002 => "mri_source_format",
        2003 => "mri_pixel_encoding", 2004 => "mri_pixel_data_offset",
        2005 => "mri_pixel_scale", 2006 => "mri_pixel_data",
        2007 => "mri_pixel_overlay_encoding", 2008 => "mri_pixel_overlay_data",
        2009 => "mri_bounding_box", 2010 => "mri_width", 2011 => "mri_width_m",
        2012 => "mri_height", 2013 => "mri_height_m", 2014 => "mri_depth",
        2015 => "mri_depth_m", 2016 => "mri_thickness", 2017 => "mri_scene_aim",
        2018 => "mri_calibration_scale", 2019 => "mri_calibration_offset",
        2020 => "mri_orig_source_path", 2021 => "mri_orig_source_format",
        2022 => "mri_orig_pixel_encoding", 2023 => "mri_orig_pixel_data_offset",
        2024 => "mri_time", 2030 => "mri_voxel_data", 2031 => "mri_voxel_encoding",
        2032 => "voxel_nchannels", 2040 => "mri_diffusion_weight",
        2041 => "mri_diffusion_param", 2100 => "mri_mrilab_setup",
        2200 => "mri_seg_region_id",
        # BEM / sphere
        3000 => "conductor_model_kind", 3001 => "sphere_origin",
        3002 => "sphere_radius", 3003 => "sphere_layers",
        3101 => "bem_surf_id", 3102 => "bem_surf_name", 3103 => "bem_surf_nnode",
        3104 => "bem_surf_ntri", 3105 => "bem_surf_nodes",
        3106 => "bem_surf_triangles", 3107 => "bem_surf_normals",
        3108 => "bem_surf_curvs", 3109 => "bem_surf_curv_values",
        3110 => "bem_pot_solution", 3111 => "bem_approx",
        3112 => "bem_coord_frame", 3113 => "bem_sigma", 3201 => "source_dipole",
        # xfit / projectors
        3300 => "beamformer_instructions", 3401 => "xfit_lead_products",
        3402 => "xfit_map_products", 3403 => "xfit_grad_map_products",
        3404 => "xfit_vol_integration", 3405 => "xfit_integration_radius",
        3406 => "xfit_conductor_model_name",
        3407 => "xfit_conductor_model_trans_name", 3408 => "xfit_cont_surf_type",
        3411 => "proj_item_kind", 3412 => "proj_item_time",
        3413 => "proj_item_ign_chs", 3414 => "proj_item_nvec",
        3415 => "proj_item_vectors", 3416 => "proj_item_definition",
        3417 => "proj_item_ch_name_list", 3501 => "xplotter_layout",
        # MNE extensions (written by MNE-Python / MNE-C)
        3502 => "mne_row_names", 3503 => "mne_col_names", 3504 => "mne_nrow",
        3505 => "mne_ncol", 3506 => "mne_coord_frame", 3507 => "mne_ch_name_list",
        3508 => "mne_file_name", 3552 => "mne_external_big_endian",
        3553 => "mne_external_little_endian", 3560 => "mne_proj_item_active",
        3561 => "mne_event_list", 3562 => "mne_hemi", 3563 => "mne_data_skip_nop",
        3564 => "mne_orig_ch_info", 3565 => "mne_event_trigger_mask",
        3566 => "mne_event_comments", 3567 => "mne_custom_ref",
        3568 => "mne_baseline_min", 3569 => "mne_baseline_max",
        3580 => "mne_ctf_comp_kind", 3581 => "mne_ctf_comp_data",
        3582 => "mne_ctf_comp_calibrated", 3585 => "mne_derivation_data",
        3612 => "mne_kit_system_id", 3800 => "mne_epochs_selection",
        3801 => "mne_epochs_drop_log", 3802 => "mne_epochs_reject_flat",
        3803 => "mne_epochs_raw_sfreq",
    ],
    :id,
    :tag,
)

# ---------------------------------------------------------------------------
# block kinds (FIFFB_*)
# ---------------------------------------------------------------------------
fiff_blocks = _fiff_table(
    [
        999 => "root",
        100 => "meas", 101 => "meas_info", 102 => "raw_data",
        103 => "processed_data", 104 => "evoked", 105 => "aspect",
        106 => "subject", 107 => "isotrak", 108 => "hpi_meas",
        109 => "hpi_result", 110 => "hpi_coil", 111 => "project",
        112 => "continuous_data", 113 => "ch_info", 114 => "void",
        115 => "events", 116 => "index", 117 => "dacq_pars", 118 => "ref",
        119 => "maxshield_raw_data", 120 => "maxshield_aspect",
        121 => "hpi_subsystem", 122 => "phantom_subsystem",
        124 => "device", 125 => "helium", 126 => "channel_info",
        200 => "structural_data", 201 => "volume_data", 202 => "volume_slice",
        203 => "scenery", 204 => "scene", 205 => "mri_seg",
        206 => "mri_seg_region",
        300 => "sphere", 310 => "bem", 311 => "bem_surf",
        312 => "conductor_model", 313 => "xfit_proj", 314 => "xfit_proj_item",
        315 => "xfit_aux",
        350 => "mne", 351 => "mne_source_space", 352 => "mne_forward_solution",
        353 => "mne_parent_mri_file", 354 => "mne_parent_meas_file",
        355 => "mne_cov", 356 => "mne_inverse_solution",
        357 => "mne_named_matrix", 358 => "mne_env", 359 => "bad_channels",
        360 => "mne_vertex_map", 361 => "mne_events", 362 => "mne_morph_map",
        363 => "mne_surface_map", 364 => "mne_surface_map_group",
        370 => "mne_ctf_comp", 371 => "mne_ctf_comp_data",
        372 => "mne_derivations", 373 => "mne_epochs", 374 => "mne_ica",
        400 => "vol_info",
        500 => "data_correction", 501 => "channels_decoupler",
        502 => "sss_info", 503 => "sss_cal_adjust", 504 => "sss_st_info",
        505 => "sss_bases", 510 => "maxshield",
        900 => "processing_history", 901 => "processing_record",
        3810 => "mne_annotations", 3811 => "mne_metadata",
    ],
    :id,
    :block,
)

# ---------------------------------------------------------------------------
# channel kinds (FIFFV_*_CH)
# ---------------------------------------------------------------------------
fiff_channel_type = _fiff_table(
    [
        1 => "meg", 2 => "eeg", 3 => "stim",
        102 => "bio", 201 => "mcg", 202 => "eog", 301 => "meg_ref", 302 => "emg",
        402 => "ecg", 502 => "misc", 602 => "resp",
        700 => "quat0", 701 => "quat1", 702 => "quat2", 703 => "quat3",
        704 => "quat4", 705 => "quat5", 706 => "quat6",
        707 => "hpi_goodness", 708 => "hpi_error", 709 => "hpi_movement",
        802 => "seeg", 803 => "dbs", 900 => "syst", 902 => "ecog", 910 => "ias",
        920 => "exci", 1000 => "dipole_wave", 1001 => "goodness_fit",
        1100 => "fnirs", 1200 => "temperature", 1300 => "galvanic",
        1400 => "eyetrack",
    ],
    :id,
    :channel_type,
)

# ---------------------------------------------------------------------------
# coil types (FIFFV_COIL_*)
# ---------------------------------------------------------------------------
fiff_coil_type = _fiff_table(
    [
        0 => "none", 1 => "eeg", 2 => "nm_122", 3 => "nm_24",
        4 => "nm_mcg_axial", 5 => "eeg_bipolar", 6 => "eeg_csd",
        200 => "dipole",
        300 => "fnirs_hbo", 301 => "fnirs_hbr", 302 => "fnirs_cw_amplitude",
        303 => "fnirs_od", 304 => "fnirs_fd_ac_amplitude", 305 => "fnirs_fd_phase",
        306 => "fnirs_td_gated_amplitude", 307 => "fnirs_td_moments_intensity",
        308 => "fnirs_td_moments_mean", 309 => "fnirs_td_moments_variance",
        400 => "eyetrack_pos", 401 => "eyetrack_pupil",
        1000 => "mcg_42",
        2000 => "point_magnetometer", 2001 => "axial_grad_5cm",
        3011 => "vv_planar_w", 3012 => "vv_planar_t1", 3013 => "vv_planar_t2",
        3014 => "vv_planar_t3", 3015 => "vv_planar_t4",
        3021 => "vv_mag_w", 3022 => "vv_mag_t1", 3023 => "vv_mag_t2",
        3024 => "vv_mag_t3", 3025 => "vv_mag_t4",
        4001 => "magnes_mag", 4002 => "magnes_grad", 4003 => "magnes_ref_mag",
        4004 => "magnes_ref_grad", 4005 => "magnes_offdiag_ref_grad",
        5001 => "ctf_grad", 5002 => "ctf_ref_mag", 5003 => "ctf_ref_grad",
        5004 => "ctf_offdiag_ref_grad",
        6001 => "kit_grad", 6002 => "kit_ref_mag",
        7001 => "baby_grad", 7002 => "baby_mag", 7003 => "baby_ref_mag",
        7004 => "baby_ref_mag2",
        7501 => "artemis123_grad", 7502 => "artemis123_ref_mag",
        7503 => "artemis123_ref_grad",
        8001 => "quspin_zfopm_mag", 8002 => "quspin_zfopm_mag2",
        8101 => "fieldline_opm_mag_gen1", 8201 => "kernel_opm_mag_gen1",
        9001 => "kriss_grad", 9101 => "compumedics_adult_grad",
        9102 => "compumedics_pediatric_grad",
    ],
    :id,
    :coil_type,
)

fiff_gantry_type = _fiff_table([0 => "fixed", 1 => "uni_axial", 2 => "free"], :id, :gantry_type)

fiff_dacq_system = _fiff_table(
    [0 => "dau", 1 => "vxi", 2 => "rpu", 3 => "orion", 4 => "triux"],
    :id,
    :dacq_system,
)

fiff_proj_item = _fiff_table(
    [
        0 => "none", 1 => "field", 2 => "dip_fix", 3 => "dip_rot",
        4 => "homog_grad", 5 => "homog_field", 10 => "eeg_avref",
    ],
    :id,
    :proj_item,
)

fiff_proj_by = _fiff_table([0 => "complement", 1 => "space"], :id, :proj_by)

fiff_aspect = _fiff_table(
    [
        100 => "average", 101 => "std_err", 102 => "single", 103 => "subaverage",
        104 => "altaverage", 105 => "sample", 106 => "power_density",
        200 => "dipole_wave", 1100 => "ifii_low", 1101 => "ifii_high",
        1102 => "gate",
    ],
    :id,
    :aspect,
)

fiff_sss_job = _fiff_table(
    [
        0 => "sss_job_nothing", 1 => "sss_job_ctc", 2 => "sss_job_filter",
        3 => "sss_job_virt", 4 => "sss_job_head_pos", 5 => "sss_job_movec_fit",
        6 => "sss_job_movec_qua", 7 => "sss_job_rec_all", 8 => "sss_job_rec_in",
        9 => "sss_job_rec_out", 10 => "sss_job_st",
    ],
    :id,
    :sss_job,
)

# ---------------------------------------------------------------------------
# element decoding
# ---------------------------------------------------------------------------

# element type of simple (non-struct) FIFF data types; nothing if not simple
function _fiff_eltype(t::Integer)
    t == 1 && return UInt8
    (t == 2 || t == 16) && return Int16   # short, dau_pack16
    t == 3 && return Int32
    t == 4 && return Float32
    t == 5 && return Float64
    t == 6 && return Int32                # julian day number
    t == 7 && return UInt16
    t == 8 && return UInt32
    t == 9 && return UInt64
    t == 11 && return Int64
    t == 20 && return Float32             # complex float: (re, im) pairs
    t == 21 && return Float64             # complex double
    return nothing
end

# decode a run of simple elements; complex types are returned as ComplexF64
function _fiff_elements(t::Integer, b::AbstractVector{UInt8})
    T = _fiff_eltype(t)
    isnothing(T) && return nothing
    v = _be(T, b)
    if t == 20 || t == 21
        return ComplexF64.(v[1:2:end], v[2:2:end])
    elseif T <: AbstractFloat
        return Float64.(v)
    else
        return Int64.(v)
    end
end

# return a scalar for one-element results, the vector otherwise
_fiff_scalar(v) = (v isa AbstractVector && length(v) == 1) ? v[1] : v

"""
    _fiff_decode(dtype, buf)

Decode a FIFF tag payload according to its data type (not its kind).

Data buffers (kind 300) are intentionally not decoded here; see `_fiff_raw_data`.
"""
function _fiff_decode(dtype::Integer, buf::AbstractVector{UInt8})
    (dtype & 0xFFFF0000) != 0 && return _fiff_matrix(dtype, buf)

    t = dtype & 0x0000FFFF
    isempty(buf) && return nothing
    # void
    t == 0 && return nothing

    # string
    t == 10 && return _fiff_string(buf)

    # julian
    if t == 6
        jdn = _be(Int32, buf[1:4])[1]
        # Julian day numbers start at noon, julian2datetime() expects that epoch
        return Date(julian2datetime(Float64(jdn)))
    end

    # old_pack
    if t == 23
        offset = _f32f64(buf[1:4])
        scale = _f32f64(buf[5:8])
        return Float64.(_be(Int16, buf[9:end])) .* scale .+ offset
    end

    # ch_info_t (96 B)
    if t == 30
        length(buf) < 96 && return nothing
        return (
            scan_no = _i32i64(buf[1:4]),
            log_no = _i32i64(buf[5:8]),
            kind = _i32i64(buf[9:12]),
            range = _f32f64(buf[13:16]),
            cal = _f32f64(buf[17:20]),
            coil_type = _i32i64(buf[21:24]),
            # r0, ex, ey, ez
            loc = Float64.(_be(Float32, buf[25:72])),
            unit = _i32i64(buf[73:76]),
            unit_mul = _i32i64(buf[77:80]),
            ch_name = _fiff_string(buf[81:96]),
        )
    end

    # id_t (20 B)
    if t == 31
        length(buf) < 20 && return nothing
        version = _i32i64(buf[1:4])
        return (
            version >> 16,
            version & 0xFFFF,
            _i32i64(buf[5:8]),
            _i32i64(buf[9:12]),
            unix2datetime(_i32i64(buf[13:16])),
            _i32i64(buf[17:20]),
        )
    end

    # dir_entry_t: not needed
    t == 32 && return nothing

    # dig_point_t (20 B)
    if t == 33
        return (
            _i32i64(buf[1:4]),
            _i32i64(buf[5:8]),
            _f32f64(buf[9:12]),
            _f32f64(buf[13:16]),
            _f32f64(buf[17:20]),
        )
    end

    # ch_pos_t (52 B)
    if t == 34
        return (_i32i64(buf[1:4]), Float64.(_be(Float32, buf[5:52]))...)
    end

    # coord_trans_t (104 B)
    if t == 35
        from = _i32i64(buf[1:4])
        to = _i32i64(buf[5:8])
        f = Float64.(_be(Float32, buf[9:104]))
        # stored row-major
        rot = permutedims(reshape(f[1:9], 3, 3))
        move = f[10:12]
        invrot = permutedims(reshape(f[13:21], 3, 3))
        invmove = f[22:24]
        return (from, to, rot, move, invrot, invmove)
    end

    # dig_string_t
    if t == 36
        kind = _i32i64(buf[1:4])
        ident = _i32i64(buf[5:8])
        np = _i32i64(buf[9:12])
        rr = Float64.(_be(Float32, buf[13:(12 + 12 * np)]))
        return (kind, ident, np, permutedims(reshape(rr, 3, np)))
    end

    v = _fiff_elements(t, buf)
    if isnothing(v)
        _warn(
            "FIFF data type $(_find_fiff_dt(t)) ($t) is not implemented; please send this file to adam.wysokinski@neuroanalyzer.org",
        )
        return nothing
    end
    return _fiff_scalar(v)
end

"""
    _fiff_matrix(dtype, buf)

Decode a FIFF matrix tag (dense, sparse CCS or sparse RCS).

Dense layout: data (row-major, C order), then dims (reversed), then ndim (Int32).

Sparse layout: data[nnz], indices[nnz], pointers[ncol+1 (CCS) | nrow+1 (RCS)], then (nnz, nrow, ncol), then ndim.

Returns a dense `Array` with the same axis order as the writer intended (e.g. `proj_item_vectors` → nvec × nchan).
"""
function _fiff_matrix(dtype::Integer, buf::AbstractVector{UInt8})
    coding = dtype & 0xFFFF0000
    t = dtype & 0x0000FFFF
    length(buf) < 4 && return nothing
    T = _fiff_eltype(t)
    if isnothing(T) || !(t in (3, 4, 5, 20, 21))
        _warn(
            "FIFF matrix of $(_find_fiff_dt(t)) is not implemented; please send this file to adam.wysokinski@neuroanalyzer.org",
        )
        return nothing
    end
    nb = (t == 20 ? 8 : t == 21 ? 16 : sizeof(T))
    ndim = _i32i64(buf[(end - 3):end])

    # dense
    if coding == _FIFFT_MATRIX
        (ndim < 1 || 4 * (ndim + 1) > length(buf)) &&
            throw(ArgumentError("Corrupted FIFF dense matrix header."))
        dims_c = reverse(Int64.(_be(Int32, buf[(end - 4 * (ndim + 1) + 1):(end - 4)])))
        n = prod(dims_c)
        nb * n + 4 * (ndim + 1) > length(buf) &&
            throw(ArgumentError("Corrupted FIFF dense matrix (size mismatch)."))
        v = _fiff_elements(t, buf[1:(nb * n)])
        ndim == 1 && return v
        # C-order data → Julia (column-major): fill reversed dims, then reverse axes
        return permutedims(reshape(v, Tuple(reverse(dims_c))), ndim:-1:1)

    elseif coding == (_FIFFT_MATRIX | _FIFFT_SPARSE_CCS) ||
           coding == (_FIFFT_MATRIX | _FIFFT_SPARSE_RCS)
        ndim == 2 || throw(ArgumentError("Only 2-D sparse FIFF matrices are supported."))
        hdr = Int64.(_be(Int32, buf[(end - 4 * (ndim + 2) + 1):(end - 4)]))
        nnz, nrow, ncol = hdr[1], hdr[2], hdr[3]
        ccs = coding == (_FIFFT_MATRIX | _FIFFT_SPARSE_CCS)
        nptr = ccs ? ncol + 1 : nrow + 1
        o = 0
        vals = _fiff_elements(t, buf[(o + 1):(o + nb * nnz)])
        o += nb * nnz
        idx = Int64.(_be(Int32, buf[(o + 1):(o + 4 * nnz)]))
        o += 4 * nnz
        ptr = Int64.(_be(Int32, buf[(o + 1):(o + 4 * nptr)]))
        (ptr[end] > nnz || any(<(0), ptr)) && throw(
            ArgumentError(
                "Sparse FIFF matrix with invalid pointers (possibly written by an old MNE-C with a byte-order bug); not supported.",
            ),
        )
        m = zeros(eltype(vals), nrow, ncol)
        if ccs
            for c = 1:ncol, k = (ptr[c] + 1):ptr[c + 1]
                m[idx[k] + 1, c] = vals[k]
            end
        else
            for r = 1:nrow, k = (ptr[r] + 1):ptr[r + 1]
                m[r, idx[k] + 1] = vals[k]
            end
        end
        return m
    end

    _warn("Unknown FIFF matrix coding 0x$(string(coding; base = 16)).")
    return nothing
end

# ---------------------------------------------------------------------------
# file traversal
# ---------------------------------------------------------------------------

function _read_fiff_tag(fid::IO)::Tuple{Int32, Int32, Int32, Vector{UInt8}, Int32}
    # header: kind, data_type, data_size, next (4 × big-endian Int32)
    raw = read(fid, 16)
    length(raw) == 16 || throw(EOFError())
    tag = ntoh.(reinterpret(Int32, raw))
    tag_kind, data_type, data_size, tag_next = tag[1], tag[2], tag[3], tag[4]
    data_size < 0 && throw(ArgumentError("Corrupted FIFF tag (negative data size)."))
    data = read(fid, data_size)
    length(data) == data_size || throw(EOFError())
    # next: 0 = sequential, -1 = last tag, >0 = absolute offset of the next tag
    tag_next > 0 && seek(fid, tag_next)
    return tag_kind, data_type, data_size, data, tag_next
end

"""
    _create_fiff_block(fid)

Read all tags sequentially and resolve the block hierarchy.

Returns:

- `Vector{Vector{UInt8}}`: raw payload of every tag
- `Matrix{Int64}`: columns = tag position, kind, data type, size, block level, block type
- `Vector{Int64}`: block instance number of every tag (unique per block occurrence)
"""
function _create_fiff_block(
    fid::IO,
)::Tuple{Vector{Vector{UInt8}}, Matrix{Int64}, Vector{Int64}}
    seek(fid, 0)

    bs = _find_fiff_tag("block_start")
    be = _find_fiff_tag("block_end")

    tag_pos = Int64[]
    tag_ids = Int64[]
    tag_type = Int64[]
    tag_size = Int64[]
    block_level = Int64[]
    block_type = Int64[]
    block_inst = Int64[]
    d = Vector{Vector{UInt8}}()

    # stack of (block type, instance number); root = 999
    stack = [(999, 0)]
    n_inst = 0
    seen = Set{Int64}()

    while !eof(fid)
        pos = position(fid)
        # protect against circular `next` pointers in corrupted files
        pos in seen && (_warn("Circular tag pointer in FIFF file; stopping."); break)
        push!(seen, pos)

        tag_kind, tag_dt, tag_sz, data, tag_next = try
            _read_fiff_tag(fid)
        catch e
            e isa EOFError || rethrow()
            _warn("FIFF file is truncated; reading stopped at byte $pos.")
            break
        end

        if tag_kind == bs
            n_inst += 1
            btype = length(data) >= 4 ? _i32i64(data[1:4]) : -1
            push!(stack, (btype, n_inst))
        end

        push!(tag_pos, pos)
        push!(tag_ids, tag_kind)
        push!(tag_type, tag_dt)
        push!(tag_size, tag_sz)
        push!(block_level, length(stack))
        push!(block_type, stack[end][1])
        push!(block_inst, stack[end][2])
        push!(d, data)

        if tag_kind == be
            length(stack) > 1 ? pop!(stack) : _warn("Unbalanced block_end tag in FIFF file.")
        end

        tag_next == -1 && break
    end

    length(stack) > 1 && _warn("FIFF file has $(length(stack) - 1) unterminated block(s).")

    return d,
    hcat(tag_pos, tag_ids, tag_type, tag_size, block_level, block_type),
    block_inst
end

function _view_fiff_block(fb::Matrix{Int64})::Nothing
    for row_idx in axes(fb, 1)
        println(
            repeat("  ", max(0, fb[row_idx, 5] - 1)) *
            "$(fb[row_idx, 2]) [$(_find_fiff_tag(fb[row_idx, 2]))]",
        )
    end
    return nothing
end

function _get_blocks(b::Matrix{Int64})::Tuple{Vector{Vector{Int64}}, Vector{Int64}}
    levels = unique(b[:, 5])
    bidx = [findall(isequal(l), b[:, 5]) for l in levels]
    btypes = [unique(b[:, 6][bidx[i]])[1] for i in eachindex(bidx)]
    return bidx, btypes
end

function _pack_fiff_blocks(
    fiff_object::Vector{Any},
    block::String,
    fields::Vector{String},
)::Dict
    block_obj = Base.filter(x -> x[3] == block, fiff_object)
    d = Dict{Symbol, Any}()
    for f in fields
        matches = Base.filter(x -> x[2] == f, block_obj)
        if length(matches) == 1
            d[Symbol(f)] = matches[1][4]
        elseif length(matches) > 1
            pad = length(string(length(matches)))
            for (n, m) in enumerate(matches)
                d[Symbol(f * "_$(lpad(n, pad, '0'))")] = m[4]
            end
        else
            d[Symbol(f)] = nothing
        end
    end
    return d
end

function _fiff_tree(fiff_object::Vector{Any})::Nothing
    for item in fiff_object
        println(item[3] * " - " * item[2])
    end
    return nothing
end
