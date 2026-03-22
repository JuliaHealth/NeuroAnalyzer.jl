# ---------------------------------------------------------------------------
# byte-reinterpretation helpers
# ntoh() operates on a scalar — the broadcast dot was unnecessary
# ---------------------------------------------------------------------------
_i16i64(x)::Int64   = Int64(ntoh(reinterpret(Int16,   x)[1]))
_i32i64(x)::Int64   = Int64(ntoh(reinterpret(Int32,   x)[1]))
_f16f64(x)::Float64 = Float64(ntoh(reinterpret(Float16, x)[1]))
_f32f64(x)::Float64 = Float64(ntoh(reinterpret(Float32, x)[1]))
_f64f64(x)::Float64 = Float64(ntoh(reinterpret(Float64, x)[1]))
_i16f64(x)::Float64 = Float64(ntoh(reinterpret(Int16,   x)[1]))
_i32f64(x)::Float64 = Float64(ntoh(reinterpret(Int32,   x)[1]))
_i8i8(x)::Int8      = ntoh(reinterpret(Int8,   x)[1])
_ui32i32(x)::Int32  = ntoh(reinterpret(UInt32, x)[1])
_i32i32(x)::Int32   = ntoh(reinterpret(Int32,  x)[1])

# ---------------------------------------------------------------------------
# FIFF lookup helpers
# ---------------------------------------------------------------------------
_find_fiff_tag(t::String)::Int64          = fiff_tags[:id][findfirst(isequal(t), fiff_tags[:tag])]
_find_fiff_tag(n::Int64)::String          = fiff_tags[:tag][findfirst(isequal(n), fiff_tags[:id])]
_find_fiff_block(t::String)::Int64        = fiff_blocks[:id][findfirst(isequal(t), fiff_blocks[:block])]
_find_fiff_block(n::Int64)::String        = fiff_blocks[:block][findfirst(isequal(n), fiff_blocks[:id])]
_find_fiff_dt(n::Int64)::String           = fiff_data_type[:name][findfirst(isequal(n & 0x00000FFF), fiff_data_type[:id])]
_find_fiff_unit(n::Int64)::String         = fiff_units[:unit][findfirst(isequal(n), fiff_units[:id])]
_find_fiff_mul(n::Int64)::String          = fiff_multipliers[:multiplier][findfirst(isequal(n), fiff_multipliers[:id])]
_find_fiff_chtype(n::Int64)::String       = fiff_channel_type[:channel_type][findfirst(isequal(n), fiff_channel_type[:id])]
_find_fiff_gantry_type(n::Int64)::String  = fiff_gantry_type[:gantry_type][findfirst(isequal(n), fiff_gantry_type[:id])]
_find_fiff_dacq_system(n::Int64)::String  = fiff_dacq_system[:dacq_system][findfirst(isequal(n), fiff_dacq_system[:id])]
_find_fiff_proj_item(n::Int64)::String    = fiff_proj_item[:proj_item][findfirst(isequal(n), fiff_proj_item[:id])]
_find_fiff_proj_by(n::Int64)::String      = fiff_proj_by[:proj_by][findfirst(isequal(n), fiff_proj_by[:id])]
_find_fiff_coiltype(n::Int64)::String     = fiff_coil_type[:coil_type][findfirst(isequal(n), fiff_coil_type[:id])]
_find_fiff_aspect(n::Int64)::String       = fiff_aspect[:aspect][findfirst(isequal(n), fiff_aspect[:id])]
_find_fiff_sss_job(n::Int64)::String      = fiff_sss_job[:sss_job][findfirst(isequal(n), fiff_sss_job[:id])]

# ---------------------------------------------------------------------------
# FIFF data-type dictionary
# ---------------------------------------------------------------------------
fiff_data_type = Dict(
    :id => [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 11, 13, 14, 16, 20, 21, 23, 30, 31, 32, 33, 34, 35],
    :name => [
        "void", "byte", "int16", "int32", "float", "double", "julian",
        "uint16", "uint32", "uint64", "string", "ascii", "int64",
        "dau_pack13", "dau_pack14", "dau_pack16", "complex_float",
        "complex_double", "old_pack", "ch_info_struct", "id_struct",
        "dir_entry_struct", "dig_point_struct", "ch_pos_struct",
        "coord_trans_struct",
    ],
    :data_type => [
        "void_t", "byte_t", "int16_t", "int32_t", "float_t", "double_t",
        "julian_t", "uint16_t", "uint32_t", "uint64_t", "byte_t", "byte_t",
        "int64_t", "dau_pack13_t", "dau_pack14_t", "dau_pack16_t",
        "complex_float_t", "complex_double_t", "old_pack_tvariable",
        "ch_info_t", "id_t", "dir_entry_t", "dig_point_t", "ch_pos_t",
        "coord_trans_t",
    ],
    :size => [1, 1, 2, 4, 4, 8, 8, 2, 4, 8, 1, 1, 8, 2, 2, 2, 8, 16, 80, 20, 16, 20, 52, 80, 80],
)

# FIFF physical units lookup table.
# Sources: MNE-Python fiff/constants.py and the FIFF standard specification.
fiff_units = Dict(
    :id => [
        -1,
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        101,
        102,
        103,
        104,
        105,
        106,
        107,
        108,
        109,
        110,
        111,
        112,
        113,
        114,
        115,
        116,
        201,
        202,
    ],
    :unit => [
        "none",
        "unitless",
        "m",        # metre
        "kg",       # kilogram
        "s",        # second
        "A",        # ampere
        "K",        # kelvin
        "mol",      # mole
        "rad",      # radian
        "sr",       # steradian
        "cd",       # candela
        "Hz",       # hertz
        "N",        # newton
        "Pa",       # pascal
        "J",        # joule
        "W",        # watt
        "C",        # coulomb
        "V",        # volt
        "F",        # farad
        "Ω",        # ohm
        "Mho",      # mho (siemens, conductance)
        "Wb",       # weber
        "T",        # tesla
        "H",        # henry
        "°C",       # degree Celsius — was: "C" (duplicate of coulomb at id=108)
        "lm",       # lumen
        "lx",       # lux
        "T/m",      # tesla per metre
        "Am",       # ampere·metre
    ],
)

# FIFF value-multiplier (SI prefix) lookup table.
# :id    — the integer exponent stored in the FIFF file (power of 10).
# :multiplier — short text identifier used in the FIFF specification.
fiff_multipliers = Dict(
    :id => [18, 15, 12, 9, 6, 3, 2, 1, 0, -1, -2, -3, -6, -9, -12, -15, -18],
    :multiplier => [
        "exa",   # 10^18  E
        "peta",  # 10^15  P
        "tera",  # 10^12  T
        "giga",  # 10^9   G
        "mega",  # 10^6   M
        "k",     # 10^3   kilo
        "h",     # 10^2   hecto
        "da",    # 10^1   deca
        "none",  # 10^0   unity
        "d",     # 10^-1  deci
        "c",     # 10^-2  centi
        "m",     # 10^-3  milli
        "μ",     # 10^-6  micro
        "n",     # 10^-9  nano
        "p",     # 10^-12 pico
        "f",     # 10^-15 femto
        "a",     # 10^-18 atto
    ],
)

# FIFF tag ID → tag name lookup table.
# Sources: MNE-Python fiff/constants.py and the FIFF standard specification.
# Tags are grouped by functional area matching the FIFF block hierarchy.
fiff_tags = Dict(
    :id => [
        # --- File structure (100–151) ---
        100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
        110, 111, 112, 113, 114, 115, 116, 117, 118, 120,
        # Note: 119 is not defined in the FIFF standard
        150, 151,

        # --- Measurement info (200–247) ---
        200, 201, 202, 203, 204, 205, 206, 207, 208, 209,
        210, 211, 212, 213, 214, 215, 216, 217, 218, 219,
        220, 221, 222, 223, 224, 225, 226, 227, 228, 229,
        230, 231,
        # Note: 232 is not defined in the FIFF standard
        233, 234, 235, 236, 237,
        240, 241, 242, 243, 244, 245, 246, 247,

        # --- Channel info (250–282) ---
        250, 251, 252, 253, 254, 255, 256, 257, 258,
        263, 264, 265, 266, 267, 268, 269, 270, 271, 272,
        273, 274, 275, 276, 277, 278, 279, 280, 281, 282,

        # --- Gantry (300–305) ---
        300, 301, 302, 303, 304, 305,

        # --- Data buffers (400–410) ---
        400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410,

        # --- Subject info (500–504) ---
        500, 501, 502, 503, 504,

        # --- Project info (600–603) ---
        600, 601, 602, 603,

        # --- Events (701–703) ---
        701, 702, 703,

        # --- SQUID / hardware (800) ---
        800,

        # --- MRI / volume (2001–2200) ---
        2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010,
        2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020,
        2021, 2022, 2023, 2024,
        2030, 2031, 2032,
        2040, 2041,
        2100,
        2200,

        # --- BEM / conductor model (3000–3201) ---
        3000, 3001, 3002, 3003,
        3101, 3102, 3103, 3104, 3105, 3106, 3107, 3108, 3109, 3110,
        3111, 3112, 3113,
        3201,

        # --- Source / xfit (3300–3417) ---
        3300,
        3401, 3402, 3403, 3404, 3405, 3406, 3407, 3408,
        3411, 3412, 3413, 3414, 3415, 3416, 3417,

        # --- Projector item (3501) ---
        3501,
    ],
    :tag => [
        # --- File structure ---
        "file_id", "dir_pointer", "dir", "block_id", "block_start",
        "block_end", "free_list", "free_block", "nop", "parent_file_id",
        "parent_block_id", "block_name", "block_version", "creator", "modifier",
        "ref_role", "ref_file_id", "ref_file_num", "ref_file_name", "ref_block_id",
        "dacq_pars", "dacq_stim",

        # --- Measurement info ---
        "nchan", "sfreq", "data_pack", "ch_info", "meas_date", "subject",
        "description", "nave", "first_sample", "last_sample", "aspect_kind",
        "ref_event", "experimenter", "dig_point", "ch_pos_vec", "hpi_slopes",
        "hpi_ncoil", "req_event", "req_limit", "lowpass", "bad_chs",
        "artef_removal", "coord_trans", "highpass", "ch_cals_vec", "hpi_bad_chs",
        "hpi_corr_coeff", "event_comment", "no_samples", "first_time",
        "subave_size", "subave_first",
        "name", "dig_string", "line_freq", "hpi_coil_freq", "signal_channel",
        "hpi_coil_moments", "hpi_fit_goodness", "hpi_fit_accept",
        "hpi_fit_good_limit", "hpi_fit_dist_limit", "hpi_coil_no",
        "hpi_coils_used", "hpi_digitization_order",

        # --- Channel info ---
        "ch_scan_no", "ch_logical_no", "ch_kind", "ch_range", "ch_cal",
        "ch_pos", "ch_unit", "ch_unit_mul", "ch_dacq_name",
        "sss_frame", "sss_job", "sss_origin", "sss_ord_in", "sss_ord_out",
        "sss_nmag", "sss_components", "sss_cal_chans", "sss_cal_corrs",
        "sss_st_corr", "sss_base_in", "sss_base_out", "sss_base_virt",
        "sss_norm", "sss_iterate", "sss_nfree", "sss_st_length",

        # --- Gantry ---
        "gantry_type", "gantry_model", "gantry_angle",
        "data_buffer", "data_skip", "epoch",

        # --- Data buffers ---
        "data_skip_samp", "data_buffer2", "time_stamp",
        "subj_id", "subj_first_name", "subj_middle_name", "subj_last_name",
        "subj_birth_day", "subj_sex", "subj_hand", "subj_weight",

        # --- Subject info ---
        "subj_height", "subj_comment", "subj_his_id",
        "proj_id", "proj_name",

        # --- Project info ---
        "proj_aim", "proj_persons", "proj_comment",
        "event_channels", "event_list",

        # --- Events ---
        "event_channel", "event_bits",
        "squid_bias",

        # --- SQUID / hardware ---
        "squid_offset", "squid_gate", "decoupler_matrix",

        # --- MRI / volume ---
        "volume_type", "mri_source_format", "mri_pixel_encoding",
        "mri_pixel_data_offset", "mri_pixel_scale", "mri_pixel_data",
        "mri_pixel_overlay_encoding", "mri_pixel_overlay_data",
        "mri_bounding_box", "mri_width", "mri_width_m", "mri_height",
        "mri_height_m", "mri_depth", "mri_depth_m", "mri_thickness",
        "mri_scene_aim", "mri_calibration_scale", "mri_calibration_offset",
        "mri_orig_source_path", "mri_orig_source_format",
        "mri_orig_pixel_encoding", "mri_orig_pixel_data_offset",
        "mri_time", "mri_voxel_data", "mri_voxel_encoding", "voxel_nchannels",
        "mri_diffusion_weight", "mri_diffusion_param", "mri_mrilab_setup",
        "mri_seg_region_id",

        # --- BEM / conductor model ---
        "conductor_model_kind", "sphere_origin", "modelsphere_coord_frame",
        "sphere_layers",
        "bem_surf_id", "bem_surf_name", "bem_surf_nnode", "bem_surf_ntri",
        "bem_surf_nodes", "bem_surf_triangles", "bem_surf_normals",
        "bem_surf_curvs", "bem_surf_curv_values", "bem_pot_solution",
        "bem_approx", "bem_coord_frame", "bem_sigma",
        "source_dipole",

        # --- Source / xfit ---
        "beamformer_instructions",
        "xfit_lead_products", "xfit_map_products", "xfit_grad_map_products",
        "xfit_vol_integration", "xfit_integration_radius",
        "xfit_conductor_model_name", "xfit_conductor_model_trans_name",
        "xfit_cont_surf_type",
        "proj_item_kind", "proj_item_time", "proj_item_ign_chs",
        "proj_item_nvec", "proj_item_vectors", "proj_item_definition",
        "proj_item_ch_name_list",

        # --- Projector item ---
        "xplotter_layout",
    ],
)

# FIFF block type ID → block name lookup table.
# Sources: MNE-Python fiff/constants.py and the FIFF standard specification.
# Block IDs define the hierarchical container structure of a FIFF file.
fiff_blocks = Dict(
    :id => [
        # --- Root ---
        999,

        # --- Measurement (100–122) ---
        100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
        110, 111, 112,
        # Note: 113 is not defined in the FIFF standard (gap is intentional)
        114, 115, 116, 117, 118, 119, 120, 121, 122,

        # --- Structural / MRI (200–206) ---
        200, 201, 202, 203, 204, 205, 206,

        # --- BEM / conductor / projection (300–359) ---
        300,
        310, 311, 312, 313, 314, 315,
        359,

        # --- Volume (400) ---
        400,

        # --- SSS / data correction (500–510) ---
        500, 501, 502, 503, 504, 505, 510,

        # --- Processing history (900–901) ---
        900, 901,
    ],
    :block => [
        # --- Root ---
        "root",

        # --- Measurement ---
        "meas",                # 100
        "meas_info",           # 101
        "raw_data",            # 102
        "processed_data",      # 103
        "evoked",              # 104
        "aspect",              # 105
        "subject",             # 106
        "isotrak",             # 107
        "hpi_meas",            # 108
        "hpi_result",          # 109
        "hpi_coil",            # 110
        "project",             # 111
        "continuous_data",     # 112
        # 113 not defined
        "void",                # 114
        "events",              # 115
        "index",               # 116
        "dacq_pars",           # 117
        "ref",                 # 118
        "maxshield_raw_data",  # 119
        "maxshield_aspect",    # 120
        "hpi_subsystem",       # 121
        "phantom_subsystem",   # 122

        # --- Structural / MRI ---
        "structural_data",     # 200
        "volume_data",         # 201
        "volume_slice",        # 202
        "scenery",             # 203
        "scene",               # 204
        "mri_seg",             # 205
        "mri_seg_region",      # 206

        # --- BEM / conductor / projection ---
        "sphere",              # 300
        "bem",                 # 310
        "bem_surf",            # 311
        "conductor_model",     # 312
        "xfit_proj",           # 313
        "xfit_proj_item",      # 314
        "xfit_aux",            # 315
        "bad_channels",        # 359

        # --- Volume ---
        "vol_info",            # 400

        # --- SSS / data correction ---
        "data_correction",     # 500
        "channels_decoupler",  # 501
        "sss_info",            # 502
        "sss_cal_adjust",      # 503
        "sss_st_info",         # 504
        "sss_bases",           # 505
        "maxshield",           # 510

        # --- Processing history ---
        "processing_history",  # 900
        "processing_record",   # 901
    ],
)

# FIFF channel type ID → channel type name lookup table.
# Sources: MNE-Python fiff/constants.py and the FIFF standard specification.
fiff_channel_type = Dict(
    :id => [
        # --- Primary recording channels (1–3) ---
        1,    # MEG
        2,    # EEG
        3,    # Stimulus / trigger

        # --- Physiological channels (102–602) ---
        102,  # Biological (generic)
        201,  # MCG (magnetocardiography)
        202,  # EOG (electro-oculography)
        301,  # MEG reference
        302,  # EMG (electromyography)
        402,  # ECG (electrocardiography)
        502,  # Miscellaneous
        602,  # Respiration

        # --- Head-position / HPI channels (700–709) ---
        700,  # Quaternion component 0 (scalar)
        701,  # Quaternion component 1
        702,  # Quaternion component 2
        703,  # Quaternion component 3
        704,  # Quaternion component 4
        705,  # Quaternion component 5
        706,  # Quaternion component 6
        707,  # HPI goodness-of-fit
        708,  # HPI error
        709,  # HPI movement indicator

        # --- System / internal channels (900–920) ---
        900,  # System status
        910,  # IAS (internal active shielding)
        920,  # Excitation

        # --- Source modelling (1000–1001) ---
        1000, # Dipole wave
        1001, # Goodness of fit
    ],
    :channel_type => [
        # --- Primary recording channels ---
        "meg",
        "eeg",
        "stim",

        # --- Physiological channels ---
        "bio",
        "mcg",
        "eog",
        "meg_ref",
        "emg",
        "ecg",
        "misc",
        "resp",

        # --- Head-position / HPI channels ---
        "quat0",
        "quat1",
        "quat2",
        "quat3",
        "quat4",
        "quat5",
        "quat6",
        "hpi_goodness",
        "hpi_error",
        "hpi_movement",

        # --- System / internal channels ---
        "syst",
        "ias",
        "exci",

        # --- Source modelling ---
        "dipole_wave",
        "goodness_fit",
    ],
)

# FIFF coil type ID → coil type name lookup table.
# Sources: MNE-Python fiff/constants.py and the FIFF standard specification.
fiff_coil_type = Dict(
    :id => [
        # --- Generic / reference (0–5) ---
        0,    # No coil (placeholder)
        1,    # EEG electrode
        2,    # Neuromag-122 planar gradiometer
        3,    # Neuromag-24 magnetometer
        4,    # Generic axial MCG gradiometer
        5,    # EEG bipolar electrode pair

        # --- Virtual / theoretical (200) ---
        200,  # Magnetic dipole (theoretical)

        # --- CTF / older systems (1000) ---
        1000, # CTF/Biomagnetic 42-channel MCG

        # --- Generic reference sensors (2000–2001) ---
        2000, # Point magnetometer
        2001, # Axial gradiometer, 5 cm baseline

        # --- Neuromag VectorView planar gradiometers (3011–3014) ---
        3011, # VectorView planar gradiometer, "W" geometry
        3012, # VectorView planar gradiometer, T1 geometry
        3013, # VectorView planar gradiometer, T2 geometry
        3014, # VectorView planar gradiometer, T3 geometry

        # --- Neuromag VectorView magnetometers (3021–3024) ---
        3021, # VectorView magnetometer, "W" geometry
        3022, # VectorView magnetometer, T1 geometry
        3023, # VectorView magnetometer, T2 geometry
        3024, # VectorView magnetometer, T3 geometry

        # --- Magnes (BTi) system (4001–4002) ---
        4001, # Magnes magnetometer
        4002, # Magnes planar gradiometer

        # --- CTF gradiometer (5001) ---
        5001, # CTF axial gradiometer
    ],
    :coil_type => [
        # --- Generic / reference ---
        "none",
        "eeg",
        "nm_122",
        "nm_24",
        "nm_mcg_axial",
        "eeg_bipolar",

        # --- Virtual / theoretical ---
        "dipole",

        # --- CTF / older systems ---
        "mcg_42",

        # --- Generic reference sensors ---
        "point_magnetometer",
        "axial_grad_5cm",

        # --- Neuromag VectorView planar gradiometers ---
        "vv_planar_w",
        "vv_planar_t1",
        "vv_planar_t2",
        "vv_planar_t3",

        # --- Neuromag VectorView magnetometers ---
        "vv_mag_w",
        "vv_mag_t1",
        "vv_mag_t2",
        "vv_mag_t3",

        # --- Magnes (BTi) system ---
        "magnes_mag",
        "magnes_grad",

        # --- CTF gradiometer ---
        "ctf_grad",
    ],
)

# FIFF gantry type ID → gantry type name lookup table.
# Describes the mechanical freedom of the MEG scanner gantry.
# Sources: MNE-Python fiff/constants.py and the FIFF standard specification.
fiff_gantry_type = Dict(
    :id         =>  [0,       1,           2     ],
    :gantry_type => ["fixed", "uni_axial", "free"],
)

# FIFF data acquisition system ID → system name lookup table.
# Identifies the hardware platform that recorded the FIFF file.
# Sources: MNE-Python fiff/constants.py and the FIFF standard specification.
fiff_dacq_system = Dict(
    :id         =>  [0,     1,     2,     3,        4     ],
    :dacq_system => ["dau", "vxi", "rpu", "orion", "triux"],
)

# FIFF SSP projection item type ID → item type name lookup table.
# Defines the type of spatial signal-space projection (SSP) vector stored
# in a projection block. Used for artefact suppression and EEG re-referencing.
# Sources: MNE-Python fiff/constants.py and the FIFF standard specification.
fiff_proj_item = Dict(
    :id        => [0,      1,       2,         3,         4,             5,             10         ],
    :proj_item => ["none", "field", "dip_fix", "dip_rot", "homog_grad",  "homog_field", "eeg_avref"],
)

# FIFF SSP projection method ID → method name lookup table.
# Defines whether the projection is applied in the signal complement space
# (nulling the projected component) or directly in the projection space.
# Sources: MNE-Python fiff/constants.py and the FIFF standard specification.
fiff_proj_by = Dict(
    :id      => [0,            1      ],
    :proj_by => ["complement", "space"],
)

# FIFF aspect type ID → aspect name lookup table.
# The "aspect" describes how averaged or processed data within an evoked
# block should be interpreted (e.g. grand average, single trial, power).
# Sources: MNE-Python fiff/constants.py and the FIFF standard specification.
fiff_aspect = Dict(
    :id => [
        # --- Standard evoked / averaged data (100–106) ---
        100,  # Grand average across trials
        101,  # Standard error of the average
        102,  # Single trial (no averaging)
        103,  # Sub-average (partial average)
        104,  # Alternating-sign average
        105,  # Single sample
        106,  # Power spectral density

        # --- Source modelling (200) ---
        200,  # Dipole wave (source estimate)

        # --- IFII gating (1100–1102) ---
        # Note: gap between 106 and 200, and between 200 and 1100, is intentional
        1100, # IFII low-frequency gate
        1101, # IFII high-frequency gate
        1102, # Gate (generic)
    ],
    :aspect => [
        # --- Standard evoked / averaged data ---
        "average",
        "std_err",
        "single",
        "subaverage",
        "altaverage",
        "sample",
        "power_density",

        # --- Source modelling ---
        "dipole_wave",

        # --- IFII gating ---
        "ifii_low",
        "ifii_high",
        "gate",
    ],
)

# FIFF SSS (Signal Space Separation) job type ID → job name lookup table.
# Identifies the type of SSS processing operation applied to MEG data.
# Sources: MNE-Python fiff/constants.py and the FIFF standard specification.
fiff_sss_job = Dict(
    :id => [
        # --- No operation (0) ---
        0,    # No SSS processing applied

        # --- Cross-talk correction (1) ---
        1,    # Cross-talk correction (CTC)

        # --- Core SSS operations (2–4) ---
        2,    # Spatial SSS filter
        3,    # Virtual sensor projection
        4,    # Head position estimation

        # --- Movement compensation (5–6) ---
        5,    # Movement compensation via least-squares fit
        6,    # Movement compensation via quaternion interpolation

        # --- Signal reconstruction (7–9) ---
        7,    # Reconstruct all components (internal + external)
        8,    # Reconstruct internal components only
        9,    # Reconstruct external components only

        # --- Spatiotemporal SSS (10) ---
        10,   # Spatiotemporal SSS (tSSS) filtering
    ],
    :sss_job => [
        # --- No operation ---
        "sss_job_nothing",

        # --- Cross-talk correction ---
        "sss_job_ctc",

        # --- Core SSS operations ---
        "sss_job_filter",
        "sss_job_virt",
        "sss_job_head_pos",

        # --- Movement compensation ---
        "sss_job_movec_fit",
        "sss_job_movec_qua",

        # --- Signal reconstruction ---
        "sss_job_rec_all",
        "sss_job_rec_in",
        "sss_job_rec_out",

        # --- Spatiotemporal SSS ---
        "sss_job_st",
    ],
)

function _fiff_matrix(fb::Int64, buf::Vector{UInt8})::Union{Vector{Float64}, Matrix{Float64}}
    df      = _find_fiff_dt(fb)
    fs_mask = fb & 0xFF000000

    if fs_mask == 0x00000000   # scalar value

        d = Float64[]
        if df == "float"
            for idx in 1:4:length(buf)
                push!(d, _f32f64(buf[idx:(idx + 3)]))
            end
        elseif df == "old_pack"
            # packed format: values stored as Int16, scaled by a and shifted by b
            a = _f32f64(buf[1:4])
            b = _f32f64(buf[5:8])
            for idx in 9:2:length(buf)
                push!(d, _i16f64(buf[idx:(idx + 1)]))
            end
            d = a .* (d .+ b)
        else
            _warn("scalar of $df is not implemented; please send this file to adam.wysokinski@neuroanalyzer.org")
        end
        return d

    elseif fs_mask == 0x40000000   # matrix value

        mc_mask = fb & 0x00FF0000

        if mc_mask == 0x00000000   # dense matrix
            # read the number of dimensions from the last 4 bytes
            n        = _i32i32(buf[(end - 3):end])
            dims_buf = buf[(end - 5 * n - 1):(end - 4)]

            # parse dimension sizes and reverse to row-major order
            dim = Int64[]
            for dim_idx in 1:4:length(dims_buf)
                push!(dim, _i32i32(dims_buf[dim_idx:(dim_idx + 3)]))
            end
            reverse!(dim)

            # data payload: everything before the dimension header
            tmp = buf[1:(end - length(dims_buf) - 4)]
            d   = Float64[]

            if df == "float"
                for idx in 1:4:length(tmp)
                    push!(d, _f32f64(tmp[idx:(idx + 3)]))
                end
            elseif df == "int32"
                for idx in 1:4:length(tmp)
                    push!(d, _i32f64(tmp[idx:(idx + 3)]))
                end
            elseif df == "old_pack"
                a = _f32f64(tmp[1:4])
                b = _f32f64(tmp[5:8])
                for idx in 9:2:length(tmp)
                    push!(d, _i16f64(tmp[idx:(idx + 1)]))
                end
                d = a .* (d .+ b)
            else
                _warn("dense matrix of $df is not implemented; please send this file to adam.wysokinski@neuroanalyzer.org")
                return Float64[]
            end
            return reshape(d, dim[1], dim[2])

        elseif mc_mask == 0x00100000   # sparse, column-compressed (CCS) matrix

            # --- parse dimension block ---
            n        = _i32i32(buf[(end - 3):end])
            dims_buf = buf[(end - 5 * n - 1):(end - 4)]
            dim      = Int64[]
            for dim_idx in 1:4:length(dims_buf)
                push!(dim, _i32i32(dims_buf[dim_idx:(dim_idx + 3)]))
            end
            reverse!(dim)

            # --- parse non-zero count (nz) — two sequential reads are intentional ---
            tmp = buf[1:(end - length(dims_buf) - 4)]
            nz  = _i32i32(tmp[(end - 3):end])
            tmp = buf[1:(end - length(dims_buf) - 8)]
            nz  = _i32i32(tmp[(end - 3):end])   # second read advances past 4 bytes
            tmp = buf[1:(end - length(dims_buf) - 12)]

            # --- parse column start indices (dim[2] Int32 values) ---
            cs_buf        = tmp[(end + 1 - dim[2] * 4):end]
            col_start_idx = Int64[]
            for idx in 1:4:length(cs_buf)
                push!(col_start_idx, _i32i32(cs_buf[idx:(idx + 3)]))
            end
            col_start_idx .+= 1   # convert 0-based → 1-based

            # --- parse row indices (nz Int32 values) ---
            tmp     = buf[1:(end - length(cs_buf) - length(dims_buf) - 12)]
            ri_buf  = tmp[(end + 1 - nz * 4):end]
            row_idx = Int64[]
            for idx in 1:4:length(ri_buf)
                push!(row_idx, _i32i32(ri_buf[idx:(idx + 3)]))
            end
            row_idx .+= 1   # convert 0-based → 1-based
            tmp = tmp[1:(end - length(ri_buf))]

            # --- parse non-zero values and fill the dense output matrix ---
            if df == "float"
                m = Float64[]
                for idx in 1:4:length(tmp)
                    push!(m, _f32f64(tmp[idx:(idx + 3)]))
                end
                d   = zeros(dim[1], dim[2])
                col = 1
                rel = 0
                nr  = reverse(diff(col_start_idx))
                @inbounds for idx in eachindex(m)
                    d[col, row_idx[idx]] = m[idx]
                    rel += 1
                    # advance to the next column when rel reaches the column's non-zero count
                    if !isempty(nr) && rel == nr[end]
                        col += 1
                        rel  = 0
                        pop!(nr)
                    end
                end
                return d

            else

                _warn("sparse CCS of $df is not implemented; please send this file to adam.wysokinski@neuroanalyzer.org")
                return Float64[]

            end

        elseif mc_mask == 0x00200000   # sparse, row-compressed (CRS) — not yet implemented

            _warn("sparse row-compressed matrix is not implemented; please send this file to adam.wysokinski@neuroanalyzer.org")
            return Float64[]

        end
    end

    # fallback: satisfies return-type constraint for unrecognised fs_mask values
    return Float64[]
end

function _read_fiff_tag(fid::IOStream)::Tuple{Int32, Int32, Int32, Vector{UInt8}, Int32}
    # read 4 × Int32 header fields: kind, data_type, data_size, next
    raw = read(fid, 4 * sizeof(Int32))
    tag = ntoh.(reinterpret(Int32, raw))
    tag_kind = tag[1]
    data_type = tag[2]
    data_size = tag[3]
    tag_next = tag[4]
    data = read(fid, data_size)
    tag_next > 0 && seek(fid, tag_next)
    return tag_kind, data_type, data_size, data, tag_next
end

function _get_fiff_block_type(
    fid::IOStream,
    tag::Tuple{Int64, Int64, Int64, Int64, Vector{UInt8}, Int64}
)::Vector{Int32}
    seek(fid, tag[1] + 16)
    buf = zeros(UInt8, tag[4])
    readbytes!(fid, buf, tag[4])
    return reinterpret(Int32, reverse(buf))
end

function _create_fiff_block(fid::IOStream)::Tuple{Vector{Vector{UInt8}}, Matrix{Int64}}
    seek(fid, 0)

    tags = Vector{Tuple{Int64, Int64, Int64, Int64, Vector{UInt8}, Int64}}()
    tag_next = nothing

    while tag_next != -1
        pos = position(fid)
        tag_kind, tag_type, tag_size, data, tag_next = _read_fiff_tag(fid)
        push!(tags, (pos, tag_kind, tag_type, tag_size, data, tag_next))
    end

    seek(fid, 0)

    n = length(tags)
    tag_pos = zeros(Int64, n)
    tag_ids = zeros(Int64, n)
    tag_type = zeros(Int64, n)
    tag_size = zeros(Int64, n)
    block_level = ones(Int64, n)
    # pre-allocate
    block_type = Vector{Int64}(undef, n)
    block_type_current = 999
    bs = _find_fiff_tag("block_start")
    be = _find_fiff_tag("block_end")
    d  = Vector{Vector{UInt8}}(undef, n)

    @inbounds for i in eachindex(tags)
        tag_pos[i] = tags[i][1]
        tag_ids[i] = tags[i][2]
        tag_type[i] = tags[i][3]
        tag_size[i] = tags[i][4]
        d[i] = tags[i][5]

        if tag_ids[i] == bs
            block_level[i:end] .+= 1
            block_type_current = _get_fiff_block_type(fid, tags[i])[]
        elseif tag_ids[i] == be
            block_level[i:end] .-= 1
        end
        block_type[i] = block_type_current
    end

    return d, hcat(tag_pos, tag_ids, tag_type, tag_size, block_level, block_type)
end

function _view_fiff_block(fb::Matrix{Int64})::Nothing
    level = 0
    for row_idx in axes(fb, 1)
        id = fb[row_idx, 2]
        # block_start → increase indent
        id == 104 && (level += 1)
        println(repeat("  ", level) * "$id [$(fb[row_idx, 5])]")
        # block_end → decrease indent
        id == 105 && (level = max(0, level - 1))
    end
    return nothing
end

function _get_blocks(b::Matrix{Int64})::Tuple{Vector{Vector{Int64}}, Vector{Int64}}
    levels = unique(b[:, 5])
    bidx = [findall(isequal(l), b[:, 5]) for l in levels]
    btypes = [unique(b[:, 6][bidx[i]])[1] for i in eachindex(bidx)]
    return bidx, btypes
end

function _pack_fiff_blocks(fiff_object::Vector{Any}, block::String, fields::Vector{String})::Dict
    block_mask = [fiff_object[i][3] for i in eachindex(fiff_object)] .== block
    block_obj  = fiff_object[block_mask]
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
