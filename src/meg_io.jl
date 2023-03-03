"""
    meg_import_fiff(file_name; detect_type)

Load FIFF (Functional Image File Format) file and return `NeuroAnalyzer.MEG` object.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `meg:MEG`

# Source

Elekta Neuromag: Functional Image File Format Description. FIFF version 1.3. March 2011
"""
function meg_import_fiff(file_name::String; detect_type::Bool=true)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    fid, fiff_block = _create_fiff_block(file_name)

    # block ids
    fiff_block_meas = 100
    fiff_block_meas_info = 101
    fiff_block_raw_data = 102
    fiff_block_processed_data = 103
    fiff_block_evoked = 104
    fiff_block_mcg_ave = 104
    fiff_block_aspect = 105
    fiff_block_subject = 106
    fiff_block_isotrak = 107
    fiff_block_hpi_meas = 108
    fiff_block_hpi_result = 109
    fiff_block_hpi_coil = 110
    fiff_block_project = 111
    fiff_block_continuous_data = 112
    fiff_block_void = 114
    fiff_block_events = 115
    fiff_block_index = 116
    fiff_block_dacq_pars = 117
    fiff_block_ref = 118
    fiff_block_maxshield_raw_data = 119
    fiff_block_maxshield_aspect = 120
    fiff_block_hpi_subsystem = 121
    fiff_block_phantom_subsystem = 122
    fiff_block_structural_data = 200
    fiff_block_volume_data = 201
    fiff_block_volume_slice = 202
    fiff_block_scenery = 203
    fiff_block_scene = 204
    fiff_block_mri_seg = 205
    fiff_block_mri_seg_region = 206
    fiff_block_sphere = 300
    fiff_block_bem = 310
    fiff_block_bem_surf = 311
    fiff_block_conductor_model = 312
    fiff_block_xfit_proj = 313
    fiff_block_xfit_proj_item = 314
    fiff_block_xfit_aux = 315
    fiff_block_bad_channels = 359
    fiff_block_vol_info = 400
    fiff_block_data_correction = 500
    fiff_block_channels_decoupler = 501
    fiff_block_sss_info = 502
    fiff_block_sss_cal_adjust = 503
    fiff_block_sss_st_info = 504
    fiff_block_sss_bases = 505
    fiff_block_maxshield = 510
    fiff_block_processing_history = 900
    fiff_block_processing_record = 901

    # # process blocks
    # for block_idx1 in unique(block_number)
    #     block_current = fiff_blocks[fiff_blocks[:, 5] .== block_idx1, :]
    #     # remove ignored blocks
    #     for block_idx2 in size(block_current, 1):-1:1
    #         if block_current[block_idx2, 7] == 999
    #             block_current = block_current[1:end .!= block_idx2, :]
    #             continue
    #         end
    #         # remove ignored tags from block tree
    #         block_current[block_idx2, 7] in [fiff_block_start, fiff_block_end, fiff_block_nop, fiff_block_root] && (block_current = block_current[1:end .!= block_idx2, :])
    #     end
    #     
    #     # skip empty blocks
    #     size(block_current, 1) == 0 && continue
    #     
    #     # continue with non-empty blocks
    #     for block_idx3 in size(block_current, 1)
    #         # current_block_type
    #         if block_current[block_idx3, 7] == fiff_block_meas
    #         end
    #     end
    # end

    # MEG data tags
    fiff_dacq_pars = 150
    fiff_dacq_stim = 151
    fiff_nchan = 200
    fiff_sfreq = 201
    fiff_data_pack = 202
    fiff_ch_info = 203
    fiff_meas_date = 204
    fiff_description = 206
    fiff_nave = 207
    fiff_first_sample = 208
    fiff_last_sample = 209
    fiff_aspect_kind = 210
    fiff_ref_event = 211
    fiff_experimenter = 212
    fiff_dig_point = 213
    fiff_hpi_slopes = 215
    fiff_hpi_ncoil = 216
    fiff_req_event = 217
    fiff_req_limit = 218
    fiff_lowpass = 219
    fiff_bad_chs = 220
    fiff_artef_removal = 221
    fiff_coord_trans = 222
    fiff_highpass = 223
    fiff_hpi_bad_chs = 225
    fiff_hpi_corr_coeff = 226
    fiff_event_comment = 227
    fiff_no_samples = 228
    fiff_first_time = 229
    fiff_subave_size = 230
    fiff_subave_first = 231
    fiff_name = 233
    fiff_dig_string = 234
    fiff_line_freq = 235
    fiff_hpi_coil_freq = 236
    fiff_signal_channel = 237
    fiff_hpi_coil_moments = 240
    fiff_hpi_fit_goodness = 241

    # read data from tags
    channel_n = round.(Int64, _read_fiff_tag(fid, fiff_blocks, fiff_nchan))
    sampling_rate = round.(Int64, _read_fiff_tag(fid, fiff_blocks, fiff_sfreq))
    data_pack = round.(Int64, _read_fiff_tag(fid, fiff_blocks, fiff_data_pack))
    date = _read_fiff_tag(fid, fiff_blocks, fiff_meas_date)
    date = unix2datetime(date[1][2])
    description = _read_fiff_tag(fid, fiff_blocks, fiff_description)
    nave = _read_fiff_tag(fid, fiff_blocks, fiff_nave)
    first_sample = _read_fiff_tag(fid, fiff_blocks, fiff_first_sample)
    last_sample = _read_fiff_tag(fid, fiff_blocks, fiff_last_sample)
    aspect_kind = _read_fiff_tag(fid, fiff_blocks, fiff_aspect_kind)
    ref_event = _read_fiff_tag(fid, fiff_blocks, fiff_ref_event)
    experimenter = _read_fiff_tag(fid, fiff_blocks, fiff_experimenter)
    dig_point = _read_fiff_tag(fid, fiff_blocks, fiff_dig_point)
    hpi_slopes = _read_fiff_tag(fid, fiff_blocks, fiff_hpi_slopes)
    hpi_ncoil = _read_fiff_tag(fid, fiff_blocks, fiff_hpi_ncoil)
    req_event = _read_fiff_tag(fid, fiff_blocks, fiff_req_event)
    req_limit = _read_fiff_tag(fid, fiff_blocks, fiff_req_limit)
    req_limit = _read_fiff_tag(fid, fiff_blocks, fiff_req_limit)
    lowpass = round.(Float64.(_read_fiff_tag(fid, fiff_blocks, fiff_lowpass)), digits=4)
    bad_chs = _read_fiff_tag(fid, fiff_blocks, fiff_bad_chs)
    artef_removal = _read_fiff_tag(fid, fiff_blocks, fiff_artef_removal)
    coord_trans = _read_fiff_tag(fid, fiff_blocks, fiff_coord_trans)
    coord_trans = _read_fiff_tag(fid, fiff_blocks, fiff_coord_trans)
    highpass = round.(Float64.(_read_fiff_tag(fid, fiff_blocks, fiff_highpass)), digits=4)
    hpi_bad_ch = _read_fiff_tag(fid, fiff_blocks, fiff_hpi_bad_chs)
    hpi_corr_coeff = _read_fiff_tag(fid, fiff_blocks, fiff_hpi_corr_coeff)
    event_comment = _read_fiff_tag(fid, fiff_blocks, fiff_event_comment)
    no_sample = _read_fiff_tag(fid, fiff_blocks, fiff_no_samples)
    first_time = _read_fiff_tag(fid, fiff_blocks, fiff_first_time)
    subave_size = _read_fiff_tag(fid, fiff_blocks, fiff_subave_size)
    subave_first = _read_fiff_tag(fid, fiff_blocks, fiff_subave_first)
    name = _read_fiff_tag(fid, fiff_blocks, fiff_name)
    dig_string = _read_fiff_tag(fid, fiff_blocks, fiff_dig_string)
    line_freq = Float64.(_read_fiff_tag(fid, fiff_blocks, fiff_line_freq))
    hpi_coil_freq = _read_fiff_tag(fid, fiff_blocks, fiff_hpi_coil_freq)
    signal_channel = _read_fiff_tag(fid, fiff_blocks, fiff_signal_channel)
    hpi_coil_moments = _read_fiff_tag(fid, fiff_blocks, fiff_hpi_coil_moments)
    hpi_fit_goodness = _read_fiff_tag(fid, fiff_blocks, fiff_hpi_fit_goodness)

    # patient data
    fiff_subj_id = 400
    fiff_subj_first_name = 401
    fiff_subj_last_name = 403
    patient = "$(_read_fiff_tag(fid, fiff_blocks, fiff_subj_id))"
    patient *= " " * _read_fiff_tag(fid, fiff_blocks, fiff_subj_first_name)
    patient *= " " * _read_fiff_tag(fid, fiff_blocks, fiff_subj_last_name)

    # ch_info_struct
    fiff_ch_info_struct = 30
    channels = findall(x -> x == fiff_ch_info_struct, tag_type)
    channels_struct = Vector{Any}()
    for channels_idx in 1:channel_n
        push!(channels_struct, _read_fiff_data(fid, fiff_blocks, channels[channels_idx]))
    end
    # identify signal type
    signal_type = ""
    channel_order = Int64.(_extract_struct(channels_struct, 2))
    l = length(string(maximum(channel_order)))
    labels = "MEG" .* lpad.(string.(channel_order), l, '0')
    range = Float64.(_extract_struct(channels_struct, 4))
    cal = Float64.(_extract_struct(channels_struct, 5))
    unit = Int64.(_extract_struct(channels_struct, 19))
    units = Vector{String}()
    for idx in 1:length(unit)
        unit[idx] == 107 && push!(units, "V")
        unit[idx] == 112 && push!(units, "T")
        unit[idx] == 201 && push!(units, "T/m")
        unit[idx] == -1 && push!(units, "")
    end
    unit_mul = Int64.(_extract_struct(channels_struct, 20))
    channel_types = _extract_struct(channels_struct, 3)
    1 in channel_types && (signal_type = "meg")
    signal_type != "meg" && throw(ArgumentError("Could not identify signal type as MEG."))
    # 2 in channel_types && (signal_type = "eeg")
    # identify channel types
    channel_type = _fiff_channel_type(channel_types)
    # identify gradiometers and magnetometers
    coils = _extract_struct(channels_struct, 6)
    gradiometers = Vector{Int64}()
    gradiometers_planar = Vector{Int64}()
    gradiometers_axial = Vector{Int64}()
    magnetometers = Vector{Int64}()
    for channel_idx in 1:channel_n
        push!(labels, string(channel_idx))
        coils[channel_idx] in [2001, 3011, 3012, 3013, 3014, 4002] && push!(gradiometers, channel_idx)
        coils[channel_idx] in [3011, 3012, 3013, 3014] && push!(gradiometers_planar, channel_idx)
        coils[channel_idx] in [2001, 5001] && push!(gradiometers_axial, channel_idx)
        coils[channel_idx] in [2000, 3021, 3022, 3023, 3024, 4001] && push!(magnetometers, channel_idx)
        coils[channel_idx] in [2001, 3011, 3012, 3013, 3014, 4002] &&(channel_type[channel_idx] = "grad")
        coils[channel_idx] in [2000, 3021, 3022, 3023, 3024, 4001] &&(channel_type[channel_idx] = "mag")
    end

    # # events
    _info("Events are not supported yet.")
    # # event channel numbers
    # fiff_event_channels = 600
    # event_channel = _read_fiff_data(fid, fiff_blocks, fiff_event_channels)
    # fiff_event_list = 601
    # # 3 integers per event: [number of samples, before, after]
    # event_list = _read_fiff_data(fid, fiff_blocks, fiff_event_list)
    # # event channel name
    # fiff_event_channel_name = 602
    # event_channel = _read_fiff_data(fid, fiff_blocks, fiff_event_channel_name)
    # # event bits array describing transition, 4 integers: [from_mask, from_state, to_mask, to_state]
    # fiff_event_bits = 603
    # event_bits = _read_fiff_data(fid, fiff_blocks, fiff_event_bits)

    # data
    # buffer containing measurement data
    fiff_data_buffer = 300
    data_buffer = _read_fiff_tag(fid, fiff_blocks, fiff_data_buffer)
    length(data_buffer) == 0 && throw(ArgumentError("Only raw data import is supported now."))
    # data skip in buffers
    fiff_data_skip = 301
    data_skip = _read_fiff_tag(fid, fiff_blocks, fiff_data_skip)
    data_skip !== nothing && _verbose("data_skip is not supported yet.")
    # # buffer containing one epoch and channel
    # fiff_epoch = 302
    # epoch = _read_fiff_tag(fid, fiff_blocks, fiff_epoch)
    # # data skip in samples
    # fiff_data_skip_smp = 303
    # data_skip_smp = _read_fiff_tag(fid, fiff_blocks, fiff_data_skip_smp)
    # # data buffer with int32 time channel
    # fiff_data_buffer2 = 304
    # data_buffer2 = _read_fiff_tag(fid, fiff_blocks, fiff_data_buffer2)
    # # data buffer with int32 time channel
    # fiff_time_stamp = 305
    # time_stamp = _read_fiff_tag(fid, fiff_blocks, fiff_time_stamp)

    # locs
    _info("Channel locations are not supported yet.")

    close(fid)

    # which data buffer to choose? 
    meg_signals = zeros(channel_n, 0)
    @inbounds for idx1 in 1:length(data_buffer)
        meg_signals = hcat(meg_signals, Float64.(reshape(data_buffer[idx1], channel_n, length(data_buffer[idx1]) ÷ channel_n)))
    end
    meg_signals = meg_signals .* range .* cal
    # for channel_idx in 1:channel_n
    #     if units[channel_idx] == "T"
    #          meg_signals[channel_idx, :, :] = meg_signals[channel_idx, :, :] .* 10^15
    #          units[channel_idx] = "fT"
    #     end
    #     if units[channel_idx] == "T/m"
    #          meg_signals[channel_idx, :, :] = meg_signals[channel_idx, :, :] .* (10^15 / 100)
    #          units[channel_idx] = "fT/cm"
    #     end
    # end
    meg_signals = reshape(meg_signals, size(meg_signals, 1), size(meg_signals, 2), 1)

    # create signal details
    duration_samples = size(meg_signals, 2)
    duration_seconds = size(meg_signals, 2) / sampling_rate
    meg_time = collect(0:(1 / sampling_rate):duration_seconds)
    meg_time = meg_time[1:end - 1]
    filesize_mb = round(filesize(file_name) / 1024^2, digits=2)

    meg_header = Dict(:signal_type => signal_type,
                      :filename => file_name,
                      :filesize_mb => filesize_mb,
                      :filetype => "FIFF",
                      :patient => patient,
                      :recording => description,
                      :recording_date => string(Dates.day(date)) * "-" * string(Dates.month(date)) * "-" * string(Dates.year(date)),
                      :recording_time => string(Dates.hour(date)) * ":" * string(Dates.minute(date)) * ":" * string(Dates.second(date)),
                      :channel_n => channel_n,
                      :channel_type => channel_type,
                      :reference => "",
                      :channel_locations => true,
                      :history => String[],
                      :components => Symbol[],
                      :coil => coils,
                      :gradiometers => gradiometers,
                      :gradiometers_planar => gradiometers_planar,
                      :gradiometers_axial => gradiometers_axial,
                      :magnetometers => magnetometers,
                      :meg_duration_samples => duration_samples,
                      :meg_duration_seconds => duration_seconds,
                      :epoch_n => 1,
                      :epoch_duration_samples => duration_samples,
                      :epoch_duration_seconds => duration_seconds,
                      :labels => labels,
                      :units => units,
                      :prefiltering => repeat(["LP: $lowpass Hz; HP: $highpass Hz"], channel_n),
                      :sampling_rate => sampling_rate,
                      :note => "",
                      :markers => false)

    meg_components = Vector{Any}()
    meg_epoch_time = meg_time
    meg_markers = DataFrame()
    meg_locs = DataFrame(:channel => Int64,
                         :labels => String[],
                         :loc_theta => Float64[],
                         :loc_radius => Float64[],
                         :loc_x => Float64[],
                         :loc_y => Float64[],
                         :loc_z => Float64[],
                         :loc_radius_sph => Float64[],
                         :loc_theta_sph => Float64[],
                         :loc_phi_sph => Float64[])

    meg = NeuroAnalyzer.MEG(meg_header, meg_time, meg_epoch_time, meg_signals[:, :, :], meg_components, meg_markers, meg_locs)
    return meg
end
