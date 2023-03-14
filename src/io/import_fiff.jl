export import_fiff

"""
    import_fiff(file_name; detect_type)

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

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    fid, fiff_block = _create_fiff_block(file_name)

    # # process blocks
    # for block_idx1 in unique(block_number)
    #     block_current = fiff_block[fiff_block[:, 5] .== block_idx1, :]
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

    # read data from tags
    ch_n = round.(Int64, _read_fiff_tag(fid, fiff_block, fiff_nchan))
    sampling_rate = round.(Int64, _read_fiff_tag(fid, fiff_block, fiff_sfreq))
    data_pack = round.(Int64, _read_fiff_tag(fid, fiff_block, fiff_data_pack))
    date = _read_fiff_tag(fid, fiff_block, fiff_meas_date)
    date = unix2datetime(date[1][2])
    description = _read_fiff_tag(fid, fiff_block, fiff_description)
    nave = _read_fiff_tag(fid, fiff_block, fiff_nave)
    first_sample = _read_fiff_tag(fid, fiff_block, fiff_first_sample)
    last_sample = _read_fiff_tag(fid, fiff_block, fiff_last_sample)
    aspect_kind = _read_fiff_tag(fid, fiff_block, fiff_aspect_kind)
    ref_event = _read_fiff_tag(fid, fiff_block, fiff_ref_event)
    experimenter = _read_fiff_tag(fid, fiff_block, fiff_experimenter)
    dig_point = _read_fiff_tag(fid, fiff_block, fiff_dig_point)
    hpi_slopes = _read_fiff_tag(fid, fiff_block, fiff_hpi_slopes)
    hpi_ncoil = _read_fiff_tag(fid, fiff_block, fiff_hpi_ncoil)
    req_event = _read_fiff_tag(fid, fiff_block, fiff_req_event)
    req_limit = _read_fiff_tag(fid, fiff_block, fiff_req_limit)
    req_limit = _read_fiff_tag(fid, fiff_block, fiff_req_limit)
    lowpass = round.(Float64.(_read_fiff_tag(fid, fiff_block, fiff_lowpass)), digits=4)
    bad_chs = _read_fiff_tag(fid, fiff_block, fiff_bad_chs)
    artef_removal = _read_fiff_tag(fid, fiff_block, fiff_artef_removal)
    coord_trans = _read_fiff_tag(fid, fiff_block, fiff_coord_trans)
    coord_trans = _read_fiff_tag(fid, fiff_block, fiff_coord_trans)
    highpass = round.(Float64.(_read_fiff_tag(fid, fiff_block, fiff_highpass)), digits=4)
    hpi_bad_ch = _read_fiff_tag(fid, fiff_block, fiff_hpi_bad_chs)
    hpi_corr_coeff = _read_fiff_tag(fid, fiff_block, fiff_hpi_corr_coeff)
    event_comment = _read_fiff_tag(fid, fiff_block, fiff_event_comment)
    no_sample = _read_fiff_tag(fid, fiff_block, fiff_no_samples)
    first_time = _read_fiff_tag(fid, fiff_block, fiff_first_time)
    subave_size = _read_fiff_tag(fid, fiff_block, fiff_subave_size)
    subave_first = _read_fiff_tag(fid, fiff_block, fiff_subave_first)
    name = _read_fiff_tag(fid, fiff_block, fiff_name)
    dig_string = _read_fiff_tag(fid, fiff_block, fiff_dig_string)
    line_freq = Float64.(_read_fiff_tag(fid, fiff_block, fiff_line_freq))
    hpi_coil_freq = _read_fiff_tag(fid, fiff_block, fiff_hpi_coil_freq)
    signal_channel = _read_fiff_tag(fid, fiff_block, fiff_signal_channel)
    hpi_coil_moments = _read_fiff_tag(fid, fiff_block, fiff_hpi_coil_moments)
    hpi_fit_goodness = _read_fiff_tag(fid, fiff_block, fiff_hpi_fit_goodness)

    # patient data
    patient = "$(_read_fiff_tag(fid, fiff_block, fiff_subj_id))"
    patient *= " " * _read_fiff_tag(fid, fiff_block, fiff_subj_first_name)
    patient *= " " * _read_fiff_tag(fid, fiff_block, fiff_subj_last_name)

    # ch_info_struct
    channels = findall(x -> x == fiff_ch_info_struct, fiff_block[:, 3])
    channels_struct = Vector{Any}()
    for channels_idx in 1:ch_n
        push!(channels_struct, _read_fiff_data(fid, fiff_block, channels[channels_idx]))
    end
    # identify signal type
    data_type = ""
    channel_order = Int64.(_extract_struct(channels_struct, 2))
    l = length(string(maximum(channel_order)))
    clabels = "MEG" .* lpad.(string.(channel_order), l, '0')
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
    1 in channel_types && (data_type = "meg")
    data_type != "meg" && throw(ArgumentError("Could not identify signal type as MEG."))
    2 in channel_types && (data_type = "eeg")
    # identify channel types
    channel_type = _fiff_channel_type(channel_types)
    # identify gradiometers and magnetometers
    coils = Int64.(_extract_struct(channels_struct, 6))
    gradiometers = Vector{Int64}()
    gradiometers_planar = Vector{Int64}()
    gradiometers_axial = Vector{Int64}()
    magnetometers = Vector{Int64}()
    for ch_idx in 1:ch_n
        push!(clabels, string(ch_idx))
        coils[ch_idx] in [2001, 3011, 3012, 3013, 3014, 4002] && push!(gradiometers, ch_idx)
        coils[ch_idx] in [3011, 3012, 3013, 3014] && push!(gradiometers_planar, ch_idx)
        coils[ch_idx] in [2001, 5001] && push!(gradiometers_axial, ch_idx)
        coils[ch_idx] in [2000, 3021, 3022, 3023, 3024, 4001] && push!(magnetometers, ch_idx)
        coils[ch_idx] in [2001, 3011, 3012, 3013, 3014, 4002] && (channel_type[ch_idx] = "grad")
        coils[ch_idx] in [2000, 3021, 3022, 3023, 3024, 4001] && (channel_type[ch_idx] = "mag")
    end

    # # events
    _info("Events are not supported yet.")
    # # event channel numbers
    # event_channel = _read_fiff_data(fid, fiff_block, fiff_event_channels)
    # # 3 integers per event: [number of samples, before, after]
    # event_list = _read_fiff_data(fid, fiff_block, fiff_event_list)
    # # event channel name
    # event_channel = _read_fiff_data(fid, fiff_block, fiff_event_ch_name)
    # # event bits array describing transition, 4 integers: [from_mask, from_state, to_mask, to_state]
    # event_bits = _read_fiff_data(fid, fiff_block, fiff_event_bits)

    # data
    # buffer containing measurement data
    data_buffer = _read_fiff_tag(fid, fiff_block, fiff_data_buffer)
    length(data_buffer) == 0 && throw(ArgumentError("Only raw data import is supported now."))
    # data skip in buffers
    data_skip = _read_fiff_tag(fid, fiff_block, fiff_data_skip)
    data_skip !== nothing && _info("data_skip is not supported yet.")
    # # buffer containing one epoch and channel
    # epoch = _read_fiff_tag(fid, fiff_block, fiff_epoch)
    # # data skip in samples
    # data_skip_smp = _read_fiff_tag(fid, fiff_block, fiff_data_skip_smp)
    # # data buffer with int32 time channel
    # data_buffer2 = _read_fiff_tag(fid, fiff_block, fiff_data_buffer2)
    # # data buffer with int32 time channel
    # time_stamp = _read_fiff_tag(fid, fiff_block, fiff_time_stamp)

    # locs
    _info("Channel locations are not supported yet.")

    close(fid)

    # which data buffer to choose? 
    data = zeros(ch_n, 0)
    @inbounds for idx1 in 1:length(data_buffer)
        data = hcat(data, Float64.(reshape(data_buffer[idx1], ch_n, length(data_buffer[idx1]) ÷ ch_n)))
    end
    data = data .* range .* cal
    # for ch_idx in 1:ch_n
    #     if units[ch_idx] == "T"
    #          data[ch_idx, :, :] = data[ch_idx, :, :] .* 10^15
    #          units[ch_idx] = "fT"
    #     end
    #     if units[ch_idx] == "T/m"
    #          data[ch_idx, :, :] = data[ch_idx, :, :] .* (10^15 / 100)
    #          units[ch_idx] = "fT/cm"
    #     end
    # end
    data = reshape(data, size(data, 1), size(data, 2), 1)

    # create signal details
    time_pts = collect(0:(1 / sampling_rate):((size(data, 2) * size(data, 3)) / sampling_rate))
    time_pts = round.(time_pts[1:end - 1], digits=3)
    epoch_time = time_pts
    file_size_mb = round(filesize(file_name) / 1024^2, digits=2)

    s = _create_subject(id="",
                        first_name="",
                        middle_name="",
                        last_name=string(patient),
                        handedness="",
                        weight=-1,
                        height=-1)

    r = _create_recording_meg(data_type=data_type,
                              file_name=file_name,
                              file_size_mb=file_size_mb,
                              file_type="FIFF",
                              recording=string(description),
                              recording_date=string(Dates.day(date)) * "-" * string(Dates.month(date)) * "-" * string(Dates.year(date)),
                              recording_time=string(Dates.hour(date)) * ":" * string(Dates.minute(date)) * ":" * string(Dates.second(date)),
                              recording_notes="",
                              channel_type=channel_type,
                              reference="",
                              clabels=clabels,
                              units=units,
                              prefiltering=repeat(["LP: $lowpass Hz; HP: $highpass Hz"], ch_n),
                              sampling_rate=sampling_rate,
                              magnetometers=magnetometers,
                              gradiometers=gradiometers,
                              gradiometers_planar=gradiometers_planar,
                              gradiometers_axial=gradiometers_axial,
                              coils=coils)
    e = _create_experiment(experiment_name="",
                           experiment_notes="",
                           experiment_design="")

    hdr = _create_header(s,
                         r,
                         e,
                         component_names=Symbol[],
                         has_markers=has_markers,
                         has_locs=false,
                         history=String[])

    components = Vector{Any}()

    markers = DataFrame()

    locs = DataFrame(:channel=>Int64,
                     :labels=>String[],
                     :loc_theta=>Float64[],
                     :loc_radius=>Float64[],
                     :loc_x=>Float64[],
                     :loc_y=>Float64[],
                     :loc_z=>Float64[],
                     :loc_radius_sph=>Float64[],
                     :loc_theta_sph=>Float64[],
                     :loc_phi_sph=>Float64[])

    return NeuroAnalyzer.NEURO(hdr, time_pts, epoch_time, data, components, markers, locs)
end
