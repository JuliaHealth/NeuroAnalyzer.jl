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

    file_name = "test/meg-test-fif.fif"

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    fid = ""
    try
        fid = open(file_name, "r")
    catch
        error("File $file_name cannot be loaded.")
    end

    # read file_id tag
    tag_kind = nothing
    tag_type = nothing
    tag_size = nothing
    tag_next = nothing
    try
        tag_kind, tag_type, tag_size, tag_next = _read_fif_tag(fid)
    catch
        error("File $file_name first tag cannot be read.")
    end

    # check file_id tag
    fiff_file_id = 100
    fiff_id_struct = 31
    fiff_next_seq = 0
    tag_kind != fiff_file_id && throw(ArgumentError("File $file_name is not a FIFF file."))
    tag_type != fiff_id_struct && throw(ArgumentError("File $file_name is not a FIFF file."))
    tag_size != 20 && throw(ArgumentError("File $file_name is not a FIFF file."))

    # read dir_pointer tag
    try
        tag_kind, tag_type, tag_size, tag_next = _read_fif_tag(fid)
    catch
        error("File $file_name dir_pointer tag cannot be read.")
    end

    # check dir_pointer tag
    fiff_dir_pointer = 101
    tag_kind != fiff_dir_pointer && throw(ArgumentError("File $file_name has no dir_pointer tag."))

    # read tags
    seek(fid, 0)
    # tags: position in file, tag_id, data_type, data_size, next
    tags = Vector{Tuple{Int64, Int64, Int64, Int64, Int64}}()
    while tag_next != -1
        current_position = position(fid)
        tag_kind, tag_type, tag_size, tag_next = _read_fif_tag(fid)
        push!(tags, (current_position, tag_kind, tag_type, tag_size, tag_next))
    end
    seek(fid, 0)

    # create list of tag IDs
    tag_ids = Vector{Int64}()
    for tag_idx in 1:length(tags)
        push!(tag_ids, tags[tag_idx][2])
    end
    tag_types = Vector{Int64}()
    for tag_idx in 1:length(tags)
        push!(tag_types, tags[tag_idx][3])
    end

    # channel_n
    fiff_nchan_id = 200
    id = _find_fif_tag(tag_ids, fiff_nchan_id)
    channel_n = _read_fif_data(fid, tags, tag_ids, id)[]

    # sr
    fiff_sfreq_id = 201
    id = _find_fif_tag(tag_ids, fiff_sfreq_id)
    sampling_rate = Int64.(_read_fif_data(fid, tags, tag_ids, id)[])

    # date
    fiff_meas_date_id = 204
    id = _find_fif_tag(tag_ids, fiff_meas_date_id)
    date = _read_fif_data(fid, tags, tag_ids, id)
    date = unix2datetime(date[2])

    # data order
    fiff_data_pack_id = 202
    id = _find_fif_tag(tag_ids, fiff_data_pack_id)
    ord = _read_fif_data(fid, tags, tag_ids, id)

    # patient data
    fiff_subj_first_name_id = 401
    id = _find_fif_tag(tag_ids, fiff_subj_first_name_id)
    patient = _read_fif_data(fid, tags, tag_ids, id)
    fiff_subj_last_name_id = 403
    id = _find_fif_tag(tag_ids, fiff_subj_last_name_id)
    patient *= " " * _read_fif_data(fid, tags, tag_ids, id)

    # ch_info_struct
    fiff_ch_info_struct_id = 30
    channels = findall(x -> x == fiff_ch_info_struct_id, tag_types)
    channels_struct = Vector{Any}()
    for channels_idx in 1:channel_n
        push!(channels_struct, _read_fif_data(fid, tags, tag_ids, channels[channels_idx]))
    end
    # identify signal type
    signal_type = ""
    channel_order = Int64.(_extract_struct(channels_struct, 2))
    channel_types = _extract_struct(channels_struct, 3)
    1 in channel_types && (signal_type = "meg")
    signal_type != "meg" && throw(ArgumentError("Could not identify signal type as MEG."))
    # 2 in channel_types && (signal_type = "eeg")
    # identify channel types
    channel_type = _fif_channel_type(channel_types)
    # identify gradiometers and magnetometers
    coils = _extract_struct(channels_struct, 6)
    gradiometers = Vector{Int64}()
    gradiometers_planar = Vector{Int64}()
    gradiometers_axial = Vector{Int64}()
    magnetometers = Vector{Int64}()
    for channel_idx in 1:channel_n
        coils[channel_idx] in [2001, 3011, 3012, 3013, 3014, 4002] && push!(gradiometers, channel_idx)
        coils[channel_idx] in [3011, 3012, 3013, 3014] && push!(gradiometers_planar, channel_idx)
        coils[channel_idx] in [2001, 5001] && push!(gradiometers_axial, channel_idx)
        coils[channel_idx] in [2000, 3021, 3022, 3023, 3024, 4001] && push!(magnetometers, channel_idx)
    end

    # events
    # event channel numbers
    fiff_event_channels_id = 600
    id = _find_fif_tag(tag_ids, fiff_event_channels_id)
    event_channel = _read_fif_data(fid, tags, tag_ids, id)
    fiff_event_list_id = 601
    id = _find_fif_tag(tag_ids, fiff_event_list_id)
    # 3 integers per event: [number of samples, before, after]
    event_list = _read_fif_data(fid, tags, tag_ids, id)
    # event channel name
    fiff_event_channel_name_id = 602
    id = _find_fif_tag(tag_ids, fiff_event_channel_name_id)
    event_channel = _read_fif_data(fid, tags, tag_ids, id)
    # event bits array describing transition, 4 integers: [from_mask, from_state, to_mask, to_state]
    fiff_event_bits_id = 603
    id = _find_fif_tag(tag_ids, fiff_event_bits_id)
    event_bits = _read_fif_data(fid, tags, tag_ids, id)

    # data
    # type of data block
    fiff_raw_data_block_id = 102
    id = _find_fif_tag(tag_ids, fiff_raw_data_block_id)
    id !== nothing && (raw_data_block = _read_fif_data(fid, tags, tag_ids, id))
    fiff_processed_data_block_id = 103
    id = _find_fif_tag(tag_ids, fiff_processed_data_block_id)
    id !== nothing && (processed_data_block = _read_fif_data(fid, tags, tag_ids, id))

    # buffer containing measurement data
    fiff_data_buffer_id = 300
    id = _find_fif_tag(tag_ids, fiff_data_buffer_id)
    data_buffer = _read_fif_data(fid, tags, tag_ids, id)
    data_buffer = Float64.(reshape(data_buffer, channel_n, length(data_buffer) ÷ channel_n))
    signals = m_sort(data_buffer, channel_order)
    signals = reshape(signals, size(data_buffer, 1), size(data_buffer, 2), 1)
    # data skip in buffers
    fiff_data_skip_id = 301
    id = _find_fif_tag(tag_ids, fiff_data_skip_id)
    data_skip = _read_fif_data(fid, tags, tag_ids, id)
    # buffer containing one epoch and channel
    fiff_epoch_id = 302
    id = _find_fif_tag(tag_ids, fiff_epoch_id)
    epoch = _read_fif_data(fid, tags, tag_ids, id)
    # data skip in samples
    fiff_data_skip_smp_id = 303
    id = _find_fif_tag(tag_ids, fiff_data_skip_smp_id)
    data_skip_smp = _read_fif_data(fid, tags, tag_ids, id)
    # data buffer with int32 time channel
    fiff_data_buffer2_id = 304
    id = _find_fif_tag(tag_ids, fiff_data_buffer2_id)
    data_buffer2 = _read_fif_data(fid, tags, tag_ids, id)
    # data buffer with int32 time channel
    fiff_time_stamp_id = 305
    id = _find_fif_tag(tag_ids, fiff_time_stamp_id)
    data_buffer2 = _read_fif_data(fid, tags, tag_ids, id)

    close(fid)

    duration_samples = size(signals, 2)
    duration_seconds = size(signals, 2) / sampling_rate
    meg_time = collect(0:(1 / sampling_rate):duration_seconds)
    meg_time = meg_time[1:end - 1]
    filesize_mb = round(filesize(file_name) / 1024^2, digits=2)

    meg_header = Dict(:signal_type => signal_type,
                      :filename => file_name,
                      :filesize_mb => filesize_mb,
                      :filetype => "FIFF",
                      :patient => string(patient),
                      :recording => string(recording),
                      :recording_date => string(Dates.day(date)) * "-" * string(Dates.month(date)) * "-" * string(Dates.year(date)),
                      :recording_time => string(Dates.hour(date)) * ":" * string(Dates.minute(date)) * ":" * string(Dates.second(date)),
                      :channel_n => channel_n,
                      :channel_type => channel_type,
                      :reference => "",
                      :channel_locations => false,
                      :history => String[],
                      :components => Symbol[],
                      :meg_duration_samples => duration_samples,
                      :meg_duration_seconds => duration_seconds,
                      :epoch_n => 1,
                      :epoch_duration_samples => duration_samples,
                      :epoch_duration_seconds => duration_seconds,
                      :labels => labels[channel_order],
                      :prefiltering => prefiltering[channel_order],
                      :sampling_rate => sampling_rate,
                      :note => "",
                      :markers => has_markers)

    meg_components = Vector{Any}()
    meg_epoch_time = meg_time
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

    meg = NeuroAnalyzer.EEG(meg_header, meg_time, meg_epoch_time, meg_signals[:, :, :], meg_components, meg_markers, meg_locs)
    return meg
end
