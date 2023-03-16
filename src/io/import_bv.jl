export import_bv
"""
    import_bv(file_name; detect_type)

Load BrainVision BVCDF file and return `NeuroAnalyzer.NEURO` object. At least two files are required: .vhdr (header) and .eeg (signal data). If available, markers are loaded from .vmrk file.

# Arguments

- `file_name::String`: name of the file to load, should point to .vhdr file.
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `::NeuroAnalyzer.NEURO`
"""
function import_bv(file_name::String; detect_type::Bool=true)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    splitext(file_name)[2] == ".vhdr" || throw(ArgumentError("file_name must specify .VHDR file."))

    vhdr = readlines(file_name)
    startswith(lowercase(replace(vhdr[1], " " => "")), "brainvision") == false && throw(ArgumentError("This is not a BrainVision .VHDR file."))

    file_type = "BrainVision"

    # delete comments
    for idx in length(vhdr):-1:1
        startswith(vhdr[idx], ';') && deleteat!(vhdr, idx)
    end

    # parse header
    eeg_file = ""
    marker_file = ""
    data_format = ""
    data_orientation = ""
    ch_n = 0
    sampling_interval = 0
    binary_format = ""
    averaged = false
    averaged_segments = 0
    averaged_points = 0
    segmentation = false
    channels_idx = 0
    locs_idx = 0
    for idx in eachindex(vhdr)
        startswith(lowercase(replace(vhdr[idx], " " => "")), "datafile=") && (eeg_file = split(vhdr[idx], '=')[2])
        replace(eeg_file, raw"$b" => split(file_name)[1])
        startswith(lowercase(replace(vhdr[idx], " " => "")), "markerfile=") && (marker_file = split(vhdr[idx], '=')[2])
        replace(marker_file, raw"$b" => split(file_name)[1])
        startswith(lowercase(replace(vhdr[idx], " " => "")), "dataformat=") && (data_format = lowercase(split(vhdr[idx], '=')[2])) # BINARY or ASCII
        startswith(lowercase(replace(vhdr[idx], " " => "")), "numberofchannels=") && (ch_n = parse(Int64, split(vhdr[idx], '=')[2])) # 32
        startswith(lowercase(replace(vhdr[idx], " " => "")), "dataorientation=") && (data_orientation = lowercase(split(vhdr[idx], '=')[2])) # MULTIPLEXED
        startswith(lowercase(replace(vhdr[idx], " " => "")), "samplinginterval=") && (sampling_interval = parse(Float64, split(vhdr[idx], '=')[2])) # 1000
        startswith(lowercase(replace(vhdr[idx], " " => "")), "binaryformat=") && (binary_format = lowercase(split(vhdr[idx], '=')[2])) # INT_16
        startswith(lowercase(replace(vhdr[idx], " " => "")), "averaged=") && (averaged = lowercase(split(vhdr[idx], '=')[2]) == "yes" ? true : false) # YES|NO
        startswith(lowercase(replace(vhdr[idx], " " => "")), "averagedsegments=") && (averaged_segments = parse(Int64, split(vhdr[idx], '=')[2]))
        startswith(lowercase(replace(vhdr[idx], " " => "")), "averageddatapoints=") && (averaged_points = parse(Int64, split(vhdr[idx], '=')[2]))
        startswith(lowercase(replace(vhdr[idx], " " => "")), "segmentation=") && (segmentation = lowercase(split(vhdr[idx], '=')[2]) == "markerbased" ? true : false) # YES|NO
        startswith(lowercase(replace(vhdr[idx], " " => "")), "[channelinfos]") && (channels_idx = idx)
        startswith(lowercase(replace(vhdr[idx], " " => "")), "[coordinates]") && (locs_idx = idx)
        startswith(lowercase(replace(vhdr[idx], " " => "")), "softwarefilters") && _info("Software filters are not supported yet.")
    end

    patient = ""
    recording = ""
    recording_date = ""
    recording_time = ""
    transducers = repeat([""], ch_n)
    units = repeat([""], ch_n)
    gain = repeat([1.0], ch_n)
    prefiltering = repeat([""], ch_n)

    clabels = repeat([""], ch_n)
    for idx in 1:ch_n
        tmp = split(split(vhdr[idx + channels_idx], '=')[2], ',')
        # channel label
        clabels[idx] = replace(split(split(vhdr[idx + channels_idx], '=')[2], ',')[1], "\1" => ",")
        # reference channel name
        # split(split(vhdr[idx + channels_idx], '=')[2], ',')[2]
        # resolution in units
        length(tmp) >= 3 && (gain[idx] = parse(Float64, split(split(vhdr[idx + channels_idx], '=')[2], ',')[3]))
        # units name, e.g. μV
        length(tmp) >= 4 && (units[idx] = split(split(vhdr[idx + channels_idx], '=')[2], ',')[4])
    end
    clabels = _clean_labels(clabels)
    if detect_type == true
        channel_type = _set_channel_types(clabels)
    else
        channel_type = repeat(["???"], ch_n)
    end
    channel_order = _sort_channels(copy(channel_type))

    # read locs
    loc_theta = zeros(ch_n)
    loc_radius = zeros(ch_n)
    loc_x = zeros(ch_n)
    loc_y = zeros(ch_n)
    loc_z = zeros(ch_n)
    loc_radius_sph = zeros(ch_n)
    loc_theta_sph = zeros(ch_n)
    loc_phi_sph = zeros(ch_n)
    if locs_idx != 0
        channel_locations = true
        for idx in 1:ch_n
            loc_radius_sph[idx] = parse(Float64, split(vhdr[locs_idx + idx], '=')[1])
            loc_theta_sph[idx] = parse(Float64, split(vhdr[locs_idx + idx], '=')[2])
            loc_phi_sph[idx] = parse(Float64, split(vhdr[locs_idx + idx], '=')[3])
            loc_theta[idx] = loc_theta_sph[idx]
            loc_radius[idx] = loc_radius_sph[idx]
            loc_x[idx], loc_y[idx], loc_z[idx] = sph2cart(loc_radius_sph[idx], loc_theta_sph[idx], loc_phi_sph[idx])
        end
    else
        channel_locations = false
    end

    # read markers
    if marker_file != ""
        has_markers = true
        if file_name != basename(file_name)
            marker_file = dirname(file_name) * "/" * marker_file
        else
            marker_file = marker_file
        end
        isfile(marker_file) || throw(ArgumentError("File $marker_file cannot be loaded."))
        vmrk = readlines(marker_file)
        # delete comments
        for idx in length(vmrk):-1:1
            startswith(vmrk[idx], ';') && deleteat!(vmrk, idx)
        end
        startswith(lowercase(replace(vmrk[1], " " => "")), "brainvision") == false && throw(ArgumentError("This is not a BrainVision .VMRK file."))
        markers_idx = 0
        for idx in eachindex(vmrk)
            startswith(lowercase(replace(vmrk[idx], " " => "")), "[markerinfos]") && (markers_idx = idx)
        end
        markers = repeat([""], length(vmrk) - markers_idx)
        for idx in eachindex(markers)
            markers[idx] = vmrk[markers_idx + idx]
        end
        # remove non-markers
        for idx in length(markers):-1:1
            startswith(lowercase(markers[idx]), "mk") == false && deleteat!(markers, idx)
        end
        m_id = repeat([""], length(markers))
        m_desc = repeat([""], length(markers))
        m_pos = zeros(Int64, length(markers))
        m_len = zeros(Int64, length(markers))
        m_ch = zeros(Int64, length(markers))
        for idx in eachindex(markers)
            m_id[idx] = replace(split(split(markers[idx], '=')[2], ',')[1], "\1" => ",")
            m_desc[idx] = replace(split(split(markers[idx], '=')[2], ',')[2], "\1" => ",")
            m_pos[idx] = parse(Int64, split(split(markers[idx], '=')[2], ',')[3])
            m_len[idx] = parse(Int64, split(split(markers[idx], '=')[2], ',')[4])
            # 0 = marker is related to all channels
            m_ch[idx] = parse(Int64, split(split(markers[idx], '=')[2], ',')[5])
        end
        markers = DataFrame(:id=>m_id, :start=>m_pos, :length=>m_len, :description=>m_desc, :channel=>m_ch)
    else
        has_markers = false
        markers = DataFrame(:id=>String[], :start=>Int64[], :length=>Int64[], :description=>String[], :channel=>Int64[])
    end

    # sampling_interval in μs to sampling rate in Hz
    sampling_rate = round(Int64, 1 / (sampling_interval / 10^6))

    # read data
    if file_name != basename(file_name)
        eeg_file = dirname(file_name) * "/" * eeg_file
    else
        eeg_file = eeg_file
    end

    isfile(eeg_file) || throw(ArgumentError("File $eeg_file cannot be loaded."))
    if data_format == "binary"
        if binary_format == "int_16"
            bytes = 2
        elseif binary_format == "ieee_float_32"
            bytes = 4
        else
            @error("Only Float32 and Int16 BVCDF binary formats are supported.")
        end

        fid = ""
        try
            fid = open(file_name, "r")
        catch
            error("File $file_name cannot be loaded.")
        end

        signal = zeros(filesize(eeg_file) ÷ bytes)
        for idx in 1:(filesize(eeg_file) ÷ bytes)
            buf = zeros(UInt8, bytes)
            readbytes!(fid, buf, bytes)
            if bytes == 4
                signal[idx] = Float64(reinterpret(Float32, buf)[1])
            else
                signal[idx] = Float64(reinterpret(Int16, buf)[1])
            end
        end
        close(fid)
        # split signal into channels
        if data_orientation == "multiplexed"
            data = zeros(ch_n, length(signal) ÷ ch_n, 1)
            idx2 = 1
            for idx1 in 1:ch_n:length(signal)
                data[:, idx2, 1] = signal[idx1:(idx1 + (ch_n - 1))]
                idx2 += 1
            end
        else
            @error "Only MULTIPLEXED data orientation is supported."
        end
    else
        @error "ASCII format is not supported yet."
    end

    time_pts = collect(0:(1 / sampling_rate):((size(data, 2) * size(data, 3)) / sampling_rate))
    time_pts = round.(time_pts[1:end - 1], digits=3)
    epoch_time = time_pts
    file_size_mb = round(filesize(eeg_file) / 1024^2, digits=2)

    data_type = "eeg"

    s = _create_subject(id="",
                        first_name="",
                        middle_name="",
                        last_name=string(patient),
                        handedness="",
                        weight=-1,
                        height=-1)

    r = _create_recording_eeg(data_type=data_type,
                              file_name=file_name,
                              file_size_mb=file_size_mb,
                              file_type=file_type,
                              recording=string(recording),
                              recording_date=recording_date,
                              recording_time=recording_time,
                              recording_notes="",
                              channel_type=channel_type[channel_order],
                              reference="",
                              clabels=clabels[channel_order],
                              transducers=transducers[channel_order],
                              units=units[channel_order],
                              prefiltering=prefiltering[channel_order],
                              sampling_rate=sampling_rate,
                              gain=gain[channel_order])

    e = _create_experiment(experiment_name="",
                           experiment_notes="",
                           experiment_design="")

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = String[]

    if channel_locations == false
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
    else
        locs = DataFrame(:ch_n=>1:ch_n,
                             :labels=>labels,
                             :loc_theta=>loc_theta,
                             :loc_radius=>loc_radius,
                             :loc_x=>loc_x,
                             :loc_y=>loc_y,
                             :loc_z=>loc_z,
                             :loc_radius_sph=>loc_radius_sph,
                             :loc_theta_sph=>loc_theta_sph,
                             :loc_phi_sph=>loc_phi_sph)
    end

    return NeuroAnalyzer.NEURO(hdr, time_pts, epoch_time, data[channel_order, :, :], components, markers, locs, history)
    
end
