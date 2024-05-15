export import_bv
"""
    import_bv(file_name; detect_type)

Load BrainVision BVCDF file and return `NeuroAnalyzer.NEURO` object. At least two files are required: .vhdr (header) and .eeg (signal data). If available, markers are loaded from .vmrk file.

# Arguments

- `file_name::String`: name of the file to load, should point to .vhdr file.
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function import_bv(file_name::String; detect_type::Bool=true)

    @assert isfile(file_name) "File $file_name cannot be loaded."
    @assert lowercase(splitext(file_name)[2]) in [".vhdr", ".ahdr"] "file_name must specify .VHDR/.AHDR file."

    vhdr = nothing
    try
        vhdr = readlines(file_name)
    catch
        @error "File $file_name cannot be loaded."
    end
    vhdr[1][1] == '\ufeff' && (vhdr[1] = vhdr[1][4:end])
    @assert startswith(lowercase(replace(vhdr[1], " " => "")), "brainvision") "File $file_name is not BrainVision .VHDR file."

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
    soft_filt_idx = 0

    vhdr = replace.(vhdr, " "=>"")

    any(startswith.(lowercase.(vhdr), "datafile=")) && (eeg_file = split(vhdr[startswith.(lowercase.(vhdr), "datafile=")][1], '=')[2])
    eeg_file = replace(eeg_file, raw"$b" => split(file_name)[1])
    any(startswith.(lowercase.(vhdr), "markerfile=")) && (marker_file = split(vhdr[startswith.(lowercase.(vhdr), "markerfile=")][1], '=')[2])
    marker_file = replace(marker_file, raw"$b" => split(file_name)[1])
    any(startswith.(lowercase.(vhdr), "dataformat=")) && (data_format = lowercase(split(vhdr[startswith.(lowercase.(vhdr), "dataformat=")][1], '=')[2])) # BINARY or ASCII
    any(startswith.(lowercase.(vhdr), "numberofchannels=")) && (ch_n = parse(Int64, split(vhdr[startswith.(lowercase.(vhdr), "numberofchannels=")][1], '=')[2]))
    any(startswith.(lowercase.(vhdr), "dataorientation=")) && (data_orientation = lowercase(split(vhdr[startswith.(lowercase.(vhdr), "dataorientation=")][1], '=')[2])) # MULTIPLEXED
    any(startswith.(lowercase.(vhdr), "samplinginterval=")) && (sampling_interval = parse(Float64, split(vhdr[startswith.(lowercase.(vhdr), "samplinginterval=")][1], '=')[2]))
    any(startswith.(lowercase.(vhdr), "binaryformat=")) && (binary_format = lowercase(split(vhdr[startswith.(lowercase.(vhdr), "binaryformat=")][1], '=')[2])) # IEEE_FLOAT_32 / INT_16
    any(startswith.(lowercase.(vhdr), "averaged=")) && (averaged = lowercase(split(vhdr[startswith.(lowercase.(vhdr), "averaged=")][1], '=')[2]) == "yes" ? true : false) # YES|NO
    any(startswith.(lowercase.(vhdr), "averagedsegments=")) && (averaged_segments = parse(Int64, split(vhdr[startswith.(lowercase.(vhdr), "averagedsegments=")][1], '=')[2]))
    any(startswith.(lowercase.(vhdr), "averageddatapoints=")) && (averaged_points = parse(Int64, split(vhdr[startswith.(lowercase.(vhdr), "averageddatapoints=")][1], '=')[2]))
    any(startswith.(lowercase.(vhdr), "segmentation=")) && (segmentation = lowercase(split(vhdr[startswith.(lowercase.(vhdr), "segmentation=")][1], '=')[2]) == "yes" ? true : false) # YES|NO

    for idx in eachindex(vhdr)
        startswith(lowercase(replace(vhdr[idx], " " => "")), "[channelinfos]") && (channels_idx = idx)
        startswith(lowercase(replace(vhdr[idx], " " => "")), "[coordinates]") && (locs_idx = idx)
        if startswith(lowercase(replace(vhdr[idx], " " => "")), "softwarefilters")
            (soft_filt_idx = idx)
        end
    end

    # software filters
    soft_filt = false
    if soft_filt_idx != 0
        if lowercase(vhdr[soft_filt_idx + 2]) != "disabled"
            _info("Embedded software filters are not implemented yet; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org")
        end
    end
    soft_filt_file = replace(splitext(file_name)[1], "eeg"=>"channels.tsv")
    isfile(soft_filt_file) && (soft_filt = CSV.read(soft_filt_file, stringtype=String, DataFrame))

    # JSON
    js_file =""
    if splitext(file_name)[2] == ".vhdr"
        js_file = replace(file_name, "vhdr"=>"json")
    elseif splitext(file_name)[2] == ".ahdr"
        js_file = replace(file_name, "ahdr"=>"json")
    end
    e_name = ""
    e_notes = ""
    r_notes = ""
    ref = ""
    if isfile(js_file)
        js = JSON.parsefile(js_file; dicttype=Dict, inttype=Int64, use_mmap=true)
        "TaskName" in keys(js) && (e_name = js["TaskName"])
        "TaskDescription" in keys(js) && (e_notes = js["TaskDescription"])
        "ManufacturersModelName" in keys(js) && (r_notes = js["ManufacturersModelName"])
        "EEGReference" in keys(js) && (ref = js["EEGReference"])
    end

    patient = ""
    recording = ""
    recording_date = ""
    recording_time = ""
    transducers = repeat([""], ch_n)
    units = repeat([""], ch_n)
    gain = ones(ch_n)
    prefiltering = repeat([""], ch_n)

    clabels = repeat([""], ch_n)
    ch_types = repeat([""], ch_n)
    ref_chs = repeat([""], ch_n)

    for idx in 1:ch_n
        tmp = split(split(vhdr[idx + channels_idx], '=')[2], ',')
        # channel label
        clabels[idx] = replace(split(split(vhdr[idx + channels_idx], '=')[2], ',')[1], "\1" => ",")
        # reference channel name
        ref_chs = split(split(vhdr[idx + channels_idx], '=')[2], ',')[2]
        # resolution in units
        if length(tmp) >= 3
            if split(split(vhdr[idx + channels_idx], '=')[2], ',')[3] != ""
                gain[idx] = parse(Float64, split(split(vhdr[idx + channels_idx], '=')[2], ',')[3])
            end
        end
        # units name, e.g. μV
        length(tmp) >= 4 && (units[idx] = split(split(vhdr[idx + channels_idx], '=')[2], ',')[4])
    end

    if soft_filt != false
        clabels = soft_filt[!, :name]
        units = soft_filt[!, :units]
        for idx in 1:ch_n
            lowercase(soft_filt[idx, :type]) in ["eeg", "eog", "ecg", "emg"] && (ch_type[idx] = lowercase(soft_filt[idx, :type]))
        end
        prefiltering = repeat(["LP: "], ch_n) .* string.(round.(soft_filt[!, :low_cutoff], digits=4)) .* repeat([" Hz, HP: "], ch_n) .* string.(round.(soft_filt[!, :high_cutoff], digits=4)) .* " Hz"
    end
    
    clabels = _clean_labels(clabels)
    if ch_types == repeat([""], ch_n)
        if detect_type
            ch_type = _set_channel_types(clabels, "eeg")
        else
            ch_type = repeat(["eeg"], ch_n)
        end
    end
    units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]
    channel_order = _sort_channels(ch_type)

    # read locs
    loc_theta = zeros(ch_n)
    loc_radius = zeros(ch_n)
    loc_x = zeros(ch_n)
    loc_y = zeros(ch_n)
    loc_z = zeros(ch_n)
    loc_radius_sph = zeros(ch_n)
    loc_theta_sph = zeros(ch_n)
    loc_phi_sph = zeros(ch_n)
    locs = _initialize_locs()
    if locs_idx != 0
        for idx in 1:ch_n
            if occursin('=', vhdr[locs_idx + idx])
                l = split(vhdr[locs_idx + idx], '=')[2]
            else
                l = vhdr[locs_idx + idx]
            end
            loc_radius_sph[idx] = parse(Float64, split(l, ',')[1])
            loc_theta_sph[idx] = parse(Float64, split(l, ',')[2])
            loc_phi_sph[idx] = parse(Float64, split(l, ',')[3])
            loc_theta[idx] = loc_theta_sph[idx]
            loc_radius[idx] = loc_radius_sph[idx]
            loc_x[idx], loc_y[idx], loc_z[idx] = sph2cart(loc_radius_sph[idx], loc_theta_sph[idx], loc_phi_sph[idx])
        end
        locs = DataFrame(:ch_n=>1:ch_n,
                         :labels=>clabels,
                         :loc_theta=>loc_theta,
                         :loc_radius=>loc_radius,
                         :loc_x=>loc_x,
                         :loc_y=>loc_y,
                         :loc_z=>loc_z,
                         :loc_radius_sph=>loc_radius_sph,
                         :loc_theta_sph=>loc_theta_sph,
                         :loc_phi_sph=>loc_phi_sph)
    end

    # sampling_interval in μs to sampling rate in Hz
    sampling_rate = round(Int64, 1 / (sampling_interval / 10^6))

    # read markers
    if marker_file != ""
        if file_name != basename(file_name)
            marker_file = dirname(file_name) * "/" * marker_file
        else
            marker_file = marker_file
        end
        @assert isfile(marker_file) "File $marker_file cannot be loaded."
        vmrk = readlines(marker_file)
        vmrk[1][1] == '\ufeff' && (vmrk[1] = vmrk[1][4:end])
        # delete comments
        for idx in length(vmrk):-1:1
            startswith(vmrk[idx], ';') && deleteat!(vmrk, idx)
        end
        @assert startswith(lowercase(replace(vmrk[1], " " => "")), "brainvision") "File $file_name is not BrainVision .VMRK file."
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
            !startswith(lowercase(markers[idx]), "mk") && deleteat!(markers, idx)
        end
        m_id = repeat([""], length(markers))
        m_desc = repeat(["marker"], length(markers))
        m_pos = zeros(Int64, length(markers))
        m_len = zeros(Int64, length(markers))
        m_ch = zeros(Int64, length(markers))
        for idx in eachindex(markers)
            m_id[idx] = replace(split(split(markers[idx], '=')[2], ',')[1], "\1" => ",")
            replace(split(split(markers[idx], '=')[2], ',')[2], "\1" => ",") != "" && (m_desc[idx] = replace(split(split(markers[idx], '=')[2], ',')[2], "\1" => ","))
            m_pos[idx] = parse(Int64, split(split(markers[idx], '=')[2], ',')[3])
            m_len[idx] = parse(Int64, split(split(markers[idx], '=')[2], ',')[4])
            # 0 = marker is related to all channels
            m_ch[idx] = parse(Int64, split(split(markers[idx], '=')[2], ',')[5])
        end
        markers = DataFrame(:id=>m_id, 
                            :start=>(m_pos ./ sampling_rate), 
                            :length=>round.(m_len ./ sampling_rate), 
                            :description=>m_desc, 
                            :channel=>m_ch)
        if markers[!, :id] == repeat([""], nrow(markers))
            markers[!, :id] == repeat(["mrk"], nrow(markers))
        end
    elseif isfile(replace(splitext(file_name)[1], "eeg"=>"events.tsv"))
        vmrk = CSV.read(replace(splitext(file_name)[1], "eeg"=>"events.tsv"), stringtype=String, DataFrame)
        markers = DataFrame(:id=>repeat(["mrk"], nrow(vmrk)), 
                            :start=>(vmrk[!, :sample] ./ sampling_rate), 
                            :length=>round.(vmrk[!, :duration] ./ sampling_rate), 
                            :description=>(vmrk[!, :trial_type] .* "_" .* string.(vmrk[!, :value])), 
                            :channel=>m_ch)
    else
        markers = DataFrame(:id=>String[], 
                            :start=>Float64[], 
                            :length=>Float64[], 
                            :description=>String[], 
                            :channel=>Int64[])
    end

    # read data
    if file_name != basename(file_name)
        eeg_file = dirname(file_name) * "/" * eeg_file
    else
        eeg_file = eeg_file
    end

    @assert isfile(eeg_file) "File $eeg_file cannot be loaded."
    if data_format == "binary"
        if binary_format == "int_16"
            bytes = 2
        elseif binary_format == "ieee_float_32"
            bytes = 4
        else
            @error "Binary formats other than Float32 and Int16 are not supported; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org"
        end

        fid = nothing
        try
            fid = open(eeg_file, "r")
        catch
            error("File $file_name cannot be loaded.")
        end

        signal = zeros(filesize(eeg_file) ÷ bytes)
        @inbounds for idx in 1:(filesize(eeg_file) ÷ bytes)
            buf = zeros(UInt8, bytes)
            readbytes!(fid, buf, bytes)
            # buf = reverse(buf)
            if bytes == 4
                # signal[idx] = Float64(reinterpret(Float32, reverse(buf))[1])
                signal[idx] = Float64(reinterpret(Float32, buf)[1])
            else
                signal[idx] = Float64(reinterpret(Int16, buf)[1])
            end
        end
        close(fid)
        # split signal into channels
        if data_orientation == "multiplexed"
            length(signal) % ch_n != 0 && (signal = signal[1:(end - length(signal) % ch_n)])
            data = zeros(ch_n, length(signal) ÷ ch_n, 1)
            idx2 = 1
            @inbounds for idx1 in 1:ch_n:length(signal)
                data[:, idx2, 1] = @views signal[idx1:(idx1 + (ch_n - 1))]
                idx2 += 1
            end
        else
            @error "Data orientation other than MULTIPLEXED is not supported; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org"
        end
    else
        @error "ASCII format is not supported; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org"
    end

    # apply gain
    data .*= gain
    # for idx in 1:ch_n
    #     data[idx, :, :] .*= gain[idx]
    # end

    # convert nV/mV to μV
    @inbounds for idx in 1:ch_n
        units[idx] == "" && (units[idx] = "μV")
        if ch_type[idx] == "eeg"
            if lowercase(units[idx]) == "mv"
                lowercase(units[idx]) == "μV"
                data[idx, :] .*= 1000
            end
            if lowercase(units[idx]) == "nv"
                lowercase(units[idx]) == "μV"
                data[idx, :] ./= 1000
            end
        end
    end

    time_pts = round.(collect(0:1/sampling_rate:size(data, 2) * size(data, 3) / sampling_rate)[1:end-1], digits=3)
    ep_time = round.((collect(0:1/sampling_rate:size(data, 2) / sampling_rate))[1:end-1], digits=3)
    
    file_size_mb = round(filesize(eeg_file) / 1024^2, digits=2)

    data_type = "eeg"

    s = _create_subject(id="",
                        first_name="",
                        middle_name="",
                        last_name=string(patient),
                        head_circumference=-1,
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
                              recording_notes=r_notes,
                              channel_type=ch_type[channel_order],
                              reference=ref,
                              clabels=clabels[channel_order],
                              transducers=transducers[channel_order],
                              units=units[channel_order],
                              prefiltering=prefiltering[channel_order],
                              sampling_rate=sampling_rate,
                              gain=gain[channel_order])
    e = _create_experiment(name=e_name,
                           notes=e_notes,
                           design="")

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = String[]

    obj = NeuroAnalyzer.NEURO(hdr, time_pts, ep_time, data[channel_order, :, :], components, markers, locs, history)
    nrow(locs) == 0 && _initialize_locs!(obj)
    
    _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits=2)) s)")

    return obj
    
end
