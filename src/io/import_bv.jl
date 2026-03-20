export import_bv
"""
    import_bv(file_name; <keyword arguments>)

Load a BrainVision BV/BVCDF recording and return a `NeuroAnalyzer.NEURO` object.

At minimum two files are required: a `.vhdr` (or `.ahdr`) header file and the corresponding `.eeg` signal file. Event markers are loaded automatically from the `.vmrk` file when present, or from a BIDS-style `_events.tsv` sidecar.

Channel locations are read from the `[Coordinates]` section of the header when available. A BIDS JSON sidecar (`.json`) is parsed for experiment metadata.

# Arguments

- `file_name::String`: path to the `.vhdr` or `.ahdr` header file
- `detect_type::Bool=true`: infer channel type from channel label

# Returns

- `NeuroAnalyzer.NEURO`

# Throws
- `ArgumentError` if the header file, signal file, or marker file cannot be loaded, or if the binary format / data orientation is unsupported
"""
function import_bv(file_name::String; detect_type::Bool = true)::NeuroAnalyzer.NEURO

    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))
    lowercase(splitext(file_name)[2]) in (".vhdr", ".ahdr") ||
        throw(ArgumentError("file_name must be a .vhdr or .ahdr file."))

    # ------------------------------------------------------------------ #
    # read and pre-process header lines                                  #
    # ------------------------------------------------------------------ #
    vhdr = try
        readlines(file_name)
    catch
        throw(ArgumentError("File $file_name cannot be loaded."))
    end

    # strip UTF-8 BOM (EF BB BF) if present on the first line
    vhdr[1][1] == '\ufeff' && (vhdr[1] = vhdr[1][4:end])
    startswith(lowercase(replace(vhdr[1], " " => "")), "brainvision") ||
        throw(ArgumentError("$file_name is not a BrainVision .VHDR file."))

    file_type = "BrainVision"

    # remove comment lines (lines starting with ';')
    Base.filter!(l -> !startswith(l, ';'), vhdr)

    # Remove all spaces so key=value parsing is whitespace-insensitive.
    vhdr = replace.(vhdr, " " => "")

    # ------------------------------------------------------------------ #
    # helper: extract the value for a given key from the header lines    #
    # returns `nothing` when the key is absent                           #
    # ------------------------------------------------------------------ #
    function vhdr_get(key)
        matches = vhdr[startswith.(lowercase.(vhdr), lowercase(key) * "=")]
        isempty(matches) ? nothing : split(matches[1], '=')[2]
    end

    # ------------------------------------------------------------------ #
    # parse scalar header fields                                         #
    # ------------------------------------------------------------------ #
    eeg_file = something(vhdr_get("DataFile"), "")
    marker_file = something(vhdr_get("MarkerFile"), "")
    data_format = lowercase(something(vhdr_get("DataFormat"), ""))
    data_orientation = lowercase(something(vhdr_get("DataOrientation"), ""))
    binary_format = lowercase(something(vhdr_get("BinaryFormat"), ""))
    ch_n = parse(Int64, something(vhdr_get("NumberOfChannels"), "0"))
    sampling_interval = parse(Float64, something(vhdr_get("SamplingInterval"), "0"))
    averaged = something(vhdr_get("Averaged"), "no")  == "yes"
    averaged_segments = parse(Int64, something(vhdr_get("AveragedSegments"), "0"))
    averaged_points = parse(Int64, something(vhdr_get("AveragedDataPoints"), "0"))
    segmentation = something(vhdr_get("Segmentation"), "no") == "yes"

    dir = dirname(file_name)
    eeg_file = replace(eeg_file, raw"$b" => dir)
    marker_file = replace(marker_file, raw"$b" => dir)

    # ------------------------------------------------------------------ #
    # locate section indices                                             #
    # ------------------------------------------------------------------ #
    channels_idx  = findfirst(l -> startswith(lowercase(l), "[channelinfos]"),  vhdr)
    locs_idx = findfirst(l -> startswith(lowercase(l), "[coordinates]"),   vhdr)
    soft_filt_idx = findfirst(l -> startswith(lowercase(l), "softwarefilters"), vhdr)
    channels_idx = something(channels_idx, 0)
    locs_idx = something(locs_idx, 0)
    soft_filt_idx = something(soft_filt_idx, 0)

    # ------------------------------------------------------------------ #
    # software filters (informational only)                              #
    # ------------------------------------------------------------------ #
    soft_filt = false
    if soft_filt_idx != 0
        if lowercase(vhdr[soft_filt_idx + 2]) != "disabled"
            _info("Embedded software filters are not implemented yet; " *
                  "if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org")
        end
    end
    soft_filt_file = replace(splitext(file_name)[1], "eeg" => "channels.tsv")
    isfile(soft_filt_file) &&
        (soft_filt = CSV.read(soft_filt_file, stringtype = String, DataFrame))

    # ------------------------------------------------------------------ #
    # optional BIDS JSON sidecar                                         #
    # ------------------------------------------------------------------ #
    # extract once
    ext = splitext(file_name)[2]
    js_file = replace(file_name, ext => ".json")
    e_name = ""
    e_notes = ""
    r_notes = ""
    ref = ""
    if isfile(js_file)
        js = JSON.parsefile(js_file, dicttype = Dict, inttype = Int64, use_mmap = true)
        "TaskName" in keys(js) && (e_name  = js["TaskName"])
        "TaskDescription" in keys(js) && (e_notes = js["TaskDescription"])
        "ManufacturersModelName" in keys(js) && (r_notes = js["ManufacturersModelName"])
        "EEGReference" in keys(js) && (ref = js["EEGReference"])
    end

    # ------------------------------------------------------------------ #
    # per-channel metadata from [ChannelInfos]                           #
    # ------------------------------------------------------------------ #
    patient = ""
    recording = ""
    recording_date = ""
    recording_time = ""
    transducers = repeat([""], ch_n)
    units = repeat([""], ch_n)
    gain = ones(ch_n)
    prefiltering = repeat([""], ch_n)
    clabels = repeat([""], ch_n)
    ref_chs = repeat([""], ch_n)

    for idx in 1:ch_n
        fields = split(split(vhdr[idx + channels_idx], '=')[2], ',')
        clabels[idx]  = replace(fields[1], "\1" => ",")
        ref_chs[idx]  = length(fields) >= 2 ? fields[2] : ""
        if length(fields) >= 3 && fields[3] != ""
            gain[idx] = parse(Float64, fields[3])
        end
        length(fields) >= 4 && (units[idx] = fields[4])
    end

    # ------------------------------------------------------------------ #
    # apply channel-list sidecar (channels.tsv) when present             #
    # ------------------------------------------------------------------ #
    if soft_filt !== false
        clabels = soft_filt[!, :name]
        units = soft_filt[!, :unit]
        prefiltering = repeat(["LP: "], ch_n) .*
                       string.(round.(soft_filt[!, :low_cutoff];  digits = 4)) .*
                       repeat([" Hz, HP: "], ch_n) .*
                       string.(round.(soft_filt[!, :high_cutoff]; digits = 4)) .*
                       " Hz"
    end
    clabels = _clean_labels(string.(clabels))

    ch_type = if detect_type
        _set_channel_types(clabels, "eeg")
    else
        repeat(["eeg"], ch_n)
    end

    # now that ch_type is defined, apply any overrides from the sidecar
    if soft_filt !== false
        for idx in 1:ch_n
            lowercase(soft_filt[idx, :type]) in ("eeg", "eog", "ecg", "emg") &&
                (ch_type[idx] = lowercase(soft_filt[idx, :type]))
        end
    end

    units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]

    # ------------------------------------------------------------------ #
    # channel locations from [Coordinates]                               #
    # ------------------------------------------------------------------ #
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
            l = occursin('=', vhdr[locs_idx + idx]) ?
                split(vhdr[locs_idx + idx], '=')[2] : vhdr[locs_idx + idx]
            parts = split(l, ',')
            loc_radius_sph[idx] = parse(Float64, parts[1])
            loc_theta_sph[idx] = parse(Float64, parts[2])
            loc_phi_sph[idx] = parse(Float64, parts[3])
            loc_theta[idx] = loc_theta_sph[idx]
            loc_radius[idx] = loc_radius_sph[idx]
            loc_x[idx], loc_y[idx], loc_z[idx] = sph2cart(
                loc_radius_sph[idx],
                loc_theta_sph[idx],
                loc_phi_sph[idx]
            )
        end
        locs = DataFrame(
            :ch_n => 1:ch_n,
            :label => clabels,
            :loc_theta => loc_theta,
            :loc_radius => loc_radius,
            :loc_x => loc_x,
            :loc_y => loc_y,
            :loc_z => loc_z,
            :loc_radius_sph => loc_radius_sph,
            :loc_theta_sph => loc_theta_sph,
            :loc_phi_sph => loc_phi_sph,
        )
    end

    # sampling interval is in μs; convert to Hz.
    sampling_rate = round(Int64, 1 / (sampling_interval / 1e6))

    # ------------------------------------------------------------------ #
    # event markers                                                      #
    # ------------------------------------------------------------------ #
    if marker_file != ""
        # Resolve relative marker path against the header file's directory.
        isabs(marker_file) || (marker_file = joinpath(dir, marker_file))
        isfile(marker_file) ||
            throw(ArgumentError("Marker file $marker_file cannot be loaded."))

        vmrk = readlines(marker_file)
        vmrk[1][1] == '\ufeff' && (vmrk[1] = vmrk[1][4:end])
        Base.filter!(l -> !startswith(l, ';'), vmrk)
        startswith(lowercase(replace(vmrk[1], " " => "")), "brainvision") ||
            throw(ArgumentError("$marker_file is not a valid BrainVision .VMRK file."))

        markers_idx = something(
            findfirst(l -> startswith(lowercase(replace(l, " " => "")), "[markerinfos]"), vmrk),
            0)

        # keep only lines that start with "Mk" (case-insensitive marker entries).
        mk_lines = Base.filter(l -> startswith(lowercase(l), "mk"), vmrk[(markers_idx + 1):end])
        n_mk = length(mk_lines)

        m_id = repeat([""], n_mk)
        m_desc = repeat(["marker"], n_mk)
        m_pos = zeros(Int64, n_mk)
        m_len = zeros(Int64, n_mk)
        m_ch = zeros(Int64, n_mk)

        for (idx, line) in enumerate(mk_lines)
            fields = split(split(line, '=')[2], ',')
            m_id[idx] = replace(fields[1], "\1" => ",")
            fields[2] != "" && (m_desc[idx] = replace(fields[2], "\1" => ","))
            m_pos[idx] = parse(Int64, fields[3])
            m_len[idx] = parse(Int64, fields[4])
            m_ch[idx]  = parse(Int64, fields[5])
        end

        markers = DataFrame(
            :id => m_id,
            :start => m_pos ./ sampling_rate,
            :length => round.(m_len ./ sampling_rate),
            :value => m_desc,
            :channel => m_ch,
        )
        if all(==(""), markers[!, :id])
            markers[!, :id] .= "mrk"
        end

    elseif isfile(replace(splitext(file_name)[1], "eeg" => "events.tsv"))
        vmrk = CSV.read(
            replace(splitext(file_name)[1], "eeg" => "events.tsv");
            stringtype = String, DataFrame)
        markers = DataFrame(
            :id => repeat(["mrk"], DataFrames.nrow(vmrk)),
            :start => vmrk[!, :sample] ./ sampling_rate,
            :length => round.(vmrk[!, :duration] ./ sampling_rate),
            :value => vmrk[!, :trial_type] .* "_" .* string.(vmrk[!, :value]),
            :channel => zeros(Int64, DataFrames.nrow(vmrk)),
        )
    else
        markers = DataFrame(
            :id => String[], :start => Float64[],
            :length => Float64[], :value => String[], :channel => Int64[])
    end

    # ------------------------------------------------------------------ #
    # signal data                                                         #
    # ------------------------------------------------------------------ #
    isabs(eeg_file) || (eeg_file = joinpath(dir, eeg_file))
    isfile(eeg_file) ||
        throw(ArgumentError("Signal file $eeg_file cannot be loaded."))

    if data_format == "binary"
        bytes = if binary_format == "int_16"
            2
        elseif binary_format == "ieee_float_32"
            4
        else
            throw(ArgumentError(
                "Binary format \"$binary_format\" is not supported. " *
                "Please send this file to adam.wysokinski@neuroanalyzer.org"))
        end

        # read the entire signal file at once then reinterpret
        raw = read(eeg_file)
        signal = if bytes == 4
            Float64.(reinterpret(Float32, raw))
        else
            Float64.(reinterpret(Int16, raw))
        end

        if data_orientation == "multiplexed"
            # trim any trailing incomplete frame
            rem = length(signal) % ch_n
            rem != 0 && (signal = signal[1:(end - rem)])
            n_samples = length(signal) ÷ ch_n
            data = zeros(ch_n, n_samples, 1)
            @inbounds for s in 1:n_samples
                data[:, s, 1] = signal[((s-1)*ch_n + 1):(s*ch_n)]
            end
        else
            throw(ArgumentError(
                "Data orientation \"$data_orientation\" is not supported. " *
                "Please send this file to adam.wysokinski@neuroanalyzer.org"))
        end
    else
        throw(ArgumentError(
            "Data format \"$data_format\" (ASCII) is not supported. " *
            "Please send this file to adam.wysokinski@neuroanalyzer.org"))
    end

    # apply per-channel gain.
    data .*= gain

    # ------------------------------------------------------------------ #
    # unit conversion: nV / mV → μV                                      #
    # ------------------------------------------------------------------ #
    @inbounds for idx in 1:ch_n
        units[idx] == "" && (units[idx] = "μV")
        ch_type[idx] == "eeg" || continue
        if lowercase(units[idx]) == "mv"
            units[idx] = "μV";  data[idx, :] .*= 1000
        elseif lowercase(units[idx]) == "nv"
            units[idx] = "μV";  data[idx, :] ./= 1000
        end
    end

    # ------------------------------------------------------------------ #
    # time axes                                                          #
    # ------------------------------------------------------------------ #
    n_samples  = size(data, 2) * size(data, 3)
    time_pts   = round.(range(0; step = 1/sampling_rate, length = n_samples);  digits = 4)
    epoch_time = round.(range(0; step = 1/sampling_rate, length = size(data,2)); digits = 4)

    # ------------------------------------------------------------------ #
    # assemble NEURO object                                              #
    # ------------------------------------------------------------------ #
    file_size_mb = round(filesize(eeg_file) / 1024^2; digits = 2)

    s = _create_subject(
        id = "", first_name = "", middle_name = "",
        last_name = string(patient), head_circumference = -1,
        handedness = "", weight = -1, height = -1)
    r = _create_recording_eeg(
        data_type = "eeg",
        file_name = file_name,
        file_size_mb = file_size_mb,
        file_type = file_type,
        recording = string(recording),
        recording_date = recording_date,
        recording_time = replace(recording_time, '.' => ':'),
        recording_notes = r_notes,
        channel_type = ch_type,
        channel_order = _sort_channels(ch_type),
        reference = ref,
        clabels = clabels,
        transducers = transducers,
        units = units,
        prefiltering = prefiltering,
        line_frequency = 50, # TODO: make this a keyword argument
        sampling_rate = sampling_rate,
        gain = gain,
        bad_channels = zeros(Bool, size(data, 1)))
    e   = _create_experiment(name = e_name, notes = e_notes, design = "")
    hdr = _create_header(subject = s, recording = r, experiment = e)

    obj = NeuroAnalyzer.NEURO(hdr, String[], markers, locs, time_pts, epoch_time, data)
    DataFrames.nrow(locs) == 0 && _initialize_locs!(obj)

    _info("Imported: " *
        uppercase(obj.header.recording[:data_type]) *
        " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj))" *
        "; $(round(obj.time_pts[end]; digits=2)) s)")

    return obj

end
