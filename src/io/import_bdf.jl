export import_bdf

"""
    import_bdf(file_name; <keyword arguments>)

Load a BDF or BDF+ file and return a `NeuroAnalyzer.NEURO` object.

BDF is BioSemi's 24-bit extension of the EDF format. Each sample is stored as a 24-bit little-endian two's-complement integer. BDF files always carry a Status channel as their last channel; BDF+ files may additionally contain annotation channels which are parsed into event markers and removed from the signal data.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: infer channel type from channel label

# Returns

- `NeuroAnalyzer.NEURO`

# Notes

- `sampling_rate = samples_per_datarecord ÷ data_record_duration`
- `gain = (physical_max − physical_min) / (digital_max − digital_min)`
- `value = (raw − digital_min) × gain + physical_min`

# References
1. https://www.biosemi.com/faq/file_format.htm
"""
function import_bdf(file_name::String; detect_type::Bool = true)::NeuroAnalyzer.NEURO

    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))
    lowercase(splitext(file_name)[2]) == ".bdf" ||
        throw(ArgumentError("$file_name does not have a .bdf extension."))

    # ------------------------------------------------------------------ #
    # parse fixed-length global header (256 bytes)                       #
    # ------------------------------------------------------------------ #
    # all header reads share a single open/close via the `do` block
    imported_object = open(file_name, "r") do fid

        buf = zeros(UInt8, 256)
        readbytes!(fid, buf, 256)
        header = String(Char.(buf))

        # byte 1: BDF identifier (0xFF); bytes 3–9: "BIOSEMI"
        Int(header[1]) == 255 ||
            throw(ArgumentError("$file_name is not a BDF file (first byte ≠ 0xFF)."))
        strip(header[3:9]) == "BIOSEMI" ||
            throw(ArgumentError("$file_name is not a valid BDF file (missing BIOSEMI identifier)."))
        file_type = "BDF"

        patient = strip(header[10:89])
        recording = strip(header[90:169])
        recording_date = header[170:177]
        recording_time = header[178:185]
        data_offset = parse(Int, strip(header[186:192]))
        reserved = strip(header[193:236])

        reserved == "BDF+D" &&
            throw(ArgumentError(
                "BDF+D (interrupted recordings) is not supported. " *
                "Please send this file to adam.wysokinski@neuroanalyzer.org"))
        reserved == "BDF+C" && (file_type = "BDF+")

        data_records = parse(Int, strip(header[237:244]))
        data_records_duration = parse(Float64, strip(header[245:252]))
        ch_n = parse(Int, strip(header[253:256]))

        # ------------------------------------------------------------ #
        # per-channel header fields (fixed-width records)              #
        # ------------------------------------------------------------ #
        read_fields(width, parse_fn = identity) = begin
            buf = zeros(UInt8, ch_n * width)
            readbytes!(fid, buf, ch_n * width)
            s = String(Char.(buf))
            [parse_fn(strip(s[(1 + (i-1)*width):(i*width)])) for i in 1:ch_n]
        end

        clabels = read_fields(16)
        transducers = read_fields(80)
        units = read_fields(8)
        units = replace(lowercase.(units), "uv" => "μV")

        physical_minimum = read_fields(8, s -> parse(Float64, s))
        physical_maximum = read_fields(8, s -> parse(Float64, s))
        digital_minimum = read_fields(8, s -> parse(Float64, s))
        digital_maximum = read_fields(8, s -> parse(Float64, s))
        prefiltering = read_fields(80)
        samples_per_datarecord = read_fields(8, s -> parse(Int, s))

        # ------------------------------------------------------------ #
        # channel types, annotation channels, markers channel          #
        # ------------------------------------------------------------ #
        clabels = _clean_labels(string.(clabels))
        ch_type = if detect_type
            _set_channel_types(clabels, "eeg")
        else
            t = repeat(["eeg"], ch_n)
            t[clabels .== "Status"] .= "mrk"
            t
        end
        units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]

        # BDF:  last channel is always the Status (markers) channel
        # BDF+: last channel is Status + possible extra annotation channels
        annotation_channels, markers_channel = if file_type == "BDF"
            Int64[], [ch_n]
        else
            ann = sort(unique(vcat(
                ch_n,
                getindex.(findall(occursin.("annotation", lowercase.(clabels))), 1),
            )))
            ann, getindex.(findall(ch_type .== "mrk"), 1)
        end

        # BDF does not support mixed sampling rates; use first signal channel.
        signal_chs = setdiff(1:ch_n, annotation_channels)
        sampling_rate = round(Int64,
            samples_per_datarecord[signal_chs[1]] / data_records_duration)
        gain = @. (physical_maximum - physical_minimum) /
                  (digital_maximum - digital_minimum)

        # ------------------------------------------------------------ #
        # signal data (24-bit little-endian two's-complement)          #
        # ------------------------------------------------------------ #
        seek(fid, data_offset)
        data = zeros(ch_n, samples_per_datarecord[signal_chs[1]] * data_records, 1)
        annotations = String[]

        @inbounds for rec in 1:data_records, ch in 1:ch_n
            n = samples_per_datarecord[ch]
            raw24 = zeros(UInt8, n * 3)
            readbytes!(fid, raw24, n * 3)

            col_start = (rec - 1) * n + 1
            col_end =  rec * n

            if ch in annotation_channels
                push!(annotations, String(Char.(raw24)))
                # annotation rows stay zero.

            elseif ch in markers_channel
                # status channel: combine only the two lower bytes (no gain).
                sig = zeros(Float64, n)
                @inbounds for (i, byte_idx) in enumerate(1:3:(n*3))
                    sig[i] = Float64(
                        Int32(raw24[byte_idx]) | Int32(raw24[byte_idx + 1]))
                end
                data[ch, col_start:col_end, 1] = sig

            else
                # signal channel: 24-bit little-endian two's-complement.
                # decode by packing bytes into the upper 24 bits of an Int32
                # then arithmetic-right-shifting by 8 to propagate the sign.
                sig = zeros(Float64, n)
                @inbounds for (i, byte_idx) in enumerate(1:3:(n*3))
                    b1 = Int32(raw24[byte_idx])     << 8
                    b2 = Int32(raw24[byte_idx + 1]) << 16
                    b3 = -Int32(-raw24[byte_idx + 2]) << 24
                    sig[i] = Float64((b1 | b2 | b3) >> 8)
                end
                data[ch, col_start:col_end, 1] = sig
            end
        end

        # Apply per-channel gain (broadcasts over samples and the epoch dim).
        data .*= gain

        # return a named tuple to avoid wide argument lists
        (;
            patient, recording, recording_date, recording_time, data_offset, reserved,
            file_type, data_records, data_records_duration, ch_n, clabels, transducers,
            units, physical_minimum, physical_maximum, digital_minimum, digital_maximum,
            prefiltering, samples_per_datarecord, ch_type, annotation_channels,
            sampling_rate, gain, annotations, data
        )
    end
    # file closed here in all cases, including exceptions

    # unpack named tuple
    (; patient, recording, recording_date, recording_time, file_type, ch_n, clabels,
       transducers, units, prefiltering, ch_type, annotation_channels, sampling_rate,
       gain, annotations, data) = imported_object
    # reuse binding

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
    # parse annotations / strip annotation channels from signal data      #
    # ------------------------------------------------------------------ #
    markers = if isempty(annotation_channels)
        DataFrame(
            :id => String[], :start => Float64[],
            :length => Float64[], :value => String[], :channel => Int64[])
    else
        m = _a2df(annotations)
        deleteat!(ch_type, annotation_channels)
        deleteat!(transducers, annotation_channels)
        deleteat!(units, annotation_channels)
        deleteat!(prefiltering, annotation_channels)
        deleteat!(clabels, annotation_channels)
        data  = data[setdiff(1:ch_n, annotation_channels), :, :]
        ch_n -= length(annotation_channels)
        m
    end

    # ------------------------------------------------------------------ #
    # time axes                                                            #
    # ------------------------------------------------------------------ #
    n_samples  = size(data, 2) * size(data, 3)
    time_pts   = round.(range(0; step = 1/sampling_rate, length = n_samples);  digits = 4)
    epoch_time = round.(range(0; step = 1/sampling_rate, length = size(data,2)); digits = 4)

    # ------------------------------------------------------------------ #
    # assemble NEURO object                                               #
    # ------------------------------------------------------------------ #
    file_size_mb = round(filesize(file_name) / 1024^2; digits = 2)

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
        recording_notes = "",
        channel_type = ch_type,
        channel_order = _sort_channels(ch_type),
        reference = _detect_montage(clabels, ch_type, "eeg"),
        clabels = clabels,
        transducers = string.(transducers),
        units = units,
        prefiltering = string.(prefiltering),
        line_frequency = 50, # TODO: make this a keyword argument
        sampling_rate = sampling_rate,
        gain = gain,
        bad_channels = zeros(Bool, size(data, 1)))
    e   = _create_experiment(name = "", notes = "", design = "")
    hdr = _create_header(subject = s, recording = r, experiment = e)

    locs = _initialize_locs()
    obj  = NeuroAnalyzer.NEURO(hdr, String[], markers, locs, time_pts, epoch_time, data)
    _initialize_locs!(obj)

    _info("Imported: " *
        uppercase(obj.header.recording[:data_type]) *
        " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj))" *
        "; $(round(obj.time_pts[end]; digits=2)) s)")

    return obj

end
