export import_alice4

"""
    import_alice4(file_name; <keyword arguments>)

Load an EDF file exported from the Alice 4 Polysomnography System and return a `NeuroAnalyzer.NEURO` object.

Alice 4 EDF files are non-conforming in two ways that prevent `import_edf` from handling them:
- The `data_records` header field is always `-1` (unknown); the true record count is inferred from file size.
- Channels may have different sampling rates; they are upsampled to the highest rate using Fourier resampling.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: infer channel type from channel label

# Returns

- `NeuroAnalyzer.NEURO`

# Throws

- `ArgumentError` if the file does not exist, is not EDF, or is not an Alice 4 recording
"""
function import_alice4(file_name::String; detect_type::Bool = true)::NeuroAnalyzer.NEURO

    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))

    # ------------------------------------------------------------------ #
    # parse header — all reads share one open/close via the `do` block.  #
    # ------------------------------------------------------------------ #
    imported_object = open(file_name, "r") do fid

        header = _v2s(_fread(fid, 256, :s))

        # bytes 1–8: version; must be 0 for EDF
        file_type_raw = strip(header[1:8])
        parse(Int, file_type_raw) == 0 ||
            throw(ArgumentError("$file_name is not a valid EDF file (version ≠ 0)."))
        file_type = "EDF"

        patient = strip(header[9:88])
        recording = strip(header[89:168])
        occursin("Alice 4", recording) ||
            throw(ArgumentError("$file_name is not an Alice 4 EDF file."))

        recording_date = header[169:176]
        recording_time = header[177:184]
        data_offset = parse(Int, strip(header[185:192]))
        reserved = strip(header[193:236])

        reserved == "EDF+D" &&
            throw(ArgumentError(
                "EDF+D (interrupted recordings) is not supported. " *
                "Please send this file to adam.wysokinski@neuroanalyzer.org"))
        reserved == "EDF+C" && (file_type = "EDF+")

        # Alice 4 always writes -1 here; validated below before use
        data_records = parse(Int, strip(header[237:244]))
        data_records != -1 &&
            throw(ArgumentError(
                "data_records ≠ -1; this looks like a standard EDF file. " *
                "Use import_edf() instead."))

        data_records_duration = parse(Float64, strip(header[245:252]))
        ch_n = parse(Int, strip(header[253:256]))

        # ------------------------------------------------------------ #
        # per-channel header fields (fixed-width records)              #
        # ------------------------------------------------------------ #
        read_fields(width, parse_fn = identity) = begin
            buf = _v2s(_fread(fid, ch_n * width, :s))
            [parse_fn(strip(buf[(1 + (i-1)*width):(i*width)])) for i in 1:ch_n]
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
        # channel types and annotation channels                        #
        # ------------------------------------------------------------ #
        clabels = _clean_labels(string.(clabels))
        ch_type = detect_type ? _set_channel_types(clabels, "eeg") : repeat(["eeg"], ch_n)
        units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]

        annotation_channels = if file_type == "EDF"
            Int64[]
        else
            sort(getindex.(findall(occursin.("annotation", lowercase.(clabels))), 1))
        end
        # FIX: `markers_channel` was computed but never used — removed.

        # ------------------------------------------------------------ #
        # sampling rate                                                 #
        # ------------------------------------------------------------ #
        if length(unique(samples_per_datarecord)) == 1
            sampling_rate = round(Int64,
                samples_per_datarecord[1] / data_records_duration)
        else
            sampling_rate = round.(Int64,
                samples_per_datarecord ./ data_records_duration)
        end

        gain = @. (physical_maximum - physical_minimum) /
                  (digital_maximum  - digital_minimum)

        # ------------------------------------------------------------ #
        # signal data                                                   #
        # ------------------------------------------------------------ #
        seek(fid, data_offset)
        annotations = String[]

        data = if sampling_rate isa Int64
            # uniform rate. Alice 4 writes data_records = -1 in the header,
            # so we infer the true record count from file size
            bytes_per_record = sum(samples_per_datarecord) * 2
            true_data_records = (filesize(file_name) - data_offset) ÷ bytes_per_record

            d = zeros(ch_n, samples_per_datarecord[1] * true_data_records, 1)

            @inbounds for rec in 1:true_data_records, ch in 1:ch_n
                raw = zeros(UInt8, samples_per_datarecord[ch] * 2)
                readbytes!(fid, raw, samples_per_datarecord[ch] * 2)
                if ch in annotation_channels
                    push!(annotations, String(Char.(raw)))
                else
                    sig = map(ltoh, reinterpret(Int16, raw))
                    col_start = (rec - 1) * samples_per_datarecord[ch] + 1
                    col_end = rec * samples_per_datarecord[ch]
                    d[ch, col_start:col_end, 1] = sig .* gain[ch]
                end
            end
            d

        else
            # Mixed rates — read the full data block, then unpack per channel.
            max_rate = maximum(sampling_rate)

            raw_all = zeros(UInt8, filesize(file_name) - data_offset)
            readbytes!(fid, raw_all, length(raw_all); all = true)
            signal = map(ltoh, reinterpret(Int16, raw_all))

            # Infer true record count from total samples and per-record sum.
            true_data_records = length(signal) ÷ sum(samples_per_datarecord)
            d = zeros(ch_n, true_data_records * max_rate)
            data_segment = max_rate

            pos = 1

            @inbounds for rec in 1:true_data_records, ch in 1:ch_n
                n = samples_per_datarecord[ch]
                tmp = Float64.(signal[pos:(pos + n - 1)])
                pos += n

                col_start = (rec - 1) * data_segment + 1
                col_end =  rec * data_segment

                if ch in annotation_channels
                    push!(annotations, String(Char.(reinterpret(UInt8, signal[pos-n:pos-1]))))
                else
                    tmp .*= gain[ch]
                    if sampling_rate[ch] != max_rate
                        tmp = FourierTools.resample(tmp, max_rate)
                    end
                    d[ch, col_start:col_end] = tmp
                end
            end

            _info("Channels upsampled to $max_rate Hz")
            sampling_rate = max_rate
            d
        end

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
    # parse annotations / strip annotation channels                       #
    # ------------------------------------------------------------------ #
    markers = if isempty(annotation_channels)
        DataFrame(
            :id => String[], :start => Float64[],
            :length => Float64[], :value => String[], :channel => Int64[])
    else
        m = _a2df(annotations)
        deleteat!(ch_type,      annotation_channels)
        deleteat!(transducers,  annotation_channels)
        deleteat!(units,        annotation_channels)
        deleteat!(prefiltering, annotation_channels)
        deleteat!(clabels,      annotation_channels)
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
    # Assemble NEURO object                                               #
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
        units = units,
        transducers = string.(transducers),
        prefiltering = string.(prefiltering),
        line_frequency = 50, # TODO: make this a keyword argument
        sampling_rate = sampling_rate,
        gain = gain,
        bad_channels = zeros(Bool, size(data, 1)))
    e = _create_experiment(name = "", notes = "", design = "")

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
