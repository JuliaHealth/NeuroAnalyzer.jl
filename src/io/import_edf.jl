export import_edf

"""
    import_edf(file_name; <keyword arguments>)

Load an EDF or EDF+ file and return a `NeuroAnalyzer.NEURO` object.

Annotation-only channels (EDF+ TAL) are automatically detected, parsed into event markers, and then removed from the signal data. When channels have mixed sampling rates they are upsampled to the highest rate using Fourier resampling.

# Arguments

- `file_name::String`: path to the EDF/EDF+ file to load
- `detect_type::Bool=true`: infer channel type from channel label

# Returns

- `NeuroAnalyzer.NEURO`

# Notes

- `sampling_rate = samples_per_datarecord ÷ data_record_duration`
- `gain = (physical_max - physical_min) / (digital_max - digital_min)`
- `value = (raw_value - digital_min) × gain + physical_min`

# References
1. Kemp B et al. A simple format for exchange of digitized polygraphic recordings. *EEG Clin Neurophysiol*. 1992;82(5):391-3.
2. Kemp B, Olivan J. European data format 'plus' (EDF+). *Clin Neurophysiol*. 2003;114:1755-61.
3. https://www.edfplus.info/specs/
"""
function import_edf(file_name::String; detect_type::Bool = true)::NeuroAnalyzer.NEURO

    # ------------------------------------------------------------------ #
    # validate file                                                      #
    # ------------------------------------------------------------------ #
    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))
    lowercase(splitext(file_name)[2]) == ".edf" ||
        throw(ArgumentError("$file_name is not an EDF file (wrong extension)."))

    # ------------------------------------------------------------------ #
    # parse fixed-length global header (256 bytes)                       #
    # ------------------------------------------------------------------ #
    # all header reads share a single open/close via the `do` block
    imported_object = open(file_name, "r") do fid

        buf = zeros(UInt8, 256)
        readbytes!(fid, buf, 256)
        hdr = String(Char.(buf))

        # byte 1-8: version. EDF = "0", EDF+ adds a reserved field
        file_type_raw = strip(hdr[1:8])
        parse(Int, file_type_raw) == 0 ||
            throw(ArgumentError("$file_name is not a valid EDF file (version field ≠ 0)."))

        patient = strip(hdr[9:88])
        recording = strip(hdr[89:168])

        # special-case: Alice 4 EDF files do not conform to the standard
        occursin("Alice 4", recording) &&
            return import_alice4(file_name; detect_type)

        recording_date = hdr[169:176]
        recording_time = hdr[177:184]
        data_offset = parse(Int, strip(hdr[185:192]))
        reserved = strip(hdr[193:236])

        reserved == "EDF+D" &&
            throw(ArgumentError(
                "EDF+D (interrupted recordings) is not supported yet. " *
                "Please send this file to adam.wysokinski@neuroanalyzer.org"))
        file_type = reserved == "EDF+C" ? "EDF+" : "EDF"

        data_records = parse(Int, strip(hdr[237:244]))
        data_records_duration = parse(Float64, strip(hdr[245:252]))
        data_records_duration > 0 ||
            throw(ArgumentError(
                "This file contains only annotations; use import_edf_annotations()."))
        ch_n = parse(Int, strip(hdr[253:256]))

        # ------------------------------------------------------------ #
        # parse per-channel header fields                              #
        # each field is a fixed-width record repeated ch_n times.      #
        # ------------------------------------------------------------ #

        # helper: read `ch_n` fixed-width fields of `width` bytes each.
        read_fields(width) = begin
            buf = UInt8[]
            readbytes!(fid, buf, ch_n * width)
            s = String(Char.(buf))
            [strip(s[(1 + (i-1)*width):(i*width)]) for i in 1:ch_n]
        end

        clabels = read_fields(16)
        transducers = read_fields(80)
        units = read_fields(8)
        units = replace(lowercase.(units), "uv" => "μV")

        physical_minimum = parse.(Float64, read_fields(8))
        physical_maximum = parse.(Float64, read_fields(8))
        digital_minimum  = parse.(Float64, read_fields(8))
        digital_maximum  = parse.(Float64, read_fields(8))
        prefiltering = read_fields(80)
        samples_per_datarecord = parse.(Int, read_fields(8))

        # ------------------------------------------------------------ #
        # resolve channel types and annotation channels                #
        # ------------------------------------------------------------ #
        clabels = _clean_labels(string.(clabels))
        ch_type = detect_type ? _set_channel_types(clabels, "eeg") : repeat(["eeg"], ch_n)
        units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]

        annotation_channels = if file_type == "EDF"
            Int64[]
        else
            sort(getindex.(findall(occursin.("annotation", lowercase.(clabels))), 1))
        end
        # FIX: `markers_channel` was computed but never used - removed

        # ------------------------------------------------------------ #
        # determine sampling rate(s)                                   #
        # ------------------------------------------------------------ #
        signal_chs = setdiff(1:ch_n, annotation_channels)

        if length(unique(samples_per_datarecord[signal_chs])) == 1
            first_sig = signal_chs[1]
            sampling_rate = round(Int64,
                samples_per_datarecord[first_sig] / data_records_duration)
        else
            sampling_rate = round.(Int64,
                samples_per_datarecord[signal_chs] ./ data_records_duration)
        end

        gain = @. (physical_maximum - physical_minimum) /
                  (digital_maximum  - digital_minimum)

        # ------------------------------------------------------------ #
        # read signal data                                             #
        # ------------------------------------------------------------ #
        seek(fid, data_offset)
        annotations = String[]

        data = if sampling_rate isa Int64
            # uniform sampling rate across all signal channels
            d = zeros(ch_n, samples_per_datarecord[signal_chs[1]] * data_records, 1)

            @inbounds for rec in 1:data_records, ch in 1:ch_n
                raw = UInt8[]
                readbytes!(fid, raw, samples_per_datarecord[ch] * 2)
                if ch in annotation_channels
                    push!(annotations, String(Char.(raw)))
                    # leave annotation channel rows as zeros.
                else
                    col_start = (rec - 1) * samples_per_datarecord[ch] + 1
                    col_end = rec * samples_per_datarecord[ch]
                    d[ch, col_start:col_end, 1] = reinterpret(Int16, raw)
                end
            end

            # apply gain to all channels (annotation rows stay 0)
            d .*= gain
            d

        else
            # mixed sampling rates — read entire data block at once,
            # then process channel by channel, upsampling to the maximum rate
            max_rate    = maximum(sampling_rate)
            max_spdr    = maximum(samples_per_datarecord[signal_chs])

            raw_all = UInt8[]
            readbytes!(fid, raw_all, filesize(file_name) - data_offset; all = true)

            data_records = length(raw_all) ÷ 2 ÷ sum(samples_per_datarecord)

            d = zeros(ch_n, data_records * max_rate)
            # byte cursor into raw_all
            pos = 1

            @inbounds for rec in 1:data_records, ch in 1:ch_n
                n_bytes = samples_per_datarecord[ch] * 2
                chunk = raw_all[pos:(pos + n_bytes - 1)]
                pos  += n_bytes

                col_start = (rec - 1) * max_spdr + 1
                col_end   =  rec      * max_spdr

                if ch in annotation_channels
                    push!(annotations, String(Char.(chunk)))
                    # leave rows as zeros   
                else
                    tmp = Float64.(reinterpret(Int16, chunk))
                    if sampling_rate[findfirst(==(ch), signal_chs)] != max_rate
                        tmp = FourierTools.resample(tmp, max_rate)
                    end
                    d[ch, col_start:col_end] = tmp
                end
            end

            d .*= gain
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
    # file closed here in all cases (including exceptions)

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
            units[idx] = "μV"
            data[idx, :] .*= 1000
        elseif lowercase(units[idx]) == "nv"
            units[idx] = "μV"
            data[idx, :] ./= 1000
        end
    end

    # ------------------------------------------------------------------ #
    # parse annotations / strip annotation channels from signal data     #
    # ------------------------------------------------------------------ #
    markers = if isempty(annotation_channels)
        DataFrame(
            :id => String[],
            :start => Float64[],
            :length => Float64[],
            :value => String[],
            :channel => Int64[],
        )
    else
        m = _a2df(annotations)
        # remove annotation channels from every per-channel vector / array
        deleteat!(ch_type, annotation_channels)
        deleteat!(transducers, annotation_channels)
        deleteat!(units, annotation_channels)
        deleteat!(prefiltering, annotation_channels)
        deleteat!(clabels, annotation_channels)
        data = data[setdiff(1:ch_n, annotation_channels), :, :]
        ch_n -= length(annotation_channels)
        m
    end

    # ------------------------------------------------------------------ #
    # build time axes                                                    #
    # ------------------------------------------------------------------ #
    n_samples = size(data, 2) * size(data, 3)
    time_pts  = round.(
        range(0; step = 1/sampling_rate, length = n_samples);
        digits = 4)
    epoch_time = round.(
        range(0; step = 1/sampling_rate, length = size(data, 2));
        digits = 4)

    # ------------------------------------------------------------------ #
    # assemble NEURO object                                               #
    # ------------------------------------------------------------------ #
    file_size_mb = round(filesize(file_name) / 1024^2; digits = 2)

    s = _create_subject(
        id = "", first_name = "", middle_name = "",
        last_name = string(patient),
        head_circumference  = -1,
        handedness = "",
        weight = -1,
        height = -1,
    )
    r = _create_recording_eeg(
        data_type = "eeg",
        file_name = file_name,
        file_size_mb = file_size_mb,
        file_type = file_type,
        recording = string(recording),
        recording_date = recording_date,
        # replace dots with colons (Alice 4 / non-conforming EDF writers
        # sometimes use '.' as a time separator instead of ':').
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
        bad_channels = zeros(Bool, size(data, 1)),
    )
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
