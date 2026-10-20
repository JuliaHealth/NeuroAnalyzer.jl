export import_bdf

"""
    import_bdf(file_name; <keyword arguments>)

Load a BDF or BDF+ file and return a `NeuroAnalyzer.NEURO` object.

BDF is BioSemi's 24-bit extension of the EDF format. Each sample is stored as a 24-bit little-endian two's-complement integer. BDF files always carry a Status channel as their last channel; BDF+ files may additionally contain annotation channels which are parsed into event markers and removed from the signal data.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: infer channel type from channel label
- `clean_labels::Bool=true`: normalize channel labels (see `import_recording()`)

# Returns

- `NeuroAnalyzer.NEURO`

# Notes

- `sampling_rate = samples_per_datarecord ÷ data_record_duration`; channels with fewer samples per record are upsampled to the highest rate
- `gain = (physical_max − physical_min) / (digital_max − digital_min)`
- `value = (raw − digital_min) × gain + physical_min`
- the Status channel keeps the 16 trigger bits (no scaling)
- channel units are taken from the file header; channels whose NeuroAnalyzer unit is μV stored in V, mV or nV are converted to μV

# References

1. https://www.biosemi.com/faq/file_format.htm
"""
function import_bdf(file_name::String; detect_type::Bool = true, clean_labels::Bool = true)::NeuroAnalyzer.NEURO
    # validate
    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))
    detect_format(file_name) === :bdf ||
        throw(ArgumentError("$file_name is not an BDF file."))

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
            throw(
                ArgumentError(
                    "$file_name is not a valid BDF file (missing BIOSEMI identifier).",
                ),
            )
        file_type = "BDF"

        patient = strip(header[10:89])
        recording = strip(header[90:169])
        recording_date = header[170:177]
        recording_time = header[178:185]
        data_offset = parse(Int, strip(header[186:192]))
        reserved = strip(header[193:236])

        reserved == "BDF+D" &&
            throw(
                ArgumentError(
                    "BDF+D (interrupted recordings) is not supported. " *
                    "Please send this file to adam.wysokinski@neuroanalyzer.org",
                ),
            )
        reserved == "BDF+C" && (file_type = "BDF+")

        data_records = parse(Int, strip(header[237:244]))
        data_records_duration = parse(Float64, strip(header[245:252]))
        ch_n = parse(Int, strip(header[253:256]))

        # ------------------------------------------------------------ #
        # per-channel header fields (fixed-width records)              #
        # ------------------------------------------------------------ #
        # bytes are split into fields before decoding, so that non-ASCII
        # (Latin-1) characters cannot shift the fixed-width boundaries
        read_fields(width, parse_fn = identity) = begin
            buf = zeros(UInt8, ch_n * width)
            readbytes!(fid, buf, ch_n * width)
            [
                parse_fn(strip(String(Char.(buf[(1 + (i - 1) * width):(i * width)]))))
                for i = 1:ch_n
            ]
        end

        clabels = read_fields(16)
        transducers = read_fields(80)
        units = read_fields(8)

        physical_minimum = read_fields(8, s -> parse(Float64, s))
        physical_maximum = read_fields(8, s -> parse(Float64, s))
        digital_minimum = read_fields(8, s -> parse(Float64, s))
        digital_maximum = read_fields(8, s -> parse(Float64, s))
        prefiltering = read_fields(80)
        samples_per_datarecord = read_fields(8, s -> parse(Int, s))

        # ------------------------------------------------------------ #
        # channel types, annotation channels, Status channel           #
        # ------------------------------------------------------------ #
        clabels = clean_labels ? _clean_labels(clabels) : String.(strip.(clabels))
        ch_type = if detect_type
            _set_channel_types(clabels, "eeg")
        else
            t = repeat(["eeg"], ch_n)
            t[lowercase.(clabels) .== "status"] .= "mrk"
            t
        end

        # BDF+ annotation channels ("BDF Annotations") hold TAL text, not samples
        annotation_channels = if file_type == "BDF"
            Int64[]
        else
            sort(getindex.(findall(occursin.("annotation", lowercase.(clabels))), 1))
        end
        # the BioSemi Status channel carries the trigger inputs in its lower 16 bits
        status_channels = setdiff(findall(lowercase.(clabels) .== "status"), annotation_channels)
        signal_chs = setdiff(1:ch_n, annotation_channels)
        isempty(signal_chs) &&
            throw(ArgumentError("$file_name does not contain any signal channels."))

        # ------------------------------------------------------------ #
        # signal data (24-bit little-endian two's-complement)          #
        # ------------------------------------------------------------ #
        record_bytes = sum(samples_per_datarecord) * 3
        n_rec = _record_count(
            data_records, filesize(file_name) - data_offset, record_bytes, file_name,
        )
        seek(fid, data_offset)
        raw = read(fid, n_rec * record_bytes)
        sig, annotations = _read_edf_records(
            raw, samples_per_datarecord, n_rec, 3;
            annotation_channels = annotation_channels,
            trigger_channels = status_channels,
        )

        # digital → physical values (not for the Status channel, which holds trigger codes),
        # then file units → NeuroAnalyzer units
        gain = ones(ch_n)
        units, u_scale = _resolve_units(units, ch_type)
        for ch in setdiff(signal_chs, status_channels)
            gain[ch] = _phys_scale!(
                sig[ch], digital_minimum[ch], digital_maximum[ch],
                physical_minimum[ch], physical_maximum[ch],
            )
            if u_scale[ch] != 1
                sig[ch] .*= u_scale[ch]
                gain[ch] *= u_scale[ch]
            end
        end

        # channels with fewer samples per record are upsampled to the highest rate
        # (marker channels by nearest sample, so that their codes are not interpolated)
        mrk_chs = union(findall(ch_type .== "mrk"), status_channels)
        data, resampled = _merge_channels(sig, annotation_channels; nearest = mrk_chs)
        fs = maximum(samples_per_datarecord[signal_chs]) / data_records_duration
        resampled && _info("Channels upsampled to $(round(fs; digits = 4)) Hz")
        # NeuroAnalyzer stores integer sampling rates
        data, sampling_rate = _integer_rate(data, fs; nearest = mrk_chs)

        # return a named tuple to avoid wide argument lists
        return (;
            patient, recording, recording_date, recording_time, file_type, ch_n,
            clabels, transducers, units, prefiltering, ch_type, annotation_channels,
            sampling_rate, gain, annotations, data,
        )
    end
    # file closed here in all cases, including exceptions

    # unpack named tuple
    (;
        patient, recording, recording_date, recording_time, file_type, ch_n, clabels,
        transducers, units, prefiltering, ch_type, annotation_channels, sampling_rate,
        gain, annotations, data,
    ) = imported_object

    # ------------------------------------------------------------------ #
    # parse annotations / strip annotation channels from signal data      #
    # ------------------------------------------------------------------ #
    markers = if isempty(annotation_channels)
        DataFrame(
            :id => String[], :start => Float64[],
            :length => Float64[], :value => String[], :channel => Int64[],
        )
    else
        m = _a2df(annotations)
        keep = setdiff(1:ch_n, annotation_channels)
        ch_type = ch_type[keep]
        transducers = transducers[keep]
        units = units[keep]
        prefiltering = prefiltering[keep]
        clabels = clabels[keep]
        gain = gain[keep]
        data = data[keep, :, :]
        ch_n = length(keep)
        m
    end

    # ------------------------------------------------------------------ #
    # time axes                                                            #
    # ------------------------------------------------------------------ #
    n_samples = size(data, 2) * size(data, 3)
    time_pts = _time_axis(n_samples, sampling_rate)
    epoch_time =
        _time_axis(size(data, 2), sampling_rate)

    # ------------------------------------------------------------------ #
    # assemble NEURO object                                               #
    # ------------------------------------------------------------------ #
    file_size_mb = round(filesize(file_name) / 1024^2; digits = 2)

    s = _create_subject(;
        id = "", first_name = "", middle_name = "",
        last_name = string(patient), head_circumference = -1,
        handedness = "", weight = -1, height = -1,
    )
    r = _create_recording_eeg(;
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
        bad_channels = zeros(Bool, ch_n),
    )
    e = _create_experiment(; name = "", notes = "", design = "")
    hdr = _create_header(; subject = s, recording = r, experiment = e)

    locs = _initialize_locs()
    obj = NeuroAnalyzer.NEURO(hdr, String[], markers, locs, time_pts, epoch_time, data)
    _initialize_locs!(obj)

    _info(
        "Imported: " *
        uppercase(obj.header.recording[:data_type]) *
        " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj))" *
        "; $(round(obj.time_pts[end], digits = 2)) s)",
    )

    return obj
end
