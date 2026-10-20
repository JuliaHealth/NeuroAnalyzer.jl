export import_edf

"""
    import_edf(file_name; <keyword arguments>)

Load an EDF or EDF+ file and return a `NeuroAnalyzer.NEURO` object.

Annotation-only channels (EDF+ TAL) are automatically detected, parsed into event markers, and then removed from the signal data. When channels have mixed sampling rates they are upsampled to the highest rate using Fourier resampling.

Channel units are taken from the file header; channels whose NeuroAnalyzer unit is μV (EEG, EOG, EMG, …) stored in V, mV or nV are converted to μV.

# Arguments

- `file_name::String`: path to the EDF/EDF+ file to load
- `detect_type::Bool=true`: infer channel type from channel label
- `clean_labels::Bool=true`: normalize channel labels (see `import_recording()`)

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
function import_edf(file_name::String; detect_type::Bool = true, clean_labels::Bool = true)::NeuroAnalyzer.NEURO
    # validate
    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))
    detect_format(file_name) === :edf ||
        throw(ArgumentError("$file_name is not EDF file."))

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

        recording_date = hdr[169:176]
        recording_time = hdr[177:184]
        data_offset = parse(Int, strip(hdr[185:192]))
        reserved = strip(hdr[193:236])

        reserved == "EDF+D" &&
            throw(
                ArgumentError(
                    "EDF+D (interrupted recordings) is not supported yet. " *
                    "Please send this file to adam.wysokinski@neuroanalyzer.org",
                ),
            )
        file_type = reserved == "EDF+C" ? "EDF+" : "EDF"

        # -1 = unknown (e.g. Alice 4); the record count is then taken from the file size
        data_records = parse(Int, strip(hdr[237:244]))
        data_records_duration = parse(Float64, strip(hdr[245:252]))
        data_records_duration > 0 ||
            throw(
                ArgumentError(
                    "This file contains only annotations; use import_edf_annotations().",
                ),
            )
        ch_n = parse(Int, strip(hdr[253:256]))

        # ------------------------------------------------------------ #
        # parse per-channel header fields                              #
        # each field is a fixed-width record repeated ch_n times.      #
        # ------------------------------------------------------------ #

        # helper: read `ch_n` fixed-width fields of `width` bytes each.
        # bytes are split into fields before decoding, so that non-ASCII
        # (Latin-1) characters cannot shift the fixed-width boundaries
        read_fields(width) = begin
            buf = UInt8[]
            readbytes!(fid, buf, ch_n * width)
            [strip(String(Char.(buf[(1 + (i - 1) * width):(i * width)]))) for i = 1:ch_n]
        end

        clabels = read_fields(16)
        transducers = read_fields(80)
        units = read_fields(8)

        physical_minimum = parse.(Float64, read_fields(8))
        physical_maximum = parse.(Float64, read_fields(8))
        digital_minimum = parse.(Float64, read_fields(8))
        digital_maximum = parse.(Float64, read_fields(8))
        prefiltering = read_fields(80)
        samples_per_datarecord = parse.(Int, read_fields(8))

        # ------------------------------------------------------------ #
        # resolve channel types and annotation channels                #
        # ------------------------------------------------------------ #
        clabels = clean_labels ? _clean_labels(clabels) : String.(strip.(clabels))
        ch_type =
            detect_type ? _set_channel_types(clabels, "eeg") : repeat(["eeg"], ch_n)

        annotation_channels = if file_type == "EDF"
            Int64[]
        else
            sort(getindex.(findall(occursin.("annotation", lowercase.(clabels))), 1))
        end
        signal_chs = setdiff(1:ch_n, annotation_channels)
        isempty(signal_chs) &&
            throw(ArgumentError("$file_name does not contain any signal channels."))

        # ------------------------------------------------------------ #
        # read signal data                                             #
        # ------------------------------------------------------------ #
        record_bytes = sum(samples_per_datarecord) * 2
        n_rec = _record_count(
            data_records, filesize(file_name) - data_offset, record_bytes, file_name,
        )
        seek(fid, data_offset)
        raw = read(fid, n_rec * record_bytes)
        sig, annotations = _read_edf_records(
            raw, samples_per_datarecord, n_rec, 2;
            annotation_channels = annotation_channels,
        )

        # digital → physical values, then file units → NeuroAnalyzer units
        gain = ones(ch_n)
        units, u_scale = _resolve_units(units, ch_type)
        for ch in signal_chs
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
        mrk_chs = findall(ch_type .== "mrk")
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
    # file closed here in all cases (including exceptions)

    # unpack named tuple
    (;
        patient, recording, recording_date, recording_time, file_type, ch_n, clabels,
        transducers, units, prefiltering, ch_type, annotation_channels, sampling_rate,
        gain, annotations, data,
    ) = imported_object

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
    # build time axes                                                    #
    # ------------------------------------------------------------------ #
    n_samples = size(data, 2) * size(data, 3)
    time_pts = _time_axis(n_samples, sampling_rate)
    epoch_time = _time_axis(size(data, 2), sampling_rate)

    # ------------------------------------------------------------------ #
    # assemble NEURO object                                              #
    # ------------------------------------------------------------------ #
    file_size_mb = round(filesize(file_name) / 1024^2; digits = 2)

    s = _create_subject(;
        id = "", first_name = "", middle_name = "",
        last_name = string(patient),
        head_circumference = -1,
        handedness = "",
        weight = -1,
        height = -1,
    )
    r = _create_recording_eeg(;
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
