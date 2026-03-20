export import_csv

"""
    import_csv(file_name; <keyword arguments>)

Load a CSV (or gzip-compressed CSV) file and return a `NeuroAnalyzer.NEURO` object.

The file layout (time × channels or channels × time) is detected automatically from the type of the first column:

- **Time × channels**: first column is numeric (time axis); remaining columns are channels labelled by their header row.
- **Channels × time**: first column is a string (channel names); remaining column headers are parsed as time points.

Sampling rate is inferred from the step between the first two time points.

Gzip-compressed files (`.csv.gz`) are decompressed automatically by `CSV.jl`.

# Arguments

- `file_name::String`: path to the `.csv` or `.csv.gz` file
- `detect_type::Bool=true`: infer channel type from channel label

# Returns
- `NeuroAnalyzer.NEURO`

# Throws
- `ArgumentError` if the file does not exist
"""
function import_csv(file_name::String; detect_type::Bool = true)::NeuroAnalyzer.NEURO

    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))

    file_type = "CSV"
    df = CSV.read(file_name; stringtype = String, DataFrame)

    # ------------------------------------------------------------------ #
    # detect layout and extract time axis + signal matrix                 #
    # ------------------------------------------------------------------ #
    if eltype(df[:, 1]) <: Real
        # layout: rows = time points, columns = channels
        # first column is the time axis; remaining columns are channels
        time_pts_raw = Float64.(df[!, 1])
        data = Matrix(df[:, 2:end])' # → (ch_n × n_samples)
        ch_n = DataFrames.ncol(df) - 1
        clabels_tmp  = String.(names(df)[2:end])
    else
        # layout: rows = channels, columns = time points
        # first column holds channel names; remaining column headers are time
        time_pts_raw = parse.(Float64, names(df)[2:end])
        data = Matrix(df[!, 2:end]) # → (ch_n × n_samples)
        ch_n = DataFrames.nrow(df)
        clabels_tmp = string.(df[:, 1])
    end

    # add epoch dimension.
    data = reshape(data, ch_n, size(data, 2), 1)

    # ------------------------------------------------------------------ #
    # sampling rate                                                      #
    # ------------------------------------------------------------------ #
    sampling_rate = round(Int64, 1 / (time_pts_raw[2] - time_pts_raw[1]))

    # ------------------------------------------------------------------ #
    # Reconstruct a canonical time axis starting at 0                    #
    # ------------------------------------------------------------------ #
    n_samples = size(data, 2) * size(data, 3)
    t0 = time_pts_raw[1]
    time_pts = round.(range(t0; step = 1/sampling_rate, length = n_samples);  digits = 4)
    epoch_time = round.(range(t0; step = 1/sampling_rate, length = size(data,2)); digits = 4)

    # ------------------------------------------------------------------ #
    # channel metadata                                                   #
    # ------------------------------------------------------------------ #
    clabels = _clean_labels(clabels_tmp)
    ch_type = detect_type ? _set_channel_types(clabels, "eeg") : repeat(["eeg"], ch_n)
    units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]

    markers = DataFrame(
        :id => String[],
        :start => Float64[],
        :length => Float64[],
        :value => String[],
        :channel => Int64[])

    # ------------------------------------------------------------------ #
    # assemble NEURO object                                               #
    # ------------------------------------------------------------------ #
    file_size_mb = round(filesize(file_name) / 1024^2; digits = 2)

    s = _create_subject(
        id = "", first_name = "", middle_name = "", last_name = "",
        head_circumference = -1, handedness = "", weight = -1, height = -1)
    r = _create_recording_eeg(
        data_type = "eeg",
        file_name = file_name,
        file_size_mb = file_size_mb,
        file_type = file_type,
        recording = "",
        recording_date = "",
        recording_time = "",
        recording_notes = "",
        channel_type = ch_type,
        channel_order = _sort_channels(ch_type),
        reference = _detect_montage(clabels, ch_type, "eeg"),
        clabels = clabels,
        transducers = repeat([""], ch_n),
        units = units,
        prefiltering = repeat([""], ch_n),
        line_frequency = 50, # TODO: make this a keyword argument
        sampling_rate = sampling_rate,
        gain = ones(ch_n),
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