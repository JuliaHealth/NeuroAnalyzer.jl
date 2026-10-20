export import_csv

"""
    import_csv(file_name; <keyword arguments>)

Load a CSV (or gzip-compressed CSV) file and return a `NeuroAnalyzer.NEURO` object.

The file layout (time × channels or channels × time) is detected automatically from the type of the first column:

- **Time × channels**: first column is numeric (time axis); remaining columns are channels labelled by their header row.
- **Channels × time**: first column is a string (channel names); remaining column headers are parsed as time points.

The time axis must be in milliseconds. Sampling rate is inferred from the mean step between time points; the imported time axis starts at 0.

Gzip-compressed files (`.csv.gz`) are decompressed automatically by `CSV.jl`.

# Arguments

- `file_name::String`: path to the `.csv` or `.csv.gz` file
- `detect_type::Bool=true`: infer channel type from channel label
- `clean_labels::Bool=true`: normalize channel labels (see `import_recording()`)

# Returns

- `NeuroAnalyzer.NEURO`
"""
function import_csv(file_name::String; detect_type::Bool = true, clean_labels::Bool = true)::NeuroAnalyzer.NEURO
    # validate
    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))
    detect_format(file_name) in [:csv, :csv_gz] ||
        throw(ArgumentError("$file_name is not CSV file."))

    file_type = "CSV"

    df = CSV.read(
        file_name,
        DataFrame;
        stringtype = String,
    )

    # ------------------------------------------------------------------ #
    # detect layout and extract time axis + signal matrix                 #
    # ------------------------------------------------------------------ #
    if eltype(df[:, 1]) <: Real
        # layout: rows = time points, columns = channels
        # first column is the time axis; remaining columns are channels
        time_pts_raw = Float64.(df[!, 1])
        data = Matrix{Float64}(df[:, 2:end])' # → (ch_n × n_samples)
        ch_n = DataFrames.ncol(df) - 1
        clabels_tmp = String.(names(df)[2:end])
    else
        # layout: rows = channels, columns = time points
        # first column holds channel names; remaining column headers are time
        time_pts_raw = parse.(Float64, names(df)[2:end])
        data = Matrix{Float64}(df[!, 2:end]) # → (ch_n × n_samples)
        ch_n = DataFrames.nrow(df)
        clabels_tmp = string.(df[:, 1])
    end

    # add epoch dimension.
    data = Array(reshape(data, ch_n, size(data, 2), 1))

    # ------------------------------------------------------------------ #
    # sampling rate                                                      #
    # ------------------------------------------------------------------ #
    length(time_pts_raw) >= 2 ||
        throw(ArgumentError("$file_name must contain at least 2 time points."))
    # time is in ms; the mean step is more robust to rounding than the first step alone
    dt = (time_pts_raw[end] - time_pts_raw[1]) / (length(time_pts_raw) - 1)
    dt > 0 || throw(ArgumentError("Time points in $file_name must be increasing."))
    # NeuroAnalyzer stores integer sampling rates (non-integer rates are resampled)
    data, sampling_rate = _integer_rate(data, 1000 / dt)

    # ------------------------------------------------------------------ #
    # Reconstruct a canonical time axis (in seconds) starting at 0       #
    # ------------------------------------------------------------------ #
    n_samples = size(data, 2) * size(data, 3)
    time_pts = _time_axis(n_samples, sampling_rate)
    epoch_time =
        _time_axis(size(data, 2), sampling_rate)

    # ------------------------------------------------------------------ #
    # channel metadata                                                   #
    # ------------------------------------------------------------------ #
    clabels = clean_labels ? _clean_labels(clabels_tmp) : String.(strip.(clabels_tmp))
    ch_type = detect_type ? _set_channel_types(clabels, "eeg") : repeat(["eeg"], ch_n)
    units = [_ch_units(ch_type[idx]) for idx = 1:ch_n]

    markers = DataFrame(
        :id => String[],
        :start => Float64[],
        :length => Float64[],
        :value => String[],
        :channel => Int64[],
    )

    # ------------------------------------------------------------------ #
    # assemble NEURO object                                               #
    # ------------------------------------------------------------------ #
    file_size_mb = round(filesize(file_name) / 1024^2; digits = 2)
    s = _create_subject(;
        id = "", first_name = "", middle_name = "", last_name = "",
        head_circumference = -1, handedness = "", weight = -1, height = -1,
    )
    r = _create_recording_eeg(;
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
