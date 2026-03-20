export import_digitrack

"""
    import_digitrack(file_name; <keyword arguments>)

Load a Digitrack ASCII EEG file and return a `NeuroAnalyzer.NEURO` object.

The Digitrack ASCII format has a 3-line text header followed by one line per channel (label and prefiltering, tab-separated), a blank separator line, and then tab-separated signal samples (one column per channel per row).

# Arguments

- `file_name::String`: path to the Digitrack ASCII file
- `detect_type::Bool=true`: infer channel type from channel label

# Returns

- `NeuroAnalyzer.NEURO`

# Throws

- `ArgumentError` if the file does not exist or is not a Digitrack file
"""
function import_digitrack(
    file_name::String;
    detect_type::Bool = true
)::NeuroAnalyzer.NEURO

    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))

    # ------------------------------------------------------------------ #
    # parse fixed-length global header (256 bytes)                       #
    # ------------------------------------------------------------------ #
    # all header reads share a single open/close via the `do` block
    recording_date, recording_time, sampling_rate, clabels, prefiltering, data = 
        open(file_name, "r") do fid

        # line 1: start time
        buffer = readline(fid)
        occursin("Start time ", buffer) ||
            throw(ArgumentError("$file_name is not a Digitrack file."))

        buffer = replace(buffer, "Start time " => "")
        parts = split(buffer, " ")
        recording_date = parts[1]
        recording_time = parts[2]

        # line 2: sampling rate (decimal separator may be a comma)
        buffer = readline(fid)
        buffer = replace(replace(buffer, "Sampling rate " => ""), "," => ".")
        sampling_rate = round(Int64, parse(Float64, replace(buffer, " Hz" => "")))

        # line 3: skip (section separator in Digitrack format)
        readline(fid)

        # channel header block: one "label\tprefiltering\n" line per channel, terminated by a blank line
        channels = String[]
        buffer = readline(fid)
        while buffer != ""
            push!(channels, buffer)
            buffer = readline(fid)
        end

        ch_n = length(channels)
        clabels = Vector{String}(undef, ch_n)
        prefiltering = Vector{String}(undef, ch_n)

        for idx in 1:ch_n
            fields = split(channels[idx], "\t")
            clabels[idx] = fields[1]
            prefiltering[idx] = fields[2]
        end

        # signal data block: one row per time point, tab-separated values,
        # trailing tab on each line (stripped with deleteat! below)
        raw_lines = readlines(fid)
        n_samples = length(raw_lines)
        data = zeros(ch_n, n_samples, 1)

        @inbounds for idx in eachindex(raw_lines)
            signals = split(raw_lines[idx], "\t")
            # remove trailing empty field produced by the trailing tab
            !isempty(signals[end]) || deleteat!(signals, length(signals))
            data[:, idx, 1] = parse.(Float64, replace.(signals, "," => "."))
        end

        recording_date, recording_time, sampling_rate, clabels, prefiltering, data
    end
    # file closed here in all cases, including exceptions

    # ------------------------------------------------------------------ #
    # channel metadata                                                    #
    # ------------------------------------------------------------------ #
    ch_n = length(clabels)
    clabels = _clean_labels(string.(clabels))
    ch_type = detect_type ? _set_channel_types(clabels, "eeg") : repeat(["eeg"], ch_n)
    units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]

    markers = DataFrame(
        :id => String[],
        :start => Float64[],
        :length => Float64[],
        :value => String[],
        :channel => Int64[])

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
        id = "", first_name = "", middle_name = "", last_name = "",
        head_circumference = -1, handedness = "", weight = -1, height = -1)
    r = _create_recording_eeg(
        data_type = "eeg",
        file_name = file_name,
        file_size_mb = file_size_mb,
        file_type = "Digitrack",
        recording = "",
        recording_date = string(recording_date),
        recording_time = replace(string(recording_time), '.' => ':'),
        recording_notes = "",
        channel_type = ch_type,
        channel_order = _sort_channels(ch_type),
        reference = _detect_montage(clabels, ch_type, "eeg"),
        clabels = clabels,
        transducers = repeat([""], ch_n),
        units = units,
        prefiltering = string.(prefiltering),
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