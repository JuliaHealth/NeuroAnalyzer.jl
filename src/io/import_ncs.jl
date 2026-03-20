export import_ncs

"""
    import_ncs(file_name)

Load a Neuralynx Continuously Sampled Channels (CSC/NCS) file and return a
`NeuroAnalyzer.NEURO` object.

NCS files have a fixed 16 KiB text header followed by data blocks of 1044 bytes each. Each block contains a timestamp, channel number, sampling frequency, valid sample count, and 512 Int16 samples. Blocks with fewer than 512 valid samples are handled correctly.

# Arguments

- `file_name::String`: path to the `.ncs` file

# Returns

- `NeuroAnalyzer.NEURO`

# Throws

- `ArgumentError` if the file does not exist or is not an NCS file
"""
function import_ncs(file_name::String)::NeuroAnalyzer.NEURO

    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))
    lowercase(splitext(file_name)[2]) == ".ncs" ||
        throw(ArgumentError("$file_name is not an NCS file."))

    file_type = "CSC"

    # fixed layout constants (Neuralynx NCS specification)
    HEADER_SIZE = 16 * 1024 # 16 KiB ASCII header
    SAMPLES_PER_BLOCK = 512
    BLOCK_SIZE = 8 + 4 + 4 + 4 + SAMPLES_PER_BLOCK * 2 # = 1044 bytes

    # ------------------------------------------------------------------ #
    # parse ASCII header                                                 #
    # ------------------------------------------------------------------ #
    ADBitVolts = 1.0
    sampling_rate = nothing
    ADChannel = nothing
    ADGain = 1.0
    AmpGain = 1.0
    AmpLowCut = 0.0
    AmpHiCut = 0.0

    qwTimeStamp = Int64[]
    dwChannelNumber = Int64[]
    dwSampleFreq = Int64[]
    dwNumValidSamples = Int64[]
    raw_samples = Int16[]

    open(file_name, "r") do fid

        # parse 16 KiB text header
        header = split(_vint2str(_fread(fid, HEADER_SIZE, :c)), "\r\n")
        header = replace.(header, "\t-" => "", "-" => "", "\t" => " ", "  " => " ")
        header = [rstrip(l) for l in header]
        Base.filter!(l -> !isempty(l) && l[1] != '#', header)
        header = split.(header, ' ')

        for h in header
            length(h) < 2 && continue
            h[1] == "ADBitVolts" && (ADBitVolts = parse(Float64, h[2]))
            h[1] == "SamplingFrequency" && (sampling_rate = parse(Int64,   h[2]))
            h[1] == "ADChannel" && (ADChannel = parse(Int64,   h[2]))
            h[1] == "ADGain" && (ADGain = parse(Float64, h[2]))
            h[1] == "AmpGain" && (AmpGain = parse(Float64, h[2]))
            h[1] == "AmpLowCut" && (AmpLowCut = parse(Float64, h[2]))
            h[1] == "AmpHiCut" && (AmpHiCut = parse(Float64, h[2]))
        end

        isnothing(sampling_rate) &&
            throw(ArgumentError("NCS header in $file_name is missing SamplingFrequency."))

        # ------------------------------------------------------------ #
        # read data blocks                                             #
        # ------------------------------------------------------------ #
        data_size = filesize(file_name) - HEADER_SIZE
        n_blocks  = div(data_size, BLOCK_SIZE)

        for _ in 1:n_blocks
            push!(qwTimeStamp,       _fread(fid, 1, :ui64))
            push!(dwChannelNumber,   _fread(fid, 1, :ui32))
            push!(dwSampleFreq,      _fread(fid, 1, :ui32))
            n_valid = _fread(fid, 1, :ui32)
            push!(dwNumValidSamples, n_valid)

            buf = UInt8[]
            readbytes!(fid, buf, SAMPLES_PER_BLOCK * 2)
            samples = map(ltoh, reinterpret(Int16, buf))
            append!(raw_samples, samples[1:n_valid])
        end
    end
    # file closed here in all cases

    # ------------------------------------------------------------------ #
    # scale to μV                                                        #
    # ------------------------------------------------------------------ #
    gain   = ADBitVolts * 1e6   # ADC LSB → μV
    data   = reshape(Float64.(raw_samples) .* gain, 1, length(raw_samples), 1)

    ch_n  = length(unique(dwChannelNumber))
    ch_n > 1 && _warn(
        "Multi-channel NCS files are not implemented yet; " *
        "please send this file to adam.wysokinski@neuroanalyzer.org")

    # use ADChannel number in label; fall back to "Ch0" if header was absent
    clabels = [isnothing(ADChannel) ? "Ch0" : "Ch$ADChannel"]
    ch_type = repeat(["ieeg"], ch_n)
    units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]  # "μV" for ieeg

    filter_str = (AmpLowCut != 0 && AmpHiCut != 0) ?
        "HP: $AmpLowCut Hz, LP: $AmpHiCut Hz" : ""

    markers = DataFrame(
        :id => String[], :start => Float64[],
        :length => Float64[], :value => String[], :channel => Int64[])

    # -------------------------------------------------------------------- #
    # time axes (6-digit precision — NCS timestamps are microsecond-based) #
    # -------------------------------------------------------------------- #
    n_samples  = size(data, 2) * size(data, 3)
    time_pts   = round.(range(0; step = 1/sampling_rate, length = n_samples);  digits = 6)
    epoch_time = round.(range(0; step = 1/sampling_rate, length = size(data,2)); digits = 6)

    # ------------------------------------------------------------------ #
    # assemble NEURO object                                              #
    # ------------------------------------------------------------------ #
    file_size_mb = round(filesize(file_name) / 1024^2; digits = 2)

    s = _create_subject(id = "",
                        first_name = "",
                        middle_name = "",
                        last_name = "",
                        head_circumference = -1,
                        handedness = "",
                        weight = -1,
                        height = -1)
    r = _create_recording_eeg(
        data_type = "ieeg",
        file_name = file_name,
        file_size_mb = file_size_mb,
        file_type = file_type,
        recording = "",
        recording_date = "",
        recording_time = "",
        recording_notes = isempty(filter_str) ? "" : "filter: $filter_str",
        channel_type = ch_type,
        channel_order = _sort_channels(ch_type),
        reference = "",
        clabels = clabels,
        transducers = repeat([""], ch_n),
        units = units,
        prefiltering = isempty(filter_str) ? repeat([""], ch_n) : repeat([filter_str], ch_n),
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
