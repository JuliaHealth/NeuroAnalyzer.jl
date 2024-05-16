export import_ncs

"""
    import_ncs(file_name)

Load Neuralinx Continuously Sampled Channels (CSC) and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function import_ncs(file_name::String)

    @assert isfile(file_name) "File $file_name cannot be loaded."
    @assert lowercase(splitext(file_name)[2]) == ".ncs" "This is not NCS file."

    fid = nothing
    try
        fid = open(file_name, "r")
    catch
        error("File $file_name cannot be loaded.")
    end

    file_type = "CSC"

    header = split(_vint2str(_fread(fid, 16 * 1024, :c)), "\r\n")
    header = replace.(header, "\t-"=>"")
    header = replace.(header, "-"=>"")
    header = replace.(header, "\t"=>" ")
    header = replace.(header, "  "=>" ")
    for idx in 1:length(header)
        header[idx][end] == ' ' && (header[idx] = header[idx][1:(end - 1)])
    end
    for idx in length(header):-1:1
        if length(header[idx]) == 0
            deleteat!(header, idx)
        elseif header[idx][1] == '#'
            deleteat!(header, idx)
        end
    end
    header = split.(header, ' ')
    ADBitVolts = 1
    sampling_rate = nothing
    ADChannel = nothing
    ADGain = 1
    AmpGain = 1
    AmpLowCut = 0
    AmpHiCut = 0
    for idx in 1:length(header)
        header[idx][1] == "ADBitVolts" && (ADBitVolts = parse(Float64, header[idx][2]))
        header[idx][1] == "SamplingFrequency" && (sampling_rate = parse(Int64, header[idx][2]))
        header[idx][1] == "ADChannel" && (ADChannel = parse(Int64, header[idx][2]))
        header[idx][1] == "ADGain" && (ADGain = parse(Float64, header[idx][2]))
        header[idx][1] == "AmpGain" && (AmpGain = parse(Float64, header[idx][2]))
        header[idx][1] == "AmpLowCut" && (AmpLowCut = parse(Float64, header[idx][2]))
        header[idx][1] == "AmpHiCut" && (AmpHiCut = parse(Float64, header[idx][2]))
    end

    if AmpLowCut != 0 && AmpHiCut != 0
        filter = "HP: $AmpLowCut Hz, LP: $AmpHiCut Hz"
    else
        filter = ""
    end

    header_size = 16 * 1024
    block_size = 8 + 4 + 4 + 4 + 1024
    data_size = filesize(file_name) - header_size
    n_blocks = div(data_size, block_size)

    data = zeros(1, 512 * n_blocks, 1)
    qwTimeStamp = zeros(Int64, n_blocks)
    dwChannelNumber = zeros(Int64, n_blocks)
    dwSampleFreq = zeros(Int64, n_blocks)
    dwNumValidSamples = zeros(Int64, n_blocks)

    # gain = ADBitVolts * ADGain * AmpGain * 10^6 # μV
    gain = ADBitVolts * 10^3 # mV

    for idx in 1:n_blocks
        qwTimeStamp[idx] = _fread(fid, 1, :ui64)
        dwChannelNumber[idx] = _fread(fid, 1, :ui32)
        dwSampleFreq[idx] = _fread(fid, 1, :ui32)
        dwNumValidSamples[idx] = _fread(fid, 1, :ui32)
        buf = UInt8[]
        readbytes!(fid, buf, dwNumValidSamples[idx] * 2)
        data[1, (1 + (512 * (idx - 1))):(idx * 512), 1] = map(ltoh, reinterpret(Int16, buf)) .* gain
    end

    ch_n = length(unique(dwChannelNumber))
    ch_n > 1 && _warn("Multi-channel files are not implemented yet; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org")

    clabels = ["Ch$ADChannel"]
    ch_type = repeat(["ieeg"], ch_n)
    units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]
    ch_order = _sort_channels(ch_type)

    markers = DataFrame(:id=>String[],
                        :start=>Float64[],
                        :length=>Float64[],
                        :description=>String[],
                        :channel=>Int64[])

    time_pts = round.(collect(0:1/sampling_rate:size(data, 2) * size(data, 3) / sampling_rate)[1:end-1], digits=6)
    ep_time = round.((collect(0:1/sampling_rate:size(data, 2) / sampling_rate))[1:end-1], digits=6)

    file_size_mb = round(filesize(file_name) / 1024^2, digits=2)

    data_type = "ieeg"

    s = _create_subject(id="",
                        first_name="",
                        middle_name="",
                        last_name="",
                        head_circumference=-1,
                        handedness="",
                        weight=-1,
                        height=-1)
    r = _create_recording_eeg(data_type=data_type,
                              file_name=file_name,
                              file_size_mb=file_size_mb,
                              file_type=file_type,
                              recording="",
                              recording_date="",
                              recording_time="",
                              recording_notes="filter: $filter",
                              channel_type=ch_type[ch_order],
                              reference="",
                              clabels=clabels[ch_order],
                              transducers=repeat([""], ch_n),
                              units=units[ch_order],
                              prefiltering=repeat([""], ch_n),
                              sampling_rate=sampling_rate,
                              gain=ones(ch_n))
    e = _create_experiment(name="", notes="", design="")

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = String[]

    locs = _initialize_locs()
    obj = NeuroAnalyzer.NEURO(hdr, time_pts, ep_time, data[ch_order, :, :], components, markers, locs, history)
    _initialize_locs!(obj)

    _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits=2)) s)")

    return obj

end
