export import_edf

"""
    import_edf(file_name; detect_type)

Load EDF/EDF+ file and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `obj::NeuroAnalyzer.NEURO`

# Notes

- sampling_rate = n.samples ÷ data.record.duration
- gain = (physical maximum - physical minimum) ÷ (digital maximum - digital minimum)
- value = (value - digital minimum ) × gain + physical minimum

# Source

1. Kemp B, Varri A, Rosa AC, Nielsen KD, Gade J. A simple format for exchange of digitized polygraphic recordings. Electroencephalography and Clinical Neurophysiology. 1992 May;82(5):391–3. 
2. Kemp B, Olivan J. European data format ‘plus’(EDF+), an EDF alike standard format for the exchange of physiological data. Clinical Neurophysiology 2003;114:1755–61.
3. https://www.edfplus.info/specs/
"""
function import_edf(file_name::String; detect_type::Bool=true)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    splitext(file_name)[2] == ".edf" || throw(ArgumentError("This is not an EDF file."))

    file_type = ""

    fid = ""
    try
        fid = open(file_name, "r")
    catch
        error("File $file_name cannot be loaded.")
    end

    header = zeros(UInt8, 256)
    readbytes!(fid, header, 256)
    header = String(Char.(header))

    file_type = parse(Int, strip(header[1:8]))
    file_type == 0 && (file_type = "EDF")
    file_type !== "EDF" && throw(ArgumentError("File $file_name is not a EDF file."))

    patient = strip(header[9:88])
    recording = strip(header[89:168])
    # EDF exported from Alice does not conform EDF standard
    occursin("Alice 4", recording) && return import_alice4(file_name, detect_type=detect_type)
    recording_date = header[169:176]
    recording_time = header[177:184]
    data_offset = parse(Int, strip(header[185:192]))
    reserved = strip(header[193:236])
    reserved == "EDF+D" && throw(ArgumentError("EDF+D format (interrupted recordings) is not supported yet."))
    reserved == "EDF+C" && (file_type = "EDF+")
    data_records = parse(Int, strip(header[237:244]))
    data_records_duration  = parse(Float64, strip(header[245:252]))
    ch_n  = parse(Int, strip(header[253:256]))

    clabels = Vector{String}(undef, ch_n)
    transducers = Vector{String}(undef, ch_n)
    units = Vector{String}(undef, ch_n)
    physical_minimum = Vector{Float64}(undef, ch_n)
    physical_maximum = Vector{Float64}(undef, ch_n)
    digital_minimum = Vector{Float64}(undef, ch_n)
    digital_maximum = Vector{Float64}(undef, ch_n)
    prefiltering = Vector{String}(undef, ch_n)
    samples_per_datarecord = Vector{Int64}(undef, ch_n)

    header = zeros(UInt8, ch_n * 16)
    readbytes!(fid, header, ch_n * 16)
    header = String(Char.(header))
    for idx in 1:ch_n
        clabels[idx] = strip(header[1 + ((idx - 1) * 16):(idx * 16)])
    end

    header = zeros(UInt8, ch_n * 80)
    readbytes!(fid, header, ch_n * 80)
    header = String(Char.(header))
    for idx in 1:ch_n
        transducers[idx] = strip(header[1 + ((idx - 1) * 80):(idx * 80)])
    end

    header = zeros(UInt8, ch_n * 8)
    readbytes!(fid, header, ch_n * 8)
    header = String(Char.(header))
    for idx in 1:ch_n
        units[idx] = strip(header[1 + ((idx - 1) * 8):(idx * 8)])
    end
    units = replace(lowercase.(units), "uv"=>"μV")

    header = zeros(UInt8, ch_n * 8)
    readbytes!(fid, header, ch_n * 8)
    header = String(Char.(header))
    for idx in 1:ch_n
        physical_minimum[idx] = parse(Float64, strip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = zeros(UInt8, ch_n * 8)
    readbytes!(fid, header, ch_n * 8)
    header = String(Char.(header))
    for idx in 1:ch_n
        physical_maximum[idx] = parse(Float64, strip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = zeros(UInt8, ch_n * 8)
    readbytes!(fid, header, ch_n * 8)
    header = String(Char.(header))
    for idx in 1:ch_n
        digital_minimum[idx] = parse(Float64, strip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = zeros(UInt8, ch_n * 8)
    readbytes!(fid, header, ch_n * 8)
    header = String(Char.(header))
    for idx in 1:ch_n
        digital_maximum[idx] = parse(Float64, strip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = zeros(UInt8, ch_n * 80)
    readbytes!(fid, header, ch_n * 80)
    header = String(Char.(header))
    for idx in 1:ch_n
        prefiltering[idx] = strip(header[1 + ((idx - 1) * 80):(idx * 80)])
    end

    header = zeros(UInt8, ch_n * 8)
    readbytes!(fid, header, ch_n * 8)
    header = String(Char.(header))
    for idx in 1:ch_n
        samples_per_datarecord[idx] = parse(Int, strip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    close(fid)

    clabels = _clean_labels(clabels)
    if detect_type == true
        channel_type = _set_channel_types(clabels, "eeg")
    else
        channel_type = repeat(["eeg"], ch_n)
    end
    channel_order = _sort_channels(copy(channel_type))
    if file_type == "EDF"
        has_markers = false
        markers = DataFrame(:id=>String[], :start=>Int64[], :length=>Int64[], :description=>String[], :channel=>Int64[])
        markers_channel = -1
    else
        has_markers, markers_channel = _has_markers(channel_type)
        markers = repeat([""], data_records)
    end

    # we assume that all channels have the same sampling rate
    sampling_rate = round(Int64, samples_per_datarecord[1] / data_records_duration)
    gain = Vector{Float64}(undef, ch_n)
    for idx in 1:ch_n
        gain[idx] = (physical_maximum[idx] - physical_minimum[idx]) / (digital_maximum[idx] - digital_minimum[idx])
    end

    fid = ""
    try
        fid = open(file_name, "r")
    catch
        error("File $file_name cannot be loaded.")
    end

    header = zeros(UInt8, data_offset)
    readbytes!(fid, header, data_offset)
    data = zeros(ch_n, samples_per_datarecord[1] * data_records, 1)
    for idx1 in 1:data_records
        for idx2 in 1:ch_n
            signal = zeros(UInt8, samples_per_datarecord[idx2] * 2)
            readbytes!(fid, signal, samples_per_datarecord[idx2] * 2)
            if idx2 != markers_channel
                signal = map(ltoh, reinterpret(Int16, signal))
                if channel_type[idx2] == "mrk"
                    for idx3 in eachindex(signal)
                        if signal[idx3] == digital_minimum[idx2]
                            signal[idx3] = 0
                        else
                            signal[idx3] = 1
                        end
                    end
                    data[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal
                elseif channel_type[idx2] == "events"
                    data[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal
                else
                    data[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal .* gain[idx2]
#=
                    if occursin("uV", units[idx2]) 
                        data[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal .* gain[idx2]
                    elseif occursin("mV", units[idx2])
                        data[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal .* gain[idx2] ./ 1000
                    else
                        data[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal .* gain[idx2]
                    end
=#
                end
            else
                markers[idx1] = String(Char.(signal))
            end
        end
    end
    close(fid)

    if has_markers
        deleteat!(channel_order, vsearch(markers_channel, channel_order))
        data = data[setdiff(1:ch_n, markers_channel), :, :]
        deleteat!(clabels, markers_channel)
        deleteat!(transducers, markers_channel)
        deleteat!(units, markers_channel)
        deleteat!(prefiltering, markers_channel)
        deleteat!(gain, markers_channel)
        ch_n -= 1
        markers = _m2df(markers)
        if markers[!, :id] == repeat([""], nrow(markers))
            # generate unique IDs
            desc = unique(markers[!, :description])
            for idx1 in 1:nrow(markers), idx2 in 1:length(desc)
                markers[idx1, :description] == desc[idx2] && (markers[idx1, :id] = string(idx2))
            end
        end
        markers[!, :start] = t2s.(markers[!, :start], sampling_rate)
        markers[!, :length] = t2s.(markers[!, :length], sampling_rate)
    end

    time_pts = round.(collect(0:1/sampling_rate:size(data, 2) * size(data, 3) / sampling_rate)[1:end-1], digits=3)
    epoch_time = round.((collect(0:1/sampling_rate:size(data, 2) / sampling_rate))[1:end-1], digits=3)
    
    file_size_mb = round(filesize(file_name) / 1024^2, digits=2)
    
    data_type = "eeg"

    s = _create_subject(id="",
                        first_name="",
                        middle_name="",
                        last_name=string(patient),
                        handedness="",
                        weight=-1,
                        height=-1)
    r = _create_recording_eeg(data_type=data_type,
                              file_name=file_name,
                              file_size_mb=file_size_mb,
                              file_type=file_type,
                              recording=string(recording),
                              recording_date=recording_date,
                              recording_time=recording_time,
                              recording_notes="",
                              channel_type=channel_type[channel_order],
                              reference="",
                              clabels=clabels[channel_order],
                              transducers=transducers[channel_order],
                              units=units[channel_order],
                              prefiltering=prefiltering[channel_order],
                              sampling_rate=sampling_rate,
                              gain=gain[channel_order])
    e = _create_experiment(name="", notes="", design="")

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = String[]

    locs = DataFrame(:channel=>Int64,
                     :labels=>String[],
                     :loc_theta=>Float64[],
                     :loc_radius=>Float64[],
                     :loc_x=>Float64[],
                     :loc_y=>Float64[],
                     :loc_z=>Float64[],
                     :loc_radius_sph=>Float64[],
                     :loc_theta_sph=>Float64[],
                     :loc_phi_sph=>Float64[])

    obj = NeuroAnalyzer.NEURO(hdr, time_pts, epoch_time, data[channel_order, :, :], components, markers, locs, history)

    _info("Imported: < " * uppercase(obj.header.recording[:data_type]) * ", $(channel_n(obj)) × $(epoch_len(obj)) × $(epoch_n(obj)) ($(signal_len(obj) / sr(obj)) s) >")

    return obj
    
end

