export import_edf

"""
    import_edf(file_name; detect_type)

Load EDF/EDF+ file and return `NeuroAnalyzer.RECORD` object.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `eeg::RECORD`

# Notes

- sampling_rate = n.samples / data.record.duration
- gain = (physical_maximum - physical_minimum) / (digital_maximum - digital_minimum)
- value = (value - digital_minimum ) * gain + physical_minimum

# Source

1. Kemp B, Värri A, Rosa AC, Nielsen KD, Gade J. A simple format for exchange of digitized polygraphic recordings. Electroencephalography and Clinical Neurophysiology. 1992 May;82(5):391–3. 
2. Kemp B, Olivan J. European data format ‘plus’(EDF+), an EDF alike standard format for the exchange of physiological data. Clinical Neurophysiology 2003;114:1755–61.
3. https://www.edfplus.info/specs/
"""
function import_edf(file_name::String; detect_type::Bool=true)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    filetype = ""

    fid = ""
    try
        fid = open(file_name, "r")
    catch
        error("File $file_name cannot be loaded.")
    end

    header = zeros(UInt8, 256)
    readbytes!(fid, header, 256)
    header = String(Char.(header))

    filetype = parse(Int, strip(header[1:8]))
    filetype == 0 && (filetype = "EDF")
    filetype !== "EDF" && throw(ArgumentError("File $file_name is not a EDF file."))

    patient = strip(header[9:88])
    recording = strip(header[89:168])
    # EDF exported from Alice does not conform EDF standard
    occursin("Alice 4", recording) && return import_alice4(file_name, detect_type=detect_type)
    recording_date = header[169:176]
    recording_time = header[177:184]
    recording_time = replace(recording_time, '.'=>':')
    data_offset = parse(Int, strip(header[185:192]))
    reserved = strip(header[193:236])
    reserved == "EDF+D" && throw(ArgumentError("EDF+D format (interrupted recordings) is not supported yet."))
    reserved == "EDF+C" && (filetype = "EDF+")
    data_records = parse(Int, strip(header[237:244]))
    data_records_duration  = parse(Float64, strip(header[245:252]))
    channel_n  = parse(Int, strip(header[253:256]))

    labels = Vector{String}(undef, channel_n)
    transducers = Vector{String}(undef, channel_n)
    physical_dimension = Vector{String}(undef, channel_n)
    physical_minimum = Vector{Float64}(undef, channel_n)
    physical_maximum = Vector{Float64}(undef, channel_n)
    digital_minimum = Vector{Float64}(undef, channel_n)
    digital_maximum = Vector{Float64}(undef, channel_n)
    prefiltering = Vector{String}(undef, channel_n)
    samples_per_datarecord = Vector{Int64}(undef, channel_n)

    header = zeros(UInt8, channel_n * 16)
    readbytes!(fid, header, channel_n * 16)
    header = String(Char.(header))
    for idx in 1:channel_n
        labels[idx] = strip(header[1 + ((idx - 1) * 16):(idx * 16)])
    end

    header = zeros(UInt8, channel_n * 80)
    readbytes!(fid, header, channel_n * 80)
    header = String(Char.(header))
    for idx in 1:channel_n
        transducers[idx] = strip(header[1 + ((idx - 1) * 80):(idx * 80)])
    end

    header = zeros(UInt8, channel_n * 8)
    readbytes!(fid, header, channel_n * 8)
    header = String(Char.(header))
    for idx in 1:channel_n
        physical_dimension[idx] = strip(header[1 + ((idx - 1) * 8):(idx * 8)])
    end

    header = zeros(UInt8, channel_n * 8)
    readbytes!(fid, header, channel_n * 8)
    header = String(Char.(header))
    for idx in 1:channel_n
        physical_minimum[idx] = parse(Float64, strip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = zeros(UInt8, channel_n * 8)
    readbytes!(fid, header, channel_n * 8)
    header = String(Char.(header))
    for idx in 1:channel_n
        physical_maximum[idx] = parse(Float64, strip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = zeros(UInt8, channel_n * 8)
    readbytes!(fid, header, channel_n * 8)
    header = String(Char.(header))
    for idx in 1:channel_n
        digital_minimum[idx] = parse(Float64, strip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = zeros(UInt8, channel_n * 8)
    readbytes!(fid, header, channel_n * 8)
    header = String(Char.(header))
    for idx in 1:channel_n
        digital_maximum[idx] = parse(Float64, strip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = zeros(UInt8, channel_n * 80)
    readbytes!(fid, header, channel_n * 80)
    header = String(Char.(header))
    for idx in 1:channel_n
        prefiltering[idx] = strip(header[1 + ((idx - 1) * 80):(idx * 80)])
    end

    header = zeros(UInt8, channel_n * 8)
    readbytes!(fid, header, channel_n * 8)
    header = String(Char.(header))
    for idx in 1:channel_n
        samples_per_datarecord[idx] = parse(Int, strip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    close(fid)

    labels = _clean_labels(labels)
    if detect_type == true
        channel_type = _set_channel_types(labels)
    else
        channel_type = repeat(["???"], channel_n)
    end
    channel_order = _sort_channels(copy(channel_type))

    if filetype == "EDF"
        has_markers = false
        markers = DataFrame(:id=>String[], :start=>Int64[], :length=>Int64[], :description=>String[], :channel=>Int64[])
        markers_channel = -1
    else
        has_markers, markers_channel = _has_markers(channel_type)
        markers = repeat([""], data_records)
    end

    # we assume that all channels have the same sampling rate
    sampling_rate = round(Int64, samples_per_datarecord[1] / data_records_duration)
    gain = Vector{Float64}(undef, channel_n)
    for idx in 1:channel_n
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
    data = zeros(channel_n, samples_per_datarecord[1] * data_records, 1)
    for idx1 in 1:data_records
        for idx2 in 1:channel_n
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
                    if occursin("uV", physical_dimension[idx2]) 
                        data[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal .* gain[idx2]
                    elseif occursin("mV", physical_dimension[idx2])
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
        data = data[setdiff(1:channel_n, markers_channel), :, :]
        deleteat!(labels, markers_channel)
        deleteat!(transducers, markers_channel)
        deleteat!(physical_dimension, markers_channel)
        deleteat!(prefiltering, markers_channel)
        deleteat!(gain, markers_channel)
        channel_n -= 1
        markers = _m2df(markers)
        markers[!, :start] = t2s.(markers[!, :start], sampling_rate)
        markers[!, :length] = t2s.(markers[!, :length], sampling_rate)
    end

    duration_samples = size(data, 2)
    duration_seconds = size(data, 2) / sampling_rate
    time_pts = collect(0:(1 / sampling_rate):duration_seconds)
    time_pts = time_pts[1:end - 1]
    file_size_mb = round(filesize(file_name) / 1024^2, digits=2)

    hdr = _create_header(last_name=string(patient), data_type="eeg", file_name=file_name, file_size_mb=file_size_mb, file_type=filetype, recording=string(recording), recording_date=recording_date, recording_time=recording_time, recording_notes="", channel_n=channel_n, channel_type=channel_type[channel_order], reference="", duration_samples=duration_samples, duration_seconds=duration_seconds, epoch_n=1, epoch_duration_samples=duration_samples, epoch_duration_seconds=duration_seconds, labels=labels[channel_order], units=physical_dimension[channel_order], prefiltering=prefiltering[channel_order], sampling_rate=sampling_rate, gain=gain[channel_order], markers=has_markers, history=[""], components=Symbol[], locations=false)

    components = Vector{Any}()
    epoch_time = time_pts
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

    return NeuroAnalyzer.RECORD(hdr, time_pts, epoch_time, data[channel_order, :, :], components, markers, locs)
end

