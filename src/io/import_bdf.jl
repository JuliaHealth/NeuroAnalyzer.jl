export import_bdf

"""
    import_bdf(file_name; detect_type)

Load BDF/BDF+ file and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `obj::NeuroAnalyzer.NEURO`

# Notes

- sampling_rate = n.samples / data.record.duration
- gain = (physical_maximum - physical_minimum) / (digital_maximum - digital_minimum)
- value = (value - digital_minimum ) * gain + physical_minimum

# Source

https://www.biosemi.com/faq/file_format.htm
"""
function import_bdf(file_name::String; detect_type::Bool=true)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

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

    file_type = Int(header[1])
    file_type == 255 && (file_type = "BDF")
    (file_type !== "BDF" && strip(header[3:9]) !== "BIOSEMI") && throw(ArgumentError("File $file_name is not a BDF file."))

    patient = strip(header[10:89])
    recording = strip(header[90:169])
    recording_date = header[170:177]
    recording_time = header[178:185]
    data_offset = parse(Int, strip(header[186:192]))
    reserved  = strip(header[193:236])
    reserved == "BDF+D" && throw(ArgumentError("BDF+D format (interrupted recordings) is not supported yet."))
    reserved == "BDF+C" && (file_type = "BDF+")
    data_records = parse(Int, strip(header[237:244]))
    data_records_duration  = parse(Float64, strip(header[245:252]))
    ch_n  = parse(Int, strip(header[253:256]))

    labels = Vector{String}(undef, ch_n)
    transducers = Vector{String}(undef, ch_n)
    physical_dimension = Vector{String}(undef, ch_n)
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
        labels[idx] = strip(header[1 + ((idx - 1) * 16):(idx * 16)])
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
        physical_dimension[idx] = strip(header[1 + ((idx - 1) * 8):(idx * 8)])
    end

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

    labels = _clean_labels(labels)
    if detect_type == true
        channel_type = _set_channel_types(labels)
    else
        channel_type = repeat(["???"], ch_n)
    end
    channel_order = _sort_channels(copy(channel_type))
    has_markers, markers_channel = _has_markers(channel_type)

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
    markers = repeat([""], data_records)
    for idx1 in 1:data_records
        for idx2 in 1:ch_n
            signal24 = zeros(UInt8, samples_per_datarecord[idx2] * 3)
            readbytes!(fid, signal24, samples_per_datarecord[idx2] * 3)
            if idx2 != markers_channel
                signal = Vector{Float64}()
                for byte_idx in 1:3:length(signal24)
                    b1 = Int32(signal24[byte_idx]) << 8
                    b2 = Int32(signal24[byte_idx + 1]) << 16
                    b3 = -Int32(-signal24[byte_idx + 2]) << 24
                    push!(signal, Float64(((b1 | b2 | b3) >> 8) * gain[idx2]))
                end
                if channel_type[idx2] == "markers"
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
                    if occursin("uV", physical_dimension[idx2]) 
                        data[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal .* gain[idx2]
                    elseif occursin("mV", physical_dimension[idx2])
                        data[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal .* gain[idx2] ./ 1000
                    else
                        data[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal .* gain[idx2]
                    end
                end
            else
                markers[idx1] = String(Char.(signal24))
            end
        end
    end
    close(fid)
    
    if has_markers
        deleteat!(channel_order, vsearch(markers_channel, channel_order))
        data = data[setdiff(1:ch_n, markers_channel), :, :]
        deleteat!(channel_type, markers_channel)
        deleteat!(labels, markers_channel)
        deleteat!(transducers, markers_channel)
        deleteat!(physical_dimension, markers_channel)
        deleteat!(prefiltering, markers_channel)
        deleteat!(gain, markers_channel)
        ch_n -= 1
        markers = _m2df(markers)
        # convert markers time to samples
        markers[!, :start] = t2s.(markers[!, :start], sampling_rate)
        markers[!, :length] = t2s.(markers[!, :length], sampling_rate)
    else
        markers = DataFrame(:id => String[], :start => Int64[], :length => Int64[], :description => String[], :channel => Int64[])
    end

    duration_samples = size(data, 2)
    duration_seconds = size(data, 2) / sampling_rate
    time_pts = collect(0:(1 / sampling_rate):duration_seconds)
    time_pts = time_pts[1:end - 1]
    file_size_mb = round(filesize(file_name) / 1024^2, digits=2)

    hdr = Dict(:data_type => "eeg",
                      :file_name => file_name,
                      :file_size_mb => file_size_mb,
                      :file_type => file_type,
                      :patient => string(patient),
                      :recording => string(recording),
                      :recording_date => recording_date,
                      :recording_time => recording_time,
                      :ch_n => ch_n,
                      :channel_type => channel_type,
                      :reference => "",
                      :channel_locations => false,
                      :history => String[],
                      :components => Symbol[],
                      :duration_samples => duration_samples,
                      :duration_seconds => duration_seconds,
                      :epoch_n => 1,
                      :epoch_duration_samples => duration_samples,
                      :epoch_duration_seconds => duration_seconds,
                      :labels => labels[channel_order],
                      :transducers => transducers[channel_order],
                      :units => physical_dimension[channel_order],
                      :prefiltering => prefiltering[channel_order],
                      :sampling_rate => sampling_rate,
                      :gain => gain[channel_order],
                      :note => "",
                      :markers => has_markers)

    components = Vector{Any}()
    epoch_time = time_pts
    locs = DataFrame(:channel => Int64,
                     :labels => String[],
                     :loc_theta => Float64[],
                     :loc_radius => Float64[],
                     :loc_x => Float64[],
                     :loc_y => Float64[],
                     :loc_z => Float64[],
                     :loc_radius_sph => Float64[],
                     :loc_theta_sph => Float64[],
                     :loc_phi_sph => Float64[])

    eeg = NeuroAnalyzer.NEURO(hdr, time_pts, epoch_time, data[channel_order, :, :], components, markers, locs)

    return eeg
end

