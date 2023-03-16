export import_alice4

"""
    import_alice4(file_name; detect_type)

Load EDF exported from Alice 4 Polysomnography System and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `::NeuroAnalyzer.NEURO`

# Notes

- EDF files exported from Alice 4 have incorrect value of `data_records` (-1) and multiple sampling rate; channels are upsampled to the highest rate.
"""
function import_alice4(file_name::String; detect_type::Bool=true)

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

    file_type = parse(Int, strip(header[1:8]))
    file_type == 0 && (file_type = "EDF")
    file_type !== "EDF" && throw(ArgumentError("File $file_name is not a EDF file."))

    patient = strip(header[9:88])
    recording = strip(header[89:168])
    occursin("Alice 4", recording) == false && throw(ArgumentError("This is not Alice 4 EDF file."))
    recording_date = header[169:176]
    recording_time = header[177:184]
    data_offset = parse(Int, strip(header[185:192]))
    reserved = strip(header[193:236])
    reserved == "EDF+D" && throw(ArgumentError("EDF+D format (interrupted recordings) is not supported."))
    reserved == "EDF+C" && (file_type = "EDF+")
    # we get -1 here
    data_records = parse(Int, strip(header[237:244]))
    # we get 1.0 here
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
        channel_type = _set_channel_types(clabels)
    else
        channel_type = repeat(["???"], ch_n)
    end
    channel_order = _sort_channels(copy(channel_type))

    if file_type == "EDF"
        markers = DataFrame(:id=>String[], :start=>Int64[], :length=>Int64[], :description=>String[], :channel=>Int64[])
        markers_channel = -1
    else
        has_markers, markers_channel = _has_markers(channel_type)
        markers = repeat([""], data_records)
    end

    gain = Vector{Float64}(undef, ch_n)
    for idx in 1:ch_n
        gain[idx] = (physical_maximum[idx] - physical_minimum[idx]) / (digital_maximum[idx] - digital_minimum[idx])
    end

    if length(unique(samples_per_datarecord)) == 1
        sampling_rate = round(Int64, samples_per_datarecord[1] / data_records_duration)

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
    else
        sampling_rate = round.(Int64, samples_per_datarecord / data_records_duration)
        max_sampling_rate = maximum(sampling_rate)

        fid = ""
        try
            fid = open(file_name, "r")
        catch
            error("File $file_name cannot be loaded.")
        end

        header = zeros(UInt8, data_offset)
        readbytes!(fid, header, data_offset)

        data_size = filesize(file_name) - data_offset
        data = zeros(UInt8, data_size)
        readbytes!(fid, data, data_size, all=true)
        signal = map(ltoh, reinterpret(Int16, data))
        data_records = length(signal) ÷ sum(sampling_rate)        
        data = zeros(ch_n, data_records * max_sampling_rate)
        data_segment = max_sampling_rate

        @inbounds for idx1 in 1:data_records            
            for idx2 in 1:ch_n
                tmp = Vector{Float64}()
                for idx3 in 1:sampling_rate[idx2]
                    push!(tmp, popat!(signal, 1))
                end
                tmp = @. (tmp - digital_minimum[idx2]) * gain[idx2] + physical_minimum[idx2]
                if sampling_rate[idx2] == max_sampling_rate
                    data[idx2, ((idx1 - 1) * data_segment + 1):idx1 * data_segment] = tmp
                else
                    tmp_upsampled = FourierTools.resample(tmp, max_sampling_rate)
                    data[idx2, ((idx1 - 1) * data_segment + 1):idx1 * data_segment] = tmp_upsampled
                end
            end
        end

        # reject weird channels
        for idx1 in 1:ch_n
            if idx1 != markers_channel
                if channel_type[idx1] == "mrk"
                    for idx2 in 1:size(data, 2)
                        if signal[idx1, idx2] == digital_minimum[idx1]
                            signal[idx1, idx2] = 0
                        else
                            signal[idx1, idx2] = 1
                        end
                    end
                end
                # if occursin("mV", units[idx1])
                #     data ./= 1000
                # end
            else
                markers[idx1] = String(Char.(signal))
            end
        end
        _info("Channels upsampled to $max_sampling_rate Hz.")
        sampling_rate = max_sampling_rate
        close(fid)
    end

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
        markers[!, :start] = t2s.(markers[!, :start], sampling_rate)
        markers[!, :length] = t2s.(markers[!, :length], sampling_rate)
    else
        markers = DataFrame(:id=>String[], :start=>Int64[], :length=>Int64[], :description=>String[], :channel=>Int64[])
    end

    time_pts = collect(0:(1 / sampling_rate):((size(data, 2) * size(data, 3)) / sampling_rate))
    time_pts = round.(linspace(time_pts[1], time_pts[end], size(data, 2) * size(data, 3)), digits=4)
    epoch_time = round.(linspace(time_pts[1], time_pts[end], size(data, 2)), digits=4)

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
                              units=units[channel_order],
                              transducers=transducers[channel_order],
                              prefiltering=prefiltering[channel_order],
                              sampling_rate=max_sampling_rate,
                              gain=gain[channel_order])

    e = _create_experiment(experiment_name="",
                           experiment_notes="",
                           experiment_design="")

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

    return NeuroAnalyzer.NEURO(hdr, time_pts, epoch_time, data[channel_order, :, :], components, markers, locs, history)
    
end
