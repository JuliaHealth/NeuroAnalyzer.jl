export import_alice4

"""
    import_alice4(file_name; detect_type)

Load EDF exported from Alice 4 Polysomnography System and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `obj::NeuroAnalyzer.NEURO`

# Notes

- EDF files exported from Alice 4 have incorrect value of `data_records` (-1) and multiple sampling rate; channels are upsampled to the highest rate.
"""
function import_alice4(file_name::String; detect_type::Bool=true)

    @assert isfile(file_name) "File $file_name cannot be loaded."

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
    @assert file_type == "EDF" "File $file_name is not EDF file."

    patient = strip(header[9:88])
    recording = strip(header[89:168])
    @assert occursin("Alice 4", recording) "This is not Alice 4 EDF file."
    recording_date = header[169:176]
    recording_time = header[177:184]
    data_offset = parse(Int, strip(header[185:192]))
    reserved = strip(header[193:236])
    @assert reserved != "EDF+D" "EDF+D format (interrupted recordings) is not supported."
    reserved == "EDF+C" && (file_type = "EDF+")
    # we get -1 here
    data_records = parse(Int, strip(header[237:244]))
    @assert data_records == -1 "This seems to be a regular EDF file, use import_edf()."
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
        ch_type = _set_channel_types(clabels, "eeg")
    else
        ch_type = repeat(["eeg"], ch_n)
    end
    units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]

    if file_type == "EDF"
        annotation_channels = Int64[]
        markers_channel = []
    else
        annotation_channels = sort(getindex.(findall(occursin.("annotation", lowercase.(clabels))), 1))
        markers_channel = getindex.(findall(ch_type .== "mrk"), 1)
    end

    if length(unique(samples_per_datarecord)) == 1
        sampling_rate = round(Int64, samples_per_datarecord[1] / data_records_duration)
    else
        sampling_rate = round.(Int64, samples_per_datarecord / data_records_duration)
    end

    gain = @. (physical_maximum - physical_minimum) / (digital_maximum - digital_minimum)

    fid = ""
    try
        fid = open(file_name, "r")
    catch
        error("File $file_name cannot be loaded.")
    end

    if sampling_rate isa Int64
        header = zeros(UInt8, data_offset)
        readbytes!(fid, header, data_offset)
        data = zeros(ch_n, samples_per_datarecord[1] * data_records, 1)
        annotations = String[]
        @inbounds for idx1 in 1:data_records
            for idx2 in 1:ch_n
                signal = zeros(UInt8, samples_per_datarecord[idx2] * 2)
                readbytes!(fid, signal, samples_per_datarecord[idx2] * 2)
                if idx2 in annotation_channels
                    push!(annotations, String(Char.(signal24)))
                    signal = zeros(samples_per_datarecord[idx2])
                else
                    signal = map(ltoh, reinterpret(Int16, signal))
                end
                data[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal .* gain[idx2]
            end
        end
        close(fid)
    else
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
        annotations = String[]

        @inbounds for idx1 in 1:data_records            
            for idx2 in 1:ch_n
                tmp = Vector{Float64}()
                for idx3 in 1:sampling_rate[idx2]
                    push!(tmp, popat!(signal, 1))
                end
                # tmp = @. (tmp - digital_minimum[idx2]) * gain[idx2] + physical_minimum[idx2]
                tmp .*= gain[idx2]
                if sampling_rate[idx2] == max_sampling_rate
                    data[idx2, ((idx1 - 1) * data_segment + 1):idx1 * data_segment] = tmp
                else
                    tmp_upsampled = FourierTools.resample(tmp, max_sampling_rate)
                    data[idx2, ((idx1 - 1) * data_segment + 1):idx1 * data_segment] = FourierTools.resample(tmp, max_sampling_rate)
                end
            end
        end

        _info("Channels upsampled to $max_sampling_rate Hz")
        sampling_rate = max_sampling_rate
        close(fid)
    end

    # convert nV/mV to μV
    @inbounds for idx in 1:ch_n
        units[idx] == "" && (units[idx] = "μV")
        if ch_type[idx] == "eeg"
            if lowercase(units[idx]) == "mv"
                lowercase(units[idx]) == "μV"
                data[idx, :] .*= 1000
            end
            if lowercase(units[idx]) == "nv"
                lowercase(units[idx]) == "μV"
                data[idx, :] ./= 1000
            end
        end
    end

    if length(annotation_channels) == 0
        markers = DataFrame(:id=>String[], :start=>Int64[], :length=>Int64[], :description=>String[], :channel=>Int64[])
    else
        markers = _a2df(annotations)
        deleteat!(ch_type, annotation_channels)
        deleteat!(transducers, annotation_channels)
        deleteat!(units, annotation_channels)
        deleteat!(prefiltering, annotation_channels)
        deleteat!(clabels, annotation_channels)
        data = data[setdiff(collect(1:ch_n), annotation_channels), :, :]
        ch_n -= length(annotation_channels)
    end

    ch_order = _sort_channels(ch_type)

    time_pts = round.(collect(0:1/sampling_rate:size(data, 2) * size(data, 3) / sampling_rate)[1:end-1], digits=3)
    ep_time = round.((collect(0:1/sampling_rate:size(data, 2) / sampling_rate))[1:end-1], digits=3)

    file_size_mb = round(filesize(file_name) / 1024^2, digits=2)

    data_type = "eeg"

    s = _create_subject(id="",
                        first_name="",
                        middle_name="",
                        last_name=string(patient),
                        head_circumference=-1,
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
                              channel_type=ch_type[ch_order],
                              reference="",
                              clabels=clabels[ch_order],
                              units=units[ch_order],
                              transducers=transducers[ch_order],
                              prefiltering=prefiltering[ch_order],
                              sampling_rate=max_sampling_rate,
                              gain=gain[ch_order])
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

    obj = NeuroAnalyzer.NEURO(hdr, time_pts, ep_time, data[ch_order, :, :], components, markers, locs, history)

    _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(obj.time_pts[end]) s)")

    return obj
    
end
