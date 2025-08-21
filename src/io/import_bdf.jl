export import_bdf

"""
    import_bdf(file_name; <keyword arguments>)

Load BDF/BDF+ file and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on channel label

# Returns

- `obj::NeuroAnalyzer.NEURO`

# Notes

- sampling_rate = n.samples ÷ data.record.duration
- gain = (physical maximum - physical minimum) ÷ (digital maximum - digital minimum)
- value = (value - digital minimum ) × gain + physical minimum

# Source

https://www.biosemi.com/faq/file_format.htm
"""
function import_bdf(file_name::String; detect_type::Bool=true)::NeuroAnalyzer.NEURO

    @assert isfile(file_name) "File $file_name cannot be loaded."
    @assert lowercase(splitext(file_name)[2]) == ".bdf" "This is not BDF file."

    file_type = ""

    fid = nothing
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
    file_type == "BDF" && @assert strip(header[3:9]) == "BIOSEMI" "File $file_name is not BDF file."

    patient = strip(header[10:89])
    recording = strip(header[90:169])
    recording_date = header[170:177]
    recording_time = header[178:185]
    data_offset = parse(Int, strip(header[186:192]))
    reserved  = strip(header[193:236])
    @assert reserved != "BDF+D" "BDF+D format (interrupted recordings) is not supported yet; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org"
    reserved == "BDF+C" && (file_type = "BDF+")
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
    if detect_type
        ch_type = _set_channel_types(clabels, "eeg")
    else
        ch_type = repeat(["eeg"], ch_n)
        ch_type[clabels .== "Status"] .= "mrk"
    end
    units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]

    if file_type == "BDF"
        # in BDF files the last channel is always the Status channel
        annotation_channels = Int64[]
        markers_channel = ch_n
    else
        # in BDF+ files the last channel is always the Status channel + additional annotations channels are possible
        annotation_channels = sort(unique(vcat(ch_n, getindex.(findall(occursin.("annotation", lowercase.(clabels))), 1))))
        markers_channel = getindex.(findall(ch_type .== "mrk"), 1)
    end

    sampling_rate = round(Int64, samples_per_datarecord[1] / data_records_duration)
    gain = @. (physical_maximum - physical_minimum) / (digital_maximum - digital_minimum)

    fid = nothing
    try
        fid = open(file_name, "r")
    catch
        error("File $file_name cannot be loaded.")
    end
    header = zeros(UInt8, data_offset)
    readbytes!(fid, header, data_offset)
    data = zeros(ch_n, samples_per_datarecord[1] * data_records, 1)
    annotations = String[]
    @inbounds for idx1 in 1:data_records
        for idx2 in 1:ch_n
            signal24 = zeros(UInt8, samples_per_datarecord[idx2] * 3)
            readbytes!(fid, signal24, samples_per_datarecord[idx2] * 3)
            signal = Vector{Float64}()
            if idx2 in annotation_channels
                push!(annotations, String(Char.(signal24)))
                signal = zeros(samples_per_datarecord[idx2])
            else
                if idx2 in markers_channel
                    # status = Vector{Int64}()
                    for byte_idx in 1:3:length(signal24)
                        b1 = Int32(signal24[byte_idx])
                        b2 = Int32(signal24[byte_idx + 1])
                        # b3 = Int64(signal24[byte_idx + 2])
                        # push!(signal, Float64(b1 | b2) * gain[idx2])
                        push!(signal, Float64(b1 | b2))
                        # push!(status, b3 * gain[idx2])
                    end
                else
                    for byte_idx in 1:3:length(signal24)
                        b1 = Int32(signal24[byte_idx]) << 8
                        b2 = Int32(signal24[byte_idx + 1]) << 16
                        b3 = -Int32(-signal24[byte_idx + 2]) << 24
                        # push!(signal, Float64(((b1 | b2 | b3) >> 8) * gain[idx2]))
                        push!(signal, Float64(((b1 | b2 | b3) >> 8)))
                    end
                end
            end
            data[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal
        end
    end

    data .*= gain

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

    close(fid)

    if length(annotation_channels) == 0
        markers = DataFrame(:id=>String[],
                            :start=>Float64[],
                            :length=>Float64[],
                            :value=>String[],
                            :channel=>Int64[])
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
                              recording_time=replace(recording_time, '.'=>':'),
                              recording_notes="",
                              channel_type=ch_type,
                              channel_order=_sort_channels(ch_type),
                              reference=_detect_montage(clabels, ch_type, data_type),
                              clabels=clabels,
                              transducers=transducers,
                              units=units,
                              prefiltering=prefiltering,
                              line_frequency=50,
                              sampling_rate=sampling_rate,
                              gain=gain,
                              bad_channels=zeros(Bool, size(data, 1), 1))
    e = _create_experiment(name="", notes="", design="")

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = String[]

    locs = _initialize_locs()
    obj = NeuroAnalyzer.NEURO(hdr, time_pts, ep_time, data, components, markers, locs, history)
    _initialize_locs!(obj)

    _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits=3)) s)")

    return obj

end
