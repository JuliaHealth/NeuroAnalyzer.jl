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
    file_type != "EDF" && throw(ArgumentError("File $file_name is not a EDF file."))

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
    if data_records_duration == 0
        _info("This file contains only annotations, use import_edf_annotations().")
    end
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

    header = UInt8[]
    readbytes!(fid, header, ch_n * 16)
    header = String(Char.(header))
    for idx in 1:ch_n
        clabels[idx] = strip(header[1 + ((idx - 1) * 16):(idx * 16)])
    end

    header = UInt8[]
    readbytes!(fid, header, ch_n * 80)
    header = String(Char.(header))
    for idx in 1:ch_n
        transducers[idx] = strip(header[1 + ((idx - 1) * 80):(idx * 80)])
    end

    header = UInt8[]
    readbytes!(fid, header, ch_n * 8)
    header = String(Char.(header))
    for idx in 1:ch_n
        units[idx] = strip(header[1 + ((idx - 1) * 8):(idx * 8)])
    end
    units = replace(lowercase.(units), "uv"=>"μV")

    header = UInt8[]
    readbytes!(fid, header, ch_n * 8)
    header = String(Char.(header))
    for idx in 1:ch_n
        physical_minimum[idx] = parse(Float64, strip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = UInt8[]
    readbytes!(fid, header, ch_n * 8)
    header = String(Char.(header))
    for idx in 1:ch_n
        physical_maximum[idx] = parse(Float64, strip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = UInt8[]
    readbytes!(fid, header, ch_n * 8)
    header = String(Char.(header))
    for idx in 1:ch_n
        digital_minimum[idx] = parse(Float64, strip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = UInt8[]
    readbytes!(fid, header, ch_n * 8)
    header = String(Char.(header))
    for idx in 1:ch_n
        digital_maximum[idx] = parse(Float64, strip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = UInt8[]
    readbytes!(fid, header, ch_n * 80)
    header = String(Char.(header))
    for idx in 1:ch_n
        prefiltering[idx] = strip(header[1 + ((idx - 1) * 80):(idx * 80)])
    end

    header = UInt8[]
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
        units = repeat(["μV"], ch_n)
    end
    if units == repeat([""], ch_n)
        for idx in 1:ch_n
            units[idx] = _get_units(ch_type[idx])
        end
    end

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
                    push!(annotations, String(Char.(signal)))
                    signal = zeros(samples_per_datarecord[idx2])
                else
                    signal = map(ltoh, reinterpret(Int16, signal))
                    if lowercase(units[idx2]) == "mv"
                        lowercase(units[idx2]) == "μV"
                        signal ./= 1000
                    end
                    if lowercase(units[idx2]) == "nv"
                        lowercase(units[idx2]) == "μV"
                        signal .*= 1000
                    end
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
        buf = UInt8[]
        readbytes!(fid, buf, data_size, all=true)
        data_records = length(buf) ÷ 2 ÷ sum(sampling_rate)        
        data = zeros(ch_n, data_records * max_sampling_rate)
        data_segment = max_sampling_rate

        annotations = String[]

        @inbounds for idx1 in 1:data_records            
            for idx2 in 1:ch_n
                signal = UInt8[]
                for idx3 in 1:samples_per_datarecord[idx2] * 2
                    push!(signal, popat!(buf, 1))
                end
                if idx2 in annotation_channels
                    push!(annotations, String(Char.(signal)))
                    signal = zeros(data_segment)
                else
                    signal = map(ltoh, reinterpret(Int16, signal))
                    signal = @. (signal - digital_minimum[idx2]) * gain[idx2] + physical_minimum[idx2]
                    sampling_rate[idx2] != max_sampling_rate && (signal = FourierTools.resample(signal, max_sampling_rate))
                    if lowercase(units[idx2]) == "mv"
                        lowercase(units[idx2]) == "μV"
                        signal ./= 1000
                    end
                    if lowercase(units[idx2]) == "nv"
                        lowercase(units[idx2]) == "μV"
                        signal .*= 1000
                    end
                end
                data[idx2, ((idx1 - 1) * data_segment + 1):idx1 * data_segment] = signal
            end
        end

        _info("Channels upsampled to $max_sampling_rate Hz.")
        sampling_rate = max_sampling_rate
        close(fid)
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
                              transducers=transducers[ch_order],
                              units=units[ch_order],
                              prefiltering=prefiltering[ch_order],
                              sampling_rate=sampling_rate,
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

    _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(channel_n(obj)) × $(epoch_len(obj)) × $(epoch_n(obj)); $(signal_len(obj) / sr(obj)) s)")
    
    return obj
    
end

