export import_edf

"""
    import_edf(file_name; <keyword arguments>)

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

1. Kemp B, Varri A, Rosa AC, Nielsen KD, Gade J. A simple format for exchange of digitized polygraphic recordings. Electroencephalography and Clinical Neurophysiology. 1992; 82(5): 391–3
2. Kemp B, Olivan J. European data format ‘plus’(EDF+), an EDF alike standard format for the exchange of physiological data. Clinical Neurophysiology 2003; 114: 1755–61
3. https://www.edfplus.info/specs/
"""
function import_edf(file_name::String; detect_type::Bool=true)::NeuroAnalyzer.NEURO

    @assert isfile(file_name) "File $file_name cannot be loaded."
    @assert lowercase(splitext(file_name)[2]) == ".edf" "This is not EDF file."

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

    file_type = parse(Int, strip(header[1:8]))
    file_type == 0 && (file_type = "EDF")
    @assert file_type == "EDF" "File $file_name is not EDF file."

    patient = strip(header[9:88])
    recording = strip(header[89:168])
    # EDF exported from Alice does not conform EDF standard
    occursin("Alice 4", recording) && return import_alice4(file_name, detect_type=detect_type)
    recording_date = header[169:176]
    recording_time = header[177:184]
    data_offset = parse(Int, strip(header[185:192]))
    reserved = strip(header[193:236])
    @assert reserved != "EDF+D" "EDF+D format (interrupted recordings) is not supported yet; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org"
    reserved == "EDF+C" && (file_type = "EDF+")
    data_records = parse(Int, strip(header[237:244]))
    data_records_duration  = parse(Float64, strip(header[245:252]))
    @assert data_records_duration > 0 "This file contains only annotations, use import_edf_annotations()."
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
    [transducers[idx] = strip(header[1 + ((idx - 1) * 80):(idx * 80)]) for idx in 1:ch_n]

    header = UInt8[]
    readbytes!(fid, header, ch_n * 8)
    header = String(Char.(header))
    [units[idx] = strip(header[1 + ((idx - 1) * 8):(idx * 8)]) for idx in 1:ch_n]
    units = replace(lowercase.(units), "uv"=>"μV")

    header = UInt8[]
    readbytes!(fid, header, ch_n * 8)
    header = String(Char.(header))
    [physical_minimum[idx] = parse(Float64, strip(header[1 + ((idx - 1) * 8):(idx * 8)])) for idx in 1:ch_n]

    header = UInt8[]
    readbytes!(fid, header, ch_n * 8)
    header = String(Char.(header))
    [physical_maximum[idx] = parse(Float64, strip(header[1 + ((idx - 1) * 8):(idx * 8)])) for idx in 1:ch_n]

    header = UInt8[]
    readbytes!(fid, header, ch_n * 8)
    header = String(Char.(header))
    [digital_minimum[idx] = parse(Float64, strip(header[1 + ((idx - 1) * 8):(idx * 8)])) for idx in 1:ch_n]

    header = UInt8[]
    readbytes!(fid, header, ch_n * 8)
    header = String(Char.(header))
    [digital_maximum[idx] = parse(Float64, strip(header[1 + ((idx - 1) * 8):(idx * 8)])) for idx in 1:ch_n]

    header = UInt8[]
    readbytes!(fid, header, ch_n * 80)
    header = String(Char.(header))
    [prefiltering[idx] = strip(header[1 + ((idx - 1) * 80):(idx * 80)]) for idx in 1:ch_n]

    header = UInt8[]
    readbytes!(fid, header, ch_n * 8)
    header = String(Char.(header))
    [samples_per_datarecord[idx] = parse(Int, strip(header[1 + ((idx - 1) * 8):(idx * 8)])) for idx in 1:ch_n]

    close(fid)

    clabels = _clean_labels(clabels)
    ch_type = detect_type ? _set_channel_types(clabels, "eeg") : repeat(["eeg"], ch_n)
    units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]

    if file_type == "EDF"
        annotation_channels = Int64[]
        markers_channel = []
    else
        annotation_channels = sort(getindex.(findall(occursin.("annotation", lowercase.(clabels))), 1))
        markers_channel = getindex.(findall(ch_type .== "mrk"), 1)
    end

    # ignore annotations channels
    if length(unique(samples_per_datarecord[setdiff(1:ch_n, annotation_channels)])) == 1
        sampling_rate = round(Int64, samples_per_datarecord[1] / data_records_duration)
    else
        sampling_rate = round.(Int64, samples_per_datarecord[setdiff(1:ch_n, annotation_channels)] / data_records_duration)
    end

    gain = @. (physical_maximum - physical_minimum) / (digital_maximum - digital_minimum)

    fid = nothing
    try
        fid = open(file_name, "r")
    catch
        error("File $file_name cannot be loaded.")
    end

    if sampling_rate isa Int64
        seek(fid, data_offset)
        data = zeros(ch_n, samples_per_datarecord[1] * data_records, 1)
        annotations = String[]

        @inbounds for idx1 in 1:data_records
            for idx2 in 1:ch_n
                signal = zeros(UInt8, samples_per_datarecord[idx2] * 2)
                readbytes!(fid, signal, samples_per_datarecord[idx2] * 2)
                if idx2 in annotation_channels
                    push!(annotations, String(Char.(signal)))
                    data[idx2, ((idx1 - 1) * samples_per_datarecord[1] + 1):(idx1 * samples_per_datarecord[1]), 1] = zeros(samples_per_datarecord[1])
                else
                    data[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = reinterpret(Int16, signal)
                end
            end
        end

        data .*= gain

        close(fid)
    else
        # ignore annotations channels
        max_sampling_rate = maximum(sampling_rate[setdiff(1:ch_n, annotation_channels)])
        max_samples_per_datarecord = maximum(samples_per_datarecord[setdiff(1:ch_n, annotation_channels)])

        fid = nothing
        try
            fid = open(file_name, "r")
        catch
            error("File $file_name cannot be loaded.")
        end

        seek(fid, data_offset)
        data_size = filesize(file_name) - data_offset
        signal = UInt8[]
        readbytes!(fid, signal, data_size, all=true)
        data_records = length(signal) ÷ 2 ÷ sum(sampling_rate)
        data = zeros(ch_n, data_records * max_sampling_rate)
        data_segment = max_samples_per_datarecord
        annotations = String[]

        @inbounds for idx1 in 1:data_records
            for idx2 in 1:ch_n
                tmp = signal[1:samples_per_datarecord[idx2] * 2]
                deleteat!(signal, 1:samples_per_datarecord[idx2] * 2)
                if idx2 in annotation_channels
                    push!(annotations, String(Char.(tmp)))
                    tmp = zeros(data_segment)
                else
                    tmp = Float64.(reinterpret(Int16, tmp))
                    sampling_rate[idx2] != max_sampling_rate && (tmp = FourierTools.resample(tmp, max_sampling_rate))
                end
                # data[idx2, ((idx1 - 1) * data_segment + 1):idx1 * data_segment] = tmp
                data[idx2, ((idx1 - 1) * data_segment + 1):idx1 * data_segment] = tmp
            end
        end

        data .*= gain

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
