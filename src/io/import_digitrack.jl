export import_digitrack

"""
    import_digitrack(file_name; <keyword arguments>)

Load Digitrack ASCII file and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function import_digitrack(file_name::String; detect_type::Bool=true)::NeuroAnalyzer.NEURO

    @assert isfile(file_name) "File $file_name cannot be loaded."

    fid = nothing
    try
        fid = open(file_name, "r")
    catch
        @error "File $file_name cannot be loaded."
    end
    buffer = readline(fid)
    @assert occursin.("Start time ", buffer) "File $file_name is not Digitrack file."

    file_type = "Digitrack"

    patient = ""
    recording = ""
    buffer = replace(buffer, "Start time " => "")
    recording_date = split(buffer, " ")[1]
    recording_time = split(buffer, " ")[2]

    buffer = readline(fid)
    buffer = replace(buffer, "Sampling rate " => "")
    buffer = replace(buffer, "," => ".")
    sampling_rate = round(Int64, parse(Float64, replace(buffer, " Hz" => "")))

    buffer = readline(fid)

    channels = Vector{String}()
    while buffer !=""
        buffer = readline(fid)
        push!(channels, buffer)
    end
    deleteat!(channels, length(channels))
    ch_n  = length(channels)

    clabels = Vector{String}(undef, ch_n)
    prefiltering = Vector{String}(undef, ch_n)
    for idx in 1:ch_n
        clabels[idx] = split(channels[idx], "\t")[1]
        prefiltering_tmp = split(channels[idx], "\t")[2]
        prefiltering[idx] = prefiltering_tmp[1:(length(prefiltering_tmp) - 1)]
    end

    transducers = repeat([""], ch_n)
    units = repeat([""], ch_n)
    gain = repeat([-1.0], ch_n)

    clabels = _clean_labels(clabels)
    if detect_type
        ch_type = _set_channel_types(clabels, "eeg")
    else
        ch_type = repeat(["eeg"], ch_n)
    end
    units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]

    buffer = readlines(fid)

    close(fid)

    data = zeros(ch_n, length(buffer), 1)
    @inbounds for idx in eachindex(buffer)
        signals = split(buffer[idx], "\t")
        deleteat!(signals, length(signals))
        signals = replace.(signals, "," => ".")
        data[:, idx, 1] = parse.(Float64, signals)
    end

    markers = DataFrame(:id=>String[],
                        :start=>Float64[],
                        :length=>Float64[],
                        :value=>String[],
                        :channel=>Int64[])

    time_pts = round.(collect(0:1/sampling_rate:size(data, 2) * size(data, 3) / sampling_rate)[1:end-1], digits=3)
    epoch_time = round.((collect(0:1/sampling_rate:size(data, 2) / sampling_rate))[1:end-1], digits=3)

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
                              recording_date=string(recording_date),
                              recording_time=replace(string(recording_time), '.'=>':'),
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
    obj = NeuroAnalyzer.NEURO(hdr, time_pts, epoch_time, data, components, markers, locs, history)
    _initialize_locs!(obj)

    _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits=2)) s)")

    return obj

end
