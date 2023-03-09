export import_digitrack

"""
    import_digitrack(file_name; detect_type)

Load Digitrack ASCII file and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `obj::NeuroAnalyzer.NEURO`

# Notes
"""
function import_digitrack(file_name::String; detect_type::Bool=true)
 
    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    fid = ""
    try
        fid = open(file_name, "r")
    catch
        error("File $file_name cannot be loaded.")
    end

    buffer = readline(fid)
    occursin("Start time ", buffer) || throw(ArgumentError("File $file_name is not a Digitrack file."))
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

    data_records = -1
    data_records_duration  = -1

    buffer = readline(fid)

    channels = Vector{String}()
    while buffer !=""
        buffer = readline(fid)
        push!(channels, buffer)
    end
    deleteat!(channels, length(channels))
    ch_n  = length(channels)

    labels = Vector{String}(undef, ch_n)
    prefiltering = Vector{String}(undef, ch_n)
    for idx in 1:ch_n
        labels[idx] = split(channels[idx], "\t")[1]
        prefiltering[idx] = split(channels[idx], "\t")[2]
        prefiltering[idx] = prefiltering[idx][1:(length(prefiltering[idx]) - 1)]
    end

    transducers = repeat([""], ch_n)
    physical_dimension = repeat([""], ch_n)
    gain = repeat([-1.0], ch_n)
    
    labels = _clean_labels(labels)
    if detect_type == true
        channel_type = _set_channel_types(labels)
    else
        channel_type = repeat(["???"], ch_n)
    end
    channel_order = _sort_channels(copy(channel_type))
    has_markers, markers_channel = _has_markers(channel_type)

    data = readlines(fid)

    close(fid)

    data = zeros(ch_n, length(data), 1)
    Threads.@threads for idx in eachindex(data)
        signals = split(data[idx], "\t")
        deleteat!(signals, length(signals))
        signals = replace.(signals, "," => ".")
        @inbounds data[:, idx, 1] = parse.(Float64, signals)
    end

    markers = DataFrame(:id => String[], :start => Int64[], :length => Int64[], :description => String[], :channel => Int64[])
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

