export import_xdf

"""
    import_xdf(file_name)

Load Extensible Data Format (XDF) and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function import_xdf(file_name::String)

    @assert isfile(file_name) "File $file_name cannot be loaded."
    @assert lowercase(splitext(file_name)[2]) == ".xdf" "This is not XDF file."

    streams = nothing
    try
        streams = read_xdf(file_name)
    catch
        error("File $file_name cannot be loaded.")
    end

    s_names = collect(keys(streams))

    n_ch = 0
    for s_idx in s_names
        n_ch += streams[s_idx]["nchannels"]
    end
    @assert n_ch > 0 "No channels data found in the $file_name."

    stream_nchannels = zeros(length(s_names))
    stream_name = repeat([""], length(s_names))
    stream_srate = zeros(length(s_names))
    stream_type = repeat([""], length(s_names))
    for s_idx in eachindex(s_names)
        stream_nchannels[s_idx] = streams[s_names[s_idx]]["nchannels"]
        stream_name[s_idx] = streams[s_names[s_idx]]["name"]
        stream_srate[s_idx] = Float64(streams[s_names[s_idx]]["srate"])
        stream_type[s_idx] = streams[s_names[s_idx]]["type"]
    end

    time = Vector{Vector{Float64}}(undef, length(s_names))
    clock = Vector{Vector{Float64}}(undef, length(s_names))
    offset = Vector{Vector{Float64}}(undef, length(s_names))
    data = Vector{Matrix{Any}}(undef, length(s_names))
    for s_idx in eachindex(s_names)
        time[s_idx] = streams[s_names[s_idx]]["time"]
        clock[s_idx] = streams[s_names[s_idx]]["clock"]
        offset[s_idx] = streams[s_names[s_idx]]["offset"]
        data[s_idx] = streams[s_names[s_idx]]["data"]
    end

    # find EEG stream
    eeg_idx = findall(n -> n == "EEG", stream_type)
    other_idx = findall(n -> n != "EEG", stream_type)
    @assert length(eeg_idx) > 0 "EEG streams not found in the $file_name."
    @assert length(eeg_idx) == 1 _info("Importing files with > 1 EEG streams is not implemented yet; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org")

    eeg_idx = eeg_idx[1]
    sampling_rate = round(Int64, streams[s_names[eeg_idx]]["srate"])
    time_pts0 = time[eeg_idx][1]
    time_pts = round.(Float64.(time[eeg_idx]) .- time_pts0, digits=3)
    ep_time = time_pts
    eeg_data = reshape(Float64.(data[eeg_idx]'), size(data[eeg_idx], 2), :, 1)
    ch_n = size(eeg_data, 1)
    clabels = repeat([""], ch_n)
    ch_type = repeat(["eeg"], ch_n)
    units = repeat(["μV"], ch_n)
    for idx in 1:ch_n
        clabels[idx] = streams[s_names[eeg_idx]]["name"] * "-$idx"
    end

    markers = DataFrame(:id=>String[],
                        :start=>Float64[],
                        :length=>Float64[],
                        :description=>String[],
                        :channel=>Int64[])
    for idx in other_idx
        length(streams[s_names[idx]]["data"]) == 0 && break
        for data_idx in 1:streams[s_names[idx]]["nchannels"]
            append!(markers, Dict(:id=>string.(data[idx][:, data_idx]), :start=>(time[idx] .- time[idx][1]), :length=>ones(length(time[idx])), :description=>repeat(["marker"], length(time[idx])), :channel=>zeros(Int64, length(time[idx]))))
        end
    end
    sort!(markers, :start)

    file_type = "XDF"

    file_size_mb = round(filesize(file_name) / 1024^2, digits=2)

    data_type = "eeg"

    s = _create_subject(id="",
                        first_name="",
                        middle_name="",
                        last_name="",
                        head_circumference=-1,
                        handedness="",
                        weight=-1,
                        height=-1)
    r = _create_recording_eeg(data_type=data_type,
                              file_name=file_name,
                              file_size_mb=file_size_mb,
                              file_type=file_type,
                              recording="",
                              recording_date="",
                              recording_time="",
                              recording_notes="",
                              channel_type=ch_type,
                              reference="",
                              clabels=clabels,
                              transducers=repeat([""], ch_n),
                              units=units,
                              prefiltering=repeat([""], ch_n),
                              sampling_rate=sampling_rate,
                              gain=ones(ch_n))
    e = _create_experiment(name="", notes="", design="")

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = String[]

    locs = _initialize_locs()
    obj = NeuroAnalyzer.NEURO(hdr, time_pts, ep_time, eeg_data, components, markers, locs, history)
    _initialize_locs!(obj)
    
    _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits=2)) s)")

    return obj

end
