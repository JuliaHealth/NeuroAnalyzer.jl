export import_set

"""
    import_set(file_name; detect_type)

Load SET file (exported from EEGLAB) and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `obj::NeuroAnalyzer.NEURO`

# Source

1. https://eeglab.org/tutorials/ConceptsGuide/Data_Structures.html
"""
function import_set(file_name::String; detect_type::Bool=true)

    @assert isfile(file_name) "File $file_name cannot be loaded."
    @assert splitext(file_name)[2] == ".set" "This is not SET file."

    file_type = "SET"

    dataset = matread(file_name)
    if "EEG" in keys(dataset) && length(keys(dataset)) == 1
        dataset = dataset["EEG"]
    end
    
    data_src = dataset["data"]
    ch_n = Int64(dataset["nbchan"])
    ep_n = 1
    # data in .FTD file
    if data_src isa String
        fid = nothing
        try
            fid = open(data_src, "r")
        catch
            error("File $data_src cannot be loaded.")
        end
        samples_per_channel = length(dataset["times"])
        @assert filesize(data_src) == ch_n * samples_per_channel * 4 "Incorrect file size."
        data = zeros(ch_n, samples_per_channel)
        for ch_idx in 1:ch_n
            buf = UInt8[]
            readbytes!(fid, buf, samples_per_channel * 4)
            data[ch_idx, :] = Float64.(reinterpret(Float32, buf))
        end
    else
        data = dataset["data"]
    end
    # there are no epochs if signal is matrix, not array
    dataset["trials"] isa Float64 && (ep_n = Int(dataset["trials"]))
    if ndims(data) == 3
        ep_n == size(data, 2)
        _info("$ep_n epochs found")
    end
    ndims(data) == 2 && (data = data[:, :, :])

    sampling_rate = round(Int64, dataset["srate"])
    gain = ones(ch_n)
    
    # channels info
    # get channel labels
    if length(dataset["chanlocs"]) > 0 && length(dataset["chanlocs"]["labels"][:]) == ch_n
        clabels = String.(dataset["chanlocs"]["labels"][:])
    else
        clabels = String[]
        for idx in 1:ch_n
            push!(clabels, "ch_$idx")
        end
    end
    # dataset["chaninfo"]
    clabels = _clean_labels(clabels)
    if detect_type == true
        ch_type = _set_channel_types(clabels, "eeg")
        units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]
    else
        if length(dataset["chanlocs"]) > 0 && string.(dataset["chanlocs"]["type"][:]) == repeat([""], ch_n)
            ch_type = repeat(["eeg"], ch_n)
            units = repeat(["μV"], ch_n)
        else
            length(dataset["chanlocs"]) > 0 && (ch_type = lowercase.(string.(dataset["chanlocs"]["type"][:])))
            units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]
        end
    end
    channel_order = _sort_channels(ch_type)
    ref = dataset["ref"]
    ref == "common" && (ref = "CAR")

    # TODO: import locations, events and other data
    # dataset["epoch"]
    # dataset["session"]
    # dataset["dipfit"]
    # dataset["group"]
    # dataset["datfile"]
    # 
    # dataset["specdata"]
    # dataset["run"]
    # dataset["epochdescription"]
    # locs data: dataset["urchanlocs"]
    # ICA: dataset["icawinv"]
    # ICA weights: dataset["icaweights"]
    # ICA: dataset["icaact"]
    # ICA: dataset["icachansind"]
    # ICA: dataset["icasplinefile"]
    # ICA: dataset["splinefile"]
    # ICA: dataset["icasphere"]
    # dataset["specicaact"]
    # dataset["stats"]
    # dataset["reject"]
    # ignore: dataset["xmin"], dataset["xmax"], dataset["filename"], dataset["filepath"], dataset["etc"], dataset["saved"], dataset["pnts"]

    # EEGLAB metadata
    patient = dataset["subject"]
    note = dataset["comments"]
    for idx in length(note):-1:1
        note[idx] == "" && deleteat!(note, idx)
    end
    note = string(note)
    note = replace(note, '"'=>"")

    history = split(dataset["history"], "\n")
    # remove first two entries, 1st is empty, second is EEGLAB version
    length(history) > 2 && (history = history[3:end])
    for idx in length(history):-1:1
        history[idx] == "" && deleteat!(history, idx)
    end
    history = string.(history)
    length(history) == 0 && (history = "")

    # LOCS
    # locs["urchan"]
    if "chanlocs" in keys(dataset) && length(dataset["chanlocs"]) > 0
        chanlocs = dataset["chanlocs"]
        # locs["ref"]
        x = zeros(ch_n)
        y = zeros(ch_n)
        z = zeros(ch_n)
        theta = zeros(ch_n)
        radius = zeros(ch_n)
        phi_sph = zeros(ch_n)
        radius_sph = zeros(ch_n)
        theta_sph = zeros(ch_n)
        for idx in 1:ch_n
            chanlocs["X"][:][idx] isa Float64 && (x[idx] = chanlocs["X"][:][idx])
            chanlocs["Y"][:][idx] isa Float64 && (y[idx] = chanlocs["Y"][:][idx])
            chanlocs["Z"][:][idx] isa Float64 && (z[idx] = chanlocs["Z"][:][idx])
            chanlocs["theta"][:][idx] isa Float64 && (theta[idx] = chanlocs["theta"][:][idx])
            chanlocs["radius"][:][idx] isa Float64 && (radius[idx] = chanlocs["radius"][:][idx])
            chanlocs["sph_phi"][:][idx] isa Float64 && (phi_sph[idx] = chanlocs["sph_phi"][:][idx])
            chanlocs["sph_radius"][:][idx] isa Float64 && (radius_sph[idx] = chanlocs["sph_radius"][:][idx])
            chanlocs["sph_theta"][:][idx] isa Float64 && (theta_sph[idx] = chanlocs["sph_theta"][:][idx])
        end
        radius_sph == zeros(ch_n) && (radius_sph = radius)
        locs = DataFrame(:labels=>clabels,
                         :loc_theta=>theta,
                         :loc_radius=>radius,
                         :loc_x=>x,
                         :loc_y=>y,
                         :loc_z=>z,
                         :loc_radius_sph=>radius_sph,
                         :loc_theta_sph=>theta_sph,
                         :loc_phi_sph=>phi_sph)
        for idx in nrow(locs):-1:1
            (chanlocs["X"][:][idx] isa Float64 && chanlocs["Y"][:][idx] isa Float64 && chanlocs["Z"][:][idx] isa Float64) || deleteat!(locs, idx)
        end
        nrow(locs) > 0 && _info("Locs for $(nrow(locs)) channel$(_pl(nrow(locs))) found")
        if nrow(locs) > 0
            dataset["chaninfo"]["nosedir"] == "+X" && locs_swapxy!(locs)
            locs_normalize!(locs)
        end
    else
        locs = _initialize_locs()
    end

    # MARKERS
    markers = DataFrame(:id=>String[],
                        :start=>Float64[],
                        :length=>Float64[],
                        :description=>String[],
                        :channel=>Int64[])

    if "event" in keys(dataset)
        events = dataset["event"]
        if length(keys(events)) > 0
            start = Float64.(events["latency"][:]) ./ sampling_rate
            # for idx in eachindex(events["position"][:])
            #     events["position"][idx] isa Matrix{Float64} && (events["position"][idx] = 0.0)
            # end
            # pos = Int.(events["position"][:])
            len = zeros(length(start))
            desc = String.(events["type"][:])
            id = repeat(["stim"], length(start))
            markers = DataFrame(:id=>id, :start=>start, :length=>len, :description=>desc, :channel=>zeros(Int64, length(start)))
        end
    end
    # dataset["eventdescription"]
    # dataset["urevent"]
    # dataset["condition"]

    if ep_n > 1
        if length(dataset["times"][:]) > 0
            epoch_time = dataset["times"][:]
        else
            epoch_time = round.((collect(0:1/sampling_rate:size(data, 2) / sampling_rate))[1:end-1], digits=3)
        end
        time_pts = round.(collect(0:1/sampling_rate:size(data, 2) * size(data, 3) / sampling_rate)[1:end-1], digits=3)
    else
        # if length(dataset["times"][:]) > 0
        #     time_pts = dataset["times"][:]
        # end
        time_pts = round.(collect(0:1/sampling_rate:size(data, 2) * size(data, 3) / sampling_rate)[1:end-1], digits=3)
        epoch_time = round.((collect(0:1/sampling_rate:size(data, 2) / sampling_rate))[1:end-1], digits=3)
    end

    if data_src isa String
        file_size_mb = round(filesize(data_src) / 1024^2, digits=2)
    else
        file_size_mb = round(filesize(file_name) / 1024^2, digits=2)
    end
    
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
                              recording="",
                              recording_date="",
                              recording_time="",
                              recording_notes=string(dataset["setname"]),
                              channel_type=ch_type,
                              reference=ref,
                              clabels=clabels,
                              transducers=repeat([""], ch_n),
                              units=repeat(["μV"], ch_n),
                              prefiltering=repeat([""], ch_n),
                              sampling_rate=sampling_rate,
                              gain=gain)
    e = _create_experiment(name="",
                           notes=note,
                           design="")

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = history

    obj = NeuroAnalyzer.NEURO(hdr, time_pts, epoch_time, data[channel_order, :, :], components, markers, locs, history)
    nrow(locs) == 0 && _initialize_locs!(obj)

    _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits=2)) s)")

    return obj

end
