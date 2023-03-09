export import_set

"""
    import_set(file_name; detect_type)

Load SET file (exported from EEGLAB) and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function import_set(file_name::String; detect_type::Bool=true)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    file_type = "SET"

    dataset = matread(file_name)
    time_pts = dataset["times"][:]
    data = dataset["data"]

    # there are no epochs if signal is matrix, not array
    ndims(data) == 2 && (data = reshape(data, size(data, 1), size(data, 2), 1))

    ch_n = size(data, 1)
    
    # get channel labels
    if length(dataset["chanlocs"]["labels"][:]) == ch_n
        labels = String.(dataset["chanlocs"]["labels"][:])
    else
        labels = repeat([""], ch_n)
    end

    labels = _clean_labels(labels)
    if detect_type == true
        channel_type = _set_channel_types(labels)
    else
        channel_type = repeat(["???"], ch_n)
    end
    channel_order = _sort_channels(copy(channel_type))

    # TODO: import locations, events and other data
    # keys(dataset) = ["event", "icawinv", "chaninfo", "epoch", "stats", "chanlocs", "reject", "icaact", "icaweights", "ref", "eventdescription", "urchanlocs", "urevent", "nbchan", "icachansind", "specicaact", "icasplinefile", "splinefile", "condition", "dipfit", "group", "icasphere", "session", "datfile", "trials", "epochdescription", "setname", "specdata", "run"]
    # epochs data: dataset["epoch"]
    # events data: dataset["event"]
    # channel data: dataset["chaninfo"]
    # locs data: dataset["chanlocs"]
    # ICA weights: dataset["icaweights"]
    # ICA weights: dataset["icaweights"]
    # ignore: xmin, xmax, filename, filepath, etc, setname, saved, pnts

    # EEGLAB metadata
    patient = dataset["subject"]
    note = dataset["comments"]
    history = split(dataset["history"], "\n")
    # remove first two entries, 1st is empty, second is EEGLAB version
    length(history) > 2 && (history = history[3:end])

    sampling_rate = round(Int64, dataset["srate"])

    has_markers = false
    markers = DataFrame(:id => String[], :start => Int64[], :length => Int64[], :description => String[], :channel => Int64[])
    gain = ones(ch_n)
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
                      :patient => patient,
                      :recording => "",
                      :recording_date => "",
                      :recording_time => "",
                      :channel_n => ch_n,
                      :channel_type => channel_type[channel_order],
                      :reference => "",
                      :channel_locations => false,
                      :history => history,
                      :components => Symbol[],
                      :duration_samples => duration_samples,
                      :duration_seconds => duration_seconds,
                      :epoch_n => 1,
                      :epoch_duration_samples => duration_samples,
                      :epoch_duration_seconds => duration_seconds,
                      :labels => labels[channel_order],
                      :transducers => repeat([""], ch_n),
                      :units => repeat([""], ch_n),
                      :prefiltering => repeat([""], ch_n),
                      :sampling_rate => sampling_rate,
                      :gain => gain[channel_order],
                      :note => note,
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
