export import_nwb

"""
    import_nwb(file_name; <keyword arguments>)

Load EEG data from Neurodata Without Borders (NWB) file and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `obj::NeuroAnalyzer.NEURO`

# Source

1. https://www.biorxiv.org/content/10.1101/523035v1
"""
function import_nwb(file_name::String; detect_type::Bool=true)

    _wip()

    @assert isfile(file_name) "File $file_name cannot be loaded."
    @assert lowercase(splitext(file_name)[2]) == ".nwb" "This is not NWB file."

    file_type = "NWB"

    # process header
    recording = ""
    recording_reference = ""
    recording_notes = ""
    sampling_rate = nothing
    exp_name = ""
    exp_design = ""
    exp_notes = ""
    data_type = ""

    file_json = splitext(file_name)[1] * ".json"
    if isfile(file_json)
        f = open(file_json, "r")
        s = read(f, String)
        close(f)
        header = JSON.parse(s)
        k = keys(header)
        recording = "Manufacturer" in k ? header["Manufacturer"] : ""
        recording_reference = "iEEGReference" in k ? header["iEEGReference"] : ""
        recording_notes = "iEEGGround" in k ? ("iEEG ground: " * header["iEEGGround"]) : ""
        sampling_rate = "SamplingFrequency" in k ? round(Int64, header["SamplingFrequency"]) : nothing
        exp_name = "TaskName" in k ? header["TaskName"] : ""
        exp_design = "TaskDescription" in k ? header["TaskDescription"] : ""
        exp_notes = "Instructions" in k ? header["Instructions"] : ""
        "RecordingType" in k && @assert header["RecordingType"] == "continuous" "Non-continuous recordings are not supported; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org"

        # what if the files contains mixed recordings (e.g. EEG + SEEG)
        if "EEGChannelCount" in k && header["EEGChannelCount"] > 0
            data_type = "eeg"
        elseif "ECOGChannelCount" in k && header["ECOGChannelCount"] > 0
            data_type = "ecog"
        elseif "SEEGChannelCount" in k && header["SEEGChannelCount"] > 0
            data_type = "seeg"
        elseif "EOGChannelCount" in k && header["EOGChannelCount"] > 0
        elseif "EMGChannelCount" in k && header["EMGChannelCount"] > 0
        elseif "EMGChannelCount" in k && header["EMGChannelCount"] > 0
        elseif "TriggerChannelCount" in k && header["TriggerChannelCount"] > 0
        end

        # header["PowerLineFrequency"]
        # header["RecordingDuration"]
    end

    # load data
    dataset = JLD2.load(file_name)
    k = keys(dataset)

    recording_date = ""
    recording_time = ""
    if "file_create_date" in k
        t = typeof(dataset["file_create_date"]) == Vector{String} ? dataset["file_create_date"][1] : dataset["file_create_date"]
        # remove microseconds
        m = match(r".*(\.[0-9]*)\+.*", t)
        m !== nothing && (t = replace(t, m.captures[1]=>""))
        recording_date = ZonedDateTime(t, DateFormat("yyyy-mm-ddTHH:MM:SSzzzz"))
        recording_date = string(Dates.year(recording_date)) * "-" * lpad(string(Dates.month(recording_date)), 2, '0') * "-" * lpad(string(Dates.day(recording_date)), 2, '0')
        # convert to UTC time
        recording_time = astimezone(ZonedDateTime(t, DateFormat("yyyy-mm-ddTHH:MM:SSzzzz")), tz"UTC")
        recording_time = lpad(string(Dates.hour(recording_time)), 2, '0') * ":" * lpad(string(Dates.minute(recording_time)), 2, '0') * ":" * lpad(string(Dates.second(recording_time)), 2, '0')
    end

    subj_id = "identifier" in k ? dataset["identifier"] : ""
    exp_name = "session_description" in k ? dataset["session_description"] : ""

    # if "session_start_time" in k
    #     t = dataset["session_start_time"][1]
    #     # remove microseconds
    #     m = match(r".*(\.[0-9]*)\+.*", t)
    #     m !== nothing && (t = replace(t, m.captures[1]=>""))
    #     Dates.Second(ZonedDateTime(t, DateFormat("yyyy-mm-ddTHH:MM:SSzzzz"))).value
    # end

    # check if there is an offset of time points
    t_start = 0
    if "timestamps_reference_time" in k
        t = typeof(dataset["timestamps_reference_time"]) == Vector{String} ? dataset["timestamps_reference_time"][1] : dataset["timestamps_reference_time"]
        # remove microseconds
        m = match(r".*(\.[0-9]*)\+.*", t)
        m !== nothing && (t = replace(t, m.captures[1]=>""))
        t_start = Dates.Second(ZonedDateTime(t, DateFormat("yyyy-mm-ddTHH:MM:SSzzzz"))).value
    end

    eeg_regexp = r"acquisition.*EEG.*data"
    m = nothing
    for idx in k
        if match(eeg_regexp, idx) !== nothing
            m = match(eeg_regexp, idx)
            break
        end
    end
    if m !== nothing
        data = dataset[m.match]
    else
        @error "File $file_name does not contain EEG data."
    end

    timepts_regexp = r"acquisition.*EEG.*timestamps"
    m = nothing
    for idx in k
        if match(timepts_regexp, idx) !== nothing
            m = match(timepts_regexp, idx)
            break
        end
    end
    if m !== nothing
        timepts = dataset[m.match]
    else
        @error "File $file_name does not contain time points data."
    end
    sampling_rate === nothing && (sampling_rate = Int64(1 / (timepts[2] - timepts[1])))

    time_pts = round.(collect(0:1/sampling_rate:size(data, 2) * size(data, 3) / sampling_rate)[1:end-1], digits=3) .+ t_start
    ep_time = round.((collect(0:1/sampling_rate:size(data, 2) / sampling_rate))[1:end-1], digits=3) .+ t_start

    # events
    "acquisition/Stimulus/data" in k && (stim = dataset["acquisition/Stimulus/data"])
    # dataset["acquisition/Stimulus/timestamps"]
    # check if events file is available
    m = match(r".*(_[is]eeg).*", file_name)
    m !== nothing && (events_file = replace(file_name, m.captures[1]=>"_events"))
    m = match(r".*(_eeg).*", file_name)
    m !== nothing && (events_file = replace(file_name, m.captures[1]=>"_events"))
    m = match(r".*(_meg).*", file_name)
    m !== nothing && (events_file = replace(file_name, m.captures[1]=>"_events"))
    m = match(r".*(_ecog).*", file_name)
    m !== nothing && (events_file = replace(file_name, m.captures[1]=>"_events"))
    occursin(".nwb", events_file) && (events_file = replace(events_file, ".nwb"=>".tsv"))
    if isfile(events_file)
        events = CSV.read(events_file, DataFrame)
        _info("$(nrow(events)) markers found")
        event_id = events[!, :trial_type]
        event_description = events[!, :value]
        event_start_sample = events[!, :sample] .+ 1
        event_start = zeros(length(event_start_sample))
        [event_start[idx] = time_pts[event_start_sample[idx]] for idx in eachindex(event_start_sample)]
        event_length = round.(events[!, :duration], digits=3)
        event_channel = zeros(Int64, nrow(events))
        markers = DataFrame(:id=>event_id,
                            :start=>event_start,
                            :length=>event_length,
                            :description=>event_description,
                            :channel=>event_channel)
    else
        markers = DataFrame(:id=>String[],
                            :start=>Float64[],
                            :length=>Float64[],
                            :description=>String[],
                            :channel=>Int64[])
    end

    @assert isfile(file_json) "$file_json not found."
    f = open(file_json, "r")
    s = read(f, String)
    close(f)

    # dataset["specifications/core/2.4.0/namespace"]
    # dataset["specifications/core/2.4.0/nwb.base"]
    # dataset["specifications/core/2.4.0/nwb.behavior"]
    # dataset["specifications/core/2.4.0/nwb.device"]
    # dataset["specifications/core/2.4.0/nwb.ecephys"]
    # dataset["specifications/core/2.4.0/nwb.epoch"]
    # dataset["specifications/core/2.4.0/nwb.file"]
    # dataset["specifications/core/2.4.0/nwb.icephys"]
    # dataset["specifications/core/2.4.0/nwb.image"]
    # dataset["specifications/core/2.4.0/nwb.misc"]
    # dataset["specifications/core/2.4.0/nwb.ogen"]
    # dataset["specifications/core/2.4.0/nwb.ophys"]
    # dataset["specifications/core/2.4.0/nwb.retinotopy"]

    # dataset["specifications/hdmf-common/1.5.1/base"]
    # dataset["specifications/hdmf-common/1.5.1/namespace"]
    # dataset["specifications/hdmf-common/1.5.1/sparse"]
    # dataset["specifications/hdmf-common/1.5.1/table"]

    # dataset["specifications/hdmf-experimental/0.2.0/experimental"]
    # dataset["specifications/hdmf-experimental/0.2.0/namespace"]
    # dataset["specifications/hdmf-experimental/0.2.0/resources"]


    file_size_mb = round(filesize(file_name) / 1024^2, digits=2)

    s = _create_subject(id=subj_id,
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
                              recording=recording,
                              recording_date=recording_date,
                              recording_time=recording_time,
                              recording_notes=recording_notes,
                              channel_type=ch_type,
                              channel_order=_sort_channels(ch_type),
                              reference=recording_reference,
                              clabels=clabels,
                              transducers=transducers,
                              units=units,
                              prefiltering=prefiltering,
                              sampling_rate=sampling_rate,
                              gain=gain,
                              bad_channels=zeros(Bool, size(data, 1), 1))
    e = _create_experiment(name=exp_name, notes=exp_notes, design=exp_design)

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = String[]

    locs = _initialize_locs()
    obj = NeuroAnalyzer.NEURO(hdr, time_pts, ep_time, data, components, markers, locs, history)
    _initialize_locs!(obj)

    _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits=2)) s)")

    return obj

end