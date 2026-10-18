export import_nwb

"""
    import_nwb(file_name; <keyword arguments>)

Load EEG data from Neurodata Without Borders (NWB) file and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `NeuroAnalyzer.NEURO`

# References

1. https://www.biorxiv.org/content/10.1101/523035v1
"""
function import_nwb(file_name::String; detect_type::Bool = true)::NeuroAnalyzer.NEURO
    _wip()

    # validate
    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    detect_format(file_name) === :nwb ||
        throw(ArgumentError("$file_name is not an NWB file."))

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
        recording_notes =
            "iEEGGround" in k ? ("iEEG ground: " * header["iEEGGround"]) : ""
        sampling_rate =
            "SamplingFrequency" in k ? round(Int64, header["SamplingFrequency"]) : nothing
        exp_name = "TaskName" in k ? header["TaskName"] : ""
        exp_design = "TaskDescription" in k ? header["TaskDescription"] : ""
        exp_notes = "Instructions" in k ? header["Instructions"] : ""
        # only continuous recordings are supported (a missing field is not an error)
        (!("RecordingType" in k) || header["RecordingType"] == "continuous") ||
            throw(
                ArgumentError(
                    "Non-continuous recordings are not supported yet; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org",
                ),
            )

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
        t = if typeof(dataset["file_create_date"]) == Vector{String}
            dataset["file_create_date"][1]
        else
            dataset["file_create_date"]
        end
        # remove microseconds
        m = match(r".*(\.[0-9]*)\+.*", t)
        m !== nothing && (t = replace(t, m.captures[1] => ""))
        recording_date = ZonedDateTime(t, DateFormat("yyyy-mm-ddTHH:MM:SSzzzz"))
        recording_date =
            string(Dates.year(recording_date)) *
            "-" *
            lpad(string(Dates.month(recording_date)), 2, '0') *
            "-" *
            lpad(string(Dates.day(recording_date)), 2, '0')
        # convert to UTC time
        recording_time = astimezone(
            ZonedDateTime(t, DateFormat("yyyy-mm-ddTHH:MM:SSzzzz")), tz"UTC",
        )
        recording_time =
            lpad(string(Dates.hour(recording_time)), 2, '0') *
            ":" *
            lpad(string(Dates.minute(recording_time)), 2, '0') *
            ":" *
            lpad(string(Dates.second(recording_time)), 2, '0')
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

    # note: timestamps_reference_time is an absolute date; the imported time axis always starts at 0

    eeg_regexp = r"acquisition.*EEG.*data"
    m = nothing
    for idx in k
        if match(eeg_regexp, idx) !== nothing
            m = match(eeg_regexp, idx)
            break
        end
    end
    m === nothing && throw(ArgumentError("File $file_name does not contain EEG data."))
    data = Float64.(dataset[m.match])
    ndims(data) == 2 && (data = reshape(data, size(data, 1), size(data, 2), 1))

    timepts_regexp = r"acquisition.*EEG.*timestamps"
    m = nothing
    for idx in k
        if match(timepts_regexp, idx) !== nothing
            m = match(timepts_regexp, idx)
            break
        end
    end
    if sampling_rate === nothing
        m === nothing &&
            throw(ArgumentError("File $file_name does not contain time points data."))
        timepts = dataset[m.match]
        sampling_rate = round(Int64, (length(timepts) - 1) / (timepts[end] - timepts[1]))
    end

    time_pts = _time_axis(size(data, 2) * size(data, 3), sampling_rate)
    epoch_time = _time_axis(size(data, 2), sampling_rate)

    # channel metadata (not read from the NWB electrodes table yet)
    data_type == "" && (data_type = "eeg")
    ch_n = size(data, 1)
    clabels = ["ch_$idx" for idx = 1:ch_n]
    ch_type = repeat([data_type], ch_n)
    units = [_ch_units(t) for t in ch_type]
    transducers = repeat([""], ch_n)
    prefiltering = repeat([""], ch_n)
    gain = ones(ch_n)

    # events
    "acquisition/Stimulus/data" in k && (stim = dataset["acquisition/Stimulus/data"])
    # dataset["acquisition/Stimulus/timestamps"]
    # check if events file is available
    # BIDS sidecar: sub-01_task-x_ieeg.nwb → sub-01_task-x_events.tsv
    base = splitext(file_name)[1]
    m = match(r"_(ieeg|seeg|ecog|eeg|meg)$", base)
    events_file = m === nothing ? "" : base[1:(end - length(m.match))] * "_events.tsv"
    if isfile(events_file)
        events = CSV.read(events_file, DataFrame; stringtype = String, missingstring = "n/a")
        _info("$(DataFrames.nrow(events)) markers found")
        event_id = string.(coalesce.(events[!, :trial_type], ""))
        event_description = string.(coalesce.(events[!, :value], ""))
        # BIDS: onset and duration in seconds; sample is 0-based
        event_start = if "onset" in names(events)
            Float64.(events[!, :onset])
        else
            Float64.(events[!, :sample]) ./ sampling_rate
        end
        event_length = Float64.(coalesce.(events[!, :duration], 0.0))
        event_channel = zeros(Int64, DataFrames.nrow(events))
        markers = DataFrame(
            :id => event_id,
            :start => event_start,
            :length => event_length,
            :value => event_description,
            :channel => event_channel,
        )
    else
        markers = DataFrame(
            :id => String[],
            :start => Float64[],
            :length => Float64[],
            :value => String[],
            :channel => Int64[],
        )
    end

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

    file_size_mb = round(filesize(file_name) / 1024^2; digits = 2)

    s = _create_subject(;
        id = subj_id,
        first_name = "",
        middle_name = "",
        last_name = "",
        head_circumference = -1,
        handedness = "",
        weight = -1,
        height = -1,
    )
    r = _create_recording_eeg(;
        data_type = data_type,
        file_name = file_name,
        file_size_mb = file_size_mb,
        file_type = file_type,
        recording = recording,
        recording_date = recording_date,
        recording_time = replace(recording_time, '.' => ':'),
        recording_notes = recording_notes,
        channel_type = ch_type,
        channel_order = _sort_channels(ch_type),
        reference = recording_reference,
        clabels = clabels,
        transducers = transducers,
        units = units,
        prefiltering = prefiltering,
        line_frequency = 50,
        sampling_rate = sampling_rate,
        gain = gain,
        bad_channels = zeros(Bool, ch_n),
    )
    e = _create_experiment(; name = exp_name, notes = exp_notes, design = exp_design)

    hdr = _create_header(; subject = s, recording = r, experiment = e)

    history = String[]

    locs = _initialize_locs()
    obj = NeuroAnalyzer.NEURO(hdr, history, markers, locs, time_pts, epoch_time, data)
    _initialize_locs!(obj)

    _info(
        "Imported: " *
        uppercase(obj.header.recording[:data_type]) *
        " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits = 2)) s)",
    )

    return obj
end
