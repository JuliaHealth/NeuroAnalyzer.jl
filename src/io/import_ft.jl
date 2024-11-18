export import_ft

"""
    import_ft(file_name; <keyword arguments>)

Load FieldTrip file (.mat) and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load
- `type::Symbol=:eeg`: type of imported data
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function import_ft(file_name::String; type::Symbol=:eeg, detect_type::Bool=true)::NeuroAnalyzer.NEURO

    _wip()

    _check_var(type, [:eeg, :meg, :nirs], "type")
    @assert type === :eeg "Currently only EEG datatype is supported."
    data_type = String(type)

    @assert isfile(file_name) "File $file_name cannot be loaded."

    file_type = "FT"

    dataset = matread(file_name)
    @assert length(keys(dataset)) == 1 "Datasets containing > 1 object are not supported; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org"
    _info("Reading object: $(String.(keys(dataset))[1])")
    dataset = dataset[String.(keys(dataset))[1]]

    @assert "cfg" in keys(dataset) "Dataset does not contain cfg field."
    cfg = dataset["cfg"]
    @assert "hdr" in keys(dataset) "Dataset does not contain header field."
    hdr = dataset["hdr"]
    @assert "nChans" in keys(hdr) "Dataset header does not contain nChans field."
    ch_n = Int64(hdr["nChans"])
    @assert "Fs" in keys(hdr) "Dataset header does not contain Fs field."
    sampling_rate = round(Int64, hdr["Fs"])

    # get channel labels, types and units
    if length(hdr["label"][:]) > 0 && length(vec(hdr["label"])) == ch_n
        clabels = String.(hdr["label"][:])
    else
        clabels = String[]
        for idx in 1:ch_n
            push!(clabels, "ch_$idx")
        end
    end
    clabels = _clean_labels(clabels)
    if detect_type
        ch_type = _set_channel_types(clabels, "eeg")
        if "chanunit" in keys(hdr)
            units = String.(hdr["chanunit"][:])
            units = replace.(units, "uV"=>"μV")
        else
            units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]
        end
    else
        if "chantype" in keys(hdr)
            ch_type = String.(hdr["chantype"][:])
        else
            ch_type = repeat(["eeg"], ch_n)
        end
        if "chanunit" in keys(hdr)
            units = String.(hdr["chanunit"][:])
            units = replace.(units, "uV"=>"μV")
        else
            units = repeat(["μV"], ch_n)
        end
    end

    # TO DO: get referencing
    if "reref" in keys(dataset["cfg"]) && dataset["cfg"]["reref"] != "no"
        _info("Embedded referencing is not supported; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org")
    else
        ref = _detect_montage(clabels, ch_type, data_type)
    end

    @assert "trial" in keys(dataset) "Dataset does not contain trial field."
    if "trial" in keys(dataset)
        ep_n = size(dataset["trial"], 2)
    else
        ep_n = 1
    end
    ep_len = size(dataset["trial"][1], 2)
    data = zeros(ch_n, ep_len, ep_n)
    for idx in 1:ep_n
        data[:, :, idx] = dataset["trial"][idx]
    end

    if "time" in keys(dataset)
        if ep_n == 1
            epoch_time = round.(dataset["time"][1][:], digits=3)
            time_pts = dataset["time"][1][:]
            time_pts .+= abs(time_pts[1])
            time_pts = round.(time_pts, digits=3)
        else
            epoch_time = round.(dataset["time"][1][:], digits=3)
            time_pts = round.(collect(0:1/sampling_rate:size(data, 2) * size(data, 3) / sampling_rate)[1:end-1], digits=3)
        end
    else
        epoch_time = round.((collect(0:1/sampling_rate:size(data, 2) / sampling_rate))[1:end-1], digits=3)
        time_pts = round.(collect(0:1/sampling_rate:size(data, 2) * size(data, 3) / sampling_rate)[1:end-1], digits=3)
    end

    # TODO: import events and other data

    file_size_mb = round(filesize(file_name) / 1024^2, digits=2)

    # MARKERS
    markers = DataFrame(:id=>String[],
                        :start=>Float64[],
                        :length=>Float64[],
                        :description=>String[],
                        :channel=>Int64[])

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
                              recording="RID" in keys(hdr["orig"]) ? string(hdr["orig"]["RID"]) : "",
                              recording_date="",
                              recording_time="",
                              recording_notes="",
                              channel_type=ch_type,
                              channel_order=_sort_channels(ch_type),
                              reference=ref,
                              clabels=clabels,
                              transducers="Transducer" in keys(hdr["orig"]) ? String.(hdr["orig"]["Transducer"]) : repeat([""], ch_n),
                              units=units,
                              prefiltering="PreFilt" in keys(hdr["orig"]) ? String.(hdr["orig"]["PreFilt"]) : repeat([""], ch_n),
                              line_frequency=50,
                              sampling_rate=sampling_rate,
                              gain=ones(ch_n),
                              bad_channels=zeros(Bool, size(data, 1), ep_n))
    e = _create_experiment(name="",
                           notes="",
                           design="")

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = [""]

    locs = _initialize_locs()
    obj = NeuroAnalyzer.NEURO(hdr, time_pts, epoch_time, data, components, markers, locs, history)
    _initialize_locs!(obj)

    _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits=2)) s)")

    return obj

end
