export import_ft

"""
    import_ft(file_name; <keyword arguments>)

Load FieldTrip file (.mat) and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load
- `type::Symbol`: type of imported data
    - `:eeg`: EEG
    - `:meg`: MEG
    - `:nirs`: fNIRS
    - `:events`: events
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `obj::NeuroAnalyzer.NEURO` - for EEG, MEG, fNIRS data
- `markers::DataFrame` - for events
"""
function import_ft(file_name::String; type::Symbol, detect_type::Bool=true)::Union{NeuroAnalyzer.NEURO, DataFrame}

    _wip()

    _check_var(type, [:eeg, :meg, :nirs, :events], "type")
    file_type = "FT"
    data_type = String(type)

    @assert isfile(file_name) "File $file_name cannot be loaded."
    dataset = matread(file_name)

    @assert length(keys(dataset)) == 1 "Datasets containing > 1 object are not supported; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org"
    _info("Reading object: $(string.(keys(dataset))[1])")
    dataset = dataset[string.(keys(dataset))[1]]

    if data_type == "events"

        f = "type"
        @assert f in keys(dataset) "Dataset does not contain $f field."
        id = dataset[f][:]
        f = "duration"
        @assert f in keys(dataset) "Dataset does not contain $f field."
        duration = dataset[f][:]
        for idx in axes(duration, 1)
            size(duration[idx]) == (0, 0) && (duration[idx] = 0.0)
        end
        f = "offset"
        @assert f in keys(dataset) "Dataset does not contain $f field."
        offset = dataset[f][:]
        for idx in axes(offset, 1)
            size(offset[idx]) == (0, 0) && (offset[idx] = 0.0)
        end
        f = "sample"
        @assert f in keys(dataset) "Dataset does not contain $f field."
        start = dataset[f][:]
        f = "value"
        @assert f in keys(dataset) "Dataset does not contain $f field."
        value = dataset["value"][:]

        # find and remove empty records
        value_idx = findall(x -> length(x) != 0, value)
        id = id[value_idx]
        duration = duration[value_idx]
        offset = offset[value_idx]
        start = start[value_idx]
        value = value[value_idx]

        markers = DataFrame(:id=>string.(id),
                            :start=>Float64.(start),
                            :length=>Float64.(duration),
                            :value=>strip.(string.(value)),
                            :channel=>zeros(Int64, length(id)))
        _info("Imported: $(nrow(markers)) events; events start and length are in samples, use `markers_s2t()` to convert to seconds")

        return markers

    else

        @assert "cfg" in keys(dataset) "Dataset does not contain cfg field."
        cfg = dataset["cfg"]
        @assert "hdr" in keys(dataset) "Dataset does not contain hdr field."
        hdr = dataset["hdr"]
        @assert "nChans" in keys(hdr) "Dataset header does not contain nChans field."
        ch_n = Int64(hdr["nChans"])
        @assert "Fs" in keys(hdr) "Dataset header does not contain Fs field."
        sampling_rate = round(Int64, hdr["Fs"])

        # get channel labels, types and units
        if length(hdr["label"][:]) > 0 && length(vec(hdr["label"])) == ch_n
            clabels = string.(hdr["label"][:])
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
                units = string.(hdr["chanunit"][:])
                units = replace.(units, "uV"=>"μV")
            else
                units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]
            end
        else
            if "chantype" in keys(hdr)
                ch_type = string.(hdr["chantype"][:])
            else
                ch_type = repeat(["eeg"], ch_n)
            end
            if "chanunit" in keys(hdr)
                units = string.(hdr["chanunit"][:])
                units = replace.(units, "uV"=>"μV")
            else
                units = repeat(["μV"], ch_n)
            end
        end

        @assert "trial" in keys(dataset) "Dataset does not contain trial field."
        ep_n = "trial" in keys(dataset) ? size(dataset["trial"], 2) : 1
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

        # FieldTrip markers are in a separate object, must be imported and added manually
        markers = DataFrame(:id=>String[],
                            :start=>Float64[],
                            :length=>Float64[],
                            :value=>String[],
                            :channel=>Int64[])

        if data_type == "eeg"

            # TO DO: get referencing
            if "reref" in keys(dataset["cfg"]) && dataset["cfg"]["reref"] != "no"
                _info("Embedded referencing is not supported; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org")
            else
                ref = _detect_montage(clabels, ch_type, data_type)
            end

            r = _create_recording_eeg(data_type=data_type,
                                      file_name=file_name,
                                      file_size_mb=round(filesize(file_name) / 1024^2, digits=2),
                                      file_type=file_type,
                                      recording="RID" in keys(hdr["orig"]) ? string(hdr["orig"]["RID"]) : "",
                                      recording_date="",
                                      recording_time="",
                                      recording_notes="",
                                      channel_type=ch_type,
                                      channel_order=_sort_channels(ch_type),
                                      reference=ref,
                                      clabels=clabels,
                                      transducers="Transducer" in keys(hdr["orig"]) ? string.(hdr["orig"]["Transducer"]) : repeat([""], ch_n),
                                      units=units,
                                      prefiltering="PreFilt" in keys(hdr["orig"]) ? string.(hdr["orig"]["PreFilt"]) : repeat([""], ch_n),
                                      line_frequency=50,
                                      sampling_rate=sampling_rate,
                                      gain=ones(ch_n),
                                      bad_channels=zeros(Bool, size(data, 1), ep_n))

        elseif data_type == "meg"

            # Internal Active Shielding (IAS)
            ias_channels = .!(isnothing.(match.(r"IAS.*", clabels)))
            ias_labels = nothing
            ias_data = nothing
            if any(ias_channels)
                ias_labels = clabels[ias_channels]
                clabels = clabels[.!(ias_channels)]
                ias_data = data[ias_channels, :, :]
                data = data[.!(ias_channels), :, :]
            end

            r = _create_recording_meg(data_type=data_type,
                                      file_name=file_name,
                                      file_size_mb=round(filesize(file_name) / 1024^2, digits=2),
                                      file_type=file_type,
                                      recording="RID" in keys(hdr["orig"]) ? string(hdr["orig"]["RID"]) : "",
                                      recording_date="",
                                      recording_time="",
                                      recording_notes="",
                                      channel_type=ch_type,
                                      channel_order=_sort_channels(ch_type),
                                      reference=ref,
                                      clabels=clabels,
                                      transducers="Transducer" in keys(hdr["orig"]) ? string.(hdr["orig"]["Transducer"]) : repeat([""], ch_n),
                                      units=units,
                                      prefiltering="PreFilt" in keys(hdr["orig"]) ? string.(hdr["orig"]["PreFilt"]) : repeat([""], ch_n),
                                      line_frequency=50,
                                      sampling_rate=sampling_rate,
                                      gain=ones(ch_n),
                                      bad_channels=zeros(Bool, size(data, 1), ep_n))

        elseif data_type == "nirs"

            r = _create_recording_nirs(data_type=data_type,
                                       file_name=file_name,
                                       file_size_mb=round(filesize(file_name) / 1024^2, digits=2),
                                       file_type=file_type,
                                       recording="RID" in keys(hdr["orig"]) ? string(hdr["orig"]["RID"]) : "",
                                       recording_date="",
                                       recording_time="",
                                       recording_notes="",
                                       channel_type=ch_type,
                                       channel_order=_sort_channels(ch_type),
                                       reference=ref,
                                       clabels=clabels,
                                       transducers="Transducer" in keys(hdr["orig"]) ? string.(hdr["orig"]["Transducer"]) : repeat([""], ch_n),
                                       units=units,
                                       prefiltering="PreFilt" in keys(hdr["orig"]) ? string.(hdr["orig"]["PreFilt"]) : repeat([""], ch_n),
                                       line_frequency=50,
                                       sampling_rate=sampling_rate,
                                       gain=ones(ch_n),
                                       bad_channels=zeros(Bool, size(data, 1), ep_n))

        end

        s = _create_subject(id="",
                            first_name="",
                            middle_name="",
                            last_name="",
                            head_circumference=-1,
                            handedness="",
                            weight=-1,
                            height=-1)
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

end
