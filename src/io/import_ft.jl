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
- `detect_type::Bool=false`: detect channel type based on its label

# Returns

- `obj::NeuroAnalyzer.NEURO` - for EEG, MEG, fNIRS data
- `markers::DataFrame` - for events
"""
function import_ft(file_name::String; type::Symbol, detect_type::Bool=false)::Union{NeuroAnalyzer.NEURO, DataFrame}

    _wip()

    _check_var(type, [:eeg, :meg, :nirs, :events], "type")
    file_type = "FT"
    data_type = String(type)
    locs = nothing

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
        @assert "fsample" in keys(dataset) "Dataset does not contain fsample field."
        sampling_rate = round(Int64, dataset["fsample"])
        @assert "nChans" in keys(hdr) "Dataset header does not contain nChans field."
        ch_n = Int64(hdr["nChans"])

        # get channel labels, types and units
        if length(hdr["label"][:]) > 0 && length(vec(hdr["label"])) == ch_n
            clabels = string.(strip.(string.(hdr["label"])))[:]
        else
            clabels = String[]
            for idx in 1:ch_n
                push!(clabels, "ch_$idx")
            end
        end
        clabels = _clean_labels(clabels)

        if detect_type
            ch_type = _set_channel_types(clabels)
            if "chanunit" in keys(hdr)
                units = strip.(string.(hdr["chanunit"][:]))
                units = replace.(units, "uV"=>"μV")
            else
                units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]
            end
        else
            if "chantype" in keys(hdr)
                ch_type = string.(hdr["chantype"][:])
            else
                ch_type = repeat([data_type], ch_n)
            end
            if "chanunit" in keys(hdr)
                units = string.(hdr["chanunit"][:])
                units = replace.(units, "uV"=>"μV")
            else
                units = repeat(["μV"], ch_n)
            end
        end
        ch_type = replace.(ch_type, "nirs"=>"nirs_od")
        ch_type = replace.(ch_type, "stimulus"=>"mrk")

        @assert "trial" in keys(dataset) "Dataset does not contain trial field."
        ep_n = "trial" in keys(dataset) ? size(dataset["trial"], 2) : 1
        ep_len = size(dataset["trial"][1], 2)
        data = zeros(ch_n, ep_len, ep_n)
        for idx in 1:ep_n
            data[:, :, idx] = dataset["trial"][idx]
        end

        if "time" in keys(dataset)
            if ep_n == 1
                epoch_time = dataset["time"][1][:]
                sign(epoch_time[1]) == 1.0 && (epoch_time .-= (epoch_time[1]))
                time_pts = dataset["time"][1][:]
                if sign(time_pts[1]) == 1.0
                    time_pts .-= (time_pts[1])
                else
                    time_pts .+= abs(time_pts[1])
                end
                epoch_time = round.(epoch_time, digits=3)
                time_pts = round.(time_pts, digits=3)
            else
                epoch_time = round.(dataset["time"][1][:], digits=3)
                time_pts = round.(collect(0:1/sampling_rate:size(data, 2) * size(data, 3) / sampling_rate)[1:end-1], digits=3)
            end
        else
            epoch_time = round.((collect(0:1/sampling_rate:size(data, 2) / sampling_rate))[1:end-1], digits=3)
            time_pts = round.(collect(0:1/sampling_rate:size(data, 2) * size(data, 3) / sampling_rate)[1:end-1], digits=3)
        end

        _info("FieldTrip markers are stored separately and must be imported using `import_ft(file_name, type=:events)` and added manually using `add_markers()`")
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

            # convert data to standard units
            @inbounds for ch_idx in 1:ch_n
                if units[ch_idx] == "V" && ch_type[ch_idx] in ["eeg", "emg", "eog", "ref"]
                    @views data[ch_idx, :, 1] .*= 10^6
                    units[ch_idx] = "μV"
                end
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
                                      transducers="Transducer" in keys(hdr["orig"]) ? string.(strip.(hdr["orig"]["Transducer"])) : repeat([""], ch_n),
                                      units=units,
                                      prefiltering="PreFilt" in keys(hdr["orig"]) ? string.(strip.(hdr["orig"]["PreFilt"])) : repeat([""], ch_n),
                                      line_frequency=50,
                                      sampling_rate=sampling_rate,
                                      gain=ones(ch_n),
                                      bad_channels=zeros(Bool, size(data, 1), ep_n))

        elseif data_type == "meg"

            @assert "grad" in keys(dataset) "Dataset does not contain grad field."
            @assert "chantype" in keys(dataset["grad"]) "Dataset does not contain chantype field."
            meg_channels = occursin.(r"meg", lowercase.(clabels))[:]
            coil_type = repeat([""], ch_n)
            mag_idx = occursin.(r".*mag.*", lowercase.(ch_type))
            coil_type[mag_idx] .= "mag"
            grad_idx = occursin.(r".*grad.*", lowercase.(ch_type))
            coil_type[grad_idx] .= "grad"
            pgrad_idx = occursin.(r".*planar.*", lowercase.(ch_type))
            coil_type[pgrad_idx] .= "pgrad"
            agrad_idx = occursin.(r".*axial.*", lowercase.(ch_type)) .|| occursin.(r".*ctf.*", lowercase.(ch_type))
            coil_type[agrad_idx] .= "agrad"
            grad_idx .+= pgrad_idx
            grad_idx .+= agrad_idx
            magnetometers = collect(1:length(ch_type))[mag_idx]
            gradiometers = collect(1:length(ch_type))[grad_idx]
            ch_type[magnetometers] .= "mag"
            ch_type[gradiometers] .= "grad"

            # Internal Active Shielding (IAS)
            ch_type[occursin.("ias", lowercase.(clabels))] .= "other"
            ch_type[occursin.("sti", lowercase.(clabels))] .= "mrk"
            ch_type[occursin.("sys", lowercase.(clabels))] .= "other"

            # convert data to standard units
            @inbounds for ch_idx in 1:ch_n
                if units[ch_idx] == "T"
                    @views data[ch_idx, :, 1] .*= 10^15
                    units[ch_idx] = "fT"
                elseif units[ch_idx] == "T/m"
                    @views data[ch_idx, :, 1] .*= (10^15 / 100)
                    units[ch_idx] = "fT/cm"
                elseif units[ch_idx] == "T/cm"
                    @views data[ch_idx, :, 1] .*= 10^15
                    units[ch_idx] = "fT/cm"
                elseif units[ch_idx] == "V"
                    @views data[ch_idx, :, 1] .*= 10^6
                    units[ch_idx] = "μV"
                end
            end

            lp = "lowpass" in keys(hdr["orig"]) ? string(round(hdr["orig"]["lowpass"][1], digits=1)) : "?"
            hp = "highpass" in keys(hdr["orig"]) ? string(round(hdr["orig"]["highpass"][1], digits=1)) : "?"
            r = _create_recording_meg(data_type=data_type,
                                      file_name=file_name,
                                      file_size_mb=round(filesize(file_name) / 1024^2, digits=2),
                                      file_type="FT",
                                      recording="dataformat" in keys(dataset["cfg"]) ? string(dataset["cfg"]["dataformat"]) : "",
                                      recording_date="",
                                      recording_time="",
                                      recording_notes="",
                                      channel_type=ch_type,
                                      channel_order=_sort_channels(ch_type),
                                      reference="",
                                      clabels=clabels,
                                      units=units,
                                      prefiltering=repeat(["LP: $lp Hz; HP: $hp Hz"], ch_n),
                                      line_frequency=50,
                                      sampling_rate=sampling_rate,
                                      magnetometers=magnetometers,
                                      gradiometers=gradiometers,
                                      coil_type=coil_type,
                                      bad_channels=zeros(Bool, size(data, 1), ep_n))

        elseif data_type == "nirs"


            clabels = replace.(clabels, ".0"=>"")
            clabels = replace.(clabels, " ["=>" ")
            clabels = replace.(clabels, "nm]"=>"")

            @assert "opto" in keys(dataset) "Dataset does not contain opto field."
            opto = dataset["opto"]
            wavelengths = opto["wavelength"][:]
            wavelength_index = Int64[]
            for idx1 in eachindex(ch_type)
                if ch_type[idx1] == "nirs"
                    for idx2 in eachindex(wavelengths)
                        occursin(string(round(Int64, wavelengths[idx2])), clabels[idx1]) && push!(wavelength_index, idx2)
                    end
                end
            end

            nirs_channels_idx = ch_type .== "nirs_od"
            nirs_channels = count(nirs_channels_idx)
            # source and detector names
            source_index = opto["optotype"] .== "transmitter"
            detector_index = opto["optotype"] .== "receiver"
            src_labels = string.(opto["optolabel"][:][source_index])
            det_labels = string.(opto["optolabel"][:][detector_index])
            opt_labels = string.(opto["optolabel"][:])

            # collect channel pairs
            opt_pairs = zeros(Int64, nirs_channels, 2)
            pairs = opto["label"][:]
            for idx in eachindex(pairs)
                p = match(r"(S\d+)-(D\d+)(.*)", pairs[idx])
                if !isnothing(p)
                    opt_pairs[idx, 1] = findfirst(isequal(p[1]), src_labels)
                    opt_pairs[idx, 2] = findfirst(isequal(p[2]), det_labels)
                end
            end

            # optode locations
            x = opto["optopos"][:, 1]
            y = opto["optopos"][:, 2]
            z = opto["optopos"][:, 3]
            global locs = DataFrame(:label=>opt_labels, :loc_radius=>zeros(length(opt_labels)), :loc_theta=>zeros(length(opt_labels)), :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>zeros(length(opt_labels)), :loc_theta_sph=>zeros(length(opt_labels)), :loc_phi_sph=>zeros(length(opt_labels)))
            locs_normalize!(locs)
            locs_cart2sph!(locs)
            locs_cart2pol!(locs)

            r = _create_recording_nirs(data_type=data_type,
                                       file_name=file_name,
                                       file_size_mb=round(filesize(file_name) / 1024^2, digits=2),
                                       file_type=file_type,
                                       recording="",
                                       recording_date="",
                                       recording_time="",
                                       recording_notes="",
                                       wavelengths=wavelengths,
                                       wavelength_index=wavelength_index,
                                       optode_pairs=opt_pairs,
                                       channel_type=ch_type,
                                       channel_order=_sort_channels(ch_type),
                                       clabels=clabels,
                                       units=units,
                                       src_labels=src_labels,
                                       det_labels=det_labels,
                                       opt_labels=opt_labels,
                                       sampling_rate=sampling_rate,
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

        if data_type == "meg" || data_type == "eeg"
            locs = _initialize_locs()
        end
        obj = NeuroAnalyzer.NEURO(hdr, time_pts, epoch_time, data, components, markers, locs, history)
        if data_type == "meg"
            _initialize_locs!(obj)
            l = import_locs_csv(joinpath(NeuroAnalyzer.res_path, "meg_306flattened.csv"))
            add_locs!(obj, locs=l)
        elseif data_type == "eeg"
            _initialize_locs!(obj)
        end

        _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits=2)) s)")

        return obj

    end

end
