export import_ft

"""
    import_ft(file_name; <keyword arguments>)

Load a FieldTrip MATLAB (`.mat`) file and return a `NeuroAnalyzer.NEURO` object or a markers `DataFrame`.

FieldTrip stores EEG, MEG, fNIRS data and event tables in separate `.mat` files. The `type` argument selects which data type to import. Markers are always stored separately in FieldTrip and must be imported with `type = :events`, then added with `add_markers()`.

# Arguments

- `file_name::String`: path to the `.mat` file
- `type::Symbol`: data type to import; one of `:eeg`, `:meg`, `:nirs`, `:events`
- `detect_type::Bool=false`: infer channel type from channel label

# Returns

- `NeuroAnalyzer.NEURO` for EEG, MEG, and fNIRS data
- `DataFrame` for event tables

# Throws

- `ArgumentError` if the file does not exist, the dataset is malformed, or required fields are missing
"""
function import_ft(
    file_name::String;
    type::Symbol,
    detect_type::Bool = false
)::Union{NeuroAnalyzer.NEURO, DataFrame}

    _check_var(type, [:eeg, :meg, :nirs, :events], "type")
    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))

    file_type = "FT"
    data_type = String(type)

    dataset = matread(file_name)

    length(keys(dataset)) == 1 ||
        throw(ArgumentError(
            "Files with > 1 dataset are not supported; " *
            "please send this file to adam.wysokinski@neuroanalyzer.org"))
    _info("Reading object: $(string.(keys(dataset))[1])")
    dataset = dataset[string.(keys(dataset))[1]]

    # ------------------------------------------------------------------ #
    # events branch                                                       #
    # ------------------------------------------------------------------ #
    if data_type == "events"

        for f in ("type", "duration", "offset", "sample", "value")
            f in keys(dataset) ||
                throw(ArgumentError("Dataset does not contain '$f' field."))
        end

        id       = dataset["type"][:]
        duration = dataset["duration"][:]
        offset   = dataset["offset"][:]
        start    = dataset["sample"][:]
        value    = dataset["value"][:]

        # replace empty MATLAB cells (0×0 matrices) with 0.0
        for idx in eachindex(duration)
            size(duration[idx]) == (0, 0) && (duration[idx] = 0.0)
            size(offset[idx])   == (0, 0) && (offset[idx]   = 0.0)
        end

        # remove records with empty value fields
        keep     = findall(x -> length(x) != 0, value)
        id = id[keep]
        duration = duration[keep]
        offset = offset[keep]
        start = start[keep]
        value = value[keep]

        markers = DataFrame(
            :id      => string.(id),
            :start   => Float64.(start),
            :length  => Float64.(duration),
            :value   => strip.(string.(value)),
            :channel => zeros(Int64, length(id)))
        _info("Imported: $(DataFrames.nrow(markers)) events; " *
              "start and length are in samples — use `markers_s2t()` to convert.")
        return markers

    end

    # ------------------------------------------------------------------ #
    # signal branches (EEG / MEG / NIRS)                                 #
    # ------------------------------------------------------------------ #
    "cfg" in keys(dataset) ||
        throw(ArgumentError("Dataset does not contain 'cfg' field."))
    "hdr" in keys(dataset) ||
        throw(ArgumentError("Dataset does not contain 'hdr' field."))
    "fsample" in keys(dataset) ||
        throw(ArgumentError("Dataset does not contain 'fsample' field."))
    "nChans" in keys(dataset["hdr"]) ||
        throw(ArgumentError("Dataset header does not contain 'nChans' field."))

    cfg = dataset["cfg"]
    hdr = dataset["hdr"]
    sampling_rate = round(Int64, dataset["fsample"])
    ch_n = Int64(hdr["nChans"])

    # ------------------------------------------------------------------ #
    # channel labels, types, units                                       #
    # ------------------------------------------------------------------ #
    clabels = if length(hdr["label"][:]) > 0 && length(vec(hdr["label"])) == ch_n
        string.(strip.(string.(hdr["label"])))[:]
    else
        ["ch_$idx" for idx in 1:ch_n]
    end

    ch_type, units = if detect_type
        # FIX: original omitted the required second argument to _set_channel_types
        ct = _set_channel_types(clabels, data_type)
        u  = "chanunit" in keys(hdr) ?
            replace.(strip.(string.(hdr["chanunit"][:])), "uV" => "μV") :
            [_ch_units(ct[i]) for i in 1:ch_n]
        ct, u
    else
        ct = "chantype" in keys(hdr) ?
            string.(hdr["chantype"][:]) : repeat([data_type], ch_n)
        u  = "chanunit" in keys(hdr) ?
            replace.(string.(hdr["chanunit"][:]), "uV" => "μV") :
            repeat(["μV"], ch_n)
        ct, u
    end

    # normalize channel-type labels to NeuroAnalyzer conventions
    ch_type = replace.(ch_type,
        "nirs" => "nirs_od",
        "aux" => "other",
        "stimulus" => "mrk",
        "analog trigger" => "mrk",
        "digital trigger" => "mrk",
        "unknown" => "other")

    # ------------------------------------------------------------------ #
    # signal data (trials)                                                #
    # ------------------------------------------------------------------ #
    "trial" in keys(dataset) ||
        throw(ArgumentError("Dataset does not contain 'trial' field."))
    ep_n   = size(dataset["trial"], 2)
    ep_len = size(dataset["trial"][1], 2)
    data   = zeros(ch_n, ep_len, ep_n)
    for idx in 1:ep_n
        data[:, :, idx] = dataset["trial"][idx]
    end

    # ------------------------------------------------------------------ #
    # time axes                                                            #
    # ------------------------------------------------------------------ #
    if "time" in keys(dataset)
        if ep_n == 1
            raw_t = Float64.(dataset["time"][1][:])
            epoch_time = round.(raw_t .- raw_t[1]; digits = 4)
            time_pts = epoch_time
        else
            epoch_time = round.(Float64.(dataset["time"][1][:]); digits = 4)
            time_pts = round.(range(0, step = 1/sampling_rate,
                                    length = size(data,2)*size(data,3)); digits = 4)
        end
    else
        epoch_time = round.(range(0, step = 1/sampling_rate, length = size(data,2)); digits = 4)
        time_pts   = round.(range(0, step = 1/sampling_rate,
                                   length = size(data,2)*size(data,3)); digits = 4)
    end

    _info("FieldTrip markers are stored separately; import with " *
          "`import_ft(file_name, type=:events)` and add with `add_markers()`.")
    markers = DataFrame(
        :id => String[], :start => Float64[],
        :length => Float64[], :value => String[], :channel => Int64[])

    locs = _initialize_locs()

    r = nothing # will be assigned in each branch below

    if data_type == "eeg"

        clabels = _clean_eeg_labels(clabels)

        ref = if "reref" in keys(cfg) && cfg["reref"] != "no"
            _info("Embedded referencing is not supported; " *
                  "please send this file to adam.wysokinski@neuroanalyzer.org")
            "" # safe fallback
        else
            _detect_montage(clabels, ch_type, data_type)
        end

        @inbounds for ch_idx in 1:ch_n
            if units[ch_idx] == "V" && ch_type[ch_idx] in ("eeg", "emg", "eog", "ref")
                @views data[ch_idx, :, 1] .*= 1e6
                units[ch_idx] = "μV"
            end
        end

        r = _create_recording_eeg(
            data_type = data_type,
            file_name = file_name,
            file_size_mb = round(filesize(file_name) / 1024^2; digits = 2),
            file_type = file_type,
            recording = "RID" in keys(hdr["orig"]) ? string(hdr["orig"]["RID"]) : "",
            recording_date = "", recording_time = "", recording_notes = "",
            channel_type = ch_type,
            channel_order = _sort_channels(ch_type),
            reference = ref,
            clabels = clabels,
            transducers = "Transducer" in keys(hdr["orig"]) ?
                                string.(strip.(hdr["orig"]["Transducer"])) :
                                repeat([""], ch_n),
            units = units,
            prefiltering = "PreFilt" in keys(hdr["orig"]) ?
                                string.(strip.(hdr["orig"]["PreFilt"])) :
                                repeat([""], ch_n),
            line_frequency = 50, # TODO: make this a keyword argument
            sampling_rate = sampling_rate,
            gain = ones(ch_n),
            bad_channels = zeros(Bool, size(data, 1)))

    elseif data_type == "meg"

        clabels = _clean_meg_labels(clabels)

        "grad" in keys(dataset) ||
            throw(ArgumentError("Dataset does not contain 'grad' field."))
        "chantype" in keys(dataset["grad"]) ||
            throw(ArgumentError("grad does not contain 'chantype' field."))

        coil_type = repeat([""], ch_n)
        mag_idx = occursin.(r".*mag.*", lowercase.(ch_type))
        grad_idx  = occursin.(r".*grad.*", lowercase.(ch_type))
        pgrad_idx = occursin.(r".*planar.*", lowercase.(ch_type))
        agrad_idx = occursin.(r".*axial.*", lowercase.(ch_type)) .|
                    occursin.(r".*ctf.*", lowercase.(ch_type))

        combined_grad_idx = grad_idx .| pgrad_idx .| agrad_idx
        coil_type[mag_idx] .= "mag"
        coil_type[pgrad_idx] .= "pgrad"
        coil_type[agrad_idx] .= "agrad"
        coil_type[grad_idx .& .!pgrad_idx .& .!agrad_idx] .= "grad"

        magnetometers = findall(mag_idx)
        gradiometers = findall(combined_grad_idx)
        ch_type[magnetometers] .= "mag"
        ch_type[gradiometers] .= "grad"
        ch_type[occursin.("ias", lowercase.(clabels))] .= "other"
        ch_type[occursin.("sti", lowercase.(clabels))] .= "mrk"
        ch_type[occursin.("sys", lowercase.(clabels))] .= "other"

        @inbounds for ch_idx in 1:ch_n
            if units[ch_idx] == "T"
                @views data[ch_idx,:,1] .*= 1e15
                units[ch_idx] = "fT"
            elseif units[ch_idx] == "T/m"
                @views data[ch_idx,:,1] .*= (1e15/100)
                units[ch_idx] = "fT/cm"
            elseif units[ch_idx] == "T/cm"
                @views data[ch_idx,:,1] .*= 1e15
                units[ch_idx] = "fT/cm"
            elseif units[ch_idx] == "V" && ch_type[ch_idx] in ("eeg","emg","eog","ref")
                @views data[ch_idx,:,1] .*= 1e6
                units[ch_idx] = "μV"
            end
        end

        if "chanpos" in keys(dataset["grad"])
            sig_mask = ch_type .∈ Ref(["meg", "mag", "grad"])
            meg_labels = clabels[sig_mask]
            pos = dataset["grad"]["chanpos"]
            meg_locs = DataFrame(
                :label => meg_labels,
                :loc_radius => zeros(length(meg_labels)),
                :loc_theta => zeros(length(meg_labels)),
                :loc_x => pos[:, 1],
                :loc_y => pos[:, 2],
                :loc_z => pos[:, 3],
                :loc_radius_sph => zeros(length(meg_labels)),
                :loc_theta_sph => zeros(length(meg_labels)),
                :loc_phi_sph => zeros(length(meg_labels)))
            locs_normalize!(meg_locs); locs_cart2sph!(meg_locs); locs_sph2pol!(meg_locs)
            locs_scale!(meg_locs; r = 1.5)
            locs = meg_locs
        else
            locs = import_locs_csv(joinpath(NeuroAnalyzer.res_path, "meg_306flattened.csv"))
        end

        if "elec" in keys(dataset)
            eeg_mask = ch_type .== "eeg"
            eeg_labels = clabels[eeg_mask]
            epos = dataset["elec"]["chanpos"]
            eeg_locs = DataFrame(
                :label => eeg_labels,
                :loc_radius => zeros(length(eeg_labels)),
                :loc_theta => zeros(length(eeg_labels)),
                :loc_x => epos[:, 1],
                :loc_y => epos[:, 2],
                :loc_z => epos[:, 3],
                :loc_radius_sph => zeros(length(eeg_labels)),
                :loc_theta_sph => zeros(length(eeg_labels)),
                :loc_phi_sph => zeros(length(eeg_labels)))
            locs_normalize!(eeg_locs); locs_cart2sph!(eeg_locs); locs_sph2pol!(eeg_locs)
            locs_scale!(eeg_locs; r = 1.5)
            locs = vcat(locs, eeg_locs)
        end

        ref = if "reref" in keys(cfg) && cfg["reref"] != "no"
            _info("Embedded referencing is not supported; " *
                  "please send this file to adam.wysokinski@neuroanalyzer.org")
            ""
        else
            _detect_montage(clabels, ch_type, data_type)
        end

        # SSP projectors
        ssp_labels = String[]
        ssp_channels = zeros(Bool, ch_n)
        ssp_data = Matrix{Float64}(undef, 0, 0)
        if "projs" in keys(hdr["orig"])
            projs = hdr["orig"]["projs"]
            if length(projs["desc"]) != 0
                ssp_labels = string.(projs["desc"][:])
                ssp_ch_names = _clean_meg_labels(string.(projs["data"][1]["col_names"][:]))
                ssp_sort_idx = [findfirst(isequal(l), clabels) for l in ssp_ch_names]
                ssp_channels[Base.filter(!isnothing, ssp_sort_idx)] .= true
                ssp_data = zeros(length(ssp_labels), projs["data"][1]["ncol"])
                for idx in axes(ssp_data, 1)
                    ssp_data[idx, :] = projs["data"][idx]["data"][:]
                end
            end
        end

        lp = "lowpass" in keys(hdr["orig"]) ?
                string(round(hdr["orig"]["lowpass"][1];  digits=1)) : "?"
        hp = "highpass" in keys(hdr["orig"]) ?
                string(round(hdr["orig"]["highpass"][1]; digits=1)) : "?"

        r = _create_recording_meg(
            data_type = data_type,
            file_name = file_name,
            file_size_mb = round(filesize(file_name) / 1024^2; digits=2),
            file_type = "FT",
            recording = "dataformat" in keys(cfg) ? string(cfg["dataformat"]) : "",
            recording_date = "", recording_time = "", recording_notes = "",
            channel_type = ch_type,
            channel_order = _sort_channels(ch_type),
            reference = ref,
            clabels = clabels,
            units = units,
            prefiltering = repeat(["LP: $lp Hz; HP: $hp Hz"], ch_n),
            line_frequency = 50, # TODO: make this a keyword argument
            sampling_rate = sampling_rate,
            magnetometers = magnetometers,
            gradiometers = gradiometers,
            coil_type = coil_type,
            bad_channels = zeros(Bool, size(data,1)),
            ssp_labels = ssp_labels,
            ssp_channels = ssp_channels,
            ssp_data = ssp_data)

    elseif data_type == "nirs"

        clabels = replace.(clabels, ".0" => "", " [" => " ", "nm]" => "")

        "opto" in keys(dataset) ||
            throw(ArgumentError("Dataset does not contain 'opto' field."))
        opto = dataset["opto"]
        wavelengths = opto["wavelength"][:]

        wavelength_index = Int64[]
        for idx1 in eachindex(ch_type), idx2 in eachindex(wavelengths)
            if ch_type[idx1] == "nirs" &&
               occursin(string(round(Int64, wavelengths[idx2])), clabels[idx1])
                push!(wavelength_index, idx2)
            end
        end

        nirs_count = count(ch_type .== "nirs_od")
        src_mask = opto["optotype"] .== "transmitter"
        det_mask = opto["optotype"] .== "receiver"
        src_labels = string.(opto["optolabel"][:][src_mask])
        det_labels = string.(opto["optolabel"][:][det_mask])
        opt_labels = string.(opto["optolabel"][:])

        opt_pairs = zeros(Int64, nirs_count, 2)
        for (idx, pair) in enumerate(opto["label"][:])
            p = match(r"(.+)-(.+) (.*)", pair)
            if !isnothing(p)
                opt_pairs[idx, 1] = findfirst(isequal(p[1]), src_labels)
                opt_pairs[idx, 2] = findfirst(isequal(p[2]), det_labels)
            end
        end

        pos  = opto["optopos"]
        locs = DataFrame(
            :label => opt_labels,
            :loc_radius => zeros(length(opt_labels)),
            :loc_theta => zeros(length(opt_labels)),
            :loc_x => pos[:, 1],
            :loc_y => pos[:, 2],
            :loc_z => pos[:, 3],
            :loc_radius_sph => zeros(length(opt_labels)),
            :loc_theta_sph => zeros(length(opt_labels)),
            :loc_phi_sph => zeros(length(opt_labels))
        )
        locs_normalize!(locs); locs_cart2sph!(locs); locs_cart2pol!(locs)

        r = _create_recording_nirs(
            data_type = data_type,
            file_name = file_name,
            file_size_mb = round(filesize(file_name) / 1024^2; digits=2),
            file_type = file_type,
            recording = "", recording_date = "", recording_time = "",
            recording_notes = "",
            wavelengths = wavelengths,
            wavelength_index = wavelength_index,
            optode_pairs = opt_pairs,
            channel_type = ch_type,
            channel_order = _sort_channels(ch_type),
            clabels = clabels,
            units = units,
            src_labels = src_labels,
            det_labels = det_labels,
            opt_labels = opt_labels,
            sampling_rate = sampling_rate,
            bad_channels = zeros(Bool, size(data, 1)))
    end

    s   = _create_subject(id = "",
                          first_name = "",
                          middle_name = "",
                          last_name = "",
                          head_circumference = -1,
                          handedness = "",
                          weight = -1,
                          height = -1)
    e   = _create_experiment(name = "", notes = "", design = "")
    hdr = _create_header(subject = s, recording = r, experiment = e)

    obj = NeuroAnalyzer.NEURO(hdr, String[], markers, locs, time_pts, epoch_time, data)
    data_type == "eeg" && _initialize_locs!(obj)

    _info("Imported: " *
        uppercase(obj.header.recording[:data_type]) *
        " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj))" *
        "; $(round(obj.time_pts[end]; digits=2)) s)")

    return obj

end