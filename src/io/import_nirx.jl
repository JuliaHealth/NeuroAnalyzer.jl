export import_nirx

"""
    import_nirx(file_name)

Load NIRX file and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load, should point to .hdr file.

# Returns

- `NeuroAnalyzer.NEURO`

# References

1. https://nirx.net/file-formats
"""
function import_nirx(file_name::String)::NeuroAnalyzer.NEURO
    # validate
    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))
    detect_format(file_name) === :nirx ||
        throw(ArgumentError("$file_name is not NIRX file."))

    hdr = try
        readlines(file_name)
    catch err
        throw(ArgumentError("File $file_name cannot be loaded: $err"))
    end
    strip(chopprefix(hdr[1], "\ufeff")) == "[GeneralInfo]" ||
        throw(ArgumentError("File $file_name is not NIRX file."))

    file_type = "NIRX"

    # parse header
    hdr = replace.(hdr, "=\"" => "=", "\"" => "")
    recording_date = ""
    recording_time = ""
    device = ""
    source = ""
    mod = ""
    apd = ""
    nirstar = ""
    subject_id = ""

    recording_date = split(hdr[startswith.(lowercase.(hdr), "date=")][1], '=')[2]
    recording_time = split(hdr[startswith.(lowercase.(hdr), "time=")][1], '=')[2]
    device = split(hdr[startswith.(lowercase.(hdr), "device=")][1], '=')[2]
    source = split(hdr[startswith.(lowercase.(hdr), "source=")][1], '=')[2]
    mod = split(hdr[startswith.(lowercase.(hdr), "mod=")][1], '=')[2]
    apd = split(hdr[startswith.(lowercase.(hdr), "apd=")][1], '=')[2]
    nirstar = split(hdr[startswith.(lowercase.(hdr), "nirstar=")][1], '=')[2]
    subject_id = split(hdr[startswith.(lowercase.(hdr), "subject=")][1], '=')[2]

    sources = -1
    detectors = -1
    shortbundles = -1
    shortdetindex = -1
    steps = -1
    wavelengths = -1
    trigins = -1
    trigouts = -1
    anins = -1
    sampling_rate = -1
    modamp = -1
    threshold = -1

    any(startswith.(lowercase.(hdr), "sources=")) &&
        (sources = split(hdr[startswith.(lowercase.(hdr), "sources=")][1], '=')[2])
    sources = parse(Int64, sources)
    any(startswith.(lowercase.(hdr), "detectors=")) &&
        (detectors = split(hdr[startswith.(lowercase.(hdr), "detectors=")][1], '=')[2])
    detectors = parse(Int64, detectors)
    shortbundles = split(hdr[startswith.(lowercase.(hdr), "shortbundles=")][1], '=')[2]
    shortbundles = parse(Int64, shortbundles)
    if any(startswith.(lowercase.(hdr), "shortdetindex="))
        shortdetindex = split(
            hdr[startswith.(lowercase.(hdr), "shortdetindex=")][1], '=',
        )[2]
        shortdetindex = parse.(Int64, split(shortdetindex, '\t'))
    end
    any(startswith.(lowercase.(hdr), "steps=")) &&
        (steps = split(hdr[startswith.(lowercase.(hdr), "steps=")][1], '=')[2])
    steps = parse(Int64, steps)
    if any(startswith.(lowercase.(hdr), "wavelengths="))
        wavelengths = split(hdr[startswith.(lowercase.(hdr), "wavelengths=")][1], '=')[2]
        wavelengths = parse.(Float64, split(wavelengths, '\t'))
    end
    any(startswith.(lowercase.(hdr), "trigins=")) &&
        (trigins = split(hdr[startswith.(lowercase.(hdr), "trigins=")][1], '=')[2])
    trigins = parse(Int64, trigins)
    any(startswith.(lowercase.(hdr), "trigouts=")) &&
        (trigouts = split(hdr[startswith.(lowercase.(hdr), "trigouts=")][1], '=')[2])
    trigouts = parse(Int64, trigouts)
    any(startswith.(lowercase.(hdr), "anins=")) &&
        (anins = split(hdr[startswith.(lowercase.(hdr), "anins=")][1], '=')[2])
    anins = parse(Int64, anins)
    any(startswith.(lowercase.(hdr), "samplingrate=")) && (
        sampling_rate =
            split(hdr[startswith.(lowercase.(hdr), "samplingrate=")][1], '=')[2]
    )
    # exact rate of the recording (converted to an integer rate after reading the data)
    sampling_rate_exact = parse(Float64, sampling_rate)
    if any(startswith.(lowercase.(hdr), "mod amp="))
        modamp = split(hdr[startswith.(lowercase.(hdr), "mod amp=")][1], '=')[2]
        modamp = parse.(Float64, split(modamp, '\t'))
    end
    if any(startswith.(lowercase.(hdr), "threshold="))
        threshold = split(hdr[startswith.(lowercase.(hdr), "threshold=")][1], '=')[2]
        threshold = parse.(Float64, split(threshold, '\t'))
    end

    stim_type = ""
    any(startswith.(lowercase.(hdr), "stimulustype=")) &&
        (stim_type = split(hdr[startswith.(lowercase.(hdr), "stimulustype=")][1], '=')[2])
    notes = ""
    any(startswith.(lowercase.(hdr), "notes=")) &&
        (notes = split(hdr[startswith.(lowercase.(hdr), "notes=")][1], '=')[2])

    # parse .inf
    subject = ["", ""]
    age = -1
    gender = ""
    study_type1 = ""
    study_type2 = ""
    study_type3 = ""
    if isfile(splitext(file_name)[1] * ".inf")
        inf = readlines(splitext(file_name)[1] * ".inf")
        inf = replace.(inf, "=\"" => "=", "\"" => "")
        subject = split(inf[findfirst(startswith.(lowercase.(inf), "name="))], "=")[2]
        subject = split(subject, "\\0")
        age = parse(
            Float64, split(inf[findfirst(startswith.(lowercase.(inf), "age="))], "=")[2],
        )
        gender = split(inf[findfirst(startswith.(lowercase.(inf), "gender="))], "=")[2]
        study_type1 = split(
            inf[findfirst(startswith.(lowercase.(inf), "study type="))], "=",
        )[2]
        study_type2 = split(
            inf[findfirst(startswith.(lowercase.(inf), "experiment history="))], "=",
        )[2]
        study_type3 = split(
            inf[findfirst(startswith.(lowercase.(inf), "additional notes="))], "=",
        )[2]
    end

    # parse gains if .set is not available (informational only)
    gains = zeros(Int64, sources, detectors)
    gain_rows = if !isfile(splitext(file_name)[1] * ".set")
        gains_start = findfirst(startswith.(hdr, "Gains="))
        isnothing(gains_start) ? String[] : hdr[(gains_start + 1):(gains_start + sources)]
    else
        readlines(splitext(file_name)[1] * ".set")
    end
    for (idx, row) in enumerate(gain_rows[1:min(end, sources)])
        v = parse.(Int64, split(strip(row)))
        gains[idx, 1:min(end, length(v))] = v[1:min(end, detectors)]
    end

    # S-D-Key: "s-d:column,…" maps each column of the .wl files to a source-detector pair
    pairs = hdr[findfirst(startswith.(lowercase.(hdr), "s-d-key="))]
    pairs = Base.filter(!isempty, split(replace(lowercase(pairs), "s-d-key=" => ""), ","))
    key_pairs = zeros(Int64, length(pairs), 2)
    for (idx, pr) in enumerate(pairs)
        sd, col = split(pr, ":")
        s_idx, d_idx = parse.(Int64, split(sd, "-"))
        key_pairs[parse(Int64, col), :] = [s_idx, d_idx]
    end

    # S-D-Mask: sources × detectors matrix of active channels
    ch_mask_start = findfirst(startswith.(lowercase.(hdr), "s-d-mask="))
    masks = split.(strip.(hdr[(ch_mask_start + 1):(ch_mask_start + sources)]))
    sd_mask = [parse(Int64, masks[s_idx][d_idx]) == 1 for s_idx = 1:sources, d_idx = 1:detectors]
    # active .wl columns, in S-D-Key order
    ch_masks = [sd_mask[key_pairs[idx, 1], key_pairs[idx, 2]] for idx in axes(key_pairs, 1)]
    opt_pairs_wl = key_pairs[ch_masks, :]
    ch_n = count(ch_masks)

    # parse dark noise ???
    dark_noise = zeros(length(wavelengths), detectors)
    for wv_idx in eachindex(wavelengths)
        dn_idx = findfirst(startswith.(lowercase.(hdr), "wavelength$wv_idx="))
        isnothing(dn_idx) && continue
        dark_noise[wv_idx, :] = parse.(Float64, split(strip(hdr[dn_idx + 1])))[1:detectors]
    end

    # channel distances (informational only)
    chd_idx = findfirst(startswith.(lowercase.(hdr), "chandis="))
    channel_distance = isnothing(chd_idx) ? Float64[] :
        parse.(Float64, split(strip(replace(lowercase(hdr[chd_idx]), "chandis=" => ""))))

    # read raw light intensity channels (V): one .wl file per wavelength,
    # rows = time points, columns = all source-detector pairs (S-D-Key order)
    read_wl(idx) = Matrix{Float64}(
        CSV.read(
            splitext(file_name)[1] * ".wl$idx",
            DataFrame;
            header = false,
            delim = ' ',
            ignorerepeated = true,
        ),
    )'[ch_masks, :]
    nirs_int = reduce(vcat, [read_wl(idx) for idx in eachindex(wavelengths)])
    wavelength_index = reduce(vcat, [repeat([idx], ch_n) for idx in eachindex(wavelengths)])
    opt_pairs = reduce(vcat, [opt_pairs_wl for _ in eachindex(wavelengths)])
    data_type_label = repeat(["nirs_int"], size(nirs_int, 1))
    data_unit = repeat(["V"], size(nirs_int, 1))
    ch_n = size(nirs_int, 1)
    data = reshape(nirs_int, ch_n, size(nirs_int, 2), 1)
    # NeuroAnalyzer stores integer sampling rates (e.g. 12.5 Hz → upsampled to 25 Hz)
    data, sampling_rate = _integer_rate(data, sampling_rate_exact)

    time_pts = _time_axis(size(data, 2), sampling_rate)
    epoch_time = time_pts

    # events: [Markers] Events= table in the header ("time  code  frame"), otherwise the .evt file
    # ("frame  bit1 … bit8", the trigger code is the binary value with bit1 as the lowest bit);
    # frames are 0-based
    stim_frame = Int64[]
    stim_id = String[]
    events_start = findfirst(startswith.(hdr, "Events="))
    event_rows = if !isnothing(events_start)
        rows = String[]
        for l in hdr[(events_start + 1):end]
            startswith(l, "#") && break
            isempty(strip(l)) || push!(rows, l)
        end
        rows
    else
        String[]
    end
    if !isempty(event_rows)
        for l in event_rows
            v = split(strip(l))
            push!(stim_id, string(round(Int64, parse(Float64, v[2]))))
            push!(stim_frame, round(Int64, parse(Float64, v[3])))
        end
    elseif isfile(splitext(file_name)[1] * ".evt")
        for l in readlines(splitext(file_name)[1] * ".evt")
            isempty(strip(l)) && continue
            v = parse.(Int64, split(strip(l)))
            push!(stim_frame, v[1])
            push!(stim_id, string(sum(v[i + 1] * 2^(i - 1) for i = 1:(length(v) - 1))))
        end
    end

    markers = DataFrame(
        # value: the condition (trigger) code, so that epochs can be selected by condition
        :id => repeat(["stim"], length(stim_id)),
        :start => stim_frame ./ sampling_rate_exact,
        :length => zeros(length(stim_id)),
        :value => stim_id,
        :channel => zeros(Int64, length(stim_id)),
    )

    # read probes data
    probes = matread(splitext(file_name)[1] * "_probeInfo.mat")
    # head_model = probes["probeInfo"]["headmodel"]
    # probes["probeInfo"]["probes"]["labels_o"]
    det_labels = probes["probeInfo"]["probes"]["labels_d"][:]
    src_labels = probes["probeInfo"]["probes"]["labels_s"][:]
    opt_labels = string.(vcat(src_labels, det_labels))

    clabels = [
        string(src_labels[opt_pairs[idx, 1]]) *
        "_" *
        string(det_labels[opt_pairs[idx, 2]]) *
        " " *
        string(wavelengths[wavelength_index[idx]]) for idx = 1:ch_n
    ]
    clabels = replace.(clabels, ".0" => "")

    # probes["probeInfo"]["probes"]["coords_c2"]
    # probes["probeInfo"]["probes"]["coords_c3"]
    # probes["probeInfo"]["probes"]["coords_o2"]
    # probes["probeInfo"]["probes"]["coords_o3"]
    src_pos2d = probes["probeInfo"]["probes"]["coords_s2"]'
    src_pos3d = probes["probeInfo"]["probes"]["coords_s3"]'
    detector_pos2d = probes["probeInfo"]["probes"]["coords_d2"]'
    detector_pos3d = probes["probeInfo"]["probes"]["coords_d3"]'

    # probes["probeInfo"]["probes"]["nChannel0"]
    # src_n = Int64(probes["probeInfo"]["probes"]["nSource0"])
    # det_n = Int64(probes["probeInfo"]["probes"]["nDetector0"])

    # probes["probeInfo"]["probes"]["normals_c"]
    # probes["probeInfo"]["probes"]["normals_o"]
    # probes["probeInfo"]["probes"]["normals_s"]
    # probes["probeInfo"]["probes"]["normals_d"]

    # TPL file ???
    # AVG file ???

    # locations: 3D coordinates when available, otherwise the 2D layout
    pos2d = Float64.(hcat(src_pos2d, detector_pos2d))
    pos3d = Float64.(hcat(src_pos3d, detector_pos3d))
    if size(pos3d, 1) >= 3 && any(!iszero, pos3d)
        x, y, z = pos3d[1, :], pos3d[2, :], pos3d[3, :]
        # normalize all axes jointly (keeps the head geometry)
        x, y, z = _locs_norm(x, y, z)
    elseif size(pos2d, 1) >= 2 && any(!iszero, pos2d)
        _warn("The data only contains 2D location information for the optode positions.")
        x, y = _locs_norm(pos2d[1, :], pos2d[2, :])
        z = zeros(length(opt_labels))
    else
        _warn(
            "The data does not contain 3D nor 2D location information for the optode positions.",
        )
        x = zeros(length(opt_labels))
        y = zeros(length(opt_labels))
        z = zeros(length(opt_labels))
    end
    radius = zeros(length(opt_labels))
    theta = zeros(length(opt_labels))
    radius_sph = zeros(length(opt_labels))
    theta_sph = zeros(length(opt_labels))
    phi_sph = zeros(length(opt_labels))
    locs = DataFrame(
        :label => opt_labels,
        :loc_radius => radius,
        :loc_theta => theta,
        :loc_x => x,
        :loc_y => y,
        :loc_z => z,
        :loc_radius_sph => radius_sph,
        :loc_theta_sph => theta_sph,
        :loc_phi_sph => phi_sph,
    )
    locs_cart2sph!(locs)
    locs_cart2pol!(locs)

    file_size_mb = round(filesize(file_name) / 1024^2; digits = 2)

    s = _create_subject(;
        id = string(subject_id[1]),
        first_name = string(subject[1]),
        middle_name = "",
        last_name = length(subject) >= 2 ? string(subject[2]) : "",
        head_circumference = -1,
        handedness = "",
        weight = -1,
        height = -1,
    )
    r = _create_recording_nirs(;
        data_type = "nirs",
        file_name = file_name,
        file_size_mb = file_size_mb,
        file_type = file_type,
        recording = string(device[1]),
        recording_date = string(recording_date[1]),
        recording_time = replace(string(recording_time[1]), '.' => ':'),
        recording_notes = "NIRStar: $nirstar",
        wavelengths = wavelengths,
        wavelength_index = wavelength_index,
        optode_pairs = opt_pairs,
        channel_type = data_type_label,
        channel_order = _sort_channels(data_type_label),
        clabels = clabels,
        units = data_unit,
        src_labels = string.(src_labels),
        det_labels = string.(det_labels),
        opt_labels = opt_labels,
        sampling_rate = sampling_rate,
        bad_channels = zeros(Bool, size(data, 1)),
    )
    e = _create_experiment(;
        name = string(study_type1),
        notes = string(study_type2),
        design = string(study_type3),
    )

    hdr = _create_header(; subject = s, recording = r, experiment = e)

    history = String[]

    obj = NeuroAnalyzer.NEURO(hdr, history, markers, locs, time_pts, epoch_time, data)

    _info(
        "Imported: " *
        uppercase(obj.header.recording[:data_type]) *
        " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits = 2)) s)",
    )

    return obj
end
