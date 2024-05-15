export import_nirx

"""
    import_nirx(file_name)

Load NIRX file and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load, should point to .hdr file.

# Returns

- `obj::NeuroAnalyzer.NEURO`

# Source

https://nirx.net/file-formats
"""
function import_nirx(file_name::String)

    @assert isfile(file_name) "File $file_name cannot be loaded."
    @assert lowercase(splitext(file_name)[2]) == ".hdr" "This is not NIRX .hdr file."

    hdr = nothing
    try
        hdr = readlines(file_name)
    catch
        @error "File $file_name cannot be loaded."
    end
    @assert hdr[1] == "[GeneralInfo]" "File $file_name is not NIRX file."

    file_type = "NIRX"

    # parse header
    hdr = replace.(hdr, "=\""=>"=", "\""=>"")
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

    any(startswith.(lowercase.(hdr), "sources=")) && (sources = split(hdr[startswith.(lowercase.(hdr), "sources=")][1], '=')[2])
    sources = parse(Int64, sources)
    any(startswith.(lowercase.(hdr), "detectors=")) && (detectors = split(hdr[startswith.(lowercase.(hdr), "detectors=")][1], '=')[2])
    detectors = parse(Int64, detectors)
    shortbundles = split(hdr[startswith.(lowercase.(hdr), "shortbundles=")][1], '=')[2]
    shortbundles = parse(Int64, shortbundles)
    if any(startswith.(lowercase.(hdr), "shortdetindex="))
        shortdetindex = split(hdr[startswith.(lowercase.(hdr), "shortdetindex=")][1], '=')[2]
        shortdetindex = parse.(Int64, split(shortdetindex, '\t'))
    end
    any(startswith.(lowercase.(hdr), "steps=")) && (steps = split(hdr[startswith.(lowercase.(hdr), "steps=")][1], '=')[2])
    steps = parse(Int64, steps)
    if any(startswith.(lowercase.(hdr), "wavelengths="))
        wavelengths = split(hdr[startswith.(lowercase.(hdr), "wavelengths=")][1], '=')[2]
        wavelengths = parse.(Float64, split(wavelengths, '\t'))
    end
    any(startswith.(lowercase.(hdr), "trigins=")) && (trigins = split(hdr[startswith.(lowercase.(hdr), "trigins=")][1], '=')[2])
    trigins = parse(Int64, trigins)
    any(startswith.(lowercase.(hdr), "trigouts=")) && (trigouts = split(hdr[startswith.(lowercase.(hdr), "trigouts=")][1], '=')[2])
    trigouts = parse(Int64, trigouts)
    any(startswith.(lowercase.(hdr), "anins=")) && (anins = split(hdr[startswith.(lowercase.(hdr), "anins=")][1], '=')[2])
    anins = parse(Int64, anins)
    any(startswith.(lowercase.(hdr), "samplingrate=")) && (sampling_rate = split(hdr[startswith.(lowercase.(hdr), "samplingrate=")][1], '=')[2])
    sampling_rate = round(Int64, parse(Float64, sampling_rate))
    if any(startswith.(lowercase.(hdr), "mod amp="))
        modamp = split(hdr[startswith.(lowercase.(hdr), "mod amp=")][1], '=')[2]
        modamp = parse.(Float64, split(modamp, '\t'))
    end
    if any(startswith.(lowercase.(hdr), "threshold="))
        threshold = split(hdr[startswith.(lowercase.(hdr), "threshold=")][1], '=')[2]
        threshold = parse.(Float64, split(threshold, '\t'))
    end

    stim_type = ""
    any(startswith.(lowercase.(hdr), "stimulustype=")) && (stim_type = split(hdr[startswith.(lowercase.(hdr), "stimulustype=")][1], '=')[2])
    notes = ""
    any(startswith.(lowercase.(hdr), "notes=")) && (notes = split(hdr[startswith.(lowercase.(hdr), "notes=")][1], '=')[2])

    # parse .inf
    #subject = ""
    age = -1
    gender = ""
    study_type1 = ""
    study_type2 = ""
    study_type3 = ""
    if isfile(splitext(file_name)[1] * ".inf")
        inf = readlines(splitext(file_name)[1] * ".inf")
        inf = replace.(inf, "=\""=>"=", "\""=>"")
        subject = split(inf[findfirst(startswith.(lowercase.(inf), "name="))], "=")[2]
        subject = split(subject, "\\0")
        age = parse(Float64, split(inf[findfirst(startswith.(lowercase.(inf), "age="))], "=")[2])
        gender = split(inf[findfirst(startswith.(lowercase.(inf), "gender="))], "=")[2]
        study_type1 = split(inf[findfirst(startswith.(lowercase.(inf), "study type="))], "=")[2]
        study_type2 = split(inf[findfirst(startswith.(lowercase.(inf), "experiment history="))], "=")[2]
        study_type3 = split(inf[findfirst(startswith.(lowercase.(inf), "additional notes="))], "=")[2]
    end

    # parse gains if .set is not available
    if !isfile(splitext(file_name)[1] * ".set")
        gains_start = findfirst(startswith.(hdr, "Gains="))
        buf = hdr[gains_start + 1:gains_start + sources]
        gains = zeros(Int64, sources, detectors)
        for idx in eachindex(buf)
            gains[idx, :] = parse.(Int64, split(buf[idx], '\t'))
        end
    else
        buf = readlines(splitext(file_name)[1] * ".set")
        buf = split.(buf, ' ')
        gains = zeros(Int64, sources, detectors)
        for idx in eachindex(buf)
            gains[idx, :] = parse.(Int64, buf[idx])
        end
    end

    # parse opt_pairs
    pairs = hdr[findfirst(startswith.(lowercase.(hdr), "s-d-key="))]
    pairs = split(replace(lowercase.(pairs), "s-d-key="=>""), ",")[1:end - 1]
    ch_n = length(pairs)
    opt_pairs = zeros(Int64, ch_n, 2)
    for idx in 1:ch_n
        opt_pairs[idx, :] = [parse(Int64, split(pairs[idx], "-")[1]), parse(Int64, split(split(pairs[idx], "-")[2], ":")[1])]
    end
    ch_mask_start = findfirst(startswith.(lowercase.(hdr), "s-d-mask="))
    masks = hdr[ch_mask_start + 1:ch_mask_start + sources]
    masks = split.(masks, '\t')
    ch_masks = zeros(Int64, sources * detectors)
    idx = 1
    for idx1 in 1:sources, idx2 in 1:detectors
        ch_masks[idx] = parse(Int64, masks[idx1][idx2])
        idx += 1
    end
    ch_n = sum(ch_masks)
    ch_masks = Bool.(ch_masks)
    opt_pairs = opt_pairs[ch_masks, :]

    # parse dark noise ???
    dark_noise = zeros(length(wavelengths), detectors)
    for wv_idx in eachindex(wavelengths)
        dn = hdr[findfirst(startswith.(lowercase.(hdr), "wavelength$wv_idx=")) + 1]
        dark_noise[wv_idx, :] = parse.(Float64, split.(dn, '\t')[1:end])
    end

    chd = hdr[findfirst(startswith.(lowercase.(hdr), "chandis="))]
    chd = replace(lowercase.(chd), "chandis="=>"")
    chd = split.(chd, '\t')
    channel_distance = zeros(length(chd))
    for idx in eachindex(chd)
        channel_distance[idx] = parse(Float64, chd[idx])
    end

    # read raw light intensity channels (V)
    nirs_int = Matrix(CSV.read(splitext(file_name)[1] * ".wl1", header=false, stringtype=String, DataFrame))'[ch_masks, :]
    wavelength_index = repeat([1], ch_n)
    data_type_label = repeat(["nirs_int"], ch_n)
    data_unit = repeat(["V"], ch_n)
    for idx in 2:length(wavelengths)
        nirs_int = vcat(nirs_int, Matrix(CSV.read(splitext(file_name)[1] * ".wl$idx", header=false, stringtype=String, DataFrame))'[ch_masks, :])
        wavelength_index = vcat(wavelength_index, repeat([idx], ch_n))
        opt_pairs = vcat(opt_pairs, opt_pairs)
        data_type_label = vcat(data_type_label, repeat(["nirs_int"], ch_n))
        data_unit = vcat(data_unit, repeat(["V"], ch_n))
    end
    ch_n = size(nirs_int, 1)

    time_pts = round.(collect(0:1/sampling_rate:(size(nirs_int, 2) / sampling_rate))[1:end - 1], digits=3)
    epoch_time = round.(collect(0:1/sampling_rate:(size(nirs_int, 2) / sampling_rate))[1:end - 1], digits=3)

    # parse events if .evt is not available
    stim_onset = nothing
    stim_id = String[]
    if any(startswith.(hdr, "Events="))
        events_start = findfirst(startswith.(hdr, "Events="))
        events_end = events_start + findfirst(startswith.(hdr[events_start:end], "#"))
        buf = hdr[events_start + 1:events_end - 2]
        buf = split.(buf, '\t')
        events = zeros(Float64, length(buf), length(buf[1]))
        for idx in eachindex(buf)
            events[idx, :] = parse.(Float64, buf[idx])
        end
        stim_onset = Int.(events[:, 3])
        stim_id = string.(Int.(events[:, 2]))
    elseif isfile(splitext(file_name)[1] * ".evt")
        buf = readlines(splitext(file_name)[1] * ".evt")
        buf = split.(buf, '\t')
        events = zeros(Int64, length(buf), length(buf[1]))
        for idx in eachindex(buf)
            events[idx, :] = parse.(Int64, buf[idx])
        end
        # what are those 0s and 1s in events[] ???
        stim_onset = events[:, 1]
        stim_id = String[]
        for idx in 1:size(events, 1)
            push!(stim_id, string(findfirst(isequal(1), events[idx, 2:end])))
        end
    end

    if stim_onset !== nothing
        markers = DataFrame(:id=>stim_id, :start=>stim_onset, :length=>repeat([1], length(stim_id)), :description=>repeat(["stim"], length(stim_id)), :channel=>zeros(Int64, length(stim_id)))
    else
        markers = DataFrame(:id=>nothing, :start=>nothing, :length=>nothing, :description=>nothing, :channel=>nothing)
    end

    # read data ???
    buf = readlines(splitext(file_name)[1] * ".dat")
    buf_r = length(parse.(Float64, split(buf[1], ' ')))
    data = zeros(buf_r, length(buf))
    for idx in eachindex(buf)
        data[:, idx] = parse.(Float64, split(buf[idx], ' '))
    end

    # read probes data
    probes = matread(splitext(file_name)[1] * "_probeInfo.mat")
    # head_model = probes["probeInfo"]["headmodel"]
    opt_pairs = Int64.(probes["probeInfo"]["probes"]["index_c"])
    opt_pairs = vcat(opt_pairs, opt_pairs)
    # probes["probeInfo"]["probes"]["labels_o"]
    det_labels = probes["probeInfo"]["probes"]["labels_d"][:]
    src_labels = probes["probeInfo"]["probes"]["labels_s"][:]
    opt_labels = string.(vcat(src_labels, det_labels))

    clabels = repeat([""], ch_n)
    for idx in 1:ch_n
        clabels[idx] = src_labels[opt_pairs[idx, :][1]] * "_" * det_labels[opt_pairs[idx, :][2]] * " " * string(wavelengths[wavelength_index[idx]])
    end
    clabels = replace.(clabels, ".0"=>"")

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

    # locations
    pos2d = hcat(src_pos2d, detector_pos2d)
    pos3d = hcat(src_pos3d, detector_pos3d)
    if src_pos3d === nothing
        if src_pos2d === nothing
            _warn("The data does not contain 3D nor 2D location information for the optode positions.")
            x = zeros(length(opt_labels))
        else
            _warn("The data only contains 2D location information for the optode positions.")
            x = pos2d[1, :]
        end
    else
        x = pos3d[1, :]
    end
    if src_pos3d === nothing
        if src_pos2d === nothing            
            y = zeros(length(opt_labels))
        else
            y = pos2d[2, :]
        end
    else
        y = pos3d[2, :]
    end
    if src_pos3d === nothing
        z = zeros(length(opt_labels))
    else
        z = pos3d[3, :]
    end
    # swap x and y
    # x, y = y, x
    # normalize to a unit-sphere
    z = normalize(z, method=:n)
    x, y = _locs_norm(x, y)
    radius = zeros(length(opt_labels))
    theta = zeros(length(opt_labels))
    radius_sph = zeros(length(opt_labels))
    theta_sph = zeros(length(opt_labels))
    phi_sph = zeros(length(opt_labels))
    locs = DataFrame(:labels=>opt_labels, :loc_radius=>radius, :loc_theta=>theta, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)
    locs_cart2sph!(locs)
    locs_cart2pol!(locs)

    file_size_mb = round(filesize(file_name) / 1024^2, digits=2)

    s = _create_subject(id=string(subject_id[1]),
                        first_name=string(subject[1]),
                        middle_name="",
                        last_name=string(subject[2]),
                        head_circumference=-1,
                        handedness="",
                        weight=-1,
                        height=-1)
    r = _create_recording_nirs(data_type="nirs",
                               file_name=file_name,
                               file_size_mb=file_size_mb,
                               file_type=file_type,
                               recording=string(device[1]),
                               recording_date=string(recording_date[1]),
                               recording_time=string(recording_time[1]),
                               recording_notes="NIRStar: $nirstar",
                               wavelengths=wavelengths,
                               wavelength_index=wavelength_index,
                               optode_pairs=opt_pairs,
                               ch_type=data_type_label,
                               clabels=clabels,
                               units=data_unit,
                               src_labels=string.(src_labels),
                               det_labels=string.(det_labels),
                               opt_labels=opt_labels,
                               sampling_rate=sampling_rate)
    e = _create_experiment(name=string(study_type1), notes=string(study_type2), design=string(study_type3))

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = String[]

    obj = NeuroAnalyzer.NEURO(hdr, time_pts, epoch_time, data[:, :, :], components, markers, locs, history)

    _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits=2)) s)")

    return obj

end
