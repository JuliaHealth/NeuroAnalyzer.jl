export import_nirs

"""
    import_nirs(file_name)

Load NIRS file and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load

# Returns

- `obj::NeuroAnalyzer.NEURO`

# Source

https://github.com/BUNPC/Homer3/wiki/HOMER3-file-formats
"""
function import_nirs(file_name::String; n::Int64=0)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    splitext(file_name)[2] == ".nirs" || throw(ArgumentError("This is not a NIRS file."))

    nirs = nothing
    try
        nirs = matread(file_name)
    catch
        throw(ArgumentError("File $file_name cannot be loaded."))
    end

    file_type = "NIRS"

    "d" in keys(nirs) || _info("This does not seem to be a NIRS file.")
    "t" in keys(nirs) || _info("This does not seem to be a NIRS file.")

    # time points
    time_pts = nirs["t"][:]
    epoch_time = time_pts
    sampling_rate = round(Int64, 1 / (time_pts[2] - time_pts[1]))

    userdata = nirs["userdata"]

    # probes
    probes = nirs["SD"]
    # ch_n × 4
    # col1: source index
    # col2: detector index
    # col3: unused
    # col4: wavelength index
    tmp = Int.(probes["MeasList"])
    ch_n = size(tmp, 1)
    ch_pairs=zeros(Int64, size(tmp, 1), 2)
    for idx in 1:size(tmp, 1)
        ch_pairs[idx, 1] = tmp[idx, 1]
        ch_pairs[idx, 2] = tmp[idx, 2]
    end
    # probes["SpringList"]
    # probes["SrcMap"]
    # probes["vrnum"]
    # probes["MeasListVis"]
    # probes["SpatialUnit"]
    # probes["AnchorList"]
    # probes["MeasListAct"]
    # probes["xmax"]
    # probes["xmin"]
    # probes["ymax"]
    # probes["ymin"]
    # probes["nDummys"]
    # probes["DummyPos"]

    #list of wavelengths (in nm)
    wavelengths = probes["Lambda"][:]
    wavelength_index = Int.(probes["MeasList"])[:, 4]

    # sources and detectors labels
    opt_labels = String[]
    for idx in 1:Int(probes["nSrcs"])
        push!(opt_labels, "S$idx")
    end
    for idx in 1:Int(probes["nDets"])
        push!(opt_labels, "D$idx")
    end

    # channel labels
    clabels = repeat([""], ch_n)
    for idx in 1:ch_n
        clabels[idx] = "S" * string(Int.(probes["MeasList"])[idx, 1]) * "-D" * string(Int.(probes["MeasList"])[idx, 2])
    end

    # data contains intensity (RAW) data
    data_unit = repeat(["V"], ch_n)

    # measurements
    data = Matrix(nirs["d"]')

    # stimuli
    s = nirs["s"][:]
    markers = DataFrame(:id=>String[], :start=>Int64[], :length=>Int64[], :description=>String[], :channel=>Int64[])

    # auxiliary measurements
    aux = Matrix(nirs["aux"]')

    # locations
    src_pos3d = Matrix(probes["SrcPos"]')
    detector_pos3d = Matrix(probes["DetPos"]')
    pos3d = hcat(src_pos3d, detector_pos3d)
    x = pos3d[1, :]
    y = pos3d[2, :]
    z = pos3d[3, :]
    # swap x and y
    # x, y = y, x
    # normalize to a unit-sphere
    if z != zeros(length(opt_labels))
        x, y, z = _locnorm(x, y, z)
    else
        x, y = _locnorm(x, y)
    end
    radius = zeros(length(opt_labels))
    theta = zeros(length(opt_labels))
    radius_sph = zeros(length(opt_labels))
    theta_sph = zeros(length(opt_labels))
    phi_sph = zeros(length(opt_labels))
    locs = DataFrame(:channel=>1:length(opt_labels), :labels=>opt_labels, :loc_theta=>theta, :loc_radius=>radius, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)
    locs = locs_cart2sph(locs)
    locs = locs_cart2pol(locs)

    file_size_mb = round(filesize(file_name) / 1024^2, digits=2)
    
    s = _create_subject(id="",
                        first_name="",
                        middle_name="",
                        last_name="",
                        handedness="",
                        weight=-1,
                        height=-1)
    r = _create_recording_nirs(data_type="nirs",
                               file_name=file_name,
                               file_size_mb=file_size_mb,
                               file_type=file_type,
                               recording="",
                               recording_date="",
                               recording_time="",
                               recording_notes="",
                               wavelengths=wavelengths,
                               wavelength_index=wavelength_index,
                               channel_pairs=ch_pairs,
                               ch_type=repeat(["nirs_int"], ch_n),
                               clabels=clabels,
                               units=data_unit,
                               opt_labels=opt_labels,
                               sampling_rate=sampling_rate)
    e = _create_experiment(experiment_name="",
                           experiment_notes="",
                           experiment_design="")

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = String[]

    return NeuroAnalyzer.NEURO(hdr, time_pts, epoch_time, data[:, :, :], components, markers, locs, history)

end
