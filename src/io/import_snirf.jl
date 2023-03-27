export import_snirf

"""
    import_snirf(file_name; n)

Load Shared Near Infrared Spectroscopy Format (SNIRF) file and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load
- `n::Int64=0`: subject number to extract in case of multi-subject file

# Returns

- `obj::NeuroAnalyzer.NEURO`

# Source

https://github.com/fNIRS/snirf/blob/v1.1/snirf_specification.md
"""
function import_snirf(file_name::String; n::Int64=0)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    splitext(file_name)[2] == ".snirf" || throw(ArgumentError("This is not a SNIRF file."))

    nirs = nothing
    try
        nirs = FileIO.load(file_name)
    catch
        throw(ArgumentError("File $file_name cannot be loaded."))
    end

    file_type = "SNIRF"

    parse(Float64, nirs["formatVersion"][1]) > 1.0 && _info("SNIRF version >1.0 detected.")

    # check for multi-subject recordings
    n_id = "nirs"
    n !== 0 && any(occursin.("nirs$n" , keys(nirs))) == false && throw(ArgumentError("No data for subject $n found in the recording."))
    if any(occursin.("nirs1" , keys(nirs))) == true
        if n == 0
            throw(ArgumentError("This is a multi-subject SfNIR file. Subject number must be specified via 'n' parameter."))
        else
            n_id = "nirs$n"
        end
    end

    # read metadata
    subject_id = nirs["$n_id/metaDataTags/SubjectID"][1]
    recording_date = nirs["$n_id/metaDataTags/MeasurementDate"][1]
    recording_time = nirs["$n_id/metaDataTags/MeasurementTime"][1]
    length_unit = nirs["$n_id/metaDataTags/LengthUnit"][1]
    time_unit = nirs["$n_id/metaDataTags/TimeUnit"][1]
    frq_unit = nirs["$n_id/metaDataTags/FrequencyUnit"][1]

    # probes

    # List of wavelengths (in nm)
    wavelengths = nirs["$n_id/probe/wavelengths"]

    wavelengths_emission = nothing
    src_pos2d = nothing
    src_pos3d = nothing
    detector_pos2d = nothing
    detector_pos3d = nothing
    frequencies = nothing
    t_delay = nothing
    t_delay_width = nothing
    moment_orders = nothing
    t_cor_delay = nothing
    t_cor_delay_width = nothing
    src_labels = nothing
    detector_labels = nothing
    landmark_pos2d = nothing
    landmark_pos3d = nothing
    landmark_labels = nothing
    coord_system = nothing
    coord_system_desc = nothing
    local_index = nothing

    # List of emission wavelengths (in nm)
    k = "$n_id/probe/wavelengthsEmission"
    k in keys(nirs) && (wavelengths_emission = nirs[k])

    # Source 2-D positions in LengthUnit
    k = "$n_id/probe/sourcePos2D"
    k in keys(nirs) && (src_pos2d = nirs[k])

    # Source 3-D positions in LengthUnit
    k = "$n_id/probe/sourcePos3D"
    k in keys(nirs) && (src_pos3d = nirs[k])

    # Detector 2-D positions in LengthUnit
    k = "$n_id/probe/detectorPos2D"
    k in keys(nirs) && (detector_pos2d = nirs[k])

    # Detector 3-D positions in LengthUnit
    k = "$n_id/probe/detectorPos3D"
    k in keys(nirs) && (detector_pos3d = nirs[k])

    # Modulation frequency list
    k = "$n_id/probe/frequencies"
    k in keys(nirs) && (frequencies = nirs[k])

    # Time delays for gated time-domain data
    k = "$n_id/probe/timeDelays"
    k in keys(nirs) && (t_delay = nirs[k])

    # Time delay width for gated time-domain data
    k = "$n_id/probe/timeDelaysWidths"
    k in keys(nirs) && (t_delay_width = nirs[k])

    # Moment orders of the moment TD data
    k = "$n_id/probe/momentOrders"
    k in keys(nirs) && (moment_orders = nirs[k])

    # Time delays for DCS measurements
    k = "$n_id/probe/correlationTimeDelays"
    k in keys(nirs) && (t_cor_delay = nirs[k])

    # Time delay width for DCS measurements
    k = "$n_id/probe/correlationTimeDelayWidths"
    k in keys(nirs) && (t_cor_delay_width = nirs[k])

    # String arrays specifying source names
    k = "$n_id/probe/sourceLabels"
    k in keys(nirs) && (src_labels = nirs[k])

    # String arrays specifying detector names
    k = "$n_id/probe/detectorLabels"
    k in keys(nirs) && (detector_labels = nirs[k])

    # Anatomical landmark 2-D positions
    k = "$n_id/probe/landmarkPos2D"
    k in keys(nirs) && (landmark_pos2d = nirs[k])

    # Anatomical landmark 3-D positions
    k = "$n_id/probe/landmarkPos3D"
    k in keys(nirs) && (landmark_pos3d = nirs[k])

    # String arrays specifying landmark names
    k = "$n_id/probe/landmarkLabels"
    k in keys(nirs) && (landmark_labels = nirs[k])

    # Coordinate system used in probe description
    k = "$n_id/probe/coordinateSystem"
    k in keys(nirs) && (coord_system = nirs[k][1])

    # Description of coordinate system
    k = "$n_id/probe/coordinateSystemDescription"
    k in keys(nirs) && (coord_system_desc = nirs[k][1])

    # If source/detector index is within a module
    k = "$n_id/probe/useLocalIndex"
    k in keys(nirs) && (local_index = nirs[k][1])

    # measurements
    data_n = 0
    while true
        data_n += 1
        if "$n_id/data$data_n/dataTimeSeries" in keys(nirs)
            continue
        else
            break
        end
    end
    data_n -= 1
    data_n > 1 && _info("Multiple data SNIRF files are not supported yet.")

    d_id = "data1"
    
    data = nirs["$n_id/$d_id/dataTimeSeries"]

    time_pts = nirs["$n_id/$d_id/time"]
    if length(time_pts) > 2
        sampling_rate = 1 / (time_pts[2] - time_pts[1])
    else
        sampling_rate = 1 / time_pts[2]
        time_pts = collect(time_pts[1]:1/sampling_rate:time_pts[1]+size(data, 2)*time_pts[2])[1:(end - 1)]
    end
    epoch_time = time_pts
    
    ch_n = size(data, 1)

    source_index = Int64[]
    detector_index = Int64[]
    wavelength_index = Int64[]
    wavelength_actual = Float64
    wavelength_emission_actual = Float64[]
    data_type = Int64[]
    data_unit = String[]
    data_type_label = String[]
    data_type_index = Int64[]
    source_power = Float64[]
    detector_gain = Float64[]
    module_index = Int64[]
    src_module_index = Int64[]
    detector_module_index = Int64[]

    for ch_idx in 1:ch_n
        # Source index for a given channel
        k = "$n_id/$d_id/measurementList$ch_idx/sourceIndex"
        k in keys(nirs) && (push!(source_index, Int.(nirs[k][1])))

        # Detector index for a given channel
        k = "$n_id/$d_id/measurementList$ch_idx/detectorIndex"
        k in keys(nirs) && (push!(detector_index, Int.(nirs[k][1])))

        # Wavelength index for a given channel
        k = "$n_id/$d_id/measurementList$ch_idx/wavelengthIndex"
        k in keys(nirs) && (push!(wavelength_index, Int.(nirs[k][1])))

        # Actual wavelength for a given channel
        k = "$n_id/$d_id/measurementList$ch_idx/wavelengthActual"
        k in keys(nirs) && (push!(wavelength_actual, nirs[k][1]))

        # Actual emission wavelength for a channel
        k = "$n_id/$d_id/measurementList$ch_idx/wavelengthEmissionActual"
        k in keys(nirs) && (push!(wavelength_emission_actual, nirs[k][1]))

        # Data type for a given channel
        k = "$n_id/$d_id/measurementList$ch_idx/dataType"
        k in keys(nirs) && (push!(data_type, Int.(nirs[k][1])))

        # SI unit for a given channel
        k = "$n_id/$d_id/measurementList$ch_idx/dataUnit"
        k in keys(nirs) && (push!(data_unit, nirs[k][1]))
        data_unit == String[] && (data_unit = repeat(["V"], ch_n))

        # Data type name for a given channel
        k = "$n_id/$d_id/measurementList$ch_idx/dataTypeLabel"
        k in keys(nirs) && (push!(data_type_label, nirs[k][1]))
        # assume its raw data (intensity) if there is no data type
        data_type_label == String[] && (data_type_label = repeat(["nirs_int"], ch_n))
        # Change in optical density
        data_type_label = replace(lowercase.(data_type_label), "dod" => "nirs_od")
        data_type_label = replace(lowercase.(data_type_label), "dmean" => "nirs_dmean")
        data_type_label = replace(lowercase.(data_type_label), "dvar" => "nirs_dvar")
        data_type_label = replace(lowercase.(data_type_label), "dskew" => "nirs_dskew")
        # Absorption coefficient
        data_type_label = replace(lowercase.(data_type_label), "mua" => "nirs_mua")
        # Scattering coefficient
        data_type_label = replace(lowercase.(data_type_label), "musp" => "nirs_musp")
        # Oxygenated hemoglobin (oxyhemoglobin) concentration
        data_type_label = replace(lowercase.(data_type_label), "hbo" => "nirs_hbo")
        # Deoxygenated hemoglobin (deoxyhemoglobin) concentration
        data_type_label = replace(lowercase.(data_type_label), "hbr" => "nirs_hbr")
        # Total hemoglobin concentration
        data_type_label = replace(lowercase.(data_type_label), "hbt" => "nirs_hbt")
        # Water content
        data_type_label = replace(lowercase.(data_type_label), "h2o" => "nirs_h2o")
        # Lipid concentration
        data_type_label = replace(lowercase.(data_type_label), "lipid" => "nirs_lipid")
        # Hemodynamic response function for blood flow index (BFi)
        data_type_label = replace(lowercase.(data_type_label), "bfi" => "nirs_bfi")
        # Hemodynamic response function for change in optical density
        data_type_label = replace(lowercase.(data_type_label), "hrf_dod" => "nirs_hrf_dod")
        data_type_label = replace(lowercase.(data_type_label), "hrf_dmean" => "nirs_hrf_dmean")
        data_type_label = replace(lowercase.(data_type_label), "hrf_dvar" => "nirs_hrf_dvar")
        data_type_label = replace(lowercase.(data_type_label), "hrf_dskew" => "nirs_hrf_dskew")
        # Hemodynamic response function for oxyhemoglobin concentration
        data_type_label = replace(lowercase.(data_type_label), "hrf_hbo" => "nirs_hrf_hbo")
        # emodynamic response function for deoxyhemoglobin concentration
        data_type_label = replace(lowercase.(data_type_label), "hrf_hbr" => "nirs_hrf_hbr")
        # Hemodynamic response function for total hemoglobin concentration
        data_type_label = replace(lowercase.(data_type_label), "hrf_hbt" => "nirs_hrf_hbt")
        # Hemodynamic response function for blood flow index (BFi)
        data_type_label = replace(lowercase.(data_type_label), "hrf_bfi" => "nirs_hrf_bfi")

        # Data type index for a given channel
        k = "$n_id/$d_id/measurementList$ch_idx/dataTypeIndex"
        k in keys(nirs) && (push!(data_type_index, Int.(nirs[k][1])))

        # Source power for a given channel
        k = "$n_id/$d_id/measurementList$ch_idx/sourcePower"
        k in keys(nirs) && (push!(source_power, nirs[k][1]))

        # Detector gain for a given channel
        k = "$n_id/$d_id/measurementList$ch_idx/detectorGain"
        k in keys(nirs) && (push!(detector_gain, nirs[k][1]))

        # Index of the parent module (if modular)
        k = "$n_id/$d_id/measurementList$ch_idx/moduleIndex"
        k in keys(nirs) && (push!(module_index, Int.(nirs[k][1])))

        # Index of the source's parent module
        k = "$n_id/$d_id/measurementList$ch_idx/sourceModuleIndex"
        k in keys(nirs) && (push!(src_module_index, Int.(nirs[k][1])))

        # Index of the detector's parent module
        k = "$n_id/$d_id/measurementList$ch_idx/detectorModuleIndex"
        k in keys(nirs) && (push!(detector_module_index, Int.(nirs[k][1])))
    end

    # currently data type is not used
    if data_type !== nothing
        tmp = String[]
        for idx in 1:length(data_type)
            data_type[idx] == 1 && push!(tmp, "Amplitude")
            data_type[idx] == 51 && push!(tmp, "Fluorescence Amplitude")
            data_type[idx] == 101 && push!(tmp, "Raw: Frequency Domain (FD): AC Amplitude")
            data_type[idx] == 102 && push!(tmp, "Raw: Frequency Domain (FD): Phase")
            data_type[idx] == 151 && push!(tmp, "Raw: Frequency Domain (FD): Fluorescence Amplitude")
            data_type[idx] == 152 && push!(tmp, "Raw: Frequency Domain (FD): Fluorescence Phase")
            data_type[idx] == 201 && push!(tmp, "Raw: Time Domain: Gated (TD Gated): Amplitude")
            data_type[idx] == 251 && push!(tmp, "Raw: Time Domain: Gated (TD Gated): Fluorescence Amplitude")
            data_type[idx] == 301 && push!(tmp, "Raw: Time Domain: Moments (TD Moments): Amplitude")
            data_type[idx] == 351 && push!(tmp, "Raw: Time Domain: Moments (TD Moments): Fluorescence Amplitude")
            data_type[idx] == 351 && push!(tmp, "Raw: Diffuse Correlation Spectroscopy (DCS): g2")
            data_type[idx] == 410 && push!(tmp, "Raw: Diffuse Correlation Spectroscopy (DCS): BFi")
            data_type[idx] == 99999 && push!(tmp, "Processed")
        end
        data_type = tmp
    end

    # collect channels
    ch_pairs = zeros(Int64, ch_n, 2)
    clabels = repeat([""], ch_n)
    for idx in 1:ch_n
        ch_pairs[idx, :] = hcat(source_index[idx], detector_index[idx])
        clabels[idx] = "S" * string(source_index[idx]) * "_D" * string(detector_index[idx]) * " " * string(wavelengths[wavelength_index[idx]])
    end
    clabels = replace.(clabels, ".0"=>"")

    # source and detector names
    if src_labels === nothing
        s = sort(unique(source_index))
        src_labels = String[]
        for idx in 1:length(s)
            push!(src_labels, "S" * string(s[idx]))
        end
    end
    if detector_labels === nothing
        d = sort(unique(detector_index))
        detector_labels = String[]
        for idx in 1:length(d)
            push!(detector_labels, "D" * string(d[idx]))
        end
    end
    opt_labels = vcat(src_labels, detector_labels)

    # stimulus measurements
    stim_n = 0
    while true
        stim_n += 1
        if "$n_id/stim$stim_n/name" in keys(nirs)
            continue
        else
            break
        end
    end
    stim_n -= 1
    stim_n > 1 && _info("Multiple stimulus SNIRF files are not supported yet.")

    s_id = "stim1"
    stim_data = nothing

    # Name of the stimulus data
    k = "$n_id/$s_id/name"
    k in keys(nirs) && (stim_name = nirs[k][1])

    # Data stream of the stimulus channel
    k = "$n_id/$s_id/data"
    k in keys(nirs) && (stim_data = nirs[k])

    # Data stream of the stimulus channel
    k = "$n_id/$s_id/dataLabels"
    k in keys(nirs) && (stim_labels = nirs[k])

    if stim_data !== nothing
        markers = DataFrame(:id=>repeat([""], size(stim_data, 2)), :start=>stim_data[1, :], :length=>stim_data[2, :], :description=>stim_name, :channel=>repeat([0], size(stim_data, 2)))
        # generate unique IDs
        desc = unique(markers[!, :description])
        for idx1 in 1:nrow(markers), idx2 in 1:length(desc)
            markers[idx1, :description] == desc[idx2] && (markers[idx1, :id] = string(idx2))
        end
    else
        markers = DataFrame(:id=>String[], :start=>Int64[], :length=>Int64[], :description=>String[], :channel=>Int64[])
    end

    # auxiliary measurements
    aux_n = 0
    while true
        aux_n += 1
        if "$n_id/aux$aux_n/name" in keys(nirs)
            continue
        else
            break
        end
    end
    aux_n -= 1
    aux_n > 1 && _info("Multiple aux SNIRF files are not supported yet.")

    a_id = "aux$aux_n"
    aux_data = nothing
    
    # Name of the auxiliary channel
    k = "$n_id/$a_id/name"
    k in keys(nirs) && (aux_name = nirs[k][1])

    # Data acquired from the auxiliary channel
    k = "$n_id/$a_id/dataTimeSeries"
    k in keys(nirs) && (aux_data = nirs[k])

    # SI unit of the auxiliary channel
    k = "$n_id/$a_id/dataUnit"
    k in keys(nirs) && (aux_unit = nirs[k])

    # Time (in TimeUnit) for auxiliary data 
    k = "$n_id/$a_id/time"
    k in keys(nirs) && (aux_time = nirs[k])

    # Time offset of auxiliary channel data
    k = "$n_id/$a_id/timeOffset"
    k in keys(nirs) && (aux_timeoffset = nirs[k])

    if aux_data !== nothing
        data = vcat(data, aux_data)
        for idx in 1:size(aux_data, 1)
            push!(clabels, "AUX$idx")
            push!(ch_type, "nirs_aux")
            data_unit = vcat(data_unit, aux_unit)
        end
    end

    # locations
    pos2d = hcat(src_pos2d, detector_pos2d)
    pos3d = hcat(src_pos3d, detector_pos3d)
    if src_pos3d === nothing
        if src_pos2d === nothing
            _info("The data does not contain 3D nor 2D location information for the optode positions.")
            x = zeros(length(opt_labels))
        else
            _info("The data only contains 2D location information for the optode positions.")
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
    
    s = _create_subject(id=subject_id,
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
                               recording_date=recording_date,
                               recording_time=recording_time,
                               recording_notes="",
                               wavelengths=wavelengths,
                               wavelength_index=wavelength_index,
                               channel_pairs=ch_pairs,
                               ch_type=data_type_label,
                               clabels=clabels,
                               units=data_unit,
                               opt_labels=opt_labels,
                               sampling_rate=round(Int64, sampling_rate))
    e = _create_experiment(name="", notes="", design="")

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = String[]

    return NeuroAnalyzer.NEURO(hdr, time_pts, epoch_time, data[:, :, :], components, markers, locs, history)

end
