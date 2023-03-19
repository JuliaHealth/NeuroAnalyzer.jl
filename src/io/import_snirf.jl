# export import_snirf

"""
    import_snirf(file_name)

Load Shared Near Infrared Spectroscopy Format (SNIRF) and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load
- `n::Int64=0`: number of subject to extract in case of multi-subject file

# Returns

- `obj::NeuroAnalyzer.NEURO`

# Source

https://github.com/fNIRS/snirf/blob/v1.1/snirf_specification.md
"""
function import_snirf(file_name::String; n::Int64=0)

    file_name = "/home/eb/Documents/Data/EEG-testing-data/SNIRF/SfNIRS/snirf_homer3/1.0.3/nirx_15_3_recording.snirf"
    file_name = "test/files/fnirs-test-snirf.snirf"

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    nirs = FileIO.load(file_name)

    for (k, v) in nirs
        println("key: $k")
        println("value: $v")
    end

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
    subject_id = nirs["$n_id/metaDataTags/SubjectID"]
    recording_date = nirs["$n_id/metaDataTags/MeasurementDate"]
    recording_time = nirs["$n_id/metaDataTags/MeasurementTime"]
    length_unit = nirs["$n_id/metaDataTags/LengthUnit"]
    time_unit = nirs["$n_id/metaDataTags/TimeUnit"]
    frq_unit = nirs["$n_id/metaDataTags/FrequencyUnit"]

    # probes

    # List of wavelengths (in nm)
    wavelengths = nirs["$n_id/probe/wavelengths"]

    # List of emission wavelengths (in nm)
    k = "$n_id/probe/wavelengthsEmission"
    k in keys(nirs) && (wavelengths_emission = nirs[k])

    # Source 2-D positions in LengthUnit
    k = "$n_id/probe/sourcePos2D"
    k in keys(nirs) && (src_pos2d = nirs[k])

    # Source 3-D positions in LengthUnit
    k = "$n_id/probe/sourcePos3D"
    k in keys(nirs) && (src_pos2d = nirs[k])

    # Detector 2-D positions in LengthUnit
    k = "$n_id/probe/detectorPos2D"
    k in keys(nirs) && (detector_pos2d = nirs[k])

    # Detector 3-D positions in LengthUnit
    k = "$n_id/probe/detectorPos3D"
    k in keys(nirs) && (detector_pos2d = nirs[k])

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
    k in keys(nirs) && (coord_system = nirs[k])[1]

    # Description of coordinate system
    k = "$n_id/probe/coordinateSystemDescription"
    k in keys(nirs) && (coord_system = nirs[k])[1]

    # If source/detector index is within a module
    k = "$n_id/probe/useLocalIndex"
    k in keys(nirs) && (coord_system = nirs[k])[1]

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

    for data_idx in 1:data_n
        d_id = "data$data_idx"
        
        data = nirs["$n_id/$d_id/dataTimeSeries"]

        time_pts = nirs["$n_id/$d_id/time"]
        if length(time_pts) > 2
            sampling_rate = 1 / time_pts[2] - time_pts[1]
        else
            sampling_rate = 1 / time_pts[2]
            time_pts = collect(time_pts[1]:1/sampling_rate:time_pts[1]+size(data, 2)*time_pts[2])[1:(end - 1)]
        end

        channel_n = size(data, 1)

        for ch_idx in 1:ch_n
            # Source index for a given channel
            k = "$n_id/$d_id/measurementList$ch_idx/sourceIndex"
            k in keys(nirs) && (source_index = Int.(nirs[k])[1])

            # Detector index for a given channel
            k = "$n_id/$d_id/measurementList$ch_idx/detectorIndex"
            k in keys(nirs) && (detector_index = Int.(nirs[k])[1])

            # Wavelength index for a given channel
            k = "$n_id/$d_id/measurementList$ch_idx/wavelengthIndex"
            k in keys(nirs) && (wavelength_index = Int.(nirs[k])[1])

            # Actual wavelength for a given channel
            k = "$n_id/$d_id/measurementList$ch_idx/wavelengthActual"
            k in keys(nirs) && (wavelength_actual = nirs[k])[1]

            # Actual emission wavelength for a channel
            k = "$n_id/$d_id/measurementList$ch_idx/wavelengthEmissionActual"
            k in keys(nirs) && (wavelength_emission_actual = nirs[k])[1]

            # Data type for a given channel
            k = "$n_id/$d_id/measurementList$ch_idx/dataType"
            k in keys(nirs) && (data_type = Int.(nirs[k])[1])

            # SI unit for a given channel
            k = "$n_id/$d_id/measurementList$ch_idx/dataUnit"
            k in keys(nirs) && (data_unit = nirs[k])[1]

            # Data type name for a given channel
            k = "$n_id/$d_id/measurementList$ch_idx/dataTypeLabel"
            k in keys(nirs) && (data_type_label = nirs[k])[1]

            # Data type index for a given channel
            k = "$n_id/$d_id/measurementList$ch_idx/dataTypeIndex"
            k in keys(nirs) && (data_type_index = Int.(nirs[k])[1])

            # Source power for a given channel
            k = "$n_id/$d_id/measurementList$ch_idx/sourcePower"
            k in keys(nirs) && (source_power = nirs[k])[1]

            # Detector gain for a given channel
            k = "$n_id/$d_id/measurementList$ch_idx/detectorGain"
            k in keys(nirs) && (detector_gain = nirs[k])[1]

            # Index of the parent module (if modular)
            k = "$n_id/$d_id/measurementList$ch_idx/moduleIndex"
            k in keys(nirs) && (module_index = Int.(nirs[k])[1])

            # Index of the source's parent module
            k = "$n_id/$d_id/measurementList$ch_idx/sourceModuleIndex"
            k in keys(nirs) && (src_module_index = Int.(nirs[k])[1])

            # Index of the detector's parent module
            k = "$n_id/$d_id/measurementList$ch_idx/detectorModuleIndex"
            k in keys(nirs) && (detector_module_index = Int.(nirs[k])[1])

        end
    end

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

    for stim_idx in 1:stim_n
        s_id = "stim$stim_idx"

        # Name of the stimulus data
        k = "$n_id/$s_id/name"
        k in keys(nirs) && (stim_name = nirs[k])[1]

        # Data stream of the stimulus channel
        k = "$n_id/$s_id/data"
        k in keys(nirs) && (stim_data = nirs[k])

        # Data stream of the stimulus channel
        k = "$n_id/$s_id/dataLabels"
        k in keys(nirs) && (stim_labels = nirs[k])
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

    for aux_idx in 1:aux_n
        a_id = "aux$aux_idx"

        # Name of the auxiliary channel
        k = "$n_id/$a_id/name"
        k in keys(nirs) && (aux_name = nirs[k])[1]

        # Data acquired from the auxiliary channel
        k = "$n_id/$a_id/dataTimeSeries"
        k in keys(nirs) && (aux_data = nirs[k])

        # SI unit of the auxiliary channel
        k = "$n_id/$a_id/dataUnit"
        k in keys(nirs) && (aux_data = nirs[k])

        # Time (in TimeUnit) for auxiliary data 
        k = "$n_id/$a_id/time"
        k in keys(nirs) && (aux_time = nirs[k])

        # Time offset of auxiliary channel data
        k = "$n_id/$a_id/timeOffset"
        k in keys(nirs) && (aux_timeoffset = nirs[k])
    end

end
