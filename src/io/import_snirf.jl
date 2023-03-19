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

    subject_id = nirs["$n_id/metaDataTags/SubjectID"]
    recording_date = nirs["$n_id/metaDataTags/MeasurementDate"]
    recording_time = nirs["$n_id/metaDataTags/MeasurementTime"]

    data = nirs["$n_id/data1/dataTimeSeries"]

    time_pts = nirs["$n_id/data1/time"]
    if length(time_pts) > 2
        sampling_rate = 1 / time_pts[2] - time_pts[1]
    else
        sampling_rate = 1 / time_pts[2]
        time_pts = collect(time_pts[1]:1/sampling_rate:time_pts[1]+size(data, 2)*time_pts[2])[1:(end - 1)]
    end

    channel_n = size(data, 1)

    for ch_idx in 1:ch_n
        k = "$n_id/data1/measurementList$ch_idx/dataType"
        println("key: $k")
        v = nirs[k]
        println("value: $v")
        k = "$n_id/data1/measurementList$ch_idx/dataTypeLabel"
        println("key: $k")
        v = nirs[k]
        println("value: $v")
        k = "$n_id/data1/measurementList$ch_idx/wavelengthIndex"
        println("key: $k")
        v = nirs[k]
        println("value: $v")
        k = "$n_id/data1/measurementList$ch_idx/sourceIndex"
        println("key: $k")
        v = nirs[k]
        println("value: $v")
        k = "$n_id/data1/measurementList$ch_idx/detectorGain"
        println("key: $k")
        v = nirs[k]
        println("value: $v")
        k = "$n_id/data1/measurementList$ch_idx/sourcePower"
        println("key: $k")
        v = nirs[k]
        println("value: $v")
    end

    # sources
    src_pos2d = nirs["nirs/probe/sourcePos2D"]
    src_pos3d = nirs["nirs/probe/sourcePos3D"]
    src_labels = nirs["nirs/probe/sourceLabels"]
    # probes
    wavelengths = nirs["nirs/probe/wavelengths"]
    probe_labels = nirs["nirs/probe/detectorLabels"]

end

