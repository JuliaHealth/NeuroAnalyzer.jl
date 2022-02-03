"""
    eeg_load_edf(in_file, read_annotations=true, header_only=false, clean_labels=true)

Loads EDF/EDFPlus file and returns EEG object.

# Arguments

- `in_file::String` - name of the file to load
- `read_annotations::Bool` - read annotations from EDF+ file (currently not implemented)
- `header_only::Bool` - do not read signal data
- `clean_labels::Bool` - only keep channel names in channel labels

# Notes

sampling_rate = n.samples / data.record.duration

gain = (physical_maximum - physical_minimum) / (digital_maximum - digital_minimum)

value = (value - digital_minimum ) * gain + physical_minimum

# Source

Kemp B, Värri A, Rosa AC, Nielsen KD, Gade J. A simple format for exchange of digitized polygraphic recordings. Electroencephalography and Clinical Neurophysiology. 1992 May;82(5):391–3. 
"""
function eeg_import_edf(in_file::String, read_annotations::Bool=true, header_only::Bool=false, clean_labels::Bool=true)
    fid = open(in_file)

    eeg_filetype = ""

    header = zeros(UInt8, 256)
    readbytes!(fid, header, 256)
    header = String(Char.(header))

    version = parse(Int, rstrip(header[1:8]))
    version == 0 && (eeg_filetype = "EDF")
    eeg_filetype !== "EDF" && throw(ArgumentError("File is not a EDF file."))

    patient = rstrip(header[9:88])
    recording = rstrip(header[89:168])
    recording_date = header[169:176]
    recording_time = header[177:184]
    data_offset = parse(Int, rstrip(header[185:192]))
    reserved  = header[193:236]
    data_records = parse(Int, rstrip(header[237:244]))
    data_records_duration  = parse(Float64, rstrip(header[245:252]))
    channels_no  = parse(Int, rstrip(header[253:256]))

    labels = Vector{String}(undef, channels_no)
    transducers = Vector{String}(undef, channels_no)
    physical_dimension = Vector{String}(undef, channels_no)
    physical_minimum = Vector{Float64}(undef, channels_no)
    physical_maximum = Vector{Float64}(undef, channels_no)
    digital_minimum = Vector{Float64}(undef, channels_no)
    digital_maximum = Vector{Float64}(undef, channels_no)
    prefiltering = Vector{String}(undef, channels_no)
    samples_per_datarecord = Vector{Int64}(undef, channels_no)

    header = zeros(UInt8, channels_no * 16)
    readbytes!(fid, header, channels_no * 16)
    header = String(Char.(header))
    for idx in 1:channels_no
        labels[idx] = rstrip(header[1 + ((idx - 1) * 16):(idx * 16)])
    end

    header = zeros(UInt8, channels_no * 80)
    readbytes!(fid, header, channels_no * 80)
    header = String(Char.(header))
    for idx in 1:channels_no
        transducers[idx] = rstrip(header[1 + ((idx - 1) * 80):(idx * 80)])
    end

    header = zeros(UInt8, channels_no * 8)
    readbytes!(fid, header, channels_no * 8)
    header = String(Char.(header))
    for idx in 1:channels_no
        physical_dimension[idx] = rstrip(header[1 + ((idx - 1) * 8):(idx * 8)])
    end

    header = zeros(UInt8, channels_no * 8)
    readbytes!(fid, header, channels_no * 8)
    header = String(Char.(header))
    for idx in 1:channels_no
        physical_minimum[idx] = parse(Float64, rstrip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = zeros(UInt8, channels_no * 8)
    readbytes!(fid, header, channels_no * 8)
    header = String(Char.(header))
    for idx in 1:channels_no
        physical_maximum[idx] = parse(Float64, rstrip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = zeros(UInt8, channels_no * 8)
    readbytes!(fid, header, channels_no * 8)
    header = String(Char.(header))
    for idx in 1:channels_no
        digital_minimum[idx] = parse(Float64, rstrip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = zeros(UInt8, channels_no * 8)
    readbytes!(fid, header, channels_no * 8)
    header = String(Char.(header))
    for idx in 1:channels_no
        digital_maximum[idx] = parse(Float64, rstrip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = zeros(UInt8, channels_no * 80)
    readbytes!(fid, header, channels_no * 80)
    header = String(Char.(header))
    for idx in 1:channels_no
        prefiltering[idx] = rstrip(header[1 + ((idx - 1) * 80):(idx * 80)])
    end

    header = zeros(UInt8, channels_no * 8)
    readbytes!(fid, header, channels_no * 8)
    header = String(Char.(header))
    for idx in 1:channels_no
        samples_per_datarecord[idx] = parse(Int, rstrip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    close(fid)

    sampling_rate = Vector{Float64}(undef, channels_no)
    gain = Vector{Float64}(undef, channels_no)
    for idx in 1:channels_no
        sampling_rate[idx] = samples_per_datarecord[idx] / data_records_duration
        gain[idx] = (physical_maximum[idx] - physical_minimum[idx]) / (digital_maximum[idx] - digital_minimum[idx])
    end

    if clean_labels == true
        labels = replace.(labels, "EEG " => "")
        labels = replace.(labels, "ECG " => "")
    end

    if header_only == true
        eeg = EEG(eeg_object_header, eeg_signal_header)
        return eeg
    end

    fid = open(in_file)
    header = zeros(UInt8, data_offset)
    readbytes!(fid, header, data_offset)
    eeg_signals = zeros(channels_no, samples_per_datarecord[1] * data_records, 1)

    for idx1 in 1:data_records
        for idx2 in 1:channels_no
            signal = zeros(UInt8, samples_per_datarecord[idx2] * 2);
            readbytes!(fid, signal, samples_per_datarecord[idx2] * 2);
            signal = map(ltoh, reinterpret(Int16, signal));
            eeg_signals[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = @. (signal - digital_minimum[idx2]) * gain[idx2] + physical_minimum[idx2];
        end
    end

    close(fid)

    eeg_duration_samples = size(eeg_signals, 2)
    eeg_duration_seconds = size(eeg_signals, 2) / sampling_rate[1]
    eeg_time = collect(0:(1 / sampling_rate[1]):eeg_duration_seconds)
    eeg_time = eeg_time[1:end - 1]
    sampling_rate = round.(Int64, sampling_rate)
    eeg_filesize_mb = round(filesize(in_file) / 1024^2, digits=2)

    eeg_header = Dict(:version => version,
                      :eeg_filename => in_file,
                      :eeg_filesize_mb => eeg_filesize_mb,
                      :eeg_filetype => eeg_filetype,
                      :patient => patient,
                      :recording => recording,
                      :recording_date => recording_date,
                      :recording_time => recording_time,
                      :data_records => data_records,
                      :data_records_duration => data_records_duration,
                      :channels_no => channels_no,
                      :reference => "",
                      :channel_locations => false,
                      :xlocs => [],
                      :ylocs => [],
                      :history => Vector{String}(),
                      :eeg_duration_samples => eeg_duration_samples,
                      :eeg_duration_seconds => eeg_duration_seconds,
                      :epochs_no => 1,
                      :epoch_duration_samples => eeg_duration_samples,
                      :epoch_duration_seconds => eeg_duration_seconds,
                      :labels => labels,
                      :transducers => transducers,
                      :physical_dimension => physical_dimension,
                      :physical_minimum => physical_minimum,
                      :physical_maximum => physical_maximum,
                      :digital_minimum => digital_minimum,
                      :digital_maximum => digital_maximum,
                      :prefiltering => prefiltering,
                      :samples_per_datarecord => samples_per_datarecord,
                      :sampling_rate => sampling_rate,
                      :gain => gain)

    eeg = EEG(eeg_header, eeg_time, eeg_signals)

    return eeg
end

"""
    eeg_load_electrode_positions(eeg, in_file)

Loads electrode positions from 
- CED

# Arguments

- `eeg::EEG`
- `in_file::String`
"""
function eeg_load_electrode_positions(eeg::EEG, in_file)
    isfile(in_file) || throw(ArgumentError("File $in_file cannot be loaded."))
    sensors = CSV.read(in_file, delim="\t", DataFrame)
    loc_x = zeros(length(sensors[:, :radius]))
    loc_y = zeros(length(sensors[:, :theta]))
    for idx in 1:length(sensors[:, :theta])
        loc_y[idx], loc_x[idx] = pol2cart(pi / 180 * sensors[idx, :theta], sensors[idx, :radius])
    end
    length(loc_x) != eeg.eeg_header[:channels_no] && throw(ArgumentError("Number of channels and number of positions do not match."))

    # create new dataset
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), deepcopy(eeg.eeg_signals))
    eeg_new.eeg_header[:channel_locations] = true
    eeg_new.eeg_header[:xlocs] = loc_x
    eeg_new.eeg_header[:ylocs] = loc_y

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_load_sensor_positions(EEG, $in_file)")

    return eeg_new
end


"""
    eeg_save(eeg, file_name)

Saves the `eeg` to `file_name` file (HDF5-based).

# Arguments

- `eeg::EEG`
- `file_name::String` - file name
"""
function eeg_save(eeg::EEG, file_name::String)
    try
        save_object(file_name, eeg)
    catch error
        throw(ArgumentError("File $file_name cannot be saved."))
        return false
    end

    return true
end

"""
    eeg_load(file_name)

Loads the `eeg` from `file_name` file (HDF5-based).

# Arguments

- `file_name::String` - file name
"""
function eeg_load(file_name::String)
    try
        eeg = load_object(file_name)
    catch error
        throw(SystemError("File $file_name cannot be loaded."))
        return false
    end

    return eeg
end
