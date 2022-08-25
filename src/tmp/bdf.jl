"""
    eeg_import_bdf(file_name; read_annotations, clean_labels)

Load BDF file and return and `NeuroAnalyzer.EEG` object.

# Arguments

- `file_name::String`: name of the file to load
- `clean_labels::Bool=true`: only keep channel names in channel labels

# Returns

- `eeg:EEG`

# Notes

- sampling_rate = n.samples / data.record.duration
- gain = (physical_maximum - physical_minimum) / (digital_maximum - digital_minimum)
- value = (value - digital_minimum ) * gain + physical_minimum

# Source

Kemp B, Värri A, Rosa AC, Nielsen KD, Gade J. A simple format for exchange of digitized polygraphic recordings. Electroencephalography and Clinical Neurophysiology. 1992 May;82(5):391–3. 
"""
function eeg_import_bdf(file_name::String; read_annotations::Bool=true, clean_labels::Bool=true)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    fid = open(file_name)

    eeg_filetype = ""

    header = zeros(UInt8, 256)
    readbytes!(fid, header, 256)
    header = String(Char.(header))

    version = Int(header[1])
    version == 255 && (eeg_filetype = "BDF")
    ftype = rstrip(header[3:9])
    (eeg_filetype !== "BDF" && ftype !== "BIOSEMI") && throw(ArgumentError("File is not a EDF file."))

    patient = rstrip(header[10:89])
    recording = rstrip(header[90:169])
    recording_date = header[170:177]
    recording_time = header[178:185]
    data_offset = parse(Int, strip(header[186:192]))
    resolution  = strip(header[193:236])
    data_records = parse(Int, strip(header[237:244]))
    data_records_duration  = parse(Float64, strip(header[245:252]))
    channel_n  = parse(Int, strip(header[253:256]))

    labels = Vector{String}(undef, channel_n)
    transducers = Vector{String}(undef, channel_n)
    physical_dimension = Vector{String}(undef, channel_n)
    physical_minimum = Vector{Float64}(undef, channel_n)
    physical_maximum = Vector{Float64}(undef, channel_n)
    digital_minimum = Vector{Float64}(undef, channel_n)
    digital_maximum = Vector{Float64}(undef, channel_n)
    prefiltering = Vector{String}(undef, channel_n)
    samples_per_datarecord = Vector{Int64}(undef, channel_n)

    header = zeros(UInt8, channel_n * 16)
    readbytes!(fid, header, channel_n * 16)
    header = String(Char.(header))
    for idx in 1:channel_n
        labels[idx] = strip(header[1 + ((idx - 1) * 16):(idx * 16)])
    end

    header = zeros(UInt8, channel_n * 80)
    readbytes!(fid, header, channel_n * 80)
    header = String(Char.(header))
    for idx in 1:channel_n
        transducers[idx] = strip(header[1 + ((idx - 1) * 80):(idx * 80)])
    end

    header = zeros(UInt8, channel_n * 8)
    readbytes!(fid, header, channel_n * 8)
    header = String(Char.(header))
    for idx in 1:channel_n
        physical_dimension[idx] = strip(header[1 + ((idx - 1) * 8):(idx * 8)])
    end

    header = zeros(UInt8, channel_n * 8)
    readbytes!(fid, header, channel_n * 8)
    header = String(Char.(header))
    for idx in 1:channel_n
        physical_minimum[idx] = parse(Float64, strip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = zeros(UInt8, channel_n * 8)
    readbytes!(fid, header, channel_n * 8)
    header = String(Char.(header))
    for idx in 1:channel_n
        physical_maximum[idx] = parse(Float64, strip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = zeros(UInt8, channel_n * 8)
    readbytes!(fid, header, channel_n * 8)
    header = String(Char.(header))
    for idx in 1:channel_n
        digital_minimum[idx] = parse(Float64, strip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = zeros(UInt8, channel_n * 8)
    readbytes!(fid, header, channel_n * 8)
    header = String(Char.(header))
    for idx in 1:channel_n
        digital_maximum[idx] = parse(Float64, strip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = zeros(UInt8, channel_n * 80)
    readbytes!(fid, header, channel_n * 80)
    header = String(Char.(header))
    for idx in 1:channel_n
        prefiltering[idx] = strip(header[1 + ((idx - 1) * 80):(idx * 80)])
    end

    header = zeros(UInt8, channel_n * 8)
    readbytes!(fid, header, channel_n * 8)
    header = String(Char.(header))
    for idx in 1:channel_n
        samples_per_datarecord[idx] = parse(Int, strip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    close(fid)

    sampling_rate = Vector{Float64}(undef, channel_n)
    gain = Vector{Float64}(undef, channel_n)
    for idx in 1:channel_n
        sampling_rate[idx] = samples_per_datarecord[idx] / data_records_duration
        gain[idx] = (physical_maximum[idx] - physical_minimum[idx]) / (digital_maximum[idx] - digital_minimum[idx])
    end

    channels_1020 = ["fp1", "fpz", "fp2", "af9", "af7", "af5", "af3", "af1", "afz", "af2", "af4", "af6", "af8", "af10", "f9", "f7", "f5", "f3", "f1", "fz", "f2", "f4", "f6", "f8", "f10", "ft9", "ft7", "fc5", "fc3", "fc1", "fcz", "fc2", "fc4", "fc6", "ft8", "ft10", "t9", "t7", "c5", "c3", "c1", "cz", "c2", "c4", "c6", "t8", "t10", "tp9", "tp7", "cp5", "cp3", "cp1", "cpz", "cp2", "cp4", "cp6", "tp8", "tp10", "p9", "p7", "p5", "p3", "p1", "pz", "p2", "p4", "p6", "p8", "p10", "po9", "po7", "po5", "po3", "po1", "poz", "po2", "po4", "po6", "po8", "po10", "o1", "oz", "o2", "o9", "o10", "t3", "t5", "t4", "t6"]
    channel_type = repeat(["unknown"], channel_n)
    for idx in 1:channel_n
        in(lowercase(labels[idx]), channels_1020) && (channel_type[idx] = "eeg")
        occursin("meg", lowercase(labels[idx])) && (channel_type[idx] = "meg")
        occursin("ecg", lowercase(labels[idx])) && (channel_type[idx] = "ecg")
        occursin("ekg", lowercase(labels[idx])) && (channel_type[idx] = "ecg")
        occursin("eog", lowercase(labels[idx])) && (channel_type[idx] = "eog")
        occursin("a1", lowercase(labels[idx])) && (channel_type[idx] = "ref")
        occursin("a2", lowercase(labels[idx])) && (channel_type[idx] = "ref")
        occursin("m1", lowercase(labels[idx])) && (channel_type[idx] = "ref")
        occursin("m2", lowercase(labels[idx])) && (channel_type[idx] = "ref")
        occursin("marker", lowercase(labels[idx])) && (channel_type[idx] = "marker")
        occursin("event", lowercase(labels[idx])) && (channel_type[idx] = "marker")
        occursin("annotation", lowercase(labels[idx])) && (channel_type[idx] = "annotation")
    end

    bdf_has_annotations = false
    if "annotation" in channel_type
        bdf_has_annotations = true
        annotation_channel = 0
        for channel_idx in 1:channel_n
            channel_type[channel_idx] == "annotation" && (annotation_channel = channel_idx)
        end
    end

    clean_labels == true && (labels = replace.(labels, "EEG " => ""))
    # clean_labels == true && (labels = replace.(labels, "MEG " => ""))
    # clean_labels == true && (labels = replace.(labels, "EOG " => ""))
    # clean_labels == true && (labels = replace.(labels, "ECG " => ""))
    clean_labels == true && (labels = replace.(labels, "BDF " => ""))
    clean_labels == true && (labels = replace.(labels, " " => ""))

    
    fid = open(file_name)
    header = zeros(UInt8, data_offset)
    readbytes!(fid, header, data_offset)
    eeg_signals = zeros(channel_n, samples_per_datarecord[1] * data_records, 1)
    annotations = repeat([""], data_records)

    for idx1 in 1:data_records
        for idx2 in 1:channel_n
            signal24 = zeros(UInt8, samples_per_datarecord[idx2] * 3)
            readbytes!(fid, signal24, samples_per_datarecord[idx2] * 3)
            if idx2 != annotation_channel
                signal64 = Vector{Float64}()
                for byte_idx in 1:3:length(signal24)
                    b1 = Int32(signal24[byte_idx]) << 8
                    b2 = Int32(signal24[byte_idx + 1]) << 16
                    b3 = -Int32(-signal24[byte_idx + 2]) << 24
                    push!(signal64, Float64(((b1 | b2 | b3) >> 8) * gain[idx2]))
                end
                eeg_signals[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal64
            else
                annotations[idx1] = String(Char.(signal24))
            end
        end
    end

    # last channel is always status channel
    close(fid)

    if bdf_has_annotations
        eeg_signals = @views eeg_signals[1:(end - 1), :, 1]
        channel_n -= 1
        annotations = _a2df(annotations)
    else
        annotations = DataFrame(id=[""], time=[""], annotation=[""])
    end

    eeg_duration_samples = size(eeg_signals, 2)
    eeg_duration_seconds = size(eeg_signals, 2) / sampling_rate[1]
    eeg_time = collect(0:(1 / sampling_rate[1]):eeg_duration_seconds)
    eeg_time = eeg_time[1:end - 1]
    sampling_rate = round.(Int64, sampling_rate)
    eeg_filesize_mb = round(filesize(file_name) / 1024^2, digits=2)

    eeg_header = Dict(:version => version,
                      :eeg_filename => file_name,
                      :eeg_filesize_mb => eeg_filesize_mb,
                      :eeg_filetype => eeg_filetype,
                      :patient => string(patient),
                      :recording => string(recording),
                      :recording_date => recording_date,
                      :recording_time => recording_time,
                      :data_records => data_records,
                      :data_records_duration => data_records_duration,
                      :channel_n => channel_n,
                      :channel_type => channel_type,
                      :reference => "",
                      :channel_locations => false,
                      :loc_theta => zeros(channel_n),
                      :loc_radius => zeros(channel_n),
                      :loc_x => zeros(channel_n),
                      :loc_y => zeros(channel_n),
                      :loc_z => zeros(channel_n),
                      :loc_radius_sph => zeros(channel_n),
                      :loc_theta_sph => zeros(channel_n),
                      :loc_phi_sph => zeros(channel_n),
                      :history => String[],
                      :components => Symbol[],
                      :eeg_duration_samples => eeg_duration_samples,
                      :eeg_duration_seconds => eeg_duration_seconds,
                      :epoch_n => 1,
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
                      :gain => gain,
                      :note => "")

    eeg_components = Vector{Any}()
    eeg_epochs_time = eeg_time

    eeg = EEG(eeg_header, eeg_time, eeg_epochs_time, eeg_signals, eeg_components)

    return eeg
end