"""
    eeg_import_edf(file_name; clean_labels)

Load EDF file and return and `NeuroAnalyzer.EEG` object.

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
function eeg_import_edf(file_name::String; clean_labels::Bool=true)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    eeg_filetype = ""
    fid = open(file_name)
    header = zeros(UInt8, 256)
    readbytes!(fid, header, 256)
    header = String(Char.(header))

    eeg_filetype = parse(Int, strip(header[1:8]))
    eeg_filetype == 0 && (eeg_filetype = "EDF")
    eeg_filetype !== "EDF" && throw(ArgumentError("File $file_name is not a EDF file."))

    patient = strip(header[9:88])
    recording = strip(header[89:168])
    recording_date = header[169:176]
    recording_time = header[177:184]
    data_offset = parse(Int, strip(header[185:192]))
    reserved  = strip(header[193:236])
    reserved == "EDF+D" && throw(ArgumentError("EDF+D format (interrupted recordings) is not supported yet."))
    reserved == "EDF+C" && throw(ArgumentError("Use eeg_import_edfplus() for EDF+C file."))
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

    clean_labels == true && (labels = _clean_labels(labels))
    channel_type = _set_channel_types(labels)
    has_annotations = false
    eeg_annotations = DataFrame(id=[""], time=[""], description=[""])

    fid = open(file_name)
    header = zeros(UInt8, data_offset)
    readbytes!(fid, header, data_offset)
    eeg_signals = zeros(channel_n, samples_per_datarecord[1] * data_records, 1)

    for idx1 in 1:data_records
        for idx2 in 1:channel_n
            signal = zeros(UInt8, samples_per_datarecord[idx2] * 2);
            readbytes!(fid, signal, samples_per_datarecord[idx2] * 2);
            signal = map(ltoh, reinterpret(Int16, signal));
            if channel_type[idx2] == "markers"
                for idx3 in 1:length(signal)
                    if signal[idx3] == digital_minimum[idx2]
                        signal[idx3] = 0
                    else
                        signal[idx3] = 1
                    end
                end
                eeg_signals[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal
            elseif channel_type[idx2] == "events"
                eeg_signals[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal
            else
                if occursin("uV", physical_dimension[idx2]) 
                    eeg_signals[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal .* gain[idx2]
                elseif occursin("mV", physical_dimension[idx2])
                    eeg_signals[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal .* gain[idx2] ./ 1000
                else
                    eeg_signals[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal .* gain[idx2]
                end
            end
        end
    end

    close(fid)

    eeg_duration_samples = size(eeg_signals, 2)
    eeg_duration_seconds = size(eeg_signals, 2) / sampling_rate[1]
    eeg_time = collect(0:(1 / sampling_rate[1]):eeg_duration_seconds)
    eeg_time = eeg_time[1:end - 1]
    sampling_rate = round.(Int64, sampling_rate)
    eeg_filesize_mb = round(filesize(file_name) / 1024^2, digits=2)

    eeg_header = Dict(:eeg_filename => file_name,
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
                      :note => "",
                      :annotations => has_annotations)

    eeg_components = Vector{Any}()
    eeg_epochs_time = eeg_time

    eeg = EEG(eeg_header, eeg_time, eeg_epochs_time, eeg_signals, eeg_components, eeg_annotations)

    return eeg
end

"""
    eeg_import_ced(file_name)

Load electrode positions from CED file.

# Arguments

- `file_name::String`

# Returns

- `sensors::DataFrame`
"""
function eeg_import_ced(file_name::String)

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    splitext(file_name)[2] == ".ced" || throw(ArgumentError("Not a CED file."))
    sensors = CSV.read(file_name, delim="\t", DataFrame)

    colnames = lowercase.(names(sensors))
    DataFrames.rename!(sensors, Symbol.(colnames))

    labels = lstrip.(sensors[!, "labels"])

    x = zeros(length(labels))
    y = zeros(length(labels))
    z = zeros(length(labels))
    radius = zeros(length(labels))
    theta = zeros(length(labels))
    radius_sph = zeros(length(labels))
    theta_sph = zeros(length(labels))
    phi_sph = zeros(length(labels))

    "x" in colnames && (x = Float64.(sensors[!, "x"]))
    "y" in colnames && (y = Float64.(sensors[!, "y"]))
    "z" in colnames && (z = Float64.(sensors[!, "z"]))
    "theta" in colnames && (theta = Float64.(sensors[!, "theta"]))
    "radius" in colnames && (radius = Float64.(sensors[!, "radius"]))
    "sph_radius" in colnames && (radius_sph = Float64.(sensors[!, "sph_radius"]))
    "sph_theta" in colnames && (theta_sph = Float64.(sensors[!, "sph_theta"]))
    "sph_phi" in colnames && (phi_sph = Float64.(sensors[!, "sph_phi"]))

    sensors = DataFrame(:labels => labels, :loc_theta => theta, :loc_radius => radius, :loc_x => x, :loc_y => y, :loc_z => z, :loc_radius_sph => radius_sph, :loc_theta_sph => theta_sph, :loc_phi_sph => phi_sph)

    return sensors
end

"""
    eeg_import_locs(file_name)

Load electrode positions from LOCS file.

# Arguments

- `file_name::String`

# Returns

- `sensors::DataFrame`
"""
function eeg_import_locs(file_name::String)

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    splitext(file_name)[2] == ".locs" || throw(ArgumentError("Not a LOCS file."))
    sensors = CSV.read(file_name, header=false, delim="\t", DataFrame)

    DataFrames.rename!(sensors, [:number, :theta, :radius, :labels])
    labels = lstrip.(sensors[!, "labels"])

    x = zeros(length(labels))
    y = zeros(length(labels))
    z = zeros(length(labels))
    radius = zeros(length(labels))
    theta = zeros(length(labels))
    radius_sph = zeros(length(labels))
    theta_sph = zeros(length(labels))
    phi_sph = zeros(length(labels))

    theta = Float64.(sensors[!, "theta"])
    radius = Float64.(sensors[!, "radius"])

    sensors = DataFrame(:labels => labels, :loc_theta => theta, :loc_radius => radius, :loc_x => x, :loc_y => y, :loc_z => z, :loc_radius_sph => radius_sph, :loc_theta_sph => theta_sph, :loc_phi_sph => phi_sph)

    return sensors
end

"""
    eeg_import_elc(file_name)

Load electrode positions from ELC file.

# Arguments

- `file_name::String`

# Returns

- `sensors::DataFrame`
"""
function eeg_import_elc(file_name::String)

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    splitext(file_name)[2] == ".elc" || throw(ArgumentError("Not a ELC file."))
    f = open(file_name, "r")
    elc_file = readlines(f)
    close(f)
    locs_n = 0
    locs_l = 0
    for idx in 1:length(elc_file)
        if occursin("NumberPositions", elc_file[idx]) == true
            locs_n = parse(Int64, replace(elc_file[idx], "NumberPositions=" => ""))
            locs_l = idx + 2
        end
    end
    labels = repeat([""], locs_n)

    x = zeros(length(labels))
    y = zeros(length(labels))
    z = zeros(length(labels))
    theta = zeros(length(labels))
    radius = zeros(length(labels))
    radius_sph = zeros(length(labels))
    theta_sph = zeros(length(labels))
    phi_sph = zeros(length(labels))

    idx2 = 1
    for idx1 in locs_l:(locs_l + locs_n - 1)
        l = elc_file[idx1]
        l[1] == ' ' && (l = l[2:end])
        x[idx2], x[idx2], z[idx2] = parse.(Float64, split(l, ' '))
        idx2 += 1
    end
    idx2 = 1
    for idx1 in (locs_l + 1 + locs_n):(locs_l + (2 * locs_n))
        labels[idx2] = elc_file[idx1]
        idx2 += 1
    end

    sensors = DataFrame(:labels => labels, :loc_theta => theta, :loc_radius => radius, :loc_x => x, :loc_y => y, :loc_z => z, :loc_radius_sph => radius_sph, :loc_theta_sph => theta_sph, :loc_phi_sph => phi_sph)

    return sensors
end

"""
    eeg_import_tsv(file_name)

Load electrode positions from TSV file.

# Arguments

- `file_name::String`

# Returns

- `sensors::DataFrame`
"""
function eeg_import_tsv(file_name::String)

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    splitext(file_name)[2] == ".tsv" || throw(ArgumentError("Not a TSV file."))
    sensors = CSV.read(file_name, header=true, delim="\t", DataFrame)

    colnames = lowercase.(names(sensors))
    DataFrames.rename!(sensors, Symbol.(colnames))

    "labels" in colnames && (labels = lstrip.(sensors[!, "labels"]))
    "label" in colnames && (labels = lstrip.(sensors[!, "label"]))
    "site" in colnames && (labels = lstrip.(sensors[!, "site"]))

    x = zeros(length(labels))
    y = zeros(length(labels))
    z = zeros(length(labels))
    radius = zeros(length(labels))
    theta = zeros(length(labels))
    radius_sph = zeros(length(labels))
    theta_sph = zeros(length(labels))
    phi_sph = zeros(length(labels))
    
    "x" in colnames && (x = Float64.(sensors[!, "x"]))
    "y" in colnames && (y = Float64.(sensors[!, "y"]))
    "z" in colnames && (z = Float64.(sensors[!, "z"]))
    "theta" in colnames && (theta = Float64.(sensors[!, "theta"]))
    "radius" in colnames && (radius = Float64.(sensors[!, "radius"]))
    "radius" in colnames && (radius_sph = Float64.(sensors[!, "radius"]))
    "radius_sph" in colnames && (radius_sph = sensors[!, "radius_sph"])
    "theta" in colnames && (theta_sph = Float64.(sensors[!, "theta"]))
    "theta_sph" in colnames && (theta_sph = sensors[!, "theta_sph"])
    "phi" in colnames && (phi_sph = sensors[!, "phi"])
    "phi_sph" in colnames && (phi_sph = sensors[!, "phi_sph"])

    sensors = DataFrame(:labels => labels, :loc_theta => theta, :loc_radius => radius, :loc_x => x, :loc_y => y, :loc_z => z, :loc_radius_sph => radius_sph, :loc_theta_sph => theta_sph, :loc_phi_sph => phi_sph)

    return sensors
end

"""
    eeg_import_sfp(file_name)

Load electrode positions from SFP file.

# Arguments

- `file_name::String`

# Returns

- `sensors::DataFrame`
"""
function eeg_import_sfp(file_name::String)

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    splitext(file_name)[2] == ".sfp" || throw(ArgumentError("Not a SFP file."))
    sensors = CSV.read(file_name, header=false, delim="\t", DataFrame)

    DataFrames.rename!(sensors, [:label, :x, :y, :z])

    labels = lstrip.(sensors[!, "label"])

    x = zeros(length(labels))
    y = zeros(length(labels))
    z = zeros(length(labels))
    radius = zeros(length(labels))
    theta = zeros(length(labels))
    radius_sph = zeros(length(labels))
    theta_sph = zeros(length(labels))
    phi_sph = zeros(length(labels))
    
    x = Float64.(sensors[!, "x"])
    y = Float64.(sensors[!, "y"])
    z = Float64.(sensors[!, "z"])

    sensors = DataFrame(:labels => labels, :loc_theta => theta, :loc_radius => radius, :loc_x => x, :loc_y => y, :loc_z => z, :loc_radius_sph => radius_sph, :loc_theta_sph => theta_sph, :loc_phi_sph => phi_sph)

    return sensors
end

"""
    eeg_load_electrodes(eeg; file_name)

Load electrode positions from `file_name` and return `NeuroAnalyzer.EEG` object with metadata: `:channel_locations`, `:loc_theta`, `:loc_radius`, `:loc_x`, `:loc_x`, `:loc_y`, `:loc_radius_sph`, `:loc_theta_sph`, `:loc_phi_sph`. 

Accepted formats:
- CED
- LOCS
- ELC
- TSV
- SFP

Electrode locations:
- loc_theta       planar polar angle
- loc_radius      planar polar radius
- loc_x           spherical Cartesian x
- loc_y           spherical Cartesian y
- loc_z           spherical Cartesian z
- loc_radius_sph  spherical radius
- loc_theta_sph   spherical horizontal angle
- loc_phi_sph     spherical azimuth angle

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `file_name::String`

# Returns

- `eeg:EEG`
"""
function eeg_load_electrodes(eeg::NeuroAnalyzer.EEG; file_name::String)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    length(eeg.eeg_header[:labels]) > 0 || throw(ArgumentError("EEG does not contain labels, use eeg_add_labels() first."))

    if splitext(file_name)[2] == ".ced"
        sensors = eeg_import_ced(file_name)
    elseif splitext(file_name)[2] == ".elc"
        sensors = eeg_import_elc(file_name)
    elseif splitext(file_name)[2] == ".locs"
        sensors = eeg_import_locs(file_name)
    elseif splitext(file_name)[2] == ".tsv"
        sensors = eeg_import_tsv(file_name)
    elseif splitext(file_name)[2] == ".sfp"
        sensors = eeg_import_sfp(file_name)
    else
        throw(ArgumentError("Unknown file format."))
    end

    f_labels = lowercase.(sensors[:, :labels])

    loc_theta = float.(sensors[:, :loc_theta])
    loc_radius = float.(sensors[:, :loc_radius])

    loc_radius_sph = float.(sensors[:, :loc_radius_sph])
    loc_theta_sph = float.(sensors[:, :loc_theta_sph])
    loc_phi_sph = float.(sensors[:, :loc_phi_sph])

    loc_x = float.(sensors[:, :loc_x])
    loc_y = float.(sensors[:, :loc_y])
    loc_z = float.(sensors[:, :loc_z])

    e_labels = lowercase.(eeg.eeg_header[:labels])
    no_match = setdiff(e_labels, f_labels)
    length(no_match) > 0 && throw(ArgumentError("Labels: $(uppercase.(no_match)) not found in $file_name."))

    labels_idx = zeros(Int64, length(e_labels))
    for idx1 in 1:length(e_labels)
        for idx2 in 1:length(f_labels)
            e_labels[idx1] == f_labels[idx2] && (labels_idx[idx1] = idx2)
        end
    end
    
    # create new dataset
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_header[:channel_locations] = true
    eeg_new.eeg_header[:loc_theta] = loc_theta[labels_idx]
    eeg_new.eeg_header[:loc_radius] = loc_radius[labels_idx]
    eeg_new.eeg_header[:loc_x] = loc_x[labels_idx]
    eeg_new.eeg_header[:loc_y] = loc_y[labels_idx]
    eeg_new.eeg_header[:loc_z] = loc_z[labels_idx]
    eeg_new.eeg_header[:loc_radius_sph] = loc_radius_sph[labels_idx]
    eeg_new.eeg_header[:loc_theta_sph] = loc_theta_sph[labels_idx]
    eeg_new.eeg_header[:loc_phi_sph] = loc_phi_sph[labels_idx]

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_load_electrodes(EEG, file_name=$file_name)")

    return eeg_new
end

"""
    eeg_load_electrodes!(eeg; file_name)

Load electrode positions from `file_name` and return `NeuroAnalyzer.EEG` object with metadata: `:channel_locations`, `:loc_theta`, `:loc_radius`, `:loc_x`, `:loc_x`, `:loc_y`, `:loc_radius_sph`, `:loc_theta_sph`, `:loc_phi_sph`. 

Accepted formats:
- CED
- LOCS
- ELC
- TSV
- SFP

Electrode locations:
- loc_theta       planar polar angle
- loc_radius      planar polar radius
- loc_x           spherical Cartesian x
- loc_y           spherical Cartesian y
- loc_z           spherical Cartesian z
- loc_radius_sph  spherical radius
- loc_theta_sph   spherical horizontal angle
- loc_phi_sph     spherical azimuth angle

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `file_name::String`
"""
function eeg_load_electrodes!(eeg::NeuroAnalyzer.EEG; file_name::String)
    
    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    length(eeg.eeg_header[:labels]) > 0 || throw(ArgumentError("EEG does not contain labels, use eeg_add_labels() first."))

    if splitext(file_name)[2] == ".ced"
        sensors = eeg_import_ced(file_name)
    elseif splitext(file_name)[2] == ".elc"
        sensors = eeg_import_elc(file_name)
    elseif splitext(file_name)[2] == ".locs"
        sensors = eeg_import_locs(file_name)
    elseif splitext(file_name)[2] == ".tsv"
        sensors = eeg_import_tsv(file_name)
    elseif splitext(file_name)[2] == ".sfp"
        sensors = eeg_import_sfp(file_name)
    else
        throw(ArgumentError("Unknown file format."))
    end

    f_labels = lowercase.(sensors[:, :labels])

    loc_theta = float.(sensors[:, :loc_theta])
    loc_radius = float.(sensors[:, :loc_radius])

    loc_radius_sph = float.(sensors[:, :loc_radius_sph])
    loc_theta_sph = float.(sensors[:, :loc_theta_sph])
    loc_phi_sph = float.(sensors[:, :loc_phi_sph])

    loc_x = float.(sensors[:, :loc_x])
    loc_y = float.(sensors[:, :loc_y])
    loc_z = float.(sensors[:, :loc_z])

    e_labels = lowercase.(eeg.eeg_header[:labels])
    no_match = setdiff(e_labels, f_labels)
    length(no_match) > 0 && throw(ArgumentError("Labels: $(uppercase.(no_match)) does not found in $file_name."))

    labels_idx = zeros(Int64, length(e_labels))
    for idx1 in 1:length(e_labels)
        for idx2 in 1:length(f_labels)
            e_labels[idx1] == f_labels[idx2] && (labels_idx[idx1] = idx2)
        end
    end
    
    # create new dataset
    eeg.eeg_header[:channel_locations] = true
    eeg.eeg_header[:loc_theta] = loc_theta[labels_idx]
    eeg.eeg_header[:loc_radius] = loc_radius[labels_idx]
    eeg.eeg_header[:loc_x] = loc_x[labels_idx]
    eeg.eeg_header[:loc_y] = loc_y[labels_idx]
    eeg.eeg_header[:loc_z] = loc_z[labels_idx]
    eeg.eeg_header[:loc_radius_sph] = loc_radius_sph[labels_idx]
    eeg.eeg_header[:loc_theta_sph] = loc_theta_sph[labels_idx]
    eeg.eeg_header[:loc_phi_sph] = loc_phi_sph[labels_idx]

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_load_electrodes!(EEG, $file_name)")

 end

"""
    eeg_save(eeg; file_name, overwrite)

Save the `eeg` to `file_name` file (HDF5-based).

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `file_name::String`: file name
- `overwrite::Bool=false`

# Returns

- `success::Bool`
"""
function eeg_save(eeg::NeuroAnalyzer.EEG; file_name::String, overwrite::Bool=false)

    (isfile(file_name) && overwrite == false) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))

    eeg.eeg_header[:eeg_filename] = file_name

    save_object("/tmp/$file_name", eeg)
    eeg.eeg_header[:eeg_filesize_mb] = round(filesize("/tmp/$file_name") / 1024, digits=2)
    rm("/tmp/$file_name")

    save_object(file_name, eeg)
end

"""
    eeg_load(file_name)

Load the `eeg` from `file_name` file (HDF5-based).

# Arguments

- `file_name::String`: file name

# Returns

- `eeg::NeuroAnalyzer.EEG`
"""
function eeg_load(file_name::String)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    eeg = load_object(file_name)

    return eeg
end

"""
    eeg_export_csv(eeg; file_name, header, components, annotations, overwrite)

Export EEG data as CSV.

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `file_name::String`
- `header::Bool=false`: export header
- `components::Bool=false`: export components
- `annotations::Bool=false`: export annotations
- `overwrite::Bool=false`

# Returns

- `success::Bool`
"""
function eeg_export_csv(eeg::NeuroAnalyzer.EEG; file_name::String, header::Bool=false, components::Bool=false, annotations::Bool=true, overwrite::Bool=false)

    (isfile(file_name) && overwrite == false) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))
    eeg.eeg_header[:components] == [""] && throw(ArgumentError("EEG does not contain components."))

    # DATA
    # unsplit epochs
    s_merged = reshape(eeg.eeg_signals,
                       size(eeg.eeg_signals, 1),
                       size(eeg.eeg_signals, 2) * size(eeg.eeg_signals, 3))
    s = s_merged[:, :, 1]'
    s = hcat(eeg.eeg_time, s)
    l = vcat("time", eeg_labels(eeg))
    CSV.write(file_name, DataFrame(s, l))

    # HEADER
    if header
        file_name = replace(file_name, ".csv" => "_header.csv")
        (isfile(file_name) && overwrite == false) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))
        f = open(file_name, "w")
        for (key, value) in eeg.eeg_header
            println(f, key, ": ", value)
        end
        close(f)
    end

    # COMPONENTS
    if components
        file_name = replace(file_name, ".csv" => "_components.csv")
        (isfile(file_name) && overwrite == false) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))
        f = open(file_name, "w")
        for idx in 1:length(eeg.eeg_header[:components])
            println(f, "component: $(eeg.eeg_header[:components][idx])")
            println(f, eeg.eeg_components[idx])
            println(f, "---")
        end
        close(f)
    end

    # ANNOTATIONS
    if annotations
        file_name = replace(file_name, ".csv" => "_annotations.csv")
        (isfile(file_name) && overwrite == false) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))
        CSV.write(file_name, eeg.eeg_annotations)
    end
end

"""
    eeg_save_electrodes(eeg; file_name, overwrite)

Export EEG channel locations data, format is based on `file_name` extension (.ced, .locs or .tsv)

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `file_name::String`
- `overwrite::Bool=false`

# Returns

- `success::Bool`
"""
function eeg_save_electrodes(eeg::NeuroAnalyzer.EEG; file_name::String, overwrite::Bool=false)

    (isfile(file_name) && overwrite == false) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before exporting."))

    labels = eeg_labels(eeg)
    channels = collect(1:eeg_channel_n(eeg))
    theta = eeg.eeg_header[:loc_theta]
    radius = eeg.eeg_header[:loc_radius]
    x = eeg.eeg_header[:loc_x]
    y = eeg.eeg_header[:loc_y]
    z = eeg.eeg_header[:loc_z]
    radius_sph = eeg.eeg_header[:loc_radius_sph]
    theta_sph = eeg.eeg_header[:loc_theta_sph]
    phi_sph = eeg.eeg_header[:loc_phi_sph]

    if splitext(file_name)[2] == ".ced"
        df = DataFrame(Number=channels, labels=labels, theta=theta, radius=radius, X=x, Y=y, Z=z, sph_theta=theta_sph, sph_phi=phi_sph, sph_radius=radius_sph)
        CSV.write(file_name, df, delim="\t", header=true)
    elseif splitext(file_name)[2] == ".locs"
        df = DataFrame(Number=channels, theta=theta, radius=radius, labels=labels)
        CSV.write(file_name, df, delim="\t", header=false)
    elseif splitext(file_name)[2] == ".tsv"
        df = DataFrame(labels=labels, x=x, y=y, z=z, theta=theta, radius=radius, radius_sph=radius_sph, theta_sph=theta_sph, phi_sph=phi_sph)
        CSV.write(file_name, df, delim="\t", header=true)
    else
        throw(ArgumentError("file_name format must be .ced, .locs or .tsv."))
    end
end

"""
    eeg_save_electrodes(locs; file_name, overwrite)

Export channel locations, format is based on `file_name` extension (.ced, .locs, .tsv)

# Arguments

- `locs::DataFrame`
- `file_name::String`
- `overwrite::Bool=false`

# Returns

- `success::Bool`
"""
function eeg_save_electrodes(locs::DataFrame; file_name::String, overwrite::Bool=false)

    (isfile(file_name) && overwrite == false) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))

    labels = locs[!, :labels]
    channels = collect(1:length(labels))
    theta = locs[!, :loc_theta]
    radius = locs[!, :loc_radius]
    x = locs[!, :loc_x]
    y = locs[!, :loc_y]
    z = locs[!, :loc_z]
    radius_sph = locs[!, :loc_radius_sph]
    theta_sph = locs[!, :loc_theta_sph]
    phi_sph = locs[!, :loc_phi_sph]

    if splitext(file_name)[2] == ".ced"
        df = DataFrame(Number=channels, labels=labels, theta=theta, radius=radius, X=x, Y=y, Z=z, sph_theta=theta_sph, sph_phi=phi_sph, sph_radius=radius_sph, head=true)
        CSV.write(file_name, df, delim="\t")
    elseif splitext(file_name)[2] == ".locs"
        df = DataFrame(Number=channels, theta=theta, radius=radius, labels=labels)
        CSV.write(file_name, df, delim="\t", header=false)
    elseif splitext(file_name)[2] == ".tsv"
        df = DataFrame(labels=labels, x=x, y=y, z=z, theta=theta, radius=radius, radius_sph=radius_sph, theta_sph=theta_sph, phi_sph=phi_sph)
        CSV.write(file_name, df, delim="\t", header=true)
    else
        throw(ArgumentError("file_name format must be .ced, .locs or .tsv."))
    end
end

"""
    eeg_add_electrodes(eeg; locs)

Add electrode positions from `locs`. 

Electrode locations:
- loc_theta       planar polar angle
- loc_radius      planar polar radius
- loc_x           spherical Cartesian x
- loc_y           spherical Cartesian y
- loc_z           spherical Cartesian z
- loc_radius_sph  spherical radius
- loc_theta_sph   spherical horizontal angle
- loc_phi_sph     spherical azimuth angle

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `locs::DataFrame`

# Returns

- `eeg:EEG`
"""
function eeg_add_electrodes(eeg::NeuroAnalyzer.EEG; locs::DataFrame)

    f_labels = lowercase.(locs[:, :labels])

    loc_theta = float.(locs[:, :loc_theta])
    loc_radius = float.(locs[:, :loc_radius])

    loc_radius_sph = float.(locs[:, :loc_radius_sph])
    loc_theta_sph = float.(locs[:, :loc_theta_sph])
    loc_phi_sph = float.(locs[:, :loc_phi_sph])

    loc_x = float.(locs[:, :loc_x])
    loc_y = float.(locs[:, :loc_y])
    loc_z = float.(locs[:, :loc_z])

    e_labels = lowercase.(eeg.eeg_header[:labels])
    no_match = setdiff(e_labels, f_labels)
    length(no_match) > 0 && throw(ArgumentError("Labels: $(uppercase.(no_match)) not found in locs object."))

    labels_idx = zeros(Int64, length(e_labels))
    for idx1 in 1:length(e_labels)
        for idx2 in 1:length(f_labels)
            e_labels[idx1] == f_labels[idx2] && (labels_idx[idx1] = idx2)
        end
    end
    
    # create new dataset
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_header[:channel_locations] = true
    eeg_new.eeg_header[:loc_theta] = loc_theta[labels_idx]
    eeg_new.eeg_header[:loc_radius] = loc_radius[labels_idx]
    eeg_new.eeg_header[:loc_x] = loc_x[labels_idx]
    eeg_new.eeg_header[:loc_y] = loc_y[labels_idx]
    eeg_new.eeg_header[:loc_z] = loc_z[labels_idx]
    eeg_new.eeg_header[:loc_radius_sph] = loc_radius_sph[labels_idx]
    eeg_new.eeg_header[:loc_theta_sph] = loc_theta_sph[labels_idx]
    eeg_new.eeg_header[:loc_phi_sph] = loc_phi_sph[labels_idx]

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_add_electrodes(EEG, locs)")

    return eeg_new
end

"""
    eeg_add_electrodes!(eeg; locs)

Load electrode positions from `locs` and return `NeuroAnalyzer.EEG` object with metadata: `:channel_locations`, `:loc_theta`, `:loc_radius`, `:loc_x`, `:loc_x`, `:loc_y`, `:loc_radius_sph`, `:loc_theta_sph`, `:loc_phi_sph`. 

Electrode locations:
- loc_theta       planar polar angle
- loc_radius      planar polar radius
- loc_x           spherical Cartesian x
- loc_y           spherical Cartesian y
- loc_z           spherical Cartesian z
- loc_radius_sph  spherical radius
- loc_theta_sph   spherical horizontal angle
- loc_phi_sph     spherical azimuth angle

# Arguments

- `eeg::NeuroAnalyzer.EEG`
- `locs::DataFrame`
"""
function eeg_add_electrodes!(eeg::NeuroAnalyzer.EEG; locs::DataFrame)
    
    f_labels = lowercase.(locs[:, :labels])

    loc_theta = float.(locs[:, :loc_theta])
    loc_radius = float.(locs[:, :loc_radius])

    loc_radius_sph = float.(locs[:, :loc_radius_sph])
    loc_theta_sph = float.(locs[:, :loc_theta_sph])
    loc_phi_sph = float.(locs[:, :loc_phi_sph])

    loc_x = float.(locs[:, :loc_x])
    loc_y = float.(locs[:, :loc_y])
    loc_z = float.(locs[:, :loc_z])

    e_labels = lowercase.(eeg.eeg_header[:labels])
    no_match = setdiff(e_labels, f_labels)
    length(no_match) > 0 && throw(ArgumentError("Labels: $(uppercase.(no_match)) not found in locs object."))

    labels_idx = zeros(Int64, length(e_labels))
    for idx1 in 1:length(e_labels)
        for idx2 in 1:length(f_labels)
            e_labels[idx1] == f_labels[idx2] && (labels_idx[idx1] = idx2)
        end
    end
    
    # create new dataset
    eeg.eeg_header[:channel_locations] = true
    eeg.eeg_header[:loc_theta] = loc_theta[labels_idx]
    eeg.eeg_header[:loc_radius] = loc_radius[labels_idx]
    eeg.eeg_header[:loc_x] = loc_x[labels_idx]
    eeg.eeg_header[:loc_y] = loc_y[labels_idx]
    eeg.eeg_header[:loc_z] = loc_z[labels_idx]
    eeg.eeg_header[:loc_radius_sph] = loc_radius_sph[labels_idx]
    eeg.eeg_header[:loc_theta_sph] = loc_theta_sph[labels_idx]
    eeg.eeg_header[:loc_phi_sph] = loc_phi_sph[labels_idx]

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_add_electrodes!(EEG, locs)")

 end

"""
    eeg_import_edfplus(file_name; clean_labels)

Load EDF/EDFPlus file and return and `NeuroAnalyzer.EEG` object.

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

1. Kemp B, Olivan J. European data format ‘plus’ (EDF+), an EDF alike standard format for the exchange of physiological data. Clinical Neurophysiology 2003;114:1755–61.
2. https://www.edfplus.info/specs/
"""
function eeg_import_edfplus(file_name::String; clean_labels::Bool=true)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    splitext(file_name)[2] == ".edf" || throw(ArgumentError("File $file_name is not a EDF file."))

    fid = open(file_name)

    eeg_filetype = ""

    header = zeros(UInt8, 256)
    readbytes!(fid, header, 256)
    header = String(Char.(header))

    version = parse(Int, strip(header[1:8]))
    version == 0 && (eeg_filetype = "EDF")
    eeg_filetype !== "EDF" && throw(ArgumentError("File is not a EDF file."))

    patient = strip(header[9:88])
    recording = strip(header[89:168])
    recording_date = header[169:176]
    recording_time = header[177:184]
    data_offset = parse(Int, strip(header[185:192]))
    reserved = strip(header[193:236])
    reserved == "EDF+D" && throw(ArgumentError("EDF+D format (interrupted recordings) is not supported yet."))
    (reserved == "EDF+C" && eeg_filetype == "EDF") && (eeg_filetype = "EDF+")
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

    clean_labels == true && (labels = _clean_labels(labels))
    channel_type = _set_channel_types(labels)
    has_annotations, annotations_channel = _has_annotations(channel_type)

    fid = open(file_name)
    header = zeros(UInt8, data_offset)
    readbytes!(fid, header, data_offset)
    eeg_signals = zeros(channel_n, samples_per_datarecord[1] * data_records, 1)
    annotations = repeat([""], data_records)
    for idx1 in 1:data_records
        for idx2 in 1:channel_n
            signal = zeros(UInt8, samples_per_datarecord[idx2] * 2)
            readbytes!(fid, signal, samples_per_datarecord[idx2] * 2)
            if idx2 != annotations_channel
                signal = map(ltoh, reinterpret(Int16, signal))
                # eeg_signals[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = @. (signal - digital_minimum[idx2]) * gain[idx2] + physical_minimum[idx2]
                if channel_type[idx2] == "markers"
                    for idx3 in 1:length(signal)
                        if signal[idx3] == digital_minimum[idx2]
                            signal[idx3] = 0
                        else
                            signal[idx3] = 1
                        end
                    end
                    eeg_signals[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal
                elseif channel_type[idx2] == "events"
                    eeg_signals[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal
                else
                    if occursin("uV", physical_dimension[idx2]) 
                        eeg_signals[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal .* gain[idx2]
                    elseif occursin("mV", physical_dimension[idx2])
                        eeg_signals[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal .* gain[idx2] ./ 1000
                    else
                        eeg_signals[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal .* gain[idx2]
                    end
                end
            else
                annotations[idx1] = String(Char.(signal))
            end
        end
    end
    close(fid)

    if has_annotations
        eeg_signals = eeg_signals[1:(end - 1), :, :]
        deleteat!(labels, channel_n)
        deleteat!(channel_type, channel_n)
        deleteat!(transducers, channel_n)
        deleteat!(physical_dimension, channel_n)
        deleteat!(physical_minimum, channel_n)
        deleteat!(physical_maximum, channel_n)
        deleteat!(digital_minimum, channel_n)
        deleteat!(digital_maximum, channel_n)
        deleteat!(prefiltering, channel_n)
        deleteat!(samples_per_datarecord, channel_n)
        deleteat!(gain, channel_n)
        deleteat!(sampling_rate, channel_n)
        channel_n -= 1
        eeg_annotations = _a2df(annotations)
    else
        eeg_annotations = DataFrame(id=[""], time=[""], description=[""])
    end

    eeg_duration_samples = size(eeg_signals, 2)
    eeg_duration_seconds = size(eeg_signals, 2) / sampling_rate[1]
    eeg_time = collect(0:(1 / sampling_rate[1]):eeg_duration_seconds)
    eeg_time = eeg_time[1:end - 1]
    sampling_rate = round.(Int64, sampling_rate)
    eeg_filesize_mb = round(filesize(file_name) / 1024^2, digits=2)

    eeg_header = Dict(:eeg_filename => file_name,
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
                      :note => "",
                      :annotations => has_annotations)

    eeg_components = Vector{Any}()
    eeg_epochs_time = eeg_time

    eeg = EEG(eeg_header, eeg_time, eeg_epochs_time, eeg_signals, eeg_components, eeg_annotations)

    return eeg
end

"""
    eeg_import_bdf(file_name; clean_labels)

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

https://www.biosemi.com/faq/file_format.htm
"""
function eeg_import_bdf(file_name::String; clean_labels::Bool=true)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    splitext(file_name)[2] == ".bdf" || throw(ArgumentError("File $file_name is not a BDF file."))

    fid = open(file_name)

    eeg_filetype = ""

    header = zeros(UInt8, 256)
    readbytes!(fid, header, 256)
    header = String(Char.(header))

    version = Int(header[1])
    version == 255 && (eeg_filetype = "BDF")
    ftype = strip(header[3:9])
    (eeg_filetype !== "BDF" && ftype !== "BIOSEMI") && throw(ArgumentError("File is not a BDF file."))

    patient = strip(header[10:89])
    recording = strip(header[90:169])
    recording_date = header[170:177]
    recording_time = header[178:185]
    data_offset = parse(Int, strip(header[186:192]))
    reserved  = strip(header[193:236])
    reserved == "BDF+D" && throw(ArgumentError("BDF+D format (interrupted recordings) is not supported yet."))
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

    clean_labels == true && (labels = _clean_labels(labels))
    channel_type = _set_channel_types(labels)
    has_annotations, annotations_channel = _has_annotations(channel_type)

    fid = open(file_name)
    header = zeros(UInt8, data_offset)
    readbytes!(fid, header, data_offset)
    eeg_signals = zeros(channel_n, samples_per_datarecord[1] * data_records, 1)
    annotations = repeat([""], data_records)
    for idx1 in 1:data_records
        for idx2 in 1:channel_n
            signal24 = zeros(UInt8, samples_per_datarecord[idx2] * 3)
            readbytes!(fid, signal24, samples_per_datarecord[idx2] * 3)
            if idx2 != annotations_channel
                signal = Vector{Float64}()
                for byte_idx in 1:3:length(signal24)
                    b1 = Int32(signal24[byte_idx]) << 8
                    b2 = Int32(signal24[byte_idx + 1]) << 16
                    b3 = -Int32(-signal24[byte_idx + 2]) << 24
                    push!(signal, Float64(((b1 | b2 | b3) >> 8) * gain[idx2]))
                end
                if channel_type[idx2] == "markers"
                    for idx3 in 1:length(signal)
                        if signal[idx3] == digital_minimum[idx2]
                            signal[idx3] = 0
                        else
                            signal[idx3] = 1
                        end
                    end
                    eeg_signals[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal
                elseif channel_type[idx2] == "events"
                    eeg_signals[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal
                else
                    if occursin("uV", physical_dimension[idx2]) 
                        eeg_signals[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal .* gain[idx2]
                    elseif occursin("mV", physical_dimension[idx2])
                        eeg_signals[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal .* gain[idx2] ./ 1000
                    else
                        eeg_signals[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2]), 1] = signal .* gain[idx2]
                    end
                end
            else
                annotations[idx1] = String(Char.(signal24))
            end
        end
    end
    close(fid)

    if has_annotations
        eeg_signals = eeg_signals[1:(end - 1), :, :]
        deleteat!(labels, channel_n)
        deleteat!(transducers, channel_n)
        deleteat!(physical_dimension, channel_n)
        deleteat!(physical_minimum, channel_n)
        deleteat!(physical_maximum, channel_n)
        deleteat!(digital_minimum, channel_n)
        deleteat!(digital_maximum, channel_n)
        deleteat!(prefiltering, channel_n)
        deleteat!(samples_per_datarecord, channel_n)
        deleteat!(gain, channel_n)
        deleteat!(sampling_rate, channel_n)
        channel_n -= 1
        eeg_annotations = _a2df(annotations)
    else
        eeg_annotations = DataFrame(id=[""], time=[""], description=[""])
    end

    eeg_duration_samples = size(eeg_signals, 2)
    eeg_duration_seconds = size(eeg_signals, 2) / sampling_rate[1]
    eeg_time = collect(0:(1 / sampling_rate[1]):eeg_duration_seconds)
    eeg_time = eeg_time[1:end - 1]
    sampling_rate = round.(Int64, sampling_rate)
    eeg_filesize_mb = round(filesize(file_name) / 1024^2, digits=2)

    eeg_header = Dict(:eeg_filename => file_name,
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
                      :note => "",
                      :annotations => has_annotations)

    eeg_components = Vector{Any}()
    eeg_epochs_time = eeg_time

    eeg = EEG(eeg_header, eeg_time, eeg_epochs_time, eeg_signals, eeg_components, eeg_annotations)

    return eeg
end