"""
    eeg_edit(file_name; read_annotations=true, header_only=false, clean_labels=true)

Loads EDF/EDFPlus file and returns EEG object.

# Arguments

- `file_name::String`: name of the file to load
- `read_annotations::Bool`: read annotations from EDF+ file (currently not implemented)
- `clean_labels::Bool`: only keep channel names in channel labels

# Returns

- `eeg:EEG`

# Notes

sampling_rate = n.samples / data.record.duration

gain = (physical_maximum - physical_minimum) / (digital_maximum - digital_minimum)

value = (value - digital_minimum ) * gain + physical_minimum

# Source

Kemp B, Värri A, Rosa AC, Nielsen KD, Gade J. A simple format for exchange of digitized polygraphic recordings. Electroencephalography and Clinical Neurophysiology. 1992 May;82(5):391–3. 
"""
function eeg_import_edf(file_name::String; read_annotations::Bool=true, clean_labels::Bool=true)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    fid = open(file_name)

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
    channel_n  = parse(Int, rstrip(header[253:256]))

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
        labels[idx] = rstrip(header[1 + ((idx - 1) * 16):(idx * 16)])
    end

    header = zeros(UInt8, channel_n * 80)
    readbytes!(fid, header, channel_n * 80)
    header = String(Char.(header))
    for idx in 1:channel_n
        transducers[idx] = rstrip(header[1 + ((idx - 1) * 80):(idx * 80)])
    end

    header = zeros(UInt8, channel_n * 8)
    readbytes!(fid, header, channel_n * 8)
    header = String(Char.(header))
    for idx in 1:channel_n
        physical_dimension[idx] = rstrip(header[1 + ((idx - 1) * 8):(idx * 8)])
    end

    header = zeros(UInt8, channel_n * 8)
    readbytes!(fid, header, channel_n * 8)
    header = String(Char.(header))
    for idx in 1:channel_n
        physical_minimum[idx] = parse(Float64, rstrip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = zeros(UInt8, channel_n * 8)
    readbytes!(fid, header, channel_n * 8)
    header = String(Char.(header))
    for idx in 1:channel_n
        physical_maximum[idx] = parse(Float64, rstrip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = zeros(UInt8, channel_n * 8)
    readbytes!(fid, header, channel_n * 8)
    header = String(Char.(header))
    for idx in 1:channel_n
        digital_minimum[idx] = parse(Float64, rstrip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = zeros(UInt8, channel_n * 8)
    readbytes!(fid, header, channel_n * 8)
    header = String(Char.(header))
    for idx in 1:channel_n
        digital_maximum[idx] = parse(Float64, rstrip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    header = zeros(UInt8, channel_n * 80)
    readbytes!(fid, header, channel_n * 80)
    header = String(Char.(header))
    for idx in 1:channel_n
        prefiltering[idx] = rstrip(header[1 + ((idx - 1) * 80):(idx * 80)])
    end

    header = zeros(UInt8, channel_n * 8)
    readbytes!(fid, header, channel_n * 8)
    header = String(Char.(header))
    for idx in 1:channel_n
        samples_per_datarecord[idx] = parse(Int, rstrip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    close(fid)

    sampling_rate = Vector{Float64}(undef, channel_n)
    gain = Vector{Float64}(undef, channel_n)
    for idx in 1:channel_n
        sampling_rate[idx] = samples_per_datarecord[idx] / data_records_duration
        gain[idx] = (physical_maximum[idx] - physical_minimum[idx]) / (digital_maximum[idx] - digital_minimum[idx])
    end

    if clean_labels == true
        labels = replace.(labels, "EEG " => "")
        labels = replace.(labels, "ECG " => "")
    end

    fid = open(file_name)
    header = zeros(UInt8, data_offset)
    readbytes!(fid, header, data_offset)
    eeg_signals = zeros(channel_n, samples_per_datarecord[1] * data_records, 1)

    for idx1 in 1:data_records
        for idx2 in 1:channel_n
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
                      :reference => "",
                      :channel_locations => false,
                      :loc_x_theta => Float64[],
                      :loc_y_theta => Float64[],
                      :loc_theta => Float64[],
                      :loc_phi => Float64[],
                      :loc_x => Float64[],
                      :loc_y => Float64[],
                      :loc_z => Float64[],
                      :loc_x_sph => Float64[],
                      :loc_x_sph => Float64[],
                      :loc_y_sph => Float64[],
                      :loc_radius_sph => Float64[],
                      :loc_theta_sph => Float64[],
                      :loc_phi_sph => Float64[],
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
                      :gain => gain)

    eeg_components = Vector{Any}()

    eeg = EEG(eeg_header, eeg_time, eeg_signals, eeg_components)

    return eeg
end

"""
    eeg_import_ced(file_name)

Loads electrode positions from CED file.

# Arguments

- `file_name::String`

# Returns

- `sensors::DataFrame`
"""
function eeg_import_ced(file_name)

    sensors = CSV.read(file_name, delim="\t", DataFrame)
    
    return sensors
end

"""
    eeg_import_locs(file_name)

Loads electrode positions from LOCS file.

# Arguments

- `file_name::String`

# Returns

- `sensors::DataFrame`
"""
function eeg_import_locs(file_name)

    sensors = CSV.read(file_name, header=false, delim="\t", DataFrame)
    DataFrames.rename!(sensors, [:Number, :theta, :radius, :labels])

    return sensors
end

"""
    eeg_import_elc(file_name)

Loads electrode positions from ELC file.

# Arguments

- `file_name::String`

# Returns

- `sensors::DataFrame`
"""
function eeg_import_elc(file_name)

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
    locx = zeros(locs_n)
    locy = zeros(locs_n)
    locz = zeros(locs_n)
    idx2 = 1
    for idx1 in locs_l:(locs_l + locs_n - 1)
        l = elc_file[idx1]
        l[1] == ' ' && (l = l[2:end])
        locx[idx2], locy[idx2], locz[idx2] = parse.(Float64, split(l, ' '))
        idx2 += 1
    end
    idx2 = 1
    for idx1 in (locs_l + 1 + locs_n):(locs_l + (2 * locs_n))
        labels[idx2] = elc_file[idx1]
        idx2 += 1
    end

    sensors = DataFrame(:labels => labels, :xlocs => locx, :ylocs => locy, :zlocs => locz)

    return sensors
end

"""
    eeg_load_electrodes(eeg; file_name)

Loads electrode positions from 
- CED

# Arguments

- `eeg::EEG`
- `file_name::String`

# Returns

- `eeg:EEG`
"""
function eeg_load_electrodes(eeg::EEG; file_name)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    length(eeg.eeg_header[:labels]) > 0 || throw(ArgumentError("EEG does not contain labels, use eeg_add_labels() first."))

    loc_x_theta = Vector{Float64}()
    loc_y_theta = Vector{Float64}()
    loc_theta = Vector{Float64}()
    loc_phi = Vector{Float64}()
    loc_x_sph = Vector{Float64}()
    loc_y_sph = Vector{Float64}()
    loc_z_sph = Vector{Float64}()
    loc_X = Vector{Float64}()
    loc_Y = Vector{Float64}()
    loc_Z = Vector{Float64}()
    loc_radius_sph = Vector{Float64}()
    loc_theta_sph = Vector{Float64}()
    loc_phi_sph = Vector{Float64}()


    if splitext(file_name)[2] == ".ced"
        sensors = eeg_import_ced(file_name)

        f_labels = lowercase.(sensors[:, :labels])

        loc_Y = sensors[:, :X]
        loc_X = sensors[:, :Y]
        loc_Z = sensors[:, :Z]

        loc_theta = sensors[:, :theta]
        loc_phi = sensors[:, :radius]

        for idx in 1:length(f_labels)
            y, x = pol2cart(pi / 180 * loc_theta[idx], loc_phi[idx])
            push!(loc_x_theta, y)
            push!(loc_y_theta, x)
        end

        loc_radius_sph = sensors[:, :sph_radius]
        loc_theta_sph = sensors[:, :sph_theta]
        loc_phi_sph = sensors[:, :sph_phi]

        for idx in 1:length(f_labels)
            x, y, z = sph2cart(loc_radius_sph[idx], loc_theta_sph[idx], loc_phi_sph[idx])
            push!(loc_x_sph, x)
            push!(loc_y_sph, y)
            push!(loc_z_sph, z)
        end
    end

    if splitext(file_name)[2] == ".elc"
        sensors = eeg_import_elc(file_name)

        f_labels = lowercase.(sensors[:, :labels])
        
        x = sensors[:, :xlocs]
        y = sensors[:, :ylocs]
        z = sensors[:, :zlocs]

        for idx in 1:length(f_labels)
            push!(loc_x_sph, x[idx])
            push!(loc_y_sph, y[idx])
            push!(loc_z_sph, z[idx])
        end
    end

    if splitext(file_name)[2] == ".locs"
        sensors = eeg_import_locs(file_name)

        f_labels = lowercase.(sensors[:, :labels])

        loc_theta = sensors[:, :theta]
        loc_phi = sensors[:, :radius]

        for idx in 1:length(f_labels)
            y, x = pol2cart(pi / 180 * loc_theta[idx], loc_phi[idx])
            push!(loc_x_theta, y)
            push!(loc_y_theta, x)
        end
    end

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
    eeg_new = EEG(deepcopy(eeg.eeg_header), deepcopy(eeg.eeg_time), deepcopy(eeg.eeg_signals), deepcopy(eeg.eeg_components))
    eeg_new.eeg_header[:channel_locations] = true
    length(loc_x_theta) > 0 && (eeg_new.eeg_header[:loc_x_theta] = loc_x_theta[labels_idx])
    length(loc_y_theta) > 0 && (eeg_new.eeg_header[:loc_y_theta] = loc_y_theta[labels_idx])
    length(loc_theta) > 0 && (eeg_new.eeg_header[:loc_theta] = loc_theta[labels_idx])
    length(loc_phi) > 0 && (eeg_new.eeg_header[:loc_phi] = loc_phi[labels_idx])
    length(loc_X) > 0 && (eeg_new.eeg_header[:loc_x] = loc_X[labels_idx])
    length(loc_Y) > 0 && (eeg_new.eeg_header[:loc_y] = loc_Y[labels_idx])
    length(loc_Z) > 0 && (eeg_new.eeg_header[:loc_z] = loc_Z[labels_idx])
    length(loc_x_sph) > 0 && (eeg_new.eeg_header[:loc_x_sph] = loc_x_sph[labels_idx])
    length(loc_y_sph) > 0 && (eeg_new.eeg_header[:loc_x_sph] = loc_y_sph[labels_idx])
    length(loc_z_sph) > 0 && (eeg_new.eeg_header[:loc_y_sph] = loc_z_sph[labels_idx])
    length(loc_radius_sph) > 0 && (eeg_new.eeg_header[:loc_radius_sph] = loc_radius_sph[labels_idx])
    length(loc_theta_sph) > 0 && (eeg_new.eeg_header[:loc_theta_sph] = loc_theta_sph[labels_idx])
    length(loc_phi_sph) > 0 && (eeg_new.eeg_header[:loc_phi_sph] = loc_phi_sph[labels_idx])

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_load_sensor_positions(EEG, $file_name)")

    return eeg_new
end

"""
    eeg_load_electrodes!(eeg; file_name)

Loads electrode positions from:
- CED

# Arguments

- `eeg::EEG`
- `file_name::String`
"""
function eeg_load_electrodes!(eeg::EEG; file_name)
    
    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    length(eeg.eeg_header[:labels]) > 0 || throw(ArgumentError("EEG does not contain labels, use eeg_add_labels() first."))

    loc_x_theta = Vector{Float64}()
    loc_y_theta = Vector{Float64}()
    loc_theta = Vector{Float64}()
    loc_phi = Vector{Float64}()
    loc_x_sph = Vector{Float64}()
    loc_y_sph = Vector{Float64}()
    loc_z_sph = Vector{Float64}()
    loc_X = Vector{Float64}()
    loc_Y = Vector{Float64}()
    loc_Z = Vector{Float64}()
    loc_radius_sph = Vector{Float64}()
    loc_theta_sph = Vector{Float64}()
    loc_phi_sph = Vector{Float64}()


    if splitext(file_name)[2] == ".ced"
        sensors = eeg_import_ced(file_name)

        f_labels = lowercase.(sensors[:, :labels])

        loc_Y = sensors[:, :X]
        loc_X = sensors[:, :Y]
        loc_Z = sensors[:, :Z]

        loc_theta = sensors[:, :theta]
        loc_phi = sensors[:, :radius]

        for idx in 1:length(f_labels)
            y, x = pol2cart(pi / 180 * loc_theta[idx], loc_phi[idx])
            push!(loc_x_theta, y)
            push!(loc_y_theta, x)
        end

        loc_radius_sph = sensors[:, :sph_radius]
        loc_theta_sph = sensors[:, :sph_theta]
        loc_phi_sph = sensors[:, :sph_phi]

        for idx in 1:length(f_labels)
            x, y, z = sph2cart(loc_radius_sph[idx], loc_theta_sph[idx], loc_phi_sph[idx])
            push!(loc_x_sph, x)
            push!(loc_y_sph, y)
            push!(loc_z_sph, z)
        end
    end

    if splitext(file_name)[2] == ".elc"
        sensors = eeg_import_elc(file_name)

        f_labels = lowercase.(sensors[:, :labels])
        
        x = sensors[:, :xlocs]
        y = sensors[:, :ylocs]
        z = sensors[:, :zlocs]

        for idx in 1:length(f_labels)
            push!(loc_x_sph, x[idx])
            push!(loc_y_sph, y[idx])
            push!(loc_z_sph, z[idx])
        end
    end

    if splitext(file_name)[2] == ".locs"
        sensors = eeg_import_locs(file_name)

        f_labels = lowercase.(sensors[:, :labels])

        loc_theta = sensors[:, :theta]
        loc_phi = sensors[:, :radius]

        for idx in 1:length(f_labels)
            y, x = pol2cart(pi / 180 * loc_theta[idx], loc_phi[idx])
            push!(loc_x_theta, y)
            push!(loc_y_theta, x)
        end
    end

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
    length(loc_x_theta) > 0 && (eeg.eeg_header[:loc_x_theta] = loc_x_theta[labels_idx])
    length(loc_y_theta) > 0 && (eeg.eeg_header[:loc_y_theta] = loc_y_theta[labels_idx])
    length(loc_theta) > 0 && (eeg.eeg_header[:loc_theta] = loc_theta[labels_idx])
    length(loc_phi) > 0 && (eeg.eeg_header[:loc_phi] = loc_phi[labels_idx])
    length(loc_X) > 0 && (eeg.eeg_header[:loc_x] = loc_X[labels_idx])
    length(loc_Y) > 0 && (eeg.eeg_header[:loc_y] = loc_Y[labels_idx])
    length(loc_Z) > 0 && (eeg.eeg_header[:loc_z] = loc_Z[labels_idx])
    length(loc_x_sph) > 0 && (eeg.eeg_header[:loc_x_sph] = loc_x_sph[labels_idx])
    length(loc_y_sph) > 0 && (eeg.eeg_header[:loc_x_sph] = loc_y_sph[labels_idx])
    length(loc_z_sph) > 0 && (eeg.eeg_header[:loc_y_sph] = loc_z_sph[labels_idx])
    length(loc_radius_sph) > 0 && (eeg.eeg_header[:loc_radius_sph] = loc_radius_sph[labels_idx])
    length(loc_theta_sph) > 0 && (eeg.eeg_header[:loc_theta_sph] = loc_theta_sph[labels_idx])
    length(loc_phi_sph) > 0 && (eeg.eeg_header[:loc_phi_sph] = loc_phi_sph[labels_idx])

    # add entry to :history field
    push!(eeg.eeg_header[:history], "eeg_load_sensor_positions(EEG, $file_name)")

    return
end

"""
    eeg_save(eeg; file_name, overwrite=false)

Saves the `eeg` to `file_name` file (HDF5-based).

# Arguments

- `eeg::EEG`
- `file_name::String`: file name
- `overwrite::Bool`

# Returns

- `success::Bool`
"""
function eeg_save(eeg::EEG; file_name::String, overwrite::Bool=false)

    (isfile(file_name) && overwrite == false) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))

    eeg.eeg_header[:eeg_filename] = file_name

    save_object("/tmp/$file_name", eeg)
    eeg.eeg_header[:eeg_filesize_mb] = round(filesize("/tmp/$file_name") / 1024, digits=2)
    rm("/tmp/$file_name")

    save_object(file_name, eeg)

    return true
end

"""
    eeg_load(file_name)

Loads the `eeg` from `file_name` file (HDF5-based).

# Arguments

- `file_name::String`: file name

# Returns

- `eeg::EEG`
"""
function eeg_load(file_name::String)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    eeg = load_object(file_name)

    return eeg
end

"""
    eeg_export_csv(eeg, file_name, header, overwrite)

Exports EEG data as CSV.

# Arguments

- `eeg::EEG`
- `file_name::String`
- `header::Bool`: export header
- `components::Bool`: export components
- `overwrite::Bool`

# Returns

- `success::Bool`
"""
function eeg_export_csv(eeg::EEG; file_name::String, header::Bool=false, components::Bool=false, overwrite::Bool=false)

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
    df = DataFrame(s, l)

    CSV.write(file_name, df)
    header == false && return true

    # HEADER
    file_name = replace(file_name, ".csv" => "_header.csv")
    (isfile(file_name) && overwrite == false) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))
    f = open(file_name, "w")
    for (key, value) in eeg.eeg_header
        println(f, key, ": ", value)
    end
    close(f)

    # COMPONENTS
    file_name = replace(file_name, ".csv" => "_components.csv")
    (isfile(file_name) && overwrite == false) && throw(ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."))
    f = open(file_name, "w")
    for idx in 1:length(eeg.eeg_header[:components])
        println(f, "component: $(eeg.eeg_header[:components][idx])")
        println(f, eeg.eeg_components[idx])
        println(f, "---")
    end
    close(f)

    return true
end

"""
    eeg_add_labels(eeg::EEG, labels::Vector{String})

Adds `labels` to `eeg` channels.

# Arguments

- `eeg::EEG`
- `labels::Vector{String}`

# Returns

- `eeg::EEG`
"""
function eeg_add_labels(eeg::EEG, labels::Vector{String})

    length(labels) == eeg_channel_n(eeg) || throw(ArgumentError("labels length must be $(eeg_channel_n(eeg))."))
    
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_header[:labels] = labels

    push!(eeg_new.eeg_header[:history], "eeg_add_labels(EEG, labels=$labels")
 
    return eeg_new
end

"""
    eeg_add_labels!(eeg::EEG, labels::Vector{String})

Adds `labels` to `eeg` channels.

# Arguments

- `eeg::EEG`
- `labels::Vector{String}`
"""
function eeg_add_labels!(eeg::EEG, labels::Vector{String})

    length(labels) == eeg_channel_n(eeg) || throw(ArgumentError("labels length must be $(eeg_channel_n(eeg))."))
    
    eeg.eeg_header[:labels] = labels

    push!(eeg.eeg_header[:history], "eeg_add_labels(EEG, labels=$labels")
end