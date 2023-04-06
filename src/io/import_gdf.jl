export import_gdf

"""
    import_gdf(file_name; detect_type)

Load GDF file and return `NeuroAnalyzer.NEURO` object.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: detect channel type based on its label

# Returns

- `obj::NeuroAnalyzer.NEURO`

# Notes

- sampling_rate = n.samples ÷ data.record.duration
- gain = (physical maximum - physical minimum) ÷ (digital maximum - digital minimum)
- value = (value - digital minimum ) × gain + physical minimum

# Source

1. Kemp B, Varri A, Rosa AC, Nielsen KD, Gade J. A simple format for exchange of digitized polygraphic recordings. Electroencephalography and Clinical Neurophysiology. 1992 May;82(5):391–3. 
2. Kemp B, Olivan J. European data format ‘plus’(EDF+), an EDF alike standard format for the exchange of physiological data. Clinical Neurophysiology 2003;114:1755–61.
3. https://www.edfplus.info/specs/
"""
function import_gdf(file_name::String; detect_type::Bool=true)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    splitext(file_name)[2] == ".gdf" || throw(ArgumentError("This is not an EDF file."))

    file_type = ""

    fid = ""
    try
        fid = open(file_name, "r")
    catch
        error("File $file_name cannot be loaded.")
    end

    header = UInt8[]
    readbytes!(fid, header, 256)
    # header = String(Char.(header))
    file_type = String(Char.(header[1:3]))
    file_type_ver = parse(Float64, String(Char.(header[5:8])))
    file_type != "GDF" && throw(ArgumentError("File $file_name is not a GDF file."))

    (file_type_ver == 1.25 || file_type_ver == 2.20) || _info("GDF versions other than 1.25 and 2.20 may not be supported correctly. ")

    if file_type_ver >= 1.25 && file_type_ver < 2.20

        patient = replace(strip(String(Char.(header[9:88]))), '\0'=>"")
        recording = replace(strip(String(Char.(header[89:168]))), '\0'=>"")
        recording_date = replace(strip(String(Char.(header[169:184]))), '\0'=>"")
        recording_date == "" && (recording_time = "")
        header_bytes = reinterpret(Int64, header[185:192])[1]
        equipment_id = reinterpret(Int64, header[193:200])[1]
        lab_id = reinterpret(Int64, header[201:208])[1]
        technician_id = reinterpret(Int64, header[209:216])[1]
        data_records = reinterpret(Int64, header[237:244])[1]
        sampling_rate = Int64(reinterpret(Int32, header[245:252])[2] ÷ reinterpret(Int32, header[245:252])[1])
        ch_n = reinterpret(Int32, header[253:256])[1]

        clabels = String[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 16)
            push!(clabels, replace(strip(String(Char.(buf))), '\0'=>""))
        end

        transducers = String[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 80)
            push!(transducers, replace(strip(String(Char.(buf))), '\0'=>""))
        end

        units = String[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 8)
            push!(units, replace(strip(String(Char.(buf))), '\0'=>"", '\x10'=>""))
            # push!(units, strip(String(Char.(buf))))
        end
        units = replace(lowercase.(units), "uv"=>"μV")

        physical_minimum = Float64[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 8)
            push!(physical_minimum, reinterpret(Float64, buf)[1])
        end

        physical_maximum = Float64[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 8)
            push!(physical_maximum, reinterpret(Float64, buf)[1])
        end

        digital_minimum = Int64[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 8)
            push!(digital_minimum, reinterpret(Int64, buf)[1])
        end

        digital_maximum = Int64[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 8)
            push!(digital_maximum, reinterpret(Int64, buf)[1])
        end

        prefiltering = String[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 80)
            push!(prefiltering, replace(strip(String(Char.(buf))), '\0'=>""))
        end

        samples_per_datarecord = Int32[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 4)
            push!(samples_per_datarecord, reinterpret(Int32, buf)[1])
        end

        ch_type = Int32[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 4)
            push!(ch_type, reinterpret(Int32, buf)[1])
        end
    
        close(fid)

    elseif file_type_ver >= 2.20

        patient = replace(strip(String(Char.(header[9:88]))), '\0'=>"")
        recording = replace(strip(String(Char.(header[89:168]))), '\0'=>"")
        recording_date = replace(strip(String(Char.(header[169:184]))), '\0'=>"")
        recording_date == "" && (recording_time = "")
        header_bytes = reinterpret(Int16, header[185:186])[1] * 256
        equipment_id = reinterpret(Int64, header[193:200])[1]
        lab_id = reinterpret(Int64, header[201:208])[1]
        technician_id = reinterpret(Int64, header[209:216])[1]
        data_records = reinterpret(Int64, header[237:244])[1]
        sampling_rate = Int64(reinterpret(Int32, header[245:252])[2] ÷ reinterpret(Int32, header[245:252])[1])
        ch_n = reinterpret(Int32, header[253:256])[1]

        clabels = String[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 16)
            push!(clabels, replace(strip(String(Char.(buf))), '\0'=>""))
        end

        transducers = String[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 80)
            push!(transducers, replace(strip(String(Char.(buf))), '\0'=>""))
        end

        units = String[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 8)
            push!(units, replace(strip(String(Char.(buf))), '\0'=>"", '\x10'=>""))
            # push!(units, strip(String(Char.(buf))))
        end
        units = replace(lowercase.(units), "uv"=>"μV")

        physical_minimum = Float64[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 8)
            push!(physical_minimum, reinterpret(Float64, buf)[1])
        end

        physical_maximum = Float64[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 8)
            push!(physical_maximum, reinterpret(Float64, buf)[1])
        end

        digital_minimum = Int64[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 8)
            push!(digital_minimum, reinterpret(Int64, buf)[1])
        end

        digital_maximum = Int64[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 8)
            push!(digital_maximum, reinterpret(Int64, buf)[1])
        end

        prefiltering = String[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 80)
            push!(prefiltering, replace(strip(String(Char.(buf))), '\0'=>""))
        end

        samples_per_datarecord = Int32[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 4)
            push!(samples_per_datarecord, reinterpret(Int32, buf)[1])
        end

        ch_type = Int32[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 4)
            push!(ch_type, reinterpret(Int32, buf)[1])
        end

        close(fid)

    end

    fid = ""
    try
        fid = open(file_name, "r")
    catch
        error("File $file_name cannot be loaded.")
    end

    header = UInt8[]
    readbytes!(fid, header, header_bytes)

    signal = Float64[]
    @inbounds for idx in 1:ch_n
        buf = UInt8[]
        if ch_type[idx] == 0
            # char => uint8
            readbytes!(fid, buf, 1 * samples_per_datarecord[idx] * data_records)
            signal = vcat(signal, Float64.(reinterpret(UInt8, buf)))
        elseif ch_type[idx] == 1
            readbytes!(fid, buf, 1 * samples_per_datarecord[idx] * data_records)
            signal = vcat(signal, Float64.(reinterpret(Int8, buf)))
        elseif ch_type[idx] == 2
            readbytes!(fid, buf, 1 * samples_per_datarecord[idx] * data_records)
            signal = vcat(signal, Float64.(reinterpret(UInt8, buf)))
        elseif ch_type[idx] == 3
            readbytes!(fid, buf, 2 * samples_per_datarecord[idx] * data_records)
            signal = vcat(signal, Float64.(reinterpret(Int16, buf)))
        elseif ch_type[idx] == 4
            readbytes!(fid, buf, 2 * samples_per_datarecord[idx] * data_records)
            signal = vcat(signal, Float64.(reinterpret(UInt16, buf)))
        elseif ch_type[idx] == 5
            readbytes!(fid, buf, 4 * samples_per_datarecord[idx] * data_records)
            signal = vcat(signal, Float64.(reinterpret(Int32, buf)))
        elseif ch_type[idx] == 6
            readbytes!(fid, buf, 4 * samples_per_datarecord[idx] * data_records)
            signal = vcat(signal, Float64.(reinterpret(UInt32, buf)))
        elseif ch_type[idx] == 7
            readbytes!(fid, buf, 8 * samples_per_datarecord[idx] * data_records)
            signal = vcat(signal, Float64.(reinterpret(Int64, buf)))
        elseif ch_type[idx] == 16
            readbytes!(fid, buf, 4 * samples_per_datarecord[idx] * data_records)
            signal = vcat(signal, Float64.(reinterpret(Float32, buf)))
        elseif ch_type[idx] == 17
            readbytes!(fid, buf, 8 * samples_per_datarecord[idx] * data_records)
            signal = vcat(signal, Float64.(reinterpret(Float64, buf)))
        else
            @error "Unknown channel type: $(ch_type[idx])."
        end
    end

    close(fid)

    # split signal into channels
    data = zeros(ch_n, data_records * samples_per_datarecord[1])
    t_idx = 1
    @inbounds for idx in 1:ch_n:length(data)
        data[:, t_idx] = @views signal[idx:(idx + ch_n - 1)]
        t_idx += 1
    end

    gain = @. (physical_maximum - physical_minimum) / (digital_maximum - digital_minimum)
    data .*= gain

    annotation_channels = Int64[]

    clabels = NeuroAnalyzer._clean_labels(clabels)
    if detect_type == true
        ch_type = NeuroAnalyzer._set_channel_types(clabels, "eeg")
    else
        ch_type = repeat(["eeg"], ch_n)
    end
    units = [_set_units(ch_type[idx]) for idx in 1:ch_n]

    if length(annotation_channels) == 0
        markers = DataFrame(:id=>String[], :start=>Int64[], :length=>Int64[], :description=>String[], :channel=>Int64[])
    else
        markers = _a2df(annotations)
        deleteat!(ch_type, annotation_channels)
        deleteat!(transducers, annotation_channels)
        deleteat!(units, annotation_channels)
        deleteat!(prefiltering, annotation_channels)
        deleteat!(clabels, annotation_channels)
        data = data[setdiff(collect(1:ch_n), annotation_channels), :, :]
        ch_n -= length(annotation_channels)
    end

    ch_order = _sort_channels(ch_type)

    time_pts = round.(collect(0:1/sampling_rate:size(data, 2) * size(data, 3) / sampling_rate)[1:end-1], digits=3)
    ep_time = round.((collect(0:1/sampling_rate:size(data, 2) / sampling_rate))[1:end-1], digits=3)
    
    file_size_mb = round(filesize(file_name) / 1024^2, digits=2)
    
    data_type = "eeg"

    s = _create_subject(id="",
                        first_name="",
                        middle_name="",
                        last_name=patient,
                        handedness="",
                        weight=-1,
                        height=-1)
    r = _create_recording_eeg(data_type=data_type,
                              file_name=file_name,
                              file_size_mb=file_size_mb,
                              file_type=file_type,
                              recording=recording,
                              recording_date=recording_date,
                              recording_time=recording_time,
                              recording_notes="",
                              channel_type=ch_type[ch_order],
                              reference="",
                              clabels=clabels[ch_order],
                              transducers=transducers[ch_order],
                              units=units[ch_order],
                              prefiltering=prefiltering[ch_order],
                              sampling_rate=sampling_rate,
                              gain=gain[ch_order])
    e = _create_experiment(name="", notes="", design="")

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = String[]

    locs = DataFrame(:channel=>Int64,
                     :labels=>String[],
                     :loc_theta=>Float64[],
                     :loc_radius=>Float64[],
                     :loc_x=>Float64[],
                     :loc_y=>Float64[],
                     :loc_z=>Float64[],
                     :loc_radius_sph=>Float64[],
                     :loc_theta_sph=>Float64[],
                     :loc_phi_sph=>Float64[])

    obj = NeuroAnalyzer.NEURO(hdr, time_pts, ep_time, data[ch_order, :, :], components, markers, locs, history)

    _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(channel_n(obj)) × $(epoch_len(obj)) × $(epoch_n(obj)); $(signal_len(obj) / sr(obj)) s)")

    return obj
    
end

