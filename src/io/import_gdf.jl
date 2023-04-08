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

1. Schlögl A, Filz O, Ramoser H, Pfurtscheller G. GDF - A General Dataformat for Biosignals Version 1.25. 1998
2. Schlögl, A. GDF - A General Dataformat for Biosignals Version 2.12. 2009
3. Schlögl, A. GDF - A General Dataformat for Biosignals Version 2.51. 2013
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
    file_type = String(Char.(header[1:3]))
    file_type_ver = parse(Float64, String(Char.(header[5:8])))
    file_type != "GDF" && throw(ArgumentError("File $file_name is not a GDF file."))

    (file_type_ver == 1.25 || file_type_ver == 2.20) || _info("GDF versions other than 1.25 and 2.20 may not be supported correctly. ")

    if file_type_ver < 2.00

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

        gdf_type = Int32[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 4)
            push!(gdf_type, reinterpret(Int32, buf)[1])
        end

    else

        patient = replace(strip(String(Char.(header[9:66]))), '\0'=>"")
        patient_weight = header[86] == 0 ? -1 : header[86]
        patient_height = header[87] == 0 ? -1 : header[86]
        bitstring(header[88])[(end - 1):end] == "00" && (patient_gender = "")
        bitstring(header[88])[(end - 1):end] == "01" && (patient_gender = "M")
        bitstring(header[88])[(end - 1):end] == "10" && (patient_gender = "W")
        recording = replace(strip(String(Char.(header[89:152]))), '\0'=>"")
        recording_date = replace(strip(String(Char.(header[169:176]))), '\0'=>"")
        recording_date == "" && (recording_time = "")
        header_bytes = reinterpret(Int16, header[185:186])[1] * 256
        etv_hdr = reinterpret(Int16, header[185:186])[1]
        equipment_id = reinterpret(Int64, header[193:200])[1]
        ref_elec_x, ref_elec_y, ref_elec_z = reinterpret(Float32, header[213:224])
        ground_elec_x, ground_elec_y, ground_elec_z = reinterpret(Float32, header[225:236])
        data_records = reinterpret(Int64, header[237:244])[1]
        data_records == -1 && @error "Number of data records cannot be -1."
        sampling_rate = Int64(reinterpret(Int32, header[245:252])[2] ÷ reinterpret(Int32, header[245:252])[1])
        ch_n = reinterpret(Int16, header[253:254])[1]

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

        # obsolete
        units = String[]
        unit = zeros(Int64, ch_n)
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 6)
            push!(units, replace(strip(String(Char.(buf))), '\0'=>"", '\x10'=>""))
        end

        units = String[]
        units_code = UInt16[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 2)
            push!(units_code, reinterpret(UInt16, buf)[1])
        end
        for idx in 1:ch_n
            unit[idx] = parse(Int, "0b" * bitstring(units_code[idx])[1:(end - 5)] * "00000")
            unit[idx] == 512 && (push!(units, ""))
            unit[idx] == 544 && (push!(units, "%"))
            unit[idx] == 736 && (push!(units, "°"))
            unit[idx] == 768 && (push!(units, "rad"))
            unit[idx] == 2496 && (push!(units, "Hz"))
            unit[idx] == 3872 && (push!(units, "mmHg"))
            unit[idx] == 4256 && (push!(units, "V"))
            unit[idx] == 4288 && (push!(units, "Ω"))
            unit[idx] == 6048 && (push!(units, "°C"))
            unit[idx] == 3072 && (push!(units, "l/min"))
            unit[idx] == 2848 && (push!(units, "l/(min m²)"))
            unit[idx] == 4128 && (push!(units, "dyn s/cm⁵"))
            unit[idx] == 6016 && (push!(units, "dyn s/m² cm⁵"))
        end
        for idx in 1:ch_n
            dec_factor = parse(Int, "0b" * bitstring(units_code[idx])[(end - 4):end])
            dec_factor == 10 && (units[idx] = "Y" * units[idx])
            dec_factor == 9 && (units[idx] = "Z" * units[idx])
            dec_factor == 8 && (units[idx] = "E" * units[idx])
            dec_factor == 7 && (units[idx] = "P" * units[idx])
            dec_factor == 6 && (units[idx] = "T" * units[idx])
            dec_factor == 5 && (units[idx] = "G" * units[idx])
            dec_factor == 4 && (units[idx] = "M" * units[idx])
            dec_factor == 3 && (units[idx] = "k" * units[idx])
            dec_factor == 2 && (units[idx] = "h" * units[idx])
            dec_factor == 1 && (units[idx] = "da" * units[idx])
            dec_factor == 0 && (units[idx] = "" * units[idx])
            dec_factor == 16 && (units[idx] = "d" * units[idx])
            dec_factor == 17 && (units[idx] = "c" * units[idx])
            dec_factor == 18 && (units[idx] = "m" * units[idx])
            dec_factor == 19 && (units[idx] = "μ" * units[idx])
            dec_factor == 20 && (units[idx] = "n" * units[idx])
            dec_factor == 21 && (units[idx] = "p" * units[idx])
            dec_factor == 22 && (units[idx] = "f" * units[idx])
            dec_factor == 23 && (units[idx] = "a" * units[idx])
            dec_factor == 24 && (units[idx] = "z" * units[idx])
            dec_factor == 25 && (units[idx] = "y" * units[idx])
        end

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

        digital_minimum = Float64[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 8)
            push!(digital_minimum, reinterpret(Float64, buf)[1])
        end

        digital_maximum = Float64[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 8)
            push!(digital_maximum, reinterpret(Float64, buf)[1])
        end

        # obsolete
        prefiltering = String[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 64)
            push!(prefiltering, replace(strip(String(Char.(buf))), '\0'=>""))
        end

        time_offset = Float32[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 4)
            push!(time_offset, reinterpret(Int32, buf)[1])
        end
        time_offset == Float32[] && (time_offset = zeros(Int64, ch_n))

        prefiltering_lp = Float32[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 4)
            push!(prefiltering_lp, reinterpret(Float32, buf)[1])
        end
        prefiltering_lp[findall(isnan, prefiltering_lp)] .= 0

        prefiltering_hp = Float32[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 4)
            push!(prefiltering_hp, reinterpret(Float32, buf)[1])
        end
        prefiltering_hp[findall(isnan, prefiltering_hp)] .= 0

        prefiltering_bs = Float32[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 4)
            push!(prefiltering_bs, reinterpret(Float32, buf)[1])
        end
        prefiltering_bs[findall(isnan, prefiltering_bs)] .= 0

        samples_per_datarecord = Int32[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 4)
            push!(samples_per_datarecord, reinterpret(Int32, buf)[1])
        end

        gdf_type = Int32[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 4)
            push!(gdf_type, reinterpret(Int32, buf)[1])
        end

        loc_x = Float32[]
        loc_y = Float32[]
        loc_z = Float32[]
        buf = UInt8[]
        for idx in 1:ch_n
            readbytes!(fid, buf, 4)
            push!(loc_x, reinterpret(Float32, buf)[1])
            readbytes!(fid, buf, 4)
            push!(loc_y, reinterpret(Float32, buf)[1])
            readbytes!(fid, buf, 4)
            push!(loc_z, reinterpret(Float32, buf)[1])
        end
        loc_x[findall(isnan, loc_x)] .= 0
        loc_y[findall(isnan, loc_y)] .= 0
        loc_z[findall(isnan, loc_z)] .= 0
        loc_x = round.(loc_x, digits=3)
        loc_y = round.(loc_y, digits=3)
        loc_z = round.(loc_z, digits=3)

        imp = UInt8[]
        if file_type_ver >= 2.19
            imp = zeros(ch_n)
            buf = UInt8[]
            for idx in 1:ch_n
                readbytes!(fid, buf, 20)
                unit[idx] == 4256 || unit[idx] == 4288 && (imp[idx] = reinterpret(Float32, buf[1:4])[1])
            end
            imp[findall(isnan, imp)] .= 0
            imp = round.(imp, digits=3)
        else
            buf = UInt8[]
            for idx in 1:ch_n
                readbytes!(fid, buf, 1)
                push!(imp, reinterpret(Int32, buf)[1])
            end
            buf = UInt8[]
            readbytes!(fid, buf, 19 * ch_n)
            imp = @. Float64(2 ^ (imp / 8))
        end

    end

    close(fid)

    fid = ""
    try
        fid = open(file_name, "r")
    catch
        error("File $file_name cannot be loaded.")
    end
    header = UInt8[]
    readbytes!(fid, header, header_bytes)

    # check TLV header
    if file_type_ver >= 2.10 && (etv_hdr - (ch_n + 1)) * 256 > 0
        tlv_tag = UInt8[]
        tlv_len = UInt8[]
        tlv_val = UInt8[]
        readbytes!(fid, tlv_tag, 1)
        if tlv_tag != [0x00]
            readbytes!(fid, tlv_len, 3)
            b1 = Int32(tlv_len[1]) << 8
            b2 = Int32(tlv_len[2]) << 16
            b3 = -Int32(-tlv_len[3]) << 24
            tlv_len = Float64(((b1 | b2 | b3) >> 8))
            readbytes!(fid, tlv_val, 1)
        else
            readbytes!(fid, tlv_tag, 255)
        end
    end

    signal = Float64[]
    data_bytes = 0
    @inbounds for idx in 1:ch_n
        buf = UInt8[]
        if gdf_type[idx] == 0
            # char => uint8
            data_bytes += samples_per_datarecord[idx] * data_records
            readbytes!(fid, buf, 1 * samples_per_datarecord[idx] * data_records)
            signal = vcat(signal, Float64.(reinterpret(UInt8, buf)))
        elseif gdf_type[idx] == 1
            data_bytes += samples_per_datarecord[idx] * data_records
            readbytes!(fid, buf, 1 * samples_per_datarecord[idx] * data_records)
            signal = vcat(signal, Float64.(reinterpret(Int8, buf)))
        elseif gdf_type[idx] == 2
            data_bytes += samples_per_datarecord[idx] * data_records
            readbytes!(fid, buf, 1 * samples_per_datarecord[idx] * data_records)
            signal = vcat(signal, Float64.(reinterpret(UInt8, buf)))
        elseif gdf_type[idx] == 3
            data_bytes += 2 * samples_per_datarecord[idx] * data_records
            readbytes!(fid, buf, 2 * samples_per_datarecord[idx] * data_records)
            signal = vcat(signal, Float64.(reinterpret(Int16, buf)))
        elseif gdf_type[idx] == 4
            data_bytes += 2 * samples_per_datarecord[idx] * data_records
            readbytes!(fid, buf, 2 * samples_per_datarecord[idx] * data_records)
            signal = vcat(signal, Float64.(reinterpret(UInt16, buf)))
        elseif gdf_type[idx] == 5
            data_bytes += 4 * samples_per_datarecord[idx] * data_records
            readbytes!(fid, buf, 4 * samples_per_datarecord[idx] * data_records)
            signal = vcat(signal, Float64.(reinterpret(Int32, buf)))
        elseif gdf_type[idx] == 6
            data_bytes += 4 * samples_per_datarecord[idx] * data_records
            readbytes!(fid, buf, 4 * samples_per_datarecord[idx] * data_records)
            signal = vcat(signal, Float64.(reinterpret(UInt32, buf)))
        elseif gdf_type[idx] == 7
            data_bytes += 8 * samples_per_datarecord[idx] * data_records
            readbytes!(fid, buf, 8 * samples_per_datarecord[idx] * data_records)
            signal = vcat(signal, Float64.(reinterpret(Int64, buf)))
        elseif gdf_type[idx] == 16
            data_bytes += 4 * samples_per_datarecord[idx] * data_records
            readbytes!(fid, buf, 4 * samples_per_datarecord[idx] * data_records)
            signal = vcat(signal, Float64.(reinterpret(Float32, buf)))
        elseif gdf_type[idx] == 17
            data_bytes += 8 * samples_per_datarecord[idx] * data_records
            readbytes!(fid, buf, 8 * samples_per_datarecord[idx] * data_records)
            signal = vcat(signal, Float64.(reinterpret(Float64, buf)))
        else
            @error "Unknown channel type: $(gdf_type[idx])."
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

    markers = DataFrame(:id=>String[], :start=>Int64[], :length=>Int64[], :description=>String[], :channel=>Int64[])
    if file_type_ver < 2.0
        if header_bytes + data_bytes < filesize(file_name)
            fid = ""
            try
                fid = open(file_name, "r")
            catch
                error("File $file_name cannot be loaded.")
            end
            header = UInt8[]
            readbytes!(fid, header, header_bytes + data_bytes)
            etp = UInt8[]
            readbytes!(fid, etp, filesize(file_name))
            etp_mode = reinterpret(Int8, etp[1])[1]
            deleteat!(etp, 1)
            (etp_mode == 1 || etp_mode == 3) || throw(ArgumentError("ETP mode must be 1 or 3, found $etp_mode."))
            b1 = Int32(etp[1]) << 8
            b2 = Int32(etp[2]) << 16
            b3 = -Int32(-etp[3]) << 24
            etp_sr = Float64(((b1 | b2 | b3) >> 8))
            deleteat!(etp, 1:3)
            etp_number = reinterpret(Int32, etp[1:4])[1]
            deleteat!(etp, 1:4)
            start = Float64[]
            len = Float64[]
            id = String[]
            desc = String[]
            ch = Int[]
            for idx in 1:etp_number
                push!(id, "mrk")
                push!(start, Float64(reinterpret(Int32, etp[1:4])[1]))
                push!(len, 0)
                push!(ch, 0)
                deleteat!(etp, 1:4)
            end
            for idx in 1:etp_number
                event = "unknown code ($(etp[1])$(etp[2]))"
                etp[1] == 0x01 && etp[2] == 0x01 && (event = "artifact:EOG (blinks)")
                etp[1] == 0x01 && etp[2] == 0x02 && (event = "artifact:ECG")
                etp[1] == 0x01 && etp[2] == 0x03 && (event = "artifact:EMG/Muscle")
                etp[1] == 0x01 && etp[2] == 0x04 && (event = "artifact:Movement")
                etp[1] == 0x01 && etp[2] == 0x05 && (event = "artifact:Failing Electrode")
                etp[1] == 0x01 && etp[2] == 0x06 && (event = "artifact:Sweat")
                etp[1] == 0x01 && etp[2] == 0x07 && (event = "artifact:50/60 Hz mains interference")
                etp[1] == 0x01 && etp[2] == 0x08 && (event = "artifact:breathing")
                etp[1] == 0x01 && etp[2] == 0x09 && (event = "artifact:pulse")
                etp[1] == 0x01 && etp[2] == 0x0a && (event = "artifact:EOG (slow)")
                etp[1] == 0x01 && etp[2] == 0x0f && (event = "calibration" )
                etp[1] == 0x01 && etp[2] == 0x11 && (event = "EEG:Sleep spindles")
                etp[1] == 0x01 && etp[2] == 0x12 && (event = "EEG:K-complexes")
                etp[1] == 0x01 && etp[2] == 0x13 && (event = "EEG:Saw-tooth waves")
                etp[1] == 0x01 && etp[2] == 0x14 && (event = "EEG:Idling EEG - eyes open")
                etp[1] == 0x01 && etp[2] == 0x15 && (event = "EEG:Idling EEG - eyes closed")
                etp[1] == 0x01 && etp[2] == 0x16 && (event = "EEG:spike")
                etp[1] == 0x01 && etp[2] == 0x17 && (event = "EEG:seizure")
                etp[1] == 0x01 && etp[2] == 0x18 && (event = "EEG:Electrographic seizure")
                etp[1] == 0x01 && etp[2] == 0x19 && (event = "EEG:Clinical seizure")
                etp[1] == 0x01 && etp[2] == 0x1a && (event = "EEG:Subclinical seizure")
                etp[1] == 0x01 && etp[2] == 0x1b && (event = "EEG:Stimulating for seizure")
                etp[1] == 0x01 && etp[2] == 0x21 && (event = "VEP: visual EP")
                etp[1] == 0x01 && etp[2] == 0x22 && (event = "AEP: auditory EP")
                etp[1] == 0x01 && etp[2] == 0x23 && (event = "SEP: somato-sensory EP")
                etp[1] == 0x01 && etp[2] == 0x2F && (event = "TMS: transcranial magnetic stimulation ")
                etp[1] == 0x01 && etp[2] == 0x31 && (event = "SSVEP")
                etp[1] == 0x01 && etp[2] == 0x32 && (event = "SSAEP")
                etp[1] == 0x01 && etp[2] == 0x33 && (event = "SSSEP")
                etp[1] == 0x01 && etp[2] == 0x40 && (event = "response code 0")
                etp[1] == 0x01 && etp[2] == 0x41 && (event = "response code 1")
                etp[1] == 0x01 && etp[2] == 0x42 && (event = "response code 2")
                etp[1] == 0x01 && etp[2] == 0x43 && (event = "response code 3")
                etp[1] == 0x01 && etp[2] == 0x44 && (event = "Go, or response code 4")
                etp[1] == 0x01 && etp[2] == 0x45 && (event = "NoGo, or response code 5")
                etp[1] == 0x02 && etp[2] == 0x01 && (event = "Spike, action potential (fiducial point)")
                etp[1] == 0x02 && etp[2] == 0x02 && (event = "Burst ")
                etp[1] == 0x02 && etp[2] == 0x03 && (event = "maximum slope time")
                etp[1] == 0x02 && etp[2] == 0x04 && (event = "peak time of spike")
                etp[1] == 0x02 && etp[2] == 0x11 && (event = "EPSP")
                etp[1] == 0x02 && etp[2] == 0x12 && (event = "IPSP")
                etp[1] == 0x02 && etp[2] == 0x13 && (event = "EPSC")
                etp[1] == 0x02 && etp[2] == 0x14 && (event = "IPSC")
                etp[1] == 0x03 && etp[2] == 0x00 && (event = "Start of Trial")
                etp[1] == 0x03 && etp[2] == 0x01 && (event = "class1, Left hand")
                etp[1] == 0x03 && etp[2] == 0x02 && (event = "class2, Right hand")
                etp[1] == 0x03 && etp[2] == 0x03 && (event = "class3, Foot, towards Right")
                etp[1] == 0x03 && etp[2] == 0x04 && (event = "class4, Tongue")
                etp[1] == 0x03 && etp[2] == 0x05 && (event = "class5")
                etp[1] == 0x03 && etp[2] == 0x06 && (event = "class6, towards Down")
                etp[1] == 0x03 && etp[2] == 0x07 && (event = "class7")
                etp[1] == 0x03 && etp[2] == 0x08 && (event = "class8")
                etp[1] == 0x03 && etp[2] == 0x09 && (event = "class9, towards Left")
                etp[1] == 0x03 && etp[2] == 0x0A && (event = "class10")
                etp[1] == 0x03 && etp[2] == 0x0B && (event = "class11")
                etp[1] == 0x03 && etp[2] == 0x0C && (event = "class12, towards Up")
                etp[1] == 0x03 && etp[2] == 0x0D && (event = "Feedback (continuous)")
                etp[1] == 0x03 && etp[2] == 0x0E && (event = "Feedback (discrete)")
                etp[1] == 0x03 && etp[2] == 0x0F && (event = "cue unknown/undefined")
                etp[1] == 0x03 && etp[2] == 0x11 && (event = "Beep")
                etp[1] == 0x03 && etp[2] == 0x12 && (event = "Cross on screen")
                etp[1] == 0x03 && etp[2] == 0x13 && (event = "Flashing light")
                #0x031b - 0x037f reserved for ASCII characters #27-#127
                etp[1] == 0x03 && etp[2] == 0x81 && (event = "target hit, task successful, correct classification")
                etp[1] == 0x03 && etp[2] == 0x82 && (event = "target missed, task not reached, incorrect classification")
                etp[1] == 0x03 && etp[2] == 0xff && (event = "Rejection of whole trial")
                etp[1] == 0x04 && etp[2] == 0x01 && (event = "OAHE")
                etp[1] == 0x04 && etp[2] == 0x02 && (event = "RERA")
                etp[1] == 0x04 && etp[2] == 0x03 && (event = "CAHE")
                etp[1] == 0x04 && etp[2] == 0x04 && (event = "CS Breathing")
                etp[1] == 0x04 && etp[2] == 0x05 && (event = "Hypoventilation ")
                etp[1] == 0x04 && etp[2] == 0x06 && (event = "Apnea")
                etp[1] == 0x04 && etp[2] == 0x07 && (event = "Obstructive apnea")
                etp[1] == 0x04 && etp[2] == 0x08 && (event = "Central apnea")
                etp[1] == 0x04 && etp[2] == 0x09 && (event = "Mixed apnea")
                etp[1] == 0x04 && etp[2] == 0x0A && (event = "Hypopnea  ")
                etp[1] == 0x04 && etp[2] == 0x0B && (event = "Periodic Breathing  ")
                etp[1] == 0x04 && etp[2] == 0x0C && (event = "Limb movement ")
                etp[1] == 0x04 && etp[2] == 0x0D && (event = "PLMS")
                etp[1] == 0x04 && etp[2] == 0x0E && (event = "(time of) maximum inspiration ")
                etp[1] == 0x04 && etp[2] == 0x0F && (event = "Start of inspiration ")
                etp[1] == 0x04 && etp[2] == 0x10 && (event = "Sleep stage Wake")
                etp[1] == 0x04 && etp[2] == 0x11 && (event = "Sleep stage 1")
                etp[1] == 0x04 && etp[2] == 0x12 && (event = "Sleep stage 2")
                etp[1] == 0x04 && etp[2] == 0x13 && (event = "Sleep stage 3")
                etp[1] == 0x04 && etp[2] == 0x14 && (event = "Sleep stage 4")
                etp[1] == 0x04 && etp[2] == 0x15 && (event = "Sleep stage REM")
                etp[1] == 0x04 && etp[2] == 0x16 && (event = "Sleep stage ?")
                etp[1] == 0x04 && etp[2] == 0x17 && (event = "Movement time")
                etp[1] == 0x04 && etp[2] == 0x18 && (event = "Bruxism")
                etp[1] == 0x04 && etp[2] == 0x19 && (event = "RBD")
                etp[1] == 0x04 && etp[2] == 0x1A && (event = "RMD")
                etp[1] == 0x04 && etp[2] == 0x1B && (event = "Sleep stage N")
                etp[1] == 0x04 && etp[2] == 0x1C && (event = "Sleep stage N1")
                etp[1] == 0x04 && etp[2] == 0x1D && (event = "Sleep stage N2")
                etp[1] == 0x04 && etp[2] == 0x1E && (event = "Sleep stage N3")
                etp[1] == 0x04 && etp[2] == 0x20 && (event = "Lights on ")
                etp[1] == 0x84 && etp[2] == 0x20 && (event = "Lights off")
                etp[1] == 0x04 && etp[2] == 0x31 && (event = "eyes left")
                etp[1] == 0x04 && etp[2] == 0x32 && (event = "eyes right")
                etp[1] == 0x04 && etp[2] == 0x33 && (event = "eyes up")
                etp[1] == 0x04 && etp[2] == 0x34 && (event = "eyes down")
                etp[1] == 0x04 && etp[2] == 0x35 && (event = "horizontal eye movement")
                etp[1] == 0x04 && etp[2] == 0x36 && (event = "vertical eye movement")
                etp[1] == 0x04 && etp[2] == 0x37 && (event = "eye rotation (clockwise)")
                etp[1] == 0x04 && etp[2] == 0x38 && (event = "eye rotation (counterclockwise)")
                etp[1] == 0x04 && etp[2] == 0x39 && (event = "eye blinks")
                etp[1] == 0x04 && etp[2] == 0x41 && (event = "left hand movement")
                etp[1] == 0x04 && etp[2] == 0x42 && (event = "right hand movement")
                etp[1] == 0x04 && etp[2] == 0x43 && (event = "head movement")
                etp[1] == 0x04 && etp[2] == 0x44 && (event = "tongue movement")
                etp[1] == 0x04 && etp[2] == 0x45 && (event = "swallowing")
                etp[1] == 0x04 && etp[2] == 0x46 && (event = "biting, chewing, teeth griding")
                etp[1] == 0x04 && etp[2] == 0x47 && (event = "foot movement")
                etp[1] == 0x04 && etp[2] == 0x48 && (event = "foot (right) movement")
                etp[1] == 0x04 && etp[2] == 0x49 && (event = "arm movement")
                etp[1] == 0x04 && etp[2] == 0x4a && (event = "arm (right) movement")
                etp[1] == 0x05 && etp[2] == 0x01 && (event = "ecg:Fiducial point of QRS complex")
                etp[1] == 0x05 && etp[2] == 0x02 && (event = "ecg:P-wave-onset")
                etp[1] == 0x85 && etp[2] == 0x02 && (event = "ecg:P-wave-end")
                etp[1] == 0x05 && etp[2] == 0x03 && (event = "ecg:Q-wave-onset")
                etp[1] == 0x85 && etp[2] == 0x03 && (event = "ecg:Q-wave-peak")
                etp[1] == 0x05 && etp[2] == 0x04 && (event = "ecg:R-point")
                etp[1] == 0x05 && etp[2] == 0x05 && (event = "ecg:S-wave-onset")
                etp[1] == 0x85 && etp[2] == 0x05 && (event = "ecg:S-wave-end")
                etp[1] == 0x05 && etp[2] == 0x06 && (event = "ecg:T-wave-onset")
                etp[1] == 0x85 && etp[2] == 0x06 && (event = "ecg:T-wave-end")
                etp[1] == 0x05 && etp[2] == 0x07 && (event = "ecg:U-wave-onset")
                etp[1] == 0x85 && etp[2] == 0x07 && (event = "ecg:U-wave-end")
                etp[1] == 0x05 && etp[2] == 0x80 && (event = "start")
                etp[1] == 0x05 && etp[2] == 0x81 && (event = "25 Watt")
                etp[1] == 0x05 && etp[2] == 0x82 && (event = "50 Watt")
                etp[1] == 0x05 && etp[2] == 0x83 && (event = "75 Watt")
                etp[1] == 0x05 && etp[2] == 0x84 && (event = "100 Watt")
                etp[1] == 0x05 && etp[2] == 0x85 && (event = "125 Watt")
                etp[1] == 0x05 && etp[2] == 0x86 && (event = "150 Watt")
                etp[1] == 0x05 && etp[2] == 0x87 && (event = "175 Watt")
                etp[1] == 0x05 && etp[2] == 0x88 && (event = "200 Watt")
                etp[1] == 0x05 && etp[2] == 0x89 && (event = "225 Watt")
                etp[1] == 0x05 && etp[2] == 0x8a && (event = "250 Watt")
                etp[1] == 0x05 && etp[2] == 0x8b && (event = "275 Watt")
                etp[1] == 0x05 && etp[2] == 0x8c && (event = "300 Watt")
                etp[1] == 0x05 && etp[2] == 0x8d && (event = "325 Watt")
                etp[1] == 0x05 && etp[2] == 0x8e && (event = "350 Watt")
                etp[1] == 0x85 && etp[2] == 0x80 && (event = "end")
                etp[1] == 0x00 && etp[2] == 0x00 && (event = "empty event")
                etp[1] == 0x00 && etp[2] == 0x01 && (event = "condition 1")
                etp[1] == 0x00 && etp[2] == 0x02 && (event = "condition 2")
                etp[1] == 0x00 && etp[2] == 0x03 && (event = "condition 3")
                etp[1] == 0x00 && etp[2] == 0x04 && (event = "condition 4")
                etp[1] == 0x00 && etp[2] == 0x05 && (event = "condition 5")
                etp[1] == 0x00 && etp[2] == 0x06 && (event = "condition 6")
                etp[1] == 0x00 && etp[2] == 0x07 && (event = "condition 7")
                etp[1] == 0x00 && etp[2] == 0x08 && (event = "condition 8")
                etp[1] == 0x00 && etp[2] == 0x09 && (event = "condition 9")
                etp[1] == 0x00 && etp[2] == 0x0a && (event = "condition 10")
                etp[1] == 0x00 && etp[2] == 0x0b && (event = "condition 11")
                etp[1] == 0x00 && etp[2] == 0x0c && (event = "condition 12")
                etp[1] == 0x00 && etp[2] == 0x0d && (event = "condition 13")
                etp[1] == 0x00 && etp[2] == 0x0e && (event = "condition 14")
                etp[1] == 0x00 && etp[2] == 0x0f && (event = "condition 15")
                etp[1] == 0x00 && etp[2] == 0x10 && (event = "condition 16")
                etp[1] == 0x00 && etp[2] == 0x11 && (event = "condition 17")
                etp[1] == 0x00 && etp[2] == 0x12 && (event = "condition 18")
                etp[1] == 0x00 && etp[2] == 0x13 && (event = "condition 19")
                etp[1] == 0x00 && etp[2] == 0x14 && (event = "condition 20")
                etp[1] == 0x00 && etp[2] == 0x15 && (event = "condition 21")
                etp[1] == 0x00 && etp[2] == 0x16 && (event = "condition 22")
                etp[1] == 0x00 && etp[2] == 0x17 && (event = "condition 23")
                etp[1] == 0x00 && etp[2] == 0x18 && (event = "condition 24")
                etp[1] == 0x00 && etp[2] == 0x19 && (event = "condition 25")
                etp[1] == 0x00 && etp[2] == 0x1a && (event = "condition 26")
                etp[1] == 0x00 && etp[2] == 0x20 && (event = "condition 32")
                etp[1] == 0x00 && etp[2] == 0x2f && (event = "condition 47")
                etp[1] == 0x00 && etp[2] == 0x30 && (event = "condition 48")
                etp[1] == 0x00 && etp[2] == 0x31 && (event = "condition 49")
                etp[1] == 0x00 && etp[2] == 0x32 && (event = "condition 50")
                etp[1] == 0x00 && etp[2] == 0x33 && (event = "condition 51")
                etp[1] == 0x00 && etp[2] == 0x34 && (event = "condition 52")
                etp[1] == 0x00 && etp[2] == 0x35 && (event = "condition 53")
                etp[1] == 0x00 && etp[2] == 0x36 && (event = "condition 54")
                etp[1] == 0x00 && etp[2] == 0x37 && (event = "condition 55")
                etp[1] == 0x00 && etp[2] == 0x38 && (event = "condition 56")
                etp[1] == 0x00 && etp[2] == 0x39 && (event = "condition 57")
                etp[1] == 0x00 && etp[2] == 0x3a && (event = "condition 58")
                etp[1] == 0x00 && etp[2] == 0x3b && (event = "condition 59")
                etp[1] == 0x00 && etp[2] == 0x3c && (event = "condition 60")
                etp[1] == 0x00 && etp[2] == 0x3d && (event = "condition 61")
                etp[1] == 0x00 && etp[2] == 0x3e && (event = "condition 62")
                etp[1] == 0x00 && etp[2] == 0x3f && (event = "condition 63")
                etp[1] == 0x00 && etp[2] == 0x40 && (event = "condition 64")
                etp[1] == 0x00 && etp[2] == 0x41 && (event = "condition 65")
                etp[1] == 0x00 && etp[2] == 0x42 && (event = "condition 66")
                etp[1] == 0x00 && etp[2] == 0x46 && (event = "condition 70")
                etp[1] == 0x00 && etp[2] == 0x51 && (event = "condition 81")
                etp[1] == 0x00 && etp[2] == 0x52 && (event = "condition 82")
                etp[1] == 0x00 && etp[2] == 0x53 && (event = "condition 83")
                etp[1] == 0x00 && etp[2] == 0x5b && (event = "condition 91")
                etp[1] == 0x00 && etp[2] == 0x5c && (event = "condition 92")
                etp[1] == 0x00 && etp[2] == 0x5d && (event = "condition 93")
                etp[1] == 0x00 && etp[2] == 0x60 && (event = "condition 96")
                etp[1] == 0x00 && etp[2] == 0x63 && (event = "condition 99")
                etp[1] == 0x00 && etp[2] == 0x80 && (event = "condition 128")
                etp[1] == 0x00 && etp[2] == 0x81 && (event = "condition 129")
                etp[1] == 0x00 && etp[2] == 0x82 && (event = "condition 130")
                etp[1] == 0x00 && etp[2] == 0x84 && (event = "condition 131")
                etp[1] == 0x00 && etp[2] == 0x85 && (event = "condition 132")
                etp[1] == 0x00 && etp[2] == 0x86 && (event = "condition 133")
                etp[1] == 0x00 && etp[2] == 0x87 && (event = "condition 134")
                etp[1] == 0x00 && etp[2] == 0xa6 && (event = "condition 166")
                etp[1] == 0x00 && etp[2] == 0xa7 && (event = "condition 167")
                etp[1] == 0x00 && etp[2] == 0xa8 && (event = "condition 168")
                etp[1] == 0x00 && etp[2] == 0xa9 && (event = "condition 169")
                etp[1] == 0x7f && etp[2] == 0xfe && (event = "start of a new segment (after a break)")
                etp[1] == 0x7f && etp[2] == 0xff && (event = "non-equidistant sampling value")
                deleteat!(etp, 1:2)
                push!(desc, event)
            end
            if etp_mode == 3
                for idx in 1:etp_number
                    push!(ch, Int64(reinterpret(Int16, etp[1:2])[1]))
                    deleteat!(etp, 1:2)
                end
                for idx in 1:etp_number
                    push!(len, Float64(reinterpret(Float32, etp[1:4])[1]))
                    deleteat!(etp, 1:4)
                end
            end
            if etp_sr == 0
                markers = DataFrame(:id=>id, :start=>round.(start ./ sampling_rate, digits=3), :length=>round.(len ./ sampling_rate, digits=3), :description=>desc, :channel=>ch)
            else
                markers = DataFrame(:id=>id, :start=>round.(start ./ etp_sr, digits=3), :length=>round.(len ./ sampling_rate, digits=3), :description=>desc, :channel=>ch)
            end
        end
    else
        if header_bytes + data_bytes < filesize(file_name)
            fid = ""
            try
                fid = open(file_name, "r")
            catch
                error("File $file_name cannot be loaded.")
            end
            header = UInt8[]
            readbytes!(fid, header, header_bytes + data_bytes)
            etp = UInt8[]
            readbytes!(fid, etp, filesize(file_name))
            etp_mode = reinterpret(Int8, etp[1])[1]
            deleteat!(etp, 1)
            (etp_mode == 1 || etp_mode == 3) || throw(ArgumentError("ETP mode must be 1 or 3, found $etp_mode."))
            b1 = Int32(etp[1]) << 8
            b2 = Int32(etp[2]) << 16
            b3 = -Int32(-etp[3]) << 24
            etp_number = Float64(((b1 | b2 | b3) >> 8))
            deleteat!(etp, 1:3)
            etp_sr = reinterpret(Float32, etp[1:4])[1]
            deleteat!(etp, 1:4)
            start = Float64[]
            len = Float64[]
            id = String[]
            desc = String[]
            ch = Int[]
            for idx in 1:etp_number
                push!(id, "mrk")
                push!(start, Float64(reinterpret(Int32, etp[1:4])[1]))
                push!(len, 0)
                push!(ch, 0)
                deleteat!(etp, 1:4)
            end
            for idx in 1:etp_number
                event = "unknown code ($(etp[1])$(etp[2]))"
                etp[1] == 0x01 && etp[2] == 0x01 && (event = "artifact:EOG (blinks)")
                etp[1] == 0x01 && etp[2] == 0x02 && (event = "artifact:ECG")
                etp[1] == 0x01 && etp[2] == 0x03 && (event = "artifact:EMG/Muscle")
                etp[1] == 0x01 && etp[2] == 0x04 && (event = "artifact:Movement")
                etp[1] == 0x01 && etp[2] == 0x05 && (event = "artifact:Failing Electrode")
                etp[1] == 0x01 && etp[2] == 0x06 && (event = "artifact:Sweat")
                etp[1] == 0x01 && etp[2] == 0x07 && (event = "artifact:50/60 Hz mains interference")
                etp[1] == 0x01 && etp[2] == 0x08 && (event = "artifact:breathing")
                etp[1] == 0x01 && etp[2] == 0x09 && (event = "artifact:pulse")
                etp[1] == 0x01 && etp[2] == 0x0a && (event = "artifact:EOG (slow)")
                etp[1] == 0x01 && etp[2] == 0x0f && (event = "calibration" )
                etp[1] == 0x01 && etp[2] == 0x11 && (event = "EEG:Sleep spindles")
                etp[1] == 0x01 && etp[2] == 0x12 && (event = "EEG:K-complexes")
                etp[1] == 0x01 && etp[2] == 0x13 && (event = "EEG:Saw-tooth waves")
                etp[1] == 0x01 && etp[2] == 0x14 && (event = "EEG:Idling EEG - eyes open")
                etp[1] == 0x01 && etp[2] == 0x15 && (event = "EEG:Idling EEG - eyes closed")
                etp[1] == 0x01 && etp[2] == 0x16 && (event = "EEG:spike")
                etp[1] == 0x01 && etp[2] == 0x17 && (event = "EEG:seizure")
                etp[1] == 0x01 && etp[2] == 0x18 && (event = "EEG:Electrographic seizure")
                etp[1] == 0x01 && etp[2] == 0x19 && (event = "EEG:Clinical seizure")
                etp[1] == 0x01 && etp[2] == 0x1a && (event = "EEG:Subclinical seizure")
                etp[1] == 0x01 && etp[2] == 0x1b && (event = "EEG:Stimulating for seizure")
                etp[1] == 0x01 && etp[2] == 0x21 && (event = "VEP: visual EP")
                etp[1] == 0x01 && etp[2] == 0x22 && (event = "AEP: auditory EP")
                etp[1] == 0x01 && etp[2] == 0x23 && (event = "SEP: somato-sensory EP")
                etp[1] == 0x01 && etp[2] == 0x2F && (event = "TMS: transcranial magnetic stimulation ")
                etp[1] == 0x01 && etp[2] == 0x31 && (event = "SSVEP")
                etp[1] == 0x01 && etp[2] == 0x32 && (event = "SSAEP")
                etp[1] == 0x01 && etp[2] == 0x33 && (event = "SSSEP")
                etp[1] == 0x01 && etp[2] == 0x40 && (event = "response code 0")
                etp[1] == 0x01 && etp[2] == 0x41 && (event = "response code 1")
                etp[1] == 0x01 && etp[2] == 0x42 && (event = "response code 2")
                etp[1] == 0x01 && etp[2] == 0x43 && (event = "response code 3")
                etp[1] == 0x01 && etp[2] == 0x44 && (event = "Go, or response code 4")
                etp[1] == 0x01 && etp[2] == 0x45 && (event = "NoGo, or response code 5")
                etp[1] == 0x02 && etp[2] == 0x01 && (event = "Spike, action potential (fiducial point)")
                etp[1] == 0x02 && etp[2] == 0x02 && (event = "Burst ")
                etp[1] == 0x02 && etp[2] == 0x03 && (event = "maximum slope time")
                etp[1] == 0x02 && etp[2] == 0x04 && (event = "peak time of spike")
                etp[1] == 0x02 && etp[2] == 0x11 && (event = "EPSP")
                etp[1] == 0x02 && etp[2] == 0x12 && (event = "IPSP")
                etp[1] == 0x02 && etp[2] == 0x13 && (event = "EPSC")
                etp[1] == 0x02 && etp[2] == 0x14 && (event = "IPSC")
                etp[1] == 0x03 && etp[2] == 0x00 && (event = "Start of Trial")
                etp[1] == 0x03 && etp[2] == 0x01 && (event = "class1, Left hand")
                etp[1] == 0x03 && etp[2] == 0x02 && (event = "class2, Right hand")
                etp[1] == 0x03 && etp[2] == 0x03 && (event = "class3, Foot, towards Right")
                etp[1] == 0x03 && etp[2] == 0x04 && (event = "class4, Tongue")
                etp[1] == 0x03 && etp[2] == 0x05 && (event = "class5")
                etp[1] == 0x03 && etp[2] == 0x06 && (event = "class6, towards Down")
                etp[1] == 0x03 && etp[2] == 0x07 && (event = "class7")
                etp[1] == 0x03 && etp[2] == 0x08 && (event = "class8")
                etp[1] == 0x03 && etp[2] == 0x09 && (event = "class9, towards Left")
                etp[1] == 0x03 && etp[2] == 0x0A && (event = "class10")
                etp[1] == 0x03 && etp[2] == 0x0B && (event = "class11")
                etp[1] == 0x03 && etp[2] == 0x0C && (event = "class12, towards Up")
                etp[1] == 0x03 && etp[2] == 0x0D && (event = "Feedback (continuous)")
                etp[1] == 0x03 && etp[2] == 0x0E && (event = "Feedback (discrete)")
                etp[1] == 0x03 && etp[2] == 0x0F && (event = "cue unknown/undefined")
                etp[1] == 0x03 && etp[2] == 0x11 && (event = "Beep")
                etp[1] == 0x03 && etp[2] == 0x12 && (event = "Cross on screen")
                etp[1] == 0x03 && etp[2] == 0x13 && (event = "Flashing light")
                #0x031b - 0x037f reserved for ASCII characters #27-#127
                etp[1] == 0x03 && etp[2] == 0x81 && (event = "target hit, task successful, correct classification")
                etp[1] == 0x03 && etp[2] == 0x82 && (event = "target missed, task not reached, incorrect classification")
                etp[1] == 0x03 && etp[2] == 0xff && (event = "Rejection of whole trial")
                etp[1] == 0x04 && etp[2] == 0x01 && (event = "OAHE")
                etp[1] == 0x04 && etp[2] == 0x02 && (event = "RERA")
                etp[1] == 0x04 && etp[2] == 0x03 && (event = "CAHE")
                etp[1] == 0x04 && etp[2] == 0x04 && (event = "CS Breathing")
                etp[1] == 0x04 && etp[2] == 0x05 && (event = "Hypoventilation ")
                etp[1] == 0x04 && etp[2] == 0x06 && (event = "Apnea")
                etp[1] == 0x04 && etp[2] == 0x07 && (event = "Obstructive apnea")
                etp[1] == 0x04 && etp[2] == 0x08 && (event = "Central apnea")
                etp[1] == 0x04 && etp[2] == 0x09 && (event = "Mixed apnea")
                etp[1] == 0x04 && etp[2] == 0x0A && (event = "Hypopnea  ")
                etp[1] == 0x04 && etp[2] == 0x0B && (event = "Periodic Breathing  ")
                etp[1] == 0x04 && etp[2] == 0x0C && (event = "Limb movement ")
                etp[1] == 0x04 && etp[2] == 0x0D && (event = "PLMS")
                etp[1] == 0x04 && etp[2] == 0x0E && (event = "(time of) maximum inspiration ")
                etp[1] == 0x04 && etp[2] == 0x0F && (event = "Start of inspiration ")
                etp[1] == 0x04 && etp[2] == 0x10 && (event = "Sleep stage Wake")
                etp[1] == 0x04 && etp[2] == 0x11 && (event = "Sleep stage 1")
                etp[1] == 0x04 && etp[2] == 0x12 && (event = "Sleep stage 2")
                etp[1] == 0x04 && etp[2] == 0x13 && (event = "Sleep stage 3")
                etp[1] == 0x04 && etp[2] == 0x14 && (event = "Sleep stage 4")
                etp[1] == 0x04 && etp[2] == 0x15 && (event = "Sleep stage REM")
                etp[1] == 0x04 && etp[2] == 0x16 && (event = "Sleep stage ?")
                etp[1] == 0x04 && etp[2] == 0x17 && (event = "Movement time")
                etp[1] == 0x04 && etp[2] == 0x18 && (event = "Bruxism")
                etp[1] == 0x04 && etp[2] == 0x19 && (event = "RBD")
                etp[1] == 0x04 && etp[2] == 0x1A && (event = "RMD")
                etp[1] == 0x04 && etp[2] == 0x1B && (event = "Sleep stage N")
                etp[1] == 0x04 && etp[2] == 0x1C && (event = "Sleep stage N1")
                etp[1] == 0x04 && etp[2] == 0x1D && (event = "Sleep stage N2")
                etp[1] == 0x04 && etp[2] == 0x1E && (event = "Sleep stage N3")
                etp[1] == 0x04 && etp[2] == 0x20 && (event = "Lights on ")
                etp[1] == 0x84 && etp[2] == 0x20 && (event = "Lights off")
                etp[1] == 0x04 && etp[2] == 0x31 && (event = "eyes left")
                etp[1] == 0x04 && etp[2] == 0x32 && (event = "eyes right")
                etp[1] == 0x04 && etp[2] == 0x33 && (event = "eyes up")
                etp[1] == 0x04 && etp[2] == 0x34 && (event = "eyes down")
                etp[1] == 0x04 && etp[2] == 0x35 && (event = "horizontal eye movement")
                etp[1] == 0x04 && etp[2] == 0x36 && (event = "vertical eye movement")
                etp[1] == 0x04 && etp[2] == 0x37 && (event = "eye rotation (clockwise)")
                etp[1] == 0x04 && etp[2] == 0x38 && (event = "eye rotation (counterclockwise)")
                etp[1] == 0x04 && etp[2] == 0x39 && (event = "eye blinks")
                etp[1] == 0x04 && etp[2] == 0x41 && (event = "left hand movement")
                etp[1] == 0x04 && etp[2] == 0x42 && (event = "right hand movement")
                etp[1] == 0x04 && etp[2] == 0x43 && (event = "head movement")
                etp[1] == 0x04 && etp[2] == 0x44 && (event = "tongue movement")
                etp[1] == 0x04 && etp[2] == 0x45 && (event = "swallowing")
                etp[1] == 0x04 && etp[2] == 0x46 && (event = "biting, chewing, teeth griding")
                etp[1] == 0x04 && etp[2] == 0x47 && (event = "foot movement")
                etp[1] == 0x04 && etp[2] == 0x48 && (event = "foot (right) movement")
                etp[1] == 0x04 && etp[2] == 0x49 && (event = "arm movement")
                etp[1] == 0x04 && etp[2] == 0x4a && (event = "arm (right) movement")
                etp[1] == 0x05 && etp[2] == 0x01 && (event = "ecg:Fiducial point of QRS complex")
                etp[1] == 0x05 && etp[2] == 0x02 && (event = "ecg:P-wave-onset")
                etp[1] == 0x85 && etp[2] == 0x02 && (event = "ecg:P-wave-end")
                etp[1] == 0x05 && etp[2] == 0x03 && (event = "ecg:Q-wave-onset")
                etp[1] == 0x85 && etp[2] == 0x03 && (event = "ecg:Q-wave-peak")
                etp[1] == 0x05 && etp[2] == 0x04 && (event = "ecg:R-point")
                etp[1] == 0x05 && etp[2] == 0x05 && (event = "ecg:S-wave-onset")
                etp[1] == 0x85 && etp[2] == 0x05 && (event = "ecg:S-wave-end")
                etp[1] == 0x05 && etp[2] == 0x06 && (event = "ecg:T-wave-onset")
                etp[1] == 0x85 && etp[2] == 0x06 && (event = "ecg:T-wave-end")
                etp[1] == 0x05 && etp[2] == 0x07 && (event = "ecg:U-wave-onset")
                etp[1] == 0x85 && etp[2] == 0x07 && (event = "ecg:U-wave-end")
                etp[1] == 0x05 && etp[2] == 0x80 && (event = "start")
                etp[1] == 0x05 && etp[2] == 0x81 && (event = "25 Watt")
                etp[1] == 0x05 && etp[2] == 0x82 && (event = "50 Watt")
                etp[1] == 0x05 && etp[2] == 0x83 && (event = "75 Watt")
                etp[1] == 0x05 && etp[2] == 0x84 && (event = "100 Watt")
                etp[1] == 0x05 && etp[2] == 0x85 && (event = "125 Watt")
                etp[1] == 0x05 && etp[2] == 0x86 && (event = "150 Watt")
                etp[1] == 0x05 && etp[2] == 0x87 && (event = "175 Watt")
                etp[1] == 0x05 && etp[2] == 0x88 && (event = "200 Watt")
                etp[1] == 0x05 && etp[2] == 0x89 && (event = "225 Watt")
                etp[1] == 0x05 && etp[2] == 0x8a && (event = "250 Watt")
                etp[1] == 0x05 && etp[2] == 0x8b && (event = "275 Watt")
                etp[1] == 0x05 && etp[2] == 0x8c && (event = "300 Watt")
                etp[1] == 0x05 && etp[2] == 0x8d && (event = "325 Watt")
                etp[1] == 0x05 && etp[2] == 0x8e && (event = "350 Watt")
                etp[1] == 0x85 && etp[2] == 0x80 && (event = "end")
                etp[1] == 0x00 && etp[2] == 0x00 && (event = "empty event")
                etp[1] == 0x00 && etp[2] == 0x01 && (event = "condition 1")
                etp[1] == 0x00 && etp[2] == 0x02 && (event = "condition 2")
                etp[1] == 0x00 && etp[2] == 0x03 && (event = "condition 3")
                etp[1] == 0x00 && etp[2] == 0x04 && (event = "condition 4")
                etp[1] == 0x00 && etp[2] == 0x05 && (event = "condition 5")
                etp[1] == 0x00 && etp[2] == 0x06 && (event = "condition 6")
                etp[1] == 0x00 && etp[2] == 0x07 && (event = "condition 7")
                etp[1] == 0x00 && etp[2] == 0x08 && (event = "condition 8")
                etp[1] == 0x00 && etp[2] == 0x09 && (event = "condition 9")
                etp[1] == 0x00 && etp[2] == 0x0a && (event = "condition 10")
                etp[1] == 0x00 && etp[2] == 0x0b && (event = "condition 11")
                etp[1] == 0x00 && etp[2] == 0x0c && (event = "condition 12")
                etp[1] == 0x00 && etp[2] == 0x0d && (event = "condition 13")
                etp[1] == 0x00 && etp[2] == 0x0e && (event = "condition 14")
                etp[1] == 0x00 && etp[2] == 0x0f && (event = "condition 15")
                etp[1] == 0x00 && etp[2] == 0x10 && (event = "condition 16")
                etp[1] == 0x00 && etp[2] == 0x11 && (event = "condition 17")
                etp[1] == 0x00 && etp[2] == 0x12 && (event = "condition 18")
                etp[1] == 0x00 && etp[2] == 0x13 && (event = "condition 19")
                etp[1] == 0x00 && etp[2] == 0x14 && (event = "condition 20")
                etp[1] == 0x00 && etp[2] == 0x15 && (event = "condition 21")
                etp[1] == 0x00 && etp[2] == 0x16 && (event = "condition 22")
                etp[1] == 0x00 && etp[2] == 0x17 && (event = "condition 23")
                etp[1] == 0x00 && etp[2] == 0x18 && (event = "condition 24")
                etp[1] == 0x00 && etp[2] == 0x19 && (event = "condition 25")
                etp[1] == 0x00 && etp[2] == 0x1a && (event = "condition 26")
                etp[1] == 0x00 && etp[2] == 0x20 && (event = "condition 32")
                etp[1] == 0x00 && etp[2] == 0x2f && (event = "condition 47")
                etp[1] == 0x00 && etp[2] == 0x30 && (event = "condition 48")
                etp[1] == 0x00 && etp[2] == 0x31 && (event = "condition 49")
                etp[1] == 0x00 && etp[2] == 0x32 && (event = "condition 50")
                etp[1] == 0x00 && etp[2] == 0x33 && (event = "condition 51")
                etp[1] == 0x00 && etp[2] == 0x34 && (event = "condition 52")
                etp[1] == 0x00 && etp[2] == 0x35 && (event = "condition 53")
                etp[1] == 0x00 && etp[2] == 0x36 && (event = "condition 54")
                etp[1] == 0x00 && etp[2] == 0x37 && (event = "condition 55")
                etp[1] == 0x00 && etp[2] == 0x38 && (event = "condition 56")
                etp[1] == 0x00 && etp[2] == 0x39 && (event = "condition 57")
                etp[1] == 0x00 && etp[2] == 0x3a && (event = "condition 58")
                etp[1] == 0x00 && etp[2] == 0x3b && (event = "condition 59")
                etp[1] == 0x00 && etp[2] == 0x3c && (event = "condition 60")
                etp[1] == 0x00 && etp[2] == 0x3d && (event = "condition 61")
                etp[1] == 0x00 && etp[2] == 0x3e && (event = "condition 62")
                etp[1] == 0x00 && etp[2] == 0x3f && (event = "condition 63")
                etp[1] == 0x00 && etp[2] == 0x40 && (event = "condition 64")
                etp[1] == 0x00 && etp[2] == 0x41 && (event = "condition 65")
                etp[1] == 0x00 && etp[2] == 0x42 && (event = "condition 66")
                etp[1] == 0x00 && etp[2] == 0x46 && (event = "condition 70")
                etp[1] == 0x00 && etp[2] == 0x51 && (event = "condition 81")
                etp[1] == 0x00 && etp[2] == 0x52 && (event = "condition 82")
                etp[1] == 0x00 && etp[2] == 0x53 && (event = "condition 83")
                etp[1] == 0x00 && etp[2] == 0x5b && (event = "condition 91")
                etp[1] == 0x00 && etp[2] == 0x5c && (event = "condition 92")
                etp[1] == 0x00 && etp[2] == 0x5d && (event = "condition 93")
                etp[1] == 0x00 && etp[2] == 0x60 && (event = "condition 96")
                etp[1] == 0x00 && etp[2] == 0x63 && (event = "condition 99")
                etp[1] == 0x00 && etp[2] == 0x80 && (event = "condition 128")
                etp[1] == 0x00 && etp[2] == 0x81 && (event = "condition 129")
                etp[1] == 0x00 && etp[2] == 0x82 && (event = "condition 130")
                etp[1] == 0x00 && etp[2] == 0x84 && (event = "condition 131")
                etp[1] == 0x00 && etp[2] == 0x85 && (event = "condition 132")
                etp[1] == 0x00 && etp[2] == 0x86 && (event = "condition 133")
                etp[1] == 0x00 && etp[2] == 0x87 && (event = "condition 134")
                etp[1] == 0x00 && etp[2] == 0xa6 && (event = "condition 166")
                etp[1] == 0x00 && etp[2] == 0xa7 && (event = "condition 167")
                etp[1] == 0x00 && etp[2] == 0xa8 && (event = "condition 168")
                etp[1] == 0x00 && etp[2] == 0xa9 && (event = "condition 169")
                etp[1] == 0x7f && etp[2] == 0xfe && (event = "start of a new segment (after a break)")
                etp[1] == 0x7f && etp[2] == 0xff && (event = "non-equidistant sampling value")
                deleteat!(etp, 1:2)
                push!(desc, event)
            end
            if etp_mode == 3
                for idx in 1:etp_number
                    push!(ch, Int64(reinterpret(Int16, etp[1:2])[1]))
                    deleteat!(etp, 1:2)
                end
                for idx in 1:etp_number
                    push!(len, Float64(reinterpret(Float32, etp[1:4])[1]))
                    deleteat!(etp, 1:4)
                end
            end
            if etp_sr == 0
                markers = DataFrame(:id=>id, :start=>round.(start ./ sampling_rate, digits=3), :length=>round.(len ./ sampling_rate, digits=3), :description=>desc, :channel=>ch)
            else
                markers = DataFrame(:id=>id, :start=>round.(start ./ etp_sr, digits=3), :length=>round.(len ./ sampling_rate, digits=3), :description=>desc, :channel=>ch)
            end
        end
    end

    clabels = NeuroAnalyzer._clean_labels(clabels)
    if detect_type == true
        ch_type = NeuroAnalyzer._set_channel_types(clabels, "eeg")
    else
        ch_type = repeat(["eeg"], ch_n)
    end
    units = [_set_units(ch_type[idx]) for idx in 1:ch_n]
    ch_order = NeuroAnalyzer._sort_channels(ch_type)
    @show ch_order
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

