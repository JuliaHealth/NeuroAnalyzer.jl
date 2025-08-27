export import_gdf

"""
    import_gdf(file_name; <keyword arguments>)

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
function import_gdf(file_name::String; detect_type::Bool=true)::NeuroAnalyzer.NEURO

    @assert isfile(file_name) "File $file_name cannot be loaded."
    @assert lowercase(splitext(file_name)[2]) == ".gdf" "This is not GDF file."

    file_type = ""

    fid = nothing
    try
        fid = open(file_name, "r")
    catch
        @error "File $file_name cannot be loaded."
    end

    header = UInt8[]
    readbytes!(fid, header, 256)
    file_type = String(Char.(header[1:3]))
    file_type_ver = parse(Float64, String(Char.(header[5:8])))
    @assert file_type == "GDF" "File $file_name is not GDF file."

    (file_type_ver == 1.25 || file_type_ver == 2.20) || _warn("GDF versions other than 1.25 and 2.20 may not be supported correctly.")

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
        [push!(clabels, _v2s(_fread(fid, 16, :c))) for idx in 1:ch_n]

        transducers = String[]
        [push!(transducers, _v2s(_fread(fid, 80, :c))) for idx in 1:ch_n]

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
        @assert data_records != -1 "Number of data records cannot be -1."
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
        for _ in 1:ch_n
            readbytes!(fid, buf, 8)
            push!(physical_minimum, reinterpret(Float64, buf)[1])
        end

        physical_maximum = Float64[]
        buf = UInt8[]
        for _ in 1:ch_n
            readbytes!(fid, buf, 8)
            push!(physical_maximum, reinterpret(Float64, buf)[1])
        end

        digital_minimum = Float64[]
        buf = UInt8[]
        for _ in 1:ch_n
            readbytes!(fid, buf, 8)
            push!(digital_minimum, reinterpret(Float64, buf)[1])
        end

        digital_maximum = Float64[]
        buf = UInt8[]
        for _ in 1:ch_n
            readbytes!(fid, buf, 8)
            push!(digital_maximum, reinterpret(Float64, buf)[1])
        end

        # obsolete
        prefiltering = String[]
        buf = UInt8[]
        for _ in 1:ch_n
            readbytes!(fid, buf, 64)
            push!(prefiltering, replace(strip(String(Char.(buf))), '\0'=>""))
        end

        time_offset = Float32[]
        for _ in 1:ch_n
            readbytes!(fid, buf, 4)
            push!(time_offset, reinterpret(Int32, buf)[1])
        end
        time_offset == Float32[] && (time_offset = zeros(Int64, ch_n))

        prefiltering_lp = Float32[]
        buf = UInt8[]
        for _ in 1:ch_n
            readbytes!(fid, buf, 4)
            push!(prefiltering_lp, reinterpret(Float32, buf)[1])
        end
        prefiltering_lp[findall(isnan, prefiltering_lp)] .= 0

        prefiltering_hp = Float32[]
        buf = UInt8[]
        for _ in 1:ch_n
            readbytes!(fid, buf, 4)
            push!(prefiltering_hp, reinterpret(Float32, buf)[1])
        end
        prefiltering_hp[findall(isnan, prefiltering_hp)] .= 0

        prefiltering_bs = Float32[]
        buf = UInt8[]
        for _ in 1:ch_n
            readbytes!(fid, buf, 4)
            push!(prefiltering_bs, reinterpret(Float32, buf)[1])
        end
        prefiltering_bs[findall(isnan, prefiltering_bs)] .= 0

        samples_per_datarecord = Int32[]
        buf = UInt8[]
        for _ in 1:ch_n
            readbytes!(fid, buf, 4)
            push!(samples_per_datarecord, reinterpret(Int32, buf)[1])
        end

        gdf_type = Int32[]
        buf = UInt8[]
        for _ in 1:ch_n
            readbytes!(fid, buf, 4)
            push!(gdf_type, reinterpret(Int32, buf)[1])
        end

        loc_x = Float32[]
        loc_y = Float32[]
        loc_z = Float32[]
        buf = UInt8[]
        for _ in 1:ch_n
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
            for _ in 1:ch_n
                readbytes!(fid, buf, 1)
                push!(imp, reinterpret(Int32, buf)[1])
            end
            buf = UInt8[]
            readbytes!(fid, buf, 19 * ch_n)
            imp = @. Float64(2 ^ (imp / 8))
        end

    end

    close(fid)

    fid = nothing
    try
        fid = open(file_name, "r")
    catch
        error("File $file_name cannot be loaded.")
    end
    seek(fid, header_bytes)

    # check TLV header
    if file_type_ver >= 2.10 && (etv_hdr - (ch_n + 1)) * 256 > 0
        tlv_tag = UInt8[]
        tlv_len = UInt8[]
        tlv_val = UInt8[]
        readbytes!(fid, tlv_tag, 1)
        if tlv_tag != [0x00]
            @error "TLV are not supported yet; if you have such a file, please send it to adam.wysokinski@neuroanalyzer.org"
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
    data = reshape(data, size(data, 1), size(data, 2), 1)

    gain = @. (physical_maximum - physical_minimum) / (digital_maximum - digital_minimum)
    data .*= gain

    markers = DataFrame(:id=>String[],
                        :start=>Float64[],
                        :length=>Float64[],
                        :value=>String[],
                        :channel=>Int64[])

    if file_type_ver < 2.0
        if header_bytes + data_bytes < filesize(file_name)
            fid = nothing
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
            @assert (etp_mode == 1 || etp_mode == 3) "ETP mode must be 1 or 3, found $etp_mode."
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
            value = String[]
            ch = Int[]
            for idx in 1:etp_number
                push!(id, "event")
                push!(start, Float64(reinterpret(Int32, etp[1:4])[1]))
                push!(len, 0)
                push!(ch, 0)
                deleteat!(etp, 1:4)
            end
            for idx in 1:etp_number
                event = _gdf_etp(etp[1:2])
                deleteat!(etp, 1:2)
                push!(value, event)
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
                markers = DataFrame(:id=>id,
                                    :start=>round.(start ./ sampling_rate, digits=3),
                                    :length=>round.(len ./ sampling_rate, digits=3),
                                    :value=>value,
                                    :channel=>ch)
            else
                markers = DataFrame(:id=>id,
                                    :start=>round.(start ./ etp_sr, digits=3),
                                    :length=>round.(len ./ sampling_rate, digits=3),
                                    :value=>value,
                                    :channel=>ch)
            end
        end
    else
        if header_bytes + data_bytes < filesize(file_name)
            fid = nothing
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
            @assert (etp_mode == 1 || etp_mode == 3) "ETP mode must be 1 or 3, found $etp_mode."
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
            value = String[]
            ch = Int[]
            for idx in 1:etp_number
                push!(id, "event")
                push!(start, Float64(reinterpret(Int32, etp[1:4])[1]))
                push!(len, 0)
                push!(ch, 0)
                deleteat!(etp, 1:4)
            end
            for idx in 1:etp_number
                event = _gdf_etp(etp[1:2])
                deleteat!(etp, 1:2)
                push!(value, event)
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
                markers = DataFrame(:id=>id,
                                    :start=>round.(start ./ sampling_rate, digits=3),
                                    :length=>round.(len ./ sampling_rate, digits=3),
                                    :value=>value,
                                    :channel=>ch)
            else
                markers = DataFrame(:id=>id,
                                    :start=>round.(start ./ etp_sr, digits=3),
                                    :length=>round.(len ./ sampling_rate, digits=3),
                                    :value=>value,
                                    :channel=>ch)
            end
        end
    end

    clabels = _clean_labels(clabels)
    if detect_type
        ch_type = _set_channel_types(clabels, "eeg")
    else
        ch_type = repeat(["eeg"], ch_n)
    end
    units = [_ch_units(ch_type[idx]) for idx in 1:ch_n]
    time_pts = round.(collect(0:1/sampling_rate:size(data, 2) * size(data, 3) / sampling_rate)[1:end-1], digits=3)
    ep_time = round.((collect(0:1/sampling_rate:size(data, 2) / sampling_rate))[1:end-1], digits=3)

    file_size_mb = round(filesize(file_name) / 1024^2, digits=2)

    "eeg" in ch_type && (data_type = "eeg")
    "meg" in ch_type && (data_type = "meg")
    "mag" in ch_type && (data_type = "meg")
    "grad" in ch_type && (data_type = "meg")

    s = _create_subject(id="",
                        first_name="",
                        middle_name="",
                        last_name=patient,
                        head_circumference=-1,
                        handedness="",
                        weight=-1,
                        height=-1)
    r = _create_recording_eeg(data_type=data_type,
                              file_name=file_name,
                              file_size_mb=file_size_mb,
                              file_type=file_type,
                              recording=recording,
                              recording_date=recording_date,
                              recording_time=replace(recording_time, '.'=>':'),
                              recording_notes="",
                              channel_type=ch_type,
                              channel_order=_sort_channels(ch_type),
                              reference=_detect_montage(clabels, ch_type, data_type),
                              clabels=clabels,
                              transducers=transducers,
                              units=units,
                              prefiltering=prefiltering,
                              line_frequency=50,
                              sampling_rate=sampling_rate,
                              gain=gain,
                              bad_channels=zeros(Bool, size(data, 1), 1))
    e = _create_experiment(name="", notes="", design="")

    hdr = _create_header(s,
                         r,
                         e)

    components = Dict()

    history = String[]

    locs = _initialize_locs()
    obj = NeuroAnalyzer.NEURO(hdr, time_pts, ep_time, data, components, markers, locs, history)
    _initialize_locs!(obj)

    _info("Imported: " * uppercase(obj.header.recording[:data_type]) * " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits=3)) s)")

    return obj

end
