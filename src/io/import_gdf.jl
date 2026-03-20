export import_gdf

"""
    import_gdf(file_name; <keyword arguments>)

Load a GDF (General Dataformat for Biosignals) file and return a `NeuroAnalyzer.NEURO` object.

Both GDF 1.x and GDF 2.x variants are supported. Channels may have mixed data types; all are converted to `Float64` on import. The event table (ETP block) is parsed into a markers `DataFrame` when present.

# Arguments

- `file_name::String`: path to the `.gdf` file
- `detect_type::Bool=true`: infer channel type from channel label

# Returns

- `NeuroAnalyzer.NEURO`

# Throws

- `ArgumentError` if the file does not exist, is not a GDF file, or contains unsupported structural features

# Notes

- `sampling_rate = samples_per_datarecord ÷ data_record_duration`
- `gain = (physical_max − physical_min) / (digital_max − digital_min)`
- `value = (raw − digital_min) × gain + physical_min`

# References

1. Schlögl A et al. GDF v1.25. 1998.
2. Schlögl A. GDF v2.12. 2009.
3. Schlögl A. GDF v2.51. 2013.
"""
function import_gdf(
    file_name::String;
    detect_type::Bool = true
)::NeuroAnalyzer.NEURO

    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))
    lowercase(splitext(file_name)[2]) == ".gdf" ||
        throw(ArgumentError("$file_name is not a GDF file."))

    # --------------------------------------------------------------------------- #
    # all reads share open/close via `do` blocks - guaranteed close on exception  #
    # --------------------------------------------------------------------------- #

    # ------------------------------------------------------------------ #
    # global 256-byte header                                              #
    # ------------------------------------------------------------------ #
    header = open(file_name, "r") do fid
        buf = UInt8[]
        readbytes!(fid, buf, 256)
        buf
    end

    file_type = String(Char.(header[1:3]))
    file_type_ver = parse(Float64, String(Char.(header[5:8])))
    file_type == "GDF" ||
        throw(ArgumentError("$file_name is not a GDF file."))
    (file_type_ver == 1.25 || file_type_ver == 2.2) ||
        _warn("GDF versions other than 1.25 and 2.20 may not be fully supported.")

    # ------------------------------------------------------------------ #
    # version-specific per-channel header + channel data                 #
    # ------------------------------------------------------------------ #

    # Helper: read n_ch fixed-width string fields, stripping NUL bytes
    read_str_fields(fid, n, width) = [
        replace(strip(String(Char.(let b = UInt8[]; readbytes!(fid, b, width); b end))),
                '\0' => "")
        for _ in 1:n]

    # helper: read n_ch fixed-width binary fields with a given reinterpret type
    read_bin_fields(fid, n, width, T) = [
        (let b = UInt8[]; readbytes!(fid, b, width); reinterpret(T, b)[1] end)
        for _ in 1:n]

    (patient, recording, recording_date, recording_time,
     header_bytes, data_records, sampling_rate, ch_n,
     clabels, transducers, units, prefiltering,
     physical_minimum, physical_maximum, digital_minimum, digital_maximum,
     samples_per_datarecord, gdf_type,
     etv_hdr) = open(file_name, "r") do fid

        # skip the already-read global header.
        seek(fid, 256)

        if file_type_ver < 2.0
            # -------------------------------------------------------- #
            # GDF 1.x per-channel header                               #
            # -------------------------------------------------------- #
            patient = replace(strip(String(Char.(header[9:88]))),   '\0' => "")
            recording = replace(strip(String(Char.(header[89:168]))), '\0' => "")
            recording_date = replace(strip(String(Char.(header[169:184]))),'\0' => "")
            recording_time = recording_date == "" ? "" :
                replace(strip(String(Char.(header[169:184]))), '\0' => "")
            header_bytes = reinterpret(Int64, header[185:192])[1]
            data_records = reinterpret(Int64, header[237:244])[1]
            sampling_rate = Int64(reinterpret(Int32, header[245:252])[2] ÷
                                  reinterpret(Int32, header[245:252])[1])
            ch_n = reinterpret(Int32, header[253:256])[1]

            clabels = read_str_fields(fid, ch_n, 16)
            transducers = read_str_fields(fid, ch_n, 80)

            units = String[]
            for _ in 1:ch_n
                buf = UInt8[]; readbytes!(fid, buf, 8)
                push!(units, replace(strip(String(Char.(buf))), '\0' => "", '\x10' => ""))
            end
            units = replace(lowercase.(units), "uv" => "μV")

            physical_minimum = read_bin_fields(fid, ch_n, 8, Float64)
            physical_maximum = read_bin_fields(fid, ch_n, 8, Float64)
            digital_minimum  = Float64.(read_bin_fields(fid, ch_n, 8, Int64))
            digital_maximum  = Float64.(read_bin_fields(fid, ch_n, 8, Int64))
            prefiltering = read_str_fields(fid, ch_n, 80)
            samples_per_datarecord = read_bin_fields(fid, ch_n, 4, Int32)
            gdf_type = read_bin_fields(fid, ch_n, 4, Int32)
            etv_hdr = 0

        else
            # -------------------------------------------------------- #
            # GDF 2.x per-channel header                               #
            # -------------------------------------------------------- #
            patient = replace(strip(String(Char.(header[9:66]))),    '\0' => "")
            patient_weight = header[86] == 0 ? -1 : Int(header[86])
            patient_height = header[87] == 0 ? -1 : Int(header[87])
            recording = replace(strip(String(Char.(header[89:152]))), '\0' => "")
            recording_date = replace(strip(String(Char.(header[169:176]))),'\0' => "")
            recording_time = recording_date == "" ? "" : recording_date
            header_bytes = reinterpret(Int16, header[185:186])[1] * 256
            etv_hdr = reinterpret(Int16, header[185:186])[1]
            data_records = reinterpret(Int64, header[237:244])[1]
            data_records == -1 &&
                throw(ArgumentError("Number of data records cannot be -1."))
            sampling_rate = Int64(reinterpret(Int32, header[245:252])[2] ÷
                                  reinterpret(Int32, header[245:252])[1])
            ch_n = reinterpret(Int16, header[253:254])[1]

            clabels = read_str_fields(fid, ch_n, 16)
            transducers = read_str_fields(fid, ch_n, 80)

            # obsolete physical units (6 bytes each) - consumed to advance position
            for _ in 1:ch_n
                buf = UInt8[]; readbytes!(fid, buf, 6)
            end

            # active unit codes (2 bytes each, bit-encoded SI prefix + unit)
            units_code = read_bin_fields(fid, ch_n, 2, UInt16)
            unit_val = zeros(Int64, ch_n)
            units = Vector{String}(undef, ch_n)
            unit_map = Dict(
                512 => "",
                544 => "%",
                736 => "°",
                768 => "rad",
                2496 => "Hz",
                3872 => "mmHg",
                4256 => "V",
                4288 => "Ω",
                6048 => "°C",
                3072 => "l/min",
                2848 => "l/(min m²)",
                4128 => "dyn s/cm⁵",
                6016 => "dyn s/m² cm⁵")
            prefix_map = Dict(
                10 => "Y",
                9 => "Z",
                8 => "E",
                7 => "P",
                6 => "T",
                5 => "G",
                4 => "M",
                3 => "k",
                2 => "h",
                1 => "da",
                0 => "",
                16 => "d",
                17 => "c",
                18 => "m",
                19 => "μ",
                20 => "n",
                21 => "p",
                22 => "f",
                23 => "a",
                24 => "z",
                25 => "y")
            for idx in 1:ch_n
                base_code = parse(Int, "0b" * bitstring(units_code[idx])[1:(end-5)] * "00000")
                dec_factor = parse(Int, "0b" * bitstring(units_code[idx])[(end-4):end])
                unit_val[idx] = base_code
                base_unit = get(unit_map, base_code, "")
                prefix = get(prefix_map, dec_factor, "")
                units[idx] = prefix * base_unit
            end

            physical_minimum = read_bin_fields(fid, ch_n, 8, Float64)
            physical_maximum = read_bin_fields(fid, ch_n, 8, Float64)
            digital_minimum = read_bin_fields(fid, ch_n, 8, Float64)
            digital_maximum = read_bin_fields(fid, ch_n, 8, Float64)

            # obsolete prefiltering (64 bytes each)
            prefiltering = read_str_fields(fid, ch_n, 64)

            # per-channel time offsets, LP/HP/BS filter frequencies
            time_offset = read_bin_fields(fid, ch_n, 4, Float32)
            prefiltering_lp = let v = read_bin_fields(fid, ch_n, 4, Float32)
                v[isnan.(v)] .= 0; v end
            prefiltering_hp = let v = read_bin_fields(fid, ch_n, 4, Float32)
                v[isnan.(v)] .= 0; v end
            prefiltering_bs = let v = read_bin_fields(fid, ch_n, 4, Float32)
                v[isnan.(v)] .= 0; v end

            # combine LP/HP into prefiltering strings (overrides the obsolete field)
            prefiltering = ["LP: $(prefiltering_lp[i]) Hz; HP: $(prefiltering_hp[i]) Hz"
                            for i in 1:ch_n]

            samples_per_datarecord = read_bin_fields(fid, ch_n, 4, Int32)
            gdf_type = read_bin_fields(fid, ch_n, 4, Int32)

            # 3D electrode locations
            loc_x = Float32[]; loc_y = Float32[]; loc_z = Float32[]
            for _ in 1:ch_n
                buf = UInt8[]
                readbytes!(fid, buf, 4); push!(loc_x, reinterpret(Float32, buf)[1])
                readbytes!(fid, buf, 4); push!(loc_y, reinterpret(Float32, buf)[1])
                readbytes!(fid, buf, 4); push!(loc_z, reinterpret(Float32, buf)[1])
            end
            loc_x[isnan.(loc_x)] .= 0; loc_y[isnan.(loc_y)] .= 0; loc_z[isnan.(loc_z)] .= 0

            # impedances (format differs between GDF 2.19+ and earlier)
            if file_type_ver >= 2.19
                imp = zeros(Float64, ch_n)
                for idx in 1:ch_n
                    buf = UInt8[]; readbytes!(fid, buf, 20)
                    (unit_val[idx] == 4256 || unit_val[idx] == 4288) &&
                        (imp[idx] = reinterpret(Float32, buf[1:4])[1])
                end
                imp[isnan.(imp)] .= 0
            else
                # GDF 1.x encodes impedance as a single byte; scale as 2^(byte/8)
                imp = zeros(Float64, ch_n)
                for idx in 1:ch_n
                    buf = UInt8[]; readbytes!(fid, buf, 1)
                    imp[idx] = Float64(2^(buf[1] / 8))
                end
                buf = UInt8[]; readbytes!(fid, buf, 19 * ch_n) # skip padding
            end
        end

        (patient, recording, recording_date, recording_time,
         header_bytes, data_records, sampling_rate, ch_n,
         clabels, transducers, units, prefiltering,
         physical_minimum, physical_maximum, digital_minimum, digital_maximum,
         samples_per_datarecord, gdf_type, etv_hdr)

    end
    # per-channel header file closed here

    # ------------------------------------------------------------------ #
    # signal data                                                        #
    # ------------------------------------------------------------------ #
    # GDF stores data in record-interleaved order:
    #   [ch1_rec1_s1..sN, ch2_rec1_s1..sN, ..., chK_rec1_s1..sN,
    #    ch1_rec2_s1..sN, ...]
    # Read each channel block sequentially then reassemble.

    # map from GDF type code → (bytes_per_sample, Julia type)
    gdf_type_map = Dict(
        0 => (1, UInt8),
        1 => (1, Int8),
        2 => (1, UInt8),
        3 => (2, Int16),
        4 => (2, UInt16),
        5 => (4, Int32),
        6 => (4, UInt32),
        7 => (8, Int64),
        16 => (4, Float32),
        17 => (8, Float64))

    # one signal buffer per channel (n_samples = spdr × data_records)
    ch_signals = Vector{Vector{Float64}}(undef, ch_n)

    open(file_name, "r") do fid
        seek(fid, header_bytes)

        # skip TLV block if present (GDF 2.1+)
        if file_type_ver >= 2.1 && (etv_hdr - (ch_n + 1)) * 256 > 0
            tlv_tag = UInt8[]; readbytes!(fid, tlv_tag, 1)
            if tlv_tag != [0x00]
                _warn("TLV not supported; please send this file to adam.wysokinski@neuroanalyzer.org")
                tlv_len_buf = UInt8[]; readbytes!(fid, tlv_len_buf, 3)
                b1 = Int32(tlv_len_buf[1]) << 8
                b2 = Int32(tlv_len_buf[2]) << 16
                b3 = -Int32(-tlv_len_buf[3]) << 24
                tlv_len = Int(Float64(((b1 | b2 | b3) >> 8)))
                skip(fid, tlv_len)
            else
                skip(fid, 255)
            end
        end

        # read each channel's raw bytes for all records at once.
        for ch in 1:ch_n
            n_samp  = Int(samples_per_datarecord[ch]) * data_records
            type_id = gdf_type[ch]
            if !haskey(gdf_type_map, type_id)
                throw(ArgumentError("Unknown GDF data type code $(type_id) for channel $ch."))
            end
            bps, T  = gdf_type_map[type_id]
            buf     = UInt8[]; readbytes!(fid, buf, bps * n_samp)
            ch_signals[ch] = Float64.(reinterpret(T, buf))
        end
    end
    # signal file closed here

    # reassemble: GDF interleaves samples across channels within each record.
    n_samp_ch1 = Int(samples_per_datarecord[1]) * data_records
    data = zeros(ch_n, n_samp_ch1, 1)
    for ch in 1:ch_n
        data[ch, :, 1] = ch_signals[ch]
    end

    gain = @. (physical_maximum - physical_minimum) / (digital_maximum - digital_minimum)
    data .*= gain

    # ------------------------------------------------------------------ #
    # event table (ETP block)                                            #
    # ------------------------------------------------------------------ #
    markers = DataFrame(
        :id => String[],
        :start => Float64[],
        :length => Float64[],
        :value => String[],
        :channel => Int64[])

    data_bytes = sum(Int.(samples_per_datarecord) .* data_records .*
        [get(gdf_type_map, Int(t), (1,))[1] for t in gdf_type])

    if header_bytes + data_bytes < filesize(file_name)
        open(file_name, "r") do fid
            skip(fid, header_bytes + data_bytes)
            etp = UInt8[]
            readbytes!(fid, etp, filesize(file_name) - header_bytes - data_bytes; all = true)

            etp_mode = reinterpret(Int8, etp[1:1])[1]
            etp = etp[2:end]
            (etp_mode == 1 || etp_mode == 3) ||
                throw(ArgumentError("ETP mode must be 1 or 3, found $etp_mode."))

            b1 = Int32(etp[1]) << 8
            b2 = Int32(etp[2]) << 16
            b3 = -Int32(-etp[3]) << 24
            etp = etp[4:end]

            if file_type_ver < 2.0
                # GDF1: [3-byte etp_sr] then [Int32 etp_number]
                etp_sr = Float64((b1 | b2 | b3) >> 8)
                etp_number = reinterpret(Int32, etp[1:4])[1]
                etp = etp[5:end]
            else
                # GDF2: [3-byte etp_number] then [Float32 etp_sr]
                etp_number = Int((b1 | b2 | b3) >> 8)
                etp_sr = reinterpret(Float32, etp[1:4])[1]
                etp = etp[5:end]
            end

            start = Float64[]; len = Float64[]
            id = String[];  value = String[]; ch = Int[]

            for _ in 1:etp_number
                push!(id, "event")
                push!(start, Float64(reinterpret(Int32, etp[1:4])[1]))
                push!(len, 0.0);  push!(ch, 0)
                etp = etp[5:end]
            end
            for _ in 1:etp_number
                push!(value, _gdf_etp(etp[1:2]))
                etp = etp[3:end]
            end
            if etp_mode == 3
                for _ in 1:etp_number
                    push!(ch, Int64(reinterpret(Int16, etp[1:2])[1]))
                    etp = etp[3:end]
                end
                for _ in 1:etp_number
                    push!(len, Float64(reinterpret(Float32, etp[1:4])[1]))
                    etp = etp[5:end]
                end
            end

            sr_denom = etp_sr == 0 ? Float64(sampling_rate) : Float64(etp_sr)
            markers = DataFrame(
                :id => id,
                :start => round.(start ./ sr_denom, digits = 4),
                :length => round.(len ./ Float64(sampling_rate), digits = 4),
                :value => value,
                :channel => ch[1:etp_number])
        end
        # event file closed here

    end

    # ------------------------------------------------------------------ #
    # channel metadata and NEURO object assembly                         #
    # ------------------------------------------------------------------ #
    clabels = _clean_labels(string.(clabels))
    ch_type = detect_type ? _set_channel_types(clabels, "eeg") : repeat(["eeg"], ch_n)
    units   = [_ch_units(ch_type[idx]) for idx in 1:ch_n]

    n_samples  = size(data, 2) * size(data, 3)
    time_pts   = round.(range(0, step = 1/sampling_rate, length = n_samples), digits = 4)
    epoch_time = round.(range(0, step = 1/sampling_rate, length = size(data,2)), digits = 4)

    data_type = "eeg"
    "meg" in ch_type && (data_type = "meg")
    "mag" in ch_type && (data_type = "meg")
    "grad" in ch_type && (data_type = "meg")

    # ------------------------------------------------------------------ #
    # assemble NEURO object                                              #
    # ------------------------------------------------------------------ #
    file_size_mb = round(filesize(file_name) / 1024^2, digits = 2)

    s = _create_subject(
        id = "",
        first_name = "",
        middle_name = "",
        last_name = patient,
        head_circumference = -1,
        handedness = "",
        weight = -1,
        height = -1)
    r = _create_recording_eeg(
        data_type = data_type,
        file_name = file_name,
        file_size_mb = file_size_mb,
        file_type = file_type,
        recording = recording,
        recording_date = recording_date,
        recording_time = replace(recording_time, '.' => ':'),
        recording_notes = "",
        channel_type = ch_type,
        channel_order = _sort_channels(ch_type),
        reference = _detect_montage(clabels, ch_type, data_type),
        clabels = clabels,
        transducers = transducers,
        units = units,
        prefiltering = prefiltering,
        line_frequency = 50, # TODO: make this a keyword argument
        sampling_rate = sampling_rate,
        gain = gain,
        bad_channels = zeros(Bool, size(data, 1)))
    e   = _create_experiment(name = "", notes = "", design = "")
    hdr = _create_header(subject = s, recording = r, experiment = e)

    locs = _initialize_locs()
    obj  = NeuroAnalyzer.NEURO(hdr, String[], markers, locs, time_pts, epoch_time, data)
    _initialize_locs!(obj)

    _info("Imported: " *
        uppercase(obj.header.recording[:data_type]) *
        " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj))" *
        "; $(round(obj.time_pts[end], digits=2)) s)")

    return obj

end