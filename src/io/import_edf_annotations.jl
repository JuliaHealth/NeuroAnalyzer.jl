export import_edf_annotations

"""
    import_edf_annotations(file_name; <keyword arguments>)

Load annotations from an EDF+ annotation-only file and return a `DataFrame` of event markers.

This function is intended for EDF+ files whose `data_records_duration` header field is 0, meaning the file contains only TAL (Time-stamped Annotations Lists) and no signal data. For regular EDF/EDF+ files with signal data, use `import_edf()` instead — it parses annotations automatically.

# Arguments

- `file_name::String`: path to the EDF+ annotation file

# Returns
- `DataFrame` with columns `:id`, `:start`, `:length`, `:value`, `:channel`

# Throws
- `ArgumentError` if the file does not exist, is not an EDF file, or is a regular EDF/EDF+ file with signal data (`data_records_duration ≠ 0`)
"""
function import_edf_annotations(file_name::String)::DataFrame

    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))
    lowercase(splitext(file_name)[2]) == ".edf" ||
        throw(ArgumentError("$file_name is not an EDF file."))

    # ------------------------------------------------------------------ #
    # parse header — single open/close via `do` block                    #
    # ------------------------------------------------------------------ #
    markers = open(file_name, "r") do fid

        buf = zeros(UInt8, 256)
        readbytes!(fid, buf, 256)
        header = String(Char.(buf))

        # bytes 1–8: version; must be 0 for EDF
        parse(Int, strip(header[1:8])) == 0 ||
            throw(ArgumentError("$file_name is not a valid EDF file (version ≠ 0)."))
        file_type = "EDF"

        patient = strip(header[9:88])
        recording = strip(header[89:168])

        # Alice 4 files are handled by a dedicated importer.
        occursin("Alice 4", recording) &&
            return _a2df(String[]) # Alice 4 annotation-only files: return empty

        data_offset = parse(Int, strip(header[185:192]))
        reserved = strip(header[193:236])

        reserved == "EDF+D" &&
            throw(ArgumentError(
                "EDF+D (interrupted recordings) is not supported. " *
                "Please send this file to adam.wysokinski@neuroanalyzer.org"))
        reserved == "EDF+C" && (file_type = "EDF+")

        data_records = parse(Int, strip(header[237:244]))
        data_records_duration = parse(Float64, strip(header[245:252]))

        data_records_duration != 0 &&
            throw(ArgumentError(
                "$file_name is a regular $file_type file with signal data; " *
                "use import_edf() instead."))

        ch_n = parse(Int, strip(header[253:256]))

        # ------------------------------------------------------------ #
        # per-channel header fields                                    #
        # ------------------------------------------------------------ #
        read_fields(width, parse_fn = identity) = begin
            buf = zeros(UInt8, ch_n * width)
            readbytes!(fid, buf, ch_n * width)
            s = String(Char.(buf))
            [parse_fn(strip(s[(1 + (i-1)*width):(i*width)])) for i in 1:ch_n]
        end

        clabels = read_fields(16)
        _transducers = read_fields(80) # not used; consumed to advance position
        _units = read_fields(8)
        physical_minimum = read_fields(8, s -> parse(Float64, s))
        physical_maximum = read_fields(8, s -> parse(Float64, s))
        digital_minimum = read_fields(8, s -> parse(Float64, s))
        digital_maximum = read_fields(8, s -> parse(Float64, s))
        _prefiltering = read_fields(80)
        samples_per_datarecord = read_fields(8,  s -> parse(Int, s))

        # ------------------------------------------------------------ #
        # identify annotation channels                                 #
        # annotation-only EDFs have no EDF+ reserved field; treat ALL  #
        # channels as annotation channels                              #
        # ------------------------------------------------------------ #
        annotation_channels = if file_type == "EDF+"
            sort(getindex.(findall(occursin.("annotation", lowercase.(clabels))), 1))
        else
            # annotation-only plain EDF: every channel is TAL
            collect(1:ch_n)
        end

        # ------------------------------------------------------------ #
        # read raw annotation bytes                                    #
        # ------------------------------------------------------------ #
        seek(fid, data_offset)
        annotations = String[]

        for _ in 1:data_records, ch in 1:ch_n
            raw = zeros(UInt8, samples_per_datarecord[ch] * 2)
            readbytes!(fid, raw, samples_per_datarecord[ch] * 2)
            ch in annotation_channels && push!(annotations, String(Char.(raw)))
        end

        isempty(annotations) ?
            DataFrame(:id => String[], :start => Float64[],
                      :length => Float64[], :value => String[], :channel => Int64[]) :
            _a2df(annotations)
    end
    # file closed here in all cases, including exceptions

    return markers

end
