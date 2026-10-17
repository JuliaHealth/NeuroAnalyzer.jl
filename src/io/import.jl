export import_recording

# ---------------------------------------------------------------------------
# format detection
# ---------------------------------------------------------------------------

# all formats import_recording() can dispatch to
const _FORMATS = (
    :edf, :bdf, :gdf, :bv, :csv, :csv_gz, :set, :npy, :xdf, :nwb, :ncs,
    :fiff, :ft, :snirf, :nirs, :duomag,
)

# extension → format
const _EXT_FORMAT = Dict{String, Symbol}(
    ".edf"   => :edf,
    ".bdf"   => :bdf,
    ".gdf"   => :gdf,
    ".vhdr"  => :bv,
    ".ahdr"  => :bv,
    ".csv"   => :csv,
    ".set"   => :set,
    ".npy"   => :npy,
    ".xdf"   => :xdf,
    ".nwb"   => :nwb,
    ".ncs"   => :ncs,
    ".fif"   => :fiff,
    ".fiff"  => :fiff,
    ".mat"   => :ft,
    ".snirf" => :snirf,
    ".nirs"  => :nirs,
    ".ascii" => :duomag,
    ".m"     => :duomag,
)

# formats whose leading bytes identify them unambiguously: content wins over extension
const _UNIQUE_SIGNATURE = (:edf, :bdf, :gdf, :bv, :npy, :xdf, :ncs, :fiff)

# generic containers → formats that may legitimately use them (extension disambiguates)
const _CONTAINER_FORMATS = Dict{Symbol, Tuple{Vararg{Symbol}}}(
    :mat5  => (:ft, :set, :nirs),
    :mat73 => (:ft, :set, :nirs),
    :hdf5  => (:nwb, :snirf),
    :gzip  => (:csv_gz,),
    :text  => (:csv, :duomag),
)

const _SIG_DESCRIPTION = Dict{Symbol, String}(
    :mat5    => "a MATLAB v5 MAT-file",
    :mat73   => "a MATLAB v7.3 MAT-file",
    :hdf5    => "an HDF5 file",
    :gzip    => "a gzip-compressed file",
    :text    => "a plain-text file",
    :unknown => "an unrecognized binary file",
)

_hasprefix(buf::AbstractVector{UInt8}, p::AbstractVector{UInt8}) =
    length(buf) >= length(p) && view(buf, 1:length(p)) == p
_hasprefix(buf::AbstractVector{UInt8}, p::String) = _hasprefix(buf, codeunits(p))

"""
    _sniff_format(file_name)

Identify file format from its leading bytes. Returns a format symbol, a container symbol (`:mat5`, `:mat73`, `:hdf5`, `:gzip`, `:text`), `:unknown` or `:empty`.
"""
function _sniff_format(file_name::String)::Symbol
    buf = open(io -> read(io, 1024), file_name)
    isempty(buf) && return :empty
    # binary signatures
    _hasprefix(buf, UInt8[0x1f, 0x8b]) && return :gzip
    # EDF/EDF+
    _hasprefix(buf, "0       ") && return :edf
    # BDF/BDF+
    _hasprefix(buf, UInt8[0xff, codeunits("BIOSEMI")...]) && return :bdf
    # GDF 1.x / 2.x
    _hasprefix(buf, "GDF ") && return :gdf
    _hasprefix(buf, UInt8[0x93, codeunits("NUMPY")...]) && return :npy
    _hasprefix(buf, "XDF:") && return :xdf
    _hasprefix(buf, UInt8[0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a]) && return :hdf5
    _hasprefix(buf, "MATLAB 5.0 MAT-file") && return :mat5
    _hasprefix(buf, "MATLAB 7.3 MAT-file") && return :mat73            # HDF5 with 512 B userblock

    # FIFF: big-endian; first tag must be FIFF_FILE_ID (100), FIFFT_ID_STRUCT (31), size 20
    if length(buf) >= 12
        kind, typ, sz = ntoh.(reinterpret(Int32, buf[1:12]))
        (kind, typ, sz) == (100, 31, 20) && return :fiff
    end
    # text-based headers (strip UTF-8 BOM)
    txt = chopprefix(String(buf[1:min(end, 256)]), "\ufeff")
    occursin(r"^Brain ?Vision Data Exchange Header File", txt) && return :bv
    startswith(txt, "########") && occursin("Neuralynx", txt) && return :ncs
    return 0x00 in buf ? :unknown : :text
end

"""
    _ext_format(file_name)

Map file extension to format symbol (`nothing` if unknown).
"""
function _ext_format(file_name::String)::Union{Nothing, Symbol}
    base, ext = splitext(file_name)
    ext = lowercase(ext)
    if ext == ".gz"
        lowercase(splitext(base)[2]) == ".csv" && return :csv_gz
        return nothing
    end
    return get(_EXT_FORMAT, ext, nothing)
end

"""
    _hdf5_flavor(file_name)

Distinguish NWB from SNIRF by inspecting HDF5 root contents. Returns `:nwb`, `:snirf` or `nothing`.
"""
function _hdf5_flavor(file_name::String)::Union{Nothing, Symbol}
    try
        return HDF5.h5open(file_name, "r") do fid
            # NWB 2.x: root attribute; NWB 1.x: root dataset
            (haskey(HDF5.attributes(fid), "nwb_version") || haskey(fid, "nwb_version")) &&
                return :nwb
            haskey(fid, "formatVersion") && return :snirf
            return nothing
        end
    catch
        return nothing
    end
end

"""
    _sibling_hint(file_name)

Suggest the correct companion file when a data/marker file was passed instead of its header.
"""
function _sibling_hint(file_name::String)::String
    base, ext = splitext(file_name)
    ext = lowercase(ext)
    candidates = if ext in (".eeg", ".vmrk", ".dat")
        (".vhdr", ".VHDR", ".ahdr", ".AHDR")
    elseif ext == ".fdt"
        (".set", ".SET")
    else
        ()
    end
    for c in candidates
        isfile(base * c) && return " Did you mean $(base * c)?"
    end
    return ""
end

"""
    _detect_format(file_name)

Determine file format using content signature first and extension second.

- unique signature (EDF, BDF, GDF, BrainVision, NPY, XDF, NCS, FIFF): content wins; warns if extension disagrees
- HDF5: NWB/SNIRF identified by content when possible, otherwise by extension
- MATLAB / gzip / text containers: extension picks the format, but must be compatible with the container
- no conclusive signature: falls back to extension (with a warning if the format should have been recognizable)
"""
function _detect_format(file_name::String)::Symbol
    sig = _sniff_format(file_name)
    sig === :empty && throw(ArgumentError("File $file_name is empty."))
    ext_fmt = _ext_format(file_name)
 
    if sig in _UNIQUE_SIGNATURE
        ext_fmt !== nothing && ext_fmt !== sig &&
            @warn "$file_name has a :$ext_fmt extension, but its contents are :$sig. Importing as :$sig."
        return sig
    end
    if sig === :hdf5
        flavor = _hdf5_flavor(file_name)
        if flavor !== nothing
            ext_fmt !== nothing && ext_fmt !== flavor &&
                @warn "$file_name has a :$ext_fmt extension, but its contents are :$flavor. Importing as :$flavor."
            return flavor
        end
    end
    allowed = get(_CONTAINER_FORMATS, sig, ())
    ext_fmt in allowed && return ext_fmt
 
    # content is inconclusive (text or unknown binary): trust the extension
    if ext_fmt !== nothing && sig in (:text, :unknown)
        ext_fmt === :duomag ||
            @warn "Contents of $file_name do not match the expected :$ext_fmt signature (detected: $(_SIG_DESCRIPTION[sig])). Trying the :$ext_fmt importer anyway."
        return ext_fmt
    end
    desc = _SIG_DESCRIPTION[sig]
    msg = if ext_fmt === nothing
        "Cannot determine the format of $file_name (contents: $desc)."
    else
        "File $file_name has a :$ext_fmt extension, but its contents are $desc."
    end
    isempty(allowed) || (msg *= " Possible formats: $(join(":" .* string.(allowed), ", ")); pass `format` explicitly.")
    throw(ArgumentError(msg * _sibling_hint(file_name)))
end

"""
    import_recording(file_name; <keyword arguments>)

Load a recording file and return a `NeuroAnalyzer.NEURO` object.

This is a meta-function: it detects the file format and delegates to the appropriate `import_*()` function. Detection inspects the file contents (magic bytes / header signature) first and the extension second, so misnamed or extensionless files are handled where the format is identifiable. A warning is shown if the extension and contents disagree. The signal data type (EEG, MEG, fNIRS, …) is auto-detected where possible and stored in `obj.header.recording[:data_type]`.

Supported formats:

|             Format             |    Extension(s)   |   `format`  | Content detection             |
| ------------------------------ | ----------------- | ----------- | ----------------------------- |
| EDF / EDF+                     | `.edf`            | `:edf`      | yes                           |
| BDF / BDF+                     | `.bdf`            | `:bdf`      | yes                           |
| GDF                            | `.gdf`            | `:gdf`      | yes                           |
| BrainVision                    | `.vhdr`, `.ahdr`  | `:bv`       | yes                           |
| CSV                            | `.csv`            | `:csv`      | text only                     |
| CSV (gzip-compressed)          | `.csv.gz`         | `:csv_gz`   | gzip only                     |
| EEGLAB dataset                 | `.set`            | `:set`      | MAT-file only                 |
| MNE NumPy array                | `.npy`            | `:npy`      | yes                           |
| Extensible Data Format         | `.xdf`            | `:xdf`      | yes                           |
| Neurodata Without Borders      | `.nwb`            | `:nwb`      | yes                           |
| Neuralynx CSC                  | `.ncs`            | `:ncs`      | yes                           |
| Neuromag FIFF                  | `.fif`, `.fiff`   | `:fiff`     | yes                           |
| FieldTrip MATLAB               | `.mat`            | `:ft`       | MAT-file only                 |
| SNIRF                          | `.snirf`          | `:snirf`    | yes                           |
| NIRS                           | `.nirs`           | `:nirs`     | MAT-file only                 |
| DuoMAG TMS MEP                 | `.ascii`, `.m`    | `:duomag`   | no                            |

"MAT-file only" / "text only" means the container is verified, but the extension (or `format`) decides which importer is used.

# Arguments

- `file_name::String`: name of the file to load
- `detect_type::Bool=true`: infer channel type from channel label
- `type::Union{Nothing, Symbol}=nothing`: data type for FieldTrip (`.mat`) files, where auto-detection is not possible; one of `:eeg`, `:ecog`, `:ecog`, `:seeg`, `:ieeg`,`:meg`, `:nirs`, `:events`
- `sampling_rate::Union{Nothing, Int64}=nothing`: sampling rate for `.npy` files, which store only raw signal data with no embedded metadata
- `n::Int64=0`: subject index to extract from multi-subject SNIRF files
- `format::Union{Nothing, Symbol}=nothing`: force a specific format (see table), skipping detection

# Returns

- `NeuroAnalyzer.NEURO`: for EEG, MEG, and fNIRS data
- `DataFrame`: when `type = :events` (FieldTrip `.mat` files)
"""
function import_recording(
    file_name::String;
    detect_type::Bool = true,
    type::Union{Nothing, Symbol} = nothing,
    sampling_rate::Union{Nothing, Int64} = nothing,
    n::Int64 = 0,
    format::Union{Nothing, Symbol} = nothing,
)::Union{NeuroAnalyzer.NEURO, DataFrame}
    # validate
    !isnothing(type) && _check_var(type, [:eeg, :ecog, :seeg, :ieeg, :meg, :nirs, :events], "type")
    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))

    fmt = if format === nothing
        _detect_format(file_name)
    else
        format in _FORMATS || throw(
            ArgumentError(
                "Unknown format :$format. Must be one of: $(join(":" .* string.(_FORMATS), ", ")).",
            ),
        )
        format
    end

    fmt === :ft && isnothing(type) && throw(ArgumentError("For FieldTrip files type must be specified."))

    fmt === :edf && return import_edf(file_name; detect_type = detect_type)
    fmt === :bdf && return import_bdf(file_name; detect_type = detect_type)
    fmt === :gdf && return import_gdf(file_name; detect_type = detect_type)
    fmt === :bv && return import_bv(file_name; detect_type = detect_type)
    fmt in (:csv, :csv_gz) && return import_csv(file_name; detect_type = detect_type)
    fmt === :set && return import_set(file_name; detect_type = detect_type)
    # .npy stores raw signal data only - sampling rate must be supplied
    fmt === :npy && return import_npy(file_name; sampling_rate)
    fmt === :xdf && return import_xdf(file_name)
    fmt === :nwb && return import_nwb(file_name; detect_type = detect_type)
    fmt === :ncs && return import_ncs(file_name)
    fmt === :fiff && return import_fiff(file_name)
    # FieldTrip requires explicit type hint - auto-detection is not possible
    fmt === :ft && return import_ft(file_name; type)
    fmt === :snirf && return import_snirf(file_name; n)
    fmt === :nirs && return import_nirs(file_name)
    fmt === :duomag && return import_duomag(file_name)

    # unreachable: every member of _FORMATS is handled above
    throw(ErrorException("Internal error: unhandled format :$fmt."))
end