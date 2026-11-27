export save_neuro
export load_neuro

# NeuroAnalyzer NEURO file layout (HDF5, see save_neuro()); the version is increased when the
# layout changes, files of older versions are converted by _neuro_migrate!()
const _NEURO_FORMAT = "NeuroAnalyzer NEURO"
const _NEURO_FORMAT_VERSION = 1

"""
    save_neuro(obj; <keyword arguments>)

Save `NeuroAnalyzer.NEURO` object to `file_name` file (HDF5).

The file layout (readable by any HDF5 software, e.g. h5py, MATLAB `h5read()`, HDFView):

- root attributes: `format` ("NeuroAnalyzer NEURO"), `format_version` (1), `neuroanalyzer_version`
- `/data`: signals, Julia array of size channels × samples × epochs (attribute `dimensions`); HDF5 stores arrays in row-major order, so other software sees it as epochs × samples × channels; chunked (one chunk per channel and epoch, up to 2^17 samples) and compressed
- `/time_pts`, `/epoch_time`: time points (s)
- `/history`: processing history (strings)
- `/header/notes`: notes (string)
- `/header/subject`, `/header/recording`: header dictionaries, one member per field
- `/markers`, `/epochs`, `/locs`: tables, one dataset per column (group attribute `columns`: column order)

Values that are not numbers, strings or arrays of them (e.g. symbols, vectors of vectors, dictionaries, tables) are stored as groups or datasets with the `na_type` attribute describing them.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `file_name::String`: name of the file to save to (`.hdf`)
- `overwrite::Bool=false`
- `compress::Int64=1`: compression (deflate) level of `/data`, 0 (no compression, fastest) to 9; higher levels make slightly smaller files but are slower

# Return

Nothing
"""
function save_neuro(
    obj::NeuroAnalyzer.NEURO;
    file_name::String,
    overwrite::Bool = false,
    compress::Int64 = 1,
)::Nothing
    # validate
    (isfile(file_name) && !overwrite) && throw(
        ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."),
    )
    lowercase(splitext(file_name)[2]) == ".hdf" ||
        throw(ArgumentError("Filename extension must be .hdf"))
    _in(compress, (0, 9), "compress")

    obj.header.recording[:file_name] = file_name
    # the file size is stored in the file: written as 0 and updated after the file is closed
    obj.header.recording[:file_size_mb] = 0.0

    # written to a temporary file: an existing file is replaced only if saving succeeds
    tmp_file = file_name * ".tmp"
    try
        HDF5.h5open(tmp_file, "w") do f
            a = HDF5.attributes(f)
            a["format"] = _NEURO_FORMAT
            a["format_version"] = _NEURO_FORMAT_VERSION
            a["neuroanalyzer_version"] = string(pkgversion(NeuroAnalyzer))

            d = obj.data
            if compress > 0 && all(>(0), size(d))
                ds = HDF5.create_dataset(
                    f, "data", Float64, size(d);
                    chunk = (1, min(size(d, 2), 2^17), 1), shuffle = true,
                    deflate = compress,
                )
                write(ds, d)
            else
                f["data"] = d
            end
            HDF5.attributes(f["data"])["dimensions"] = ["channel", "sample", "epoch"]
            f["time_pts"] = Vector{Float64}(obj.time_pts)
            f["epoch_time"] = Vector{Float64}(obj.epoch_time)
            f["history"] = Vector{String}(obj.history)

            h = HDF5.create_group(f, "header")
            h["notes"] = obj.header.notes
            _h5_write(h, "subject", obj.header.subject; path = "header/subject")
            _h5_write(h, "recording", obj.header.recording; path = "header/recording")

            _h5_write(f, "markers", obj.markers)
            _h5_write(f, "epochs", obj.epochs)
            _h5_write(f, "locs", obj.locs)
        end
        fs = round(filesize(tmp_file) / 1024^2; digits = 2)
        HDF5.h5open(f -> write(f["header/recording/file_size_mb"], fs), tmp_file, "r+")
        obj.header.recording[:file_size_mb] = fs
        mv(tmp_file, file_name; force = true)
    finally
        rm(tmp_file; force = true)
    end

    return nothing
end

"""
    load_neuro(file_name; <keyword arguments>)

Load `NeuroAnalyzer.NEURO` object from `file_name` file (saved by `save_neuro()`). Files saved by NeuroAnalyzer before v0.26.12 (JLD2) are also loaded.

Selected channels and/or epochs may be loaded: only these signals are read from the file.

# Arguments

- `file_name::String`: name of the file to load
- `ch::Union{Nothing, String, Vector{String}, Regex}=nothing`: channel name(s) to load (default: all); NIRS: selecting an intensity or OD channel also loads the channels of the same optode pair at the other wavelengths
- `ep::Union{Nothing, Int64, Vector{Int64}, AbstractUnitRange{Int64}}=nothing`: epoch number(s) to load (default: all)

# Returns

- `NeuroAnalyzer.NEURO`
"""
function load_neuro(
    file_name::String;
    ch::Union{Nothing, String, Vector{String}, Regex} = nothing,
    ep::Union{Nothing, Int64, Vector{Int64}, AbstractUnitRange{Int64}} = nothing,
)::NeuroAnalyzer.NEURO
    # validate
    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    fmt = _neuro_format_version(file_name)
    obj = if isnothing(fmt)
        # not written by save_neuro() of v0.26.12 or later: JLD2 file
        o = _load_neuro_jld2(file_name)
        if !isnothing(ch)
            n0 = length(o.history)
            keep_channel!(o; ch = ch)
            length(o.history) > n0 && pop!(o.history)
        end
        if !isnothing(ep)
            n0 = length(o.history)
            keep_epoch!(o; ep = ep)
            length(o.history) > n0 && pop!(o.history)
        end
        o
    else
        fmt > _NEURO_FORMAT_VERSION && throw(
            ArgumentError(
                "$file_name was saved by a newer version of NeuroAnalyzer (file format version $fmt); update NeuroAnalyzer to load it.",
            ),
        )
        HDF5.h5open(f -> _load_neuro_hdf5(f, fmt; ch = ch, ep = ep), file_name, "r")
    end

    # objects saved by older versions
    _neuro_upgrade!(obj)
    (isnothing(ch) && isnothing(ep)) || _history!(obj, "load_neuro", (; ch, ep))

    _info(
        "Loaded: " *
        uppercase(obj.header.recording[:data_type]) *
        " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits = 2)) s)",
    )

    return obj
end

# file format version of a NEURO file written by save_neuro(); nothing: not such a file
# (e.g. a JLD2 file of an older version)
function _neuro_format_version(file_name::String)::Union{Nothing, Int64}
    HDF5.ishdf5(file_name) || throw(
        ArgumentError("$file_name is not a NeuroAnalyzer file (not an HDF5 file)."),
    )
    return HDF5.h5open(file_name, "r") do f
        a = HDF5.attributes(f)
        haskey(a, "format") || return nothing
        fmt = string(read(a["format"]))
        fmt == _STUDY_FORMAT && throw(
            ArgumentError(
                "$file_name contains a NeuroAnalyzer study (use load_study()).",
            ),
        )
        fmt == _NEURO_FORMAT ||
            throw(ArgumentError("$file_name is not a NeuroAnalyzer file ($fmt)."))
        return Int64(read(a["format_version"]))
    end
end

# conversion of files of older format versions (to the current layout, before reading)
function _neuro_migrate(::HDF5.File, ::Int64)::Nothing
    # format version 1 is the first version
    return nothing
end

# read the HDF5 layout (see save_neuro())
function _load_neuro_hdf5(
    f, fmt::Int64;
    ch::Union{Nothing, String, Vector{String}, Regex} = nothing,
    ep::Union{Nothing, Int64, Vector{Int64}, AbstractUnitRange{Int64}} = nothing,
)::NeuroAnalyzer.NEURO
    _neuro_migrate(f, fmt)
    for k in ("data", "time_pts", "epoch_time", "history", "header", "markers", "epochs",
              "locs")
        haskey(f, k) || throw(ArgumentError("The file is damaged (no /$k)."))
    end

    h = f["header"]
    subject = _h5_read(h["subject"])
    recording = _h5_read(h["recording"])
    header = NeuroAnalyzer.HEADER(
        Dict{Symbol, Any}(subject), Dict{Symbol, Any}(recording), read(h["notes"]),
    )
    hist = read(f["history"])
    history = isempty(hist) ? String[] : Vector{String}(hist)

    ds = f["data"]
    ch_n, ep_len, ep_n = size(ds)
    # the header and metadata are selected on an object without samples (the data are
    # read afterwards)
    obj = NeuroAnalyzer.NEURO(
        header, history, _h5_read(f["markers"]), _h5_read(f["epochs"]),
        _h5_read(f["locs"]), Vector{Float64}(read(f["time_pts"])),
        Vector{Float64}(read(f["epoch_time"])), zeros(Float64, ch_n, 0, ep_n),
    )

    # channels
    ch_idx = collect(1:ch_n)
    if !isnothing(ch)
        ch_idx = get_channel(obj; ch = ch)
        ch_idx = unique(ch_idx isa Integer ? [Int64(ch_idx)] : Vector{Int64}(ch_idx))
        isempty(ch_idx) && throw(ArgumentError("No channels selected."))
        _nirs_n_opt(obj.header.recording) > 0 &&
            (ch_idx = _nirs_pair_wl(obj, ch_idx; msg = "Also loading"))
        sort!(ch_idx)
        length(ch_idx) < ch_n && _select_channels!(obj, ch_idx)
    end

    # epochs
    ep_idx = collect(1:ep_n)
    if !isnothing(ep)
        ep_idx = sort(unique(_n2v(ep)))
        isempty(ep_idx) && throw(ArgumentError("No epochs selected."))
        all(e -> 1 <= e <= ep_n, ep_idx) ||
            throw(ArgumentError("ep must be in [1, $ep_n]."))
        if length(ep_idx) < ep_n
            obj.epochs = DataFrames.nrow(obj.epochs) == ep_n ? obj.epochs[ep_idx, :] :
                         DataFrame()
            _delete_epoch_markers!(
                obj.markers, setdiff(1:ep_n, ep_idx), ep_len / sr(obj),
            )
        end
    end

    # signals: contiguous runs of channels of each epoch (HDF5 reads ranges)
    if length(ch_idx) == ch_n && length(ep_idx) == ep_n
        obj.data = read(ds)
    else
        data = Array{Float64, 3}(undef, length(ch_idx), ep_len, length(ep_idx))
        runs = UnitRange{Int64}[]
        for c in ch_idx
            if !isempty(runs) && last(runs[end]) == c - 1
                runs[end] = first(runs[end]):c
            else
                push!(runs, c:c)
            end
        end
        for (k, e) in enumerate(ep_idx)
            pos = 1
            for r in runs
                data[pos:(pos + length(r) - 1), :, k] = ds[r, :, e]
                pos += length(r)
            end
        end
        obj.data = data
        length(ep_idx) < ep_n && ((obj.time_pts, obj.epoch_time) = _get_t(obj))
    end
    obj.epochs = _epochs_table(obj)

    return obj
end

# upgrades of objects saved by older versions of NeuroAnalyzer
function _neuro_upgrade!(obj::NeuroAnalyzer.NEURO)::Nothing
    # marker channels saved as numbers → channel labels
    _markers_upgrade!(obj)
    _epochs_upgrade!(obj)
    # recording notes (before v0.26.12: header.recording[:recording_notes]) → header notes
    rn = strip(string(pop!(obj.header.recording, :recording_notes, "")))
    isempty(rn) ||
        (obj.header.notes = isempty(obj.header.notes) ? rn : obj.header.notes * "\n" * rn)
    # :prefiltering descriptions → per-channel low-pass and high-pass frequencies
    _rec_filters!(obj.header.recording)
    return nothing
end

# NEURO objects saved with JLD2 (before v0.26.12)
function _load_neuro_jld2(file_name::String)::NeuroAnalyzer.NEURO
    # objects saved before v0.26.12 have no epoch metadata (`epochs`) and a header with an
    # `experiment` dictionary instead of `notes`: NEURO and HEADER objects are read through
    # `rconvert`, which adds the default table and converts the experiment to notes
    obj = try
        JLD2.jldopen(
        file_name, "r";
        typemap = Dict(
            "NeuroAnalyzer.NEURO" => JLD2.Upgrade(NeuroAnalyzer.NEURO),
            "NeuroAnalyzer.HEADER" => JLD2.Upgrade(NeuroAnalyzer.HEADER),
        ),
    ) do f
        k = keys(f)
        length(k) == 1 || throw(
            ArgumentError("$file_name must contain one object (found: $(length(k)))."),
        )
        f[k[1]]
    end
    catch e
        e isa ArgumentError && rethrow()
        throw(
            ArgumentError(
                "$file_name is not a NeuroAnalyzer file ($(sprint(showerror, e)))."
            ),
        )
    end
    # files saved by old versions of NeuroAnalyzer (e.g. EEG objects) are reconstructed by JLD2 as
    # anonymous types without the NEURO fields
    obj isa NeuroAnalyzer.NEURO || throw(
        ArgumentError(
            "$file_name does not contain a NEURO object (found: $(typeof(obj).name.name)); it may have been saved by an incompatible version of NeuroAnalyzer.",
        ),
    )
    return obj
end

# upgrade path of `load_neuro()`: a NEURO object read from a file as a NamedTuple of its
# fields; objects saved before v0.26.12 have no `epochs` field
function JLD2.rconvert(::Type{NeuroAnalyzer.NEURO}, nt::NamedTuple)::NeuroAnalyzer.NEURO
    haskey(nt, :epochs) && return NeuroAnalyzer.NEURO(
        nt.header, nt.history, nt.markers, nt.epochs, nt.locs, nt.time_pts,
        nt.epoch_time, nt.data,
    )
    return NeuroAnalyzer.NEURO(
        nt.header, nt.history, nt.markers, nt.locs, nt.time_pts, nt.epoch_time, nt.data,
    )
end

# upgrade path of `load_neuro()`: a header read from a file as a NamedTuple of its fields;
# headers saved before v0.26.12 have an `experiment` dictionary (name, notes, design) instead
# of `notes`
function JLD2.rconvert(::Type{NeuroAnalyzer.HEADER}, nt::NamedTuple)::NeuroAnalyzer.HEADER
    haskey(nt, :notes) && return NeuroAnalyzer.HEADER(nt.subject, nt.recording, nt.notes)
    e = haskey(nt, :experiment) ? nt.experiment : Dict()
    notes = _header_notes(
        "Experiment" => get(e, :name, ""), "" => get(e, :notes, ""),
        "Design" => get(e, :design, ""),
    )
    return NeuroAnalyzer.HEADER(nt.subject, nt.recording, notes)
end
