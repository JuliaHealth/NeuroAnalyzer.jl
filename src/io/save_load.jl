export save_neuro
export load_neuro

"""
    save_neuro(obj; <keyword arguments>)

Save `NeuroAnalyzer.NEURO` object to `file_name` file (HDF5-based).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `file_name::String`: name of the file to save to
- `overwrite::Bool=false`

# Return

Nothing
"""
function save_neuro(
    obj::NeuroAnalyzer.NEURO; file_name::String, overwrite::Bool = false,
)::Nothing
    # validate
    (isfile(file_name) && !overwrite) && throw(
        ArgumentError("File $file_name cannot be saved, to overwrite use overwrite=true."),
    )
    lowercase(splitext(file_name)[2]) == ".hdf" ||
        throw(ArgumentError("Filename extension must be .hdf"))

    obj.header.recording[:file_name] = file_name

    # the file size is stored in the object, so it is measured on a temporary copy first
    tmp_file = tempname() * ".hdf"
    try
        JLD2.save_object(tmp_file, obj)
        obj.header.recording[:file_size_mb] = round(filesize(tmp_file) / 1024^2; digits = 2)
    finally
        rm(tmp_file; force = true)
    end

    JLD2.save_object(file_name, obj)

    return nothing
end

"""
    load_neuro(file_name)

Load `NeuroAnalyzer.NEURO` object from `file_name` file (HDF5-based).

# Arguments

- `file_name::String`: name of the file to load

# Returns

- `NeuroAnalyzer.NEURO`
"""
function load_neuro(file_name::String)::NeuroAnalyzer.NEURO
    # validate
    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    # objects saved before v0.26.12 have no epoch metadata (`epochs`) and a header with an
    # `experiment` dictionary instead of `notes`: NEURO and HEADER objects are read through
    # `rconvert`, which adds the default table and converts the experiment to notes
    obj = JLD2.jldopen(
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
    # files saved by old versions of NeuroAnalyzer (e.g. EEG objects) are reconstructed by JLD2 as
    # anonymous types without the NEURO fields
    obj isa NeuroAnalyzer.NEURO || throw(
        ArgumentError(
            "$file_name does not contain a NEURO object (found: $(typeof(obj).name.name)); it may have been saved by an incompatible version of NeuroAnalyzer.",
        ),
    )

    # marker channels saved as numbers by older versions → channel labels
    _markers_upgrade!(obj)
    _epochs_upgrade!(obj)
    # recording notes (before v0.26.12: header.recording[:recording_notes]) → header notes
    rn = strip(string(pop!(obj.header.recording, :recording_notes, "")))
    isempty(rn) ||
        (obj.header.notes = isempty(obj.header.notes) ? rn : obj.header.notes * "\n" * rn)
    # :prefiltering descriptions → per-channel low-pass and high-pass frequencies
    _rec_filters!(obj.header.recording)

    _info(
        "Loaded: " *
        uppercase(obj.header.recording[:data_type]) *
        " ($(nchannels(obj)) × $(epoch_len(obj)) × $(nepochs(obj)); $(round(obj.time_pts[end], digits = 2)) s)",
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
