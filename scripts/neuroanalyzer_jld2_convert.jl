# Conversion of NeuroAnalyzer files saved with JLD2 (NeuroAnalyzer before v0.26.12) to the
# HDF5 file format of NeuroAnalyzer v0.26.12 and later.
#
# NeuroAnalyzer v0.26.12 saves NEURO objects (save_neuro()) and studies (save_study()) in an
# explicit HDF5 layout and no longer depends on JLD2; this script contains the JLD2 loaders
# removed from NeuroAnalyzer (load_neuro() and load_study() of older files).
#
# Requirements: an environment with NeuroAnalyzer (v0.26.12 or later), JLD2 and DataFrames, e.g.
#
#     julia> using Pkg; Pkg.activate(temp = true); Pkg.add(["NeuroAnalyzer", "JLD2", "DataFrames"])
#     julia> include("/tmp/neuroanalyzer_jld2_convert.jl")
#
# Usage:
#
#     obj = load_neuro_jld2("old.hdf")                 # NEURO object
#     study = load_study_jld2("old_study.hdf")         # STUDY object
#     convert_jld2("old.hdf", "new.hdf")               # NEURO or study file → new format
#     convert_jld2_dir("folder"; backup = true)        # all old .hdf files of a folder (and
#                                                      # its subfolders), converted in place
#
# Command line:
#
#     julia --project=<env> neuroanalyzer_jld2_convert.jl old.hdf new.hdf
#     julia --project=<env> neuroanalyzer_jld2_convert.jl folder

using NeuroAnalyzer
using JLD2
using DataFrames

# ---------------------------------------------------------------------------
# JLD2 upgrade paths of the NeuroAnalyzer types
# ---------------------------------------------------------------------------

# a NEURO object read from a file as a NamedTuple of its fields; objects saved before
# v0.26.12 have no `epochs` field (the default epoch table is created)
function JLD2.rconvert(::Type{NeuroAnalyzer.NEURO}, nt::NamedTuple)::NeuroAnalyzer.NEURO
    haskey(nt, :epochs) && return NeuroAnalyzer.NEURO(
        nt.header, nt.history, nt.markers, nt.epochs, nt.locs, nt.time_pts,
        nt.epoch_time, nt.data,
    )
    return NeuroAnalyzer.NEURO(
        nt.header, nt.history, nt.markers, nt.locs, nt.time_pts, nt.epoch_time, nt.data,
    )
end

# a header read from a file as a NamedTuple of its fields; headers saved before v0.26.12
# have an `experiment` dictionary (name, notes, design) instead of `notes`
function JLD2.rconvert(::Type{NeuroAnalyzer.HEADER}, nt::NamedTuple)::NeuroAnalyzer.HEADER
    haskey(nt, :notes) && return NeuroAnalyzer.HEADER(nt.subject, nt.recording, nt.notes)
    e = haskey(nt, :experiment) ? nt.experiment : Dict()
    notes = NeuroAnalyzer._header_notes(
        "Experiment" => get(e, :name, ""), "" => get(e, :notes, ""),
        "Design" => get(e, :design, ""),
    )
    return NeuroAnalyzer.HEADER(nt.subject, nt.recording, notes)
end

# ---------------------------------------------------------------------------
# loaders
# ---------------------------------------------------------------------------

"""
    load_neuro_jld2(file_name)

Load a NEURO object saved with JLD2 by NeuroAnalyzer before v0.26.12.
"""
function load_neuro_jld2(file_name::String)::NeuroAnalyzer.NEURO
    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

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
                ArgumentError(
                    "$file_name must contain one object (found: $(length(k))).",
                ),
            )
            f[k[1]]
        end
    catch e
        e isa ArgumentError && rethrow()
        throw(
            ArgumentError(
                "$file_name is not a NeuroAnalyzer JLD2 file ($(sprint(showerror, e))).",
            ),
        )
    end
    # files saved by old versions of NeuroAnalyzer (e.g. EEG objects) are reconstructed by
    # JLD2 as anonymous types without the NEURO fields
    obj isa NeuroAnalyzer.NEURO || throw(
        ArgumentError(
            "$file_name does not contain a NEURO object (found: $(typeof(obj).name.name)); it may have been saved by an incompatible version of NeuroAnalyzer.",
        ),
    )

    # marker channels saved as numbers → channel labels; marker `id` → `type`
    NeuroAnalyzer._markers_upgrade!(obj)
    # :epoch_id → epoch conditions
    NeuroAnalyzer._epochs_upgrade!(obj)
    # recording notes (header.recording[:recording_notes]) → header notes
    rn = strip(string(pop!(obj.header.recording, :recording_notes, "")))
    isempty(rn) ||
        (obj.header.notes = isempty(obj.header.notes) ? rn : obj.header.notes * "\n" * rn)
    # :prefiltering descriptions → per-channel low-pass and high-pass frequencies
    NeuroAnalyzer._rec_filters!(obj.header.recording)

    return obj
end

"""
    load_study_jld2(file_name)

Load a study saved with JLD2 by NeuroAnalyzer before v0.26.12 (`root` as saved).
"""
function load_study_jld2(file_name::String)::NeuroAnalyzer.STUDY
    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    d = try
        JLD2.load(file_name)
    catch e
        throw(
            ArgumentError(
                "$file_name is not a NeuroAnalyzer study ($(sprint(showerror, e))).",
            ),
        )
    end
    haskey(d, "neuroanalyzer_study") || throw(
        ArgumentError("$file_name does not contain a NeuroAnalyzer study."),
    )
    return NeuroAnalyzer.STUDY(
        d["name"],
        # (studies saved before v0.26.12 have no notes)
        get(d, "notes", ""),
        d["root"],
        d["design"],
        Dict{Symbol, Any}(Symbol(k) => v for (k, v) in d["measures"]),
        d["history"],
    )
end

# ---------------------------------------------------------------------------
# conversion
# ---------------------------------------------------------------------------

# :neuro, :study (JLD2 files), :new (already in the HDF5 format of v0.26.12), :other
function jld2_kind(file_name::String)::Symbol
    NeuroAnalyzer.HDF5.ishdf5(file_name) || return :other
    fmt = NeuroAnalyzer.HDF5.h5open(file_name, "r") do f
        a = NeuroAnalyzer.HDF5.attributes(f)
        haskey(a, "format") ? string(read(a["format"])) : ""
    end
    startswith(fmt, "NeuroAnalyzer") && return :new
    k = try
        JLD2.jldopen(f -> keys(f), file_name, "r")
    catch
        return :other
    end
    "neuroanalyzer_study" in k && return :study
    length(k) == 1 && return :neuro
    return :other
end

"""
    convert_jld2(src, dst; overwrite = false)

Convert a NEURO or study file saved with JLD2 (`src`) to the HDF5 format of NeuroAnalyzer v0.26.12 (`dst`, `.hdf`). Returns the kind of the file (`:neuro` or `:study`).
"""
function convert_jld2(src::String, dst::String; overwrite::Bool = false)::Symbol
    kind = jld2_kind(src)
    if kind === :neuro
        save_neuro(load_neuro_jld2(src); file_name = dst, overwrite = overwrite)
    elseif kind === :study
        save_study(load_study_jld2(src); file_name = dst, overwrite = overwrite)
    elseif kind === :new
        throw(ArgumentError("$src is already in the current NeuroAnalyzer format."))
    else
        throw(ArgumentError("$src is not a NeuroAnalyzer JLD2 file."))
    end
    return kind
end

"""
    convert_jld2_dir(dir; backup = true)

Convert all `.hdf` files saved with JLD2 in `dir` and its subfolders in place; with `backup = true` the original files are kept as `<name>.hdf.jld2`. Files already in the current format and other files are skipped. Returns a table: `file`, `kind`, `status`.
"""
function convert_jld2_dir(dir::String; backup::Bool = true)::DataFrame
    isdir(dir) || throw(ArgumentError("Folder $dir does not exist."))
    out = DataFrame(:file => String[], :kind => Symbol[], :status => String[])
    for (root, _, files) in walkdir(dir), f in files
        lowercase(splitext(f)[2]) == ".hdf" || continue
        p = joinpath(root, f)
        kind = jld2_kind(p)
        if !(kind in (:neuro, :study))
            push!(out, (p, kind, "skipped"))
            continue
        end
        tmp = p * ".converting.hdf"
        try
            convert_jld2(p, tmp; overwrite = true)
            backup && cp(p, p * ".jld2"; force = true)
            mv(tmp, p; force = true)
            push!(out, (p, kind, "converted"))
        catch e
            rm(tmp; force = true)
            push!(out, (p, kind, "error: " * sprint(showerror, e)))
        end
    end
    return out
end

# command line
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) == 2
        k = convert_jld2(ARGS[1], ARGS[2])
        println("Converted ($k): $(ARGS[1]) → $(ARGS[2])")
    elseif length(ARGS) == 1 && isdir(ARGS[1])
        show(convert_jld2_dir(ARGS[1]); allrows = true, truncate = 200)
        println()
    else
        println("Usage: julia neuroanalyzer_jld2_convert.jl old.hdf new.hdf | folder")
    end
end
