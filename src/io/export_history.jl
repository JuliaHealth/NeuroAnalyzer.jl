export export_history

"""
    export_history(obj; <keyword arguments>)

Export the processing history (`obj.history`) to a Julia script that repeats the processing: the recording is imported with `import_recording()` and each history entry becomes one call (`obj = name(obj; …)`, or `name!(obj; …)` if an in-place version exists).

Entries that cannot be repeated are written as comments, with the reason:

- values that are not stored in the history (data arrays, tables, models, other objects, interactive edits; written as `<not stored: …>`)
- entries that are not calls of a NeuroAnalyzer function of `obj` (e.g. the history imported from EEGLAB files or entries written by NeuroAnalyzer before v0.26.12)

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `file_name::String`: path of the output `.jl` file
- `source::String=obj.header.recording[:file_name]`: recording imported by the script (after `save()` this is the `.hdf` file, which already contains the processing; give the original recording here)
- `overwrite::Bool=false`: allow overwriting an existing file

# Returns

- `Nothing`

# Notes

The history contains the processing that changes the object; analyses and plots are not recorded. Channels are written as labels; the arguments are written as given in the call, so the script repeats the processing with the same settings (the filter design measured when the filter was applied is written as a comment).
"""
function export_history(
    obj::NeuroAnalyzer.NEURO;
    file_name::String,
    source::String = string(obj.header.recording[:file_name]),
    overwrite::Bool = false,
)::Nothing
    # validate
    lowercase(splitext(file_name)[2]) == ".jl" ||
        throw(ArgumentError("Filename extension must be .jl"))
    isfile(file_name) && !overwrite && throw(
        ArgumentError("File $file_name already exists; use overwrite=true to overwrite."),
    )
    isempty(obj.history) && _warn("The processing history is empty.")

    io = IOBuffer()
    println(io, "# processing history exported by NeuroAnalyzer v$(NeuroAnalyzer.VER)")
    println(io, "# $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM"))")
    println(io)
    println(io, "using NeuroAnalyzer")
    println(io)
    println(io, "obj = import_recording($(repr(source)))")
    n_skip = 0
    for h in obj.history
        l, ok = _hist_script(h)
        ok || (n_skip += 1)
        println(io, l)
    end
    println(io)
    println(io, "obj")

    write(file_name, take!(io))
    n_skip > 0 && _warn("$n_skip history entr$(n_skip == 1 ? "y" : "ies") written as comments.")

    return nothing
end

# one history entry as a line of the script; returns (line, repeatable)
function _hist_script(h::AbstractString)::Tuple{String, Bool}
    h = strip(h)
    occursin("<not stored", h) &&
        return ("# (not repeatable: values not stored in the history) " * h, false)
    ex = try
        Meta.parse(h)
    catch
        nothing
    end
    name = _hist_call(ex)
    isnothing(name) &&
        return ("# (not repeatable: not a NeuroAnalyzer call of obj) " * h, false)
    # in-place version if it exists: no copy of the object
    endswith(name, "!") && return (h, true)
    nb = Symbol(name * "!")
    isdefined(NeuroAnalyzer, nb) && return (name * "!" * h[(length(name) + 1):end], true)
    return ("obj = " * h, true)
end

# name of the NeuroAnalyzer function called with obj as the first positional argument
function _hist_call(ex)::Union{Nothing, String}
    ex isa Expr && ex.head === :call || return nothing
    f = ex.args[1]
    f isa Symbol && isdefined(NeuroAnalyzer, f) || return nothing
    pos = [a for a in ex.args[2:end] if !(a isa Expr && a.head in (:parameters, :kw))]
    !isempty(pos) && pos[1] === :obj || return nothing
    return string(f)
end
