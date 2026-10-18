# Generate the body of docs/src/index.md: one `@docs` block per source directory.
#
# Entries are taken from the docstrings that NeuroAnalyzer itself defines (Base.Docs.meta),
# so the list contains only:
# - functions that are actually loaded (files not included in NeuroAnalyzer.jl are skipped),
# - public functions: exported by NeuroAnalyzer, or NeuroAnalyzer methods added to functions of
#   Base and other packages (e.g. `Base.size(obj)`, `StatsBase.zscore(x, m, sd)`),
# - documented functions.
# Private functions (`_*`, or not exported) and plugin functions are left out. Because
# `makedocs(modules = [NeuroAnalyzer])`, an entry such as `NeuroAnalyzer.zscore` shows only the
# docstrings written in NeuroAnalyzer, never those of the package that owns the function.
#
# usage: julia --project=.. template.jl src/index.md   (appends to the file)

# output goes to a file: loading the package may print to stdout (e.g. plugin messages)
length(ARGS) == 1 || error("usage: julia --project=.. template.jl <output file>")
const OUT = ARGS[1]

using NeuroAnalyzer

const NA = NeuroAnalyzer
const SRC = joinpath(pkgdir(NA), "src")

# section title => source directory (in this order)
const SECTIONS = [
    "NeuroAnalyzer" => "na",
    "Utils" => "utils",
    "Stats" => "stats",
    "IO" => "io",
    "Edit" => "edit",
    "Process" => "process",
    "Locs" => "locs",
    "Analyze" => "analyze",
    "Model" => "model",
    "Plot" => "plots",
    # "GUI" => "gui",
    "NeuroRecorder" => "recorder",
    "NeuroStim" => "stim",
    "NeuroTester" => "tester",
]

# source directory of a docstring (nothing if it is not in src/, e.g. a plugin)
function _doc_dir(d::Base.Docs.DocStr)::Union{Nothing, String}
    path = get(d.data, :path, nothing)
    path isa AbstractString || return nothing
    rel = relpath(normpath(path), SRC)
    startswith(rel, "..") && return nothing
    parts = splitpath(rel)
    return length(parts) > 1 ? parts[1] : nothing
end

_ispublic(b::Base.Docs.Binding)::Bool =
    !startswith(string(b.var), "_") && (Base.isexported(NA, b.var) || b.mod !== NA)

entries = Dict(dir => Set{String}() for (_, dir) in SECTIONS)
for (b, multidoc) in Base.Docs.meta(NA)
    _ispublic(b) || continue
    # one section per function: that of its first docstring
    dirs = filter(!isnothing, [_doc_dir(multidoc.docs[sig]) for sig in multidoc.order])
    isempty(dirs) && continue
    haskey(entries, dirs[1]) && push!(entries[dirs[1]], "NeuroAnalyzer.$(b.var)")
end

# exported NeuroAnalyzer functions without a docstring are missing from the documentation
# (functions of other packages that are only re-exported, e.g. `sem`, have no docs here)
for n in names(NA)
    (n === :NeuroAnalyzer || startswith(string(n), "_")) && continue
    isdefined(NA, n) && getfield(NA, n) isa Function || continue
    b = Base.Docs.Binding(NA, n)
    b.mod === NA || continue
    haskey(Base.Docs.meta(NA), b) ||
        @warn "Exported function without a docstring (not in the documentation): $n"
end

open(OUT, "a") do io
    println(io, "# NeuroAnalyzer.jl documentation")
    println(io)
    println(io, "This documentation has been generated using [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/).")
    for (idx, (title, dir)) in enumerate(SECTIONS)
        isempty(entries[dir]) && continue
        idx > 1 && (println(io); println(io, "---"))
        println(io)
        println(io, "## $title")
        println(io)
        println(io, "```@docs")
        foreach(e -> println(io, e), sort(collect(entries[dir])))
        println(io, "```")
    end
end
