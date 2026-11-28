# ------------------------------------------------------------------------------
# processing history (obj.history): one entry per call, written as Julia code, e.g.
#   filter(obj; ch=["Fp1", "Fp2"], fprototype=:butterworth, ftype=:lp, cutoff=30, …)
# so that the history can be exported as a script (export_history())
# ------------------------------------------------------------------------------

# values that are stored in the history (as Julia code, `repr`)
_hist_ok(v)::Bool =
    v isa Union{Nothing, Missing, Bool, Number, AbstractString, Symbol, Regex, AbstractRange} ||
    (v isa Tuple && all(_hist_ok, v)) ||
    (v isa NamedTuple && all(_hist_ok, values(v))) ||
    (v isa AbstractArray && ndims(v) <= 2 && length(v) <= 256 && all(_hist_ok, v)) ||
    # labels (e.g. the channel order of a MEG recording)
    (v isa AbstractVector{<:AbstractString} && length(v) <= 1024) ||
    (v isa AbstractDict && length(v) <= 64 && all(x -> _hist_ok(x[1]) && _hist_ok(x[2]), v))

"""
    _hist_val(obj, k, v)

History representation of the value `v` of argument `k`: Julia code (`repr`); channel numbers of channel arguments (`ch`, `ch1`, `ch2`, `sch`) are written as channel labels of `obj`; values that are not stored (large arrays, tables, models, objects) are written as `<not stored: T>` (not valid Julia code: `export_history()` comments such calls out).
"""
function _hist_val(obj::NeuroAnalyzer.NEURO, k::Symbol, v)::String
    if k in (:ch, :ch1, :ch2, :sch) && !(v isa Bool) &&
       (v isa Integer || v isa AbstractVector{<:Integer})
        l = labels(obj)[v]
        return repr(l isa AbstractVector ? collect(String, l) : String(l))
    end
    _hist_ok(v) && return repr(v)
    return "<not stored: $(nameof(typeof(v)))>"
end

"""
    _history!(obj, name, kw; <keyword arguments>)

Add the call `name(obj; kw...)` to the processing history of `obj`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: object
- `name::String`: name of the function (the out-of-place version, e.g. `"filter"`; in-place only functions with `!`)
- `kw::NamedTuple=NamedTuple()`: all keyword arguments of the call, with their values
- `args::Tuple=()`: positional arguments after `obj`
- `from::Union{Nothing, Int64}=nothing`: number of history entries before the call; the entries added after it (by functions called internally, e.g. `filter_apply()` called by `filter()`) are replaced by this entry
- `note::String=""`: note written after the call as a comment (e.g. a measured filter response)

# Returns

- `Nothing`
"""
function _history!(
    obj::NeuroAnalyzer.NEURO,
    name::String,
    kw::NamedTuple = NamedTuple();
    args::Tuple = (),
    from::Union{Nothing, Int64} = nothing,
    note::String = "",
)::Nothing
    if !isnothing(from) && length(obj.history) > from
        deleteat!(obj.history, (from + 1):length(obj.history))
    end
    a = join(vcat("obj", [_hist_val(obj, :arg, x) for x in args]), ", ")
    k = join(["$k=$(_hist_val(obj, k, v))" for (k, v) in pairs(kw)], ", ")
    push!(
        obj.history,
        "$name($a" * (isempty(k) ? "" : "; $k") * ")" * (isempty(note) ? "" : " # $note"),
    )
    return nothing
end
