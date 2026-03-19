# ---------------------------------------------------------------------------
# range / bound predicates
# ---------------------------------------------------------------------------

"""Return `true` if `x ∈ [r[1], r[2]]` (closed interval)."""
_in(x::Real, r::Tuple{Real, Real})::Bool = x >= r[1] && x <= r[2]

"""Return `true` if `x ∈ (r[1], r[2])` (open interval)."""
_bin(x::Real, r::Tuple{Real, Real})::Bool = x > r[1] && x < r[2]

"""Assert that `x ∈ [r[1], r[2]]`; throw with `v` in the message if not."""
function _in(x::Real, r::Tuple{Real, Real}, v::String)::Nothing
    x >= r[1] || throw(ArgumentError("$v must be ≥ $(r[1])."))
    x <= r[2] || throw(ArgumentError("$v must be ≤ $(r[2])."))
    return nothing
end

"""Assert that `x ∈ (r[1], r[2])`; throw with `v` in the message if not."""
function _bin(x::Real, r::Tuple{Real, Real}, v::String)::Nothing
    x > r[1] || throw(ArgumentError("$v must be > $(r[1])."))
    x < r[2] || throw(ArgumentError("$v must be < $(r[2])."))
    return nothing
end

# ---------------------------------------------------------------------------
# dimensionality checks
# ---------------------------------------------------------------------------

"""Assert that `a` is 2-dimensional."""
function _chk2d(a::AbstractArray)::Nothing
    ndims(a) == 2 || throw(ArgumentError("Input array must be 2-dimensional; got $(ndims(a))."))
    return nothing
end

"""Assert that `a` is 3-dimensional."""
function _chk3d(a::AbstractArray)::Nothing
    ndims(a) == 3 || throw(ArgumentError("Input array must be 3-dimensional; got $(ndims(a))."))
    return nothing
end

"""Assert that `a` is 4-dimensional."""
function _chk4d(a::AbstractArray)::Nothing
    ndims(a) == 4 || throw(ArgumentError("Input array must be 4-dimensional; got $(ndims(a))."))
    return nothing
end

# ---------------------------------------------------------------------------
# Tuple validation
# ---------------------------------------------------------------------------

"""
Assert that tuple `t` contains two values in strict ascending order and lies within the reference range `r`. The `type` argument selects the bound type: `:in` (closed) or `:bin` (open).
"""
function _check_tuple(
    t::Tuple{Real, Real},
    r::Tuple{Real, Real},
    name::Union{Nothing, String}=nothing,
    type::Symbol=:in,
)::Nothing
    _check_var(type, [:in, :bin], "type")
    label = isnothing(name) ? "Tuple" : name
    t[1] < t[2] || throw(ArgumentError("$label must contain two strictly ascending values."))
    if r != t
        if type === :bin
            (t[1] > r[1] && t[2] < r[2]) || throw(ArgumentError("$label must be in ($(r[1]), $(r[2]))."))
        else
            (t[1] >= r[1] && t[2] <= r[2]) || throw(ArgumentError("$label must be in [$(r[1]), $(r[2])]."))
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# channel checks
# ---------------------------------------------------------------------------

"""Assert that all integer channel indices in `ch` are within `[1, size(s, 1)]`."""
function _check_channels(s::AbstractArray, ch::Union{Int64, Vector{Int64}, AbstractRange})::Nothing
    isa(ch, Int64) && (ch = [ch])
    n = size(s, 1)
    for ch_idx in ch
        (1 <= ch_idx <= n) || throw(ArgumentError("ch must be in [1, $n]; got $ch_idx."))
    end
    return nothing
end

"""Assert that all channel names in `ch` exist in the object."""
function _check_channels(obj::NeuroAnalyzer.NEURO, ch::Union{String, Vector{String}, Regex})::Nothing
    _check_channels(get_channel(obj; type="all"), ch)
    return nothing
end

"""Assert that all channel names in `ch` exist among channels of `type`."""
function _check_channels(obj::NeuroAnalyzer.NEURO, ch::Union{String, Vector{String}, Regex}, type::String)::Nothing
    _check_channels(get_channel(obj; type=type), ch)
    return nothing
end

"""Assert that all channel name(s) in `ch` are present in `ch_ref`."""
function _check_channels(ch_ref::Union{String, Vector{String}}, ch::Union{String, Vector{String}, Regex})::Nothing
    isa(ch_ref, String) && (ch_ref = [ch_ref])
    isa(ch, String) && (ch = [ch])
    length(ch) > 0 || throw(ArgumentError("ch must not be empty."))
    length(ch_ref) > 0 || throw(ArgumentError("ch_ref must not be empty."))
    for label in ch
        label in ch_ref || throw(ArgumentError("$label does not match any label in ch_ref."))
    end
    return nothing
end

# ---------------------------------------------------------------------------
# epoch / segment checks
# ---------------------------------------------------------------------------

"""Assert that all epoch indices are within `[1, nepochs(obj)]`."""
function _check_epochs(obj::NeuroAnalyzer.NEURO, epoch::Union{Int64, Vector{Int64}, AbstractRange})::Nothing
    n = nepochs(obj)
    for idx in epoch
        (1 <= idx <= n) || throw(ArgumentError("epoch must be in [1, $n]; got $idx."))
    end
    return nothing
end

"""Assert that the time segment `seg` is valid for `obj`."""
function _check_segment(obj::NeuroAnalyzer.NEURO, seg::Tuple{Real, Real})::Nothing
    _check_segment(obj, seg[1], seg[2])
    return nothing
end

"""Assert that the time segment `[from, to]` lies within `obj`'s time axis."""
function _check_segment(obj::NeuroAnalyzer.NEURO, from::Real, to::Real)::Nothing
    t0, t1 = obj.time_pts[1], obj.time_pts[end]
    to >= from || throw(ArgumentError("Segment end ($to) must be ≥ than segment start ($from)."))
    from >= t0 || throw(ArgumentError("Segment start must be ≥ $t0."))
    to <= t1 || throw(ArgumentError("Segment end must be ≤ $t1."))
    return nothing
end

"""Assert that the sample-index segment `[from, to]` is valid for `signal`."""
function _check_segment(signal::AbstractVector, from::Real, to::Real)::Nothing
    n = length(signal)
    from >= 1 || throw(ArgumentError("Segment start must be ≥ 1."))
    to >= 1 || throw(ArgumentError("Segment end must be ≥ 1."))
    to >= from || throw(ArgumentError("Segment end ($to) must be ≥ segment start ($from)."))
    from <= n || throw(ArgumentError("Segment start must be ≤ $n."))
    to <= n || throw(ArgumentError("Segment end must be ≤ $n."))
    return nothing
end

# ---------------------------------------------------------------------------
# variable / symbol validation
# ---------------------------------------------------------------------------

"""
Assert that symbol `s1` is one of the allowed values `s2`, producing a human-readable error listing all valid options.
"""
function _check_var(s1::Symbol, s2::Vector{Symbol}, var::String)::Nothing
    if s1 ∉ s2
        options = Base.join([":" * string(s) for s in s2], ", ", " or ")
        false || throw(ArgumentError("$var must be $options."))
    end
    return nothing
end

"""
Assert that string `s1` is one of the allowed values `s2`, producing a human-readable error listing all valid options.
"""
function _check_var(s1::String, s2::Vector{String}, var::String)::Nothing
    if s1 ∉ s2
        options = Base.join(s2, ", ", " or ")
        false || throw(ArgumentError("$var must be $options."))
    end
    return nothing
end

# ---------------------------------------------------------------------------
# marker checks
# ---------------------------------------------------------------------------

"""Assert that `marker` is present in the `markers` vector."""
function _check_markers(markers::Vector{String}, marker::String)::Nothing
    marker in markers || throw(ArgumentError("Marker '$marker' not found."))
    return nothing
end

"""Assert that `marker` is a known marker value in `obj`."""
function _check_markers(obj::NeuroAnalyzer.NEURO, marker::String)::Nothing
    marker in unique(obj.markers[!, :value]) || throw(ArgumentError("Marker '$marker' not found."))
    return nothing
end

# ---------------------------------------------------------------------------
# data-type check
# ---------------------------------------------------------------------------

"""Assert that `obj`'s data type matches `type` (String or Vector{String})."""
function _check_datatype(obj::NeuroAnalyzer.NEURO, type::Union{String, Vector{String}})::Nothing
    dt = datatype(obj)
    if type isa String
        dt == type || throw(ArgumentError("This function requires a $(uppercase(type)) object; got $(uppercase(dt))."))
    else
        dt in type|| throw(ArgumentError("This function requires one of $(uppercase.(type)); got $(uppercase(dt))."))
    end
    return nothing
end

# ---------------------------------------------------------------------------
# String-format validation helpers
# ---------------------------------------------------------------------------

"""Return `true` if `s` is a valid comma-separated integer vector string, e.g. `"[1,2,3]"`."""
function _check_svec(s::String)::Bool
    s = replace(s, " " => "", "[" => "", "]" => "")
    # all characters must be digits or commas
    all(c -> c in ('0':'9'..., ','), s) || return false
    parts = split(s, ",")
    # require at least 2 non-empty parts
    return length(parts) >= 2 && all(!isempty, parts)
end

"""Return `true` if `s` is a valid integer range string, e.g. `"1:10"` with start < end."""
function _check_srange(s::String)::Bool
    s = replace(s, " " => "")
    all(c -> c in ('0':'9'..., ':'), s) || return false
    parts = split(s, ":")
    length(parts) == 2 || return false
    all(!isempty, parts) || return false
    a, b = tryparse(Int64, parts[1]), tryparse(Int64, parts[2])
    return !isnothing(a) && !isnothing(b) && a < b
end

"""Return `true` if `s` is a valid integer 2-tuple string, e.g. `"(1,2)"`."""
function _check_stuplei(s::String)::Bool
    s = replace(s, " " => "")
    all(c -> c in ('0':'9'..., ',', '(', ')'), s) || return false
    return startswith(s, "(") && endswith(s, ")") &&
           length(split(s, ",")) == 2 && length(s) > 2
end

"""Return `true` if `s` is a valid float 2-tuple string, e.g. `"(1.5,2.0)"`."""
function _check_stuplef(s::String)::Bool
    s = replace(s, " " => "")
    all(c -> c in ('0':'9'..., '.', ',', '(', ')'), s) || return false
    return startswith(s, "(") && endswith(s, ")") &&
           length(split(s, ",")) == 2 && length(s) > 2
end

"""Return `true` if `s` consists entirely of digits and at most one decimal point."""
function _check_sfloat(s::String)::Bool
    all(c -> c in ('0':'9'..., '.'), s) || return false
    return count(==('.'), s) <= 1
end

"""Return `true` if `s` consists entirely of decimal digit characters."""
function _check_sint(s::String)::Bool
    return !isempty(s) && all(c -> c in '0':'9', s)
end