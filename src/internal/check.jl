_in(x::Real, r::Tuple{Real, Real})::Bool = x >= r[1] && x <= r[2]
_bin(x::Real, r::Tuple{Real, Real})::Bool = x > r[1] && x < r[2]

function _in(x::Real, r::Tuple{Real, Real}, v::String)::Nothing
    @assert x >= r[1] "$v must be ≥ $(r[1])."
    @assert x <= r[2] "$v must be ≤ $(r[2])."
    return nothing
end

function _bin(x::Real, r::Tuple{Real, Real}, v::String)::Nothing
    @assert x > r[1] "$v must be > $(r[1])."
    @assert x < r[2] "$v must be < $(r[2])."
    return nothing
end

function _chk2d(a::AbstractArray)::Nothing
    @assert ndims(a) == 2 "Input array must have 2 dimensions."
    return nothing
end

function _chk3d(a::AbstractArray)::Nothing
    @assert ndims(a) == 3 "Input array must have 3 dimensions."
    return nothing
end

function _chk4d(a::AbstractArray)::Nothing
    @assert ndims(a) == 4 "Input array must have 4 dimensions."
    return nothing
end

function _check_tuple(
    t::Tuple{Real, Real},
    r::Union{Nothing, Tuple{Real, Real}},
    name::Union{Nothing, String} = nothing,
    type::Symbol=:bin
)::Nothing
    _check_var(type, [:in, :bin], "type")
    @assert t == tuple_order(t) "$name must contain two values in ascending order."
    @assert t[1] < t[2] "$name must contain two different values in ascending order."
    isnothing(name) && (name = "Tuple")
    if r !== nothing
        if type === :bin
            @assert !(t[1] > r[1] || t[2] > r[1] || t[1] < r[2] || t[2] < r[2]) "$name must be in <$(r[1]), $(r[2])>."
        else
            @assert t[1] >= r[1] && t[2] >= r[1] && t[1] < r[2] && t[2] <= r[2] "$name must be in [$(r[1]), $(r[2])]."
        end
    end
    return nothing
end

function _check_channels(s::AbstractArray, ch::Union{Int64, Vector{Int64}, AbstractRange})::Nothing
    isa(ch, Int64) && (ch = [ch])
    [@assert !(ch_idx < 1 || ch_idx > size(s, 1)) "ch must be in [1, $(size(s, 1))]." for ch_idx in ch]
    return nothing
end

function _check_channels(obj::NeuroAnalyzer.NEURO, ch::Union{String, Vector{String}, Regex})::Nothing
    _check_channels(get_channel(obj; type = "all"), ch)
    return nothing
end

function _check_channels(obj::NeuroAnalyzer.NEURO, ch::Union{String, Vector{String}, Regex}, type::String)::Nothing
    _check_channels(get_channel(obj; type = type), ch)
    return nothing
end

function _check_channels(ch_ref::Union{String, Vector{String}}, ch::Union{String, Vector{String}, Regex})::Nothing
    isa(ch_ref, String) && (ch_ref = [ch_ref])
    isa(ch, String) && (ch = [ch])
    @assert length(ch) > 0 "ch is empty."
    @assert length(ch_ref) > 0 "ch_ref is empty."
    [@assert ch_idx in ch_ref "$ch_idx does not match labels." for ch_idx in ch]
    return nothing
end

function _check_epochs(obj::NeuroAnalyzer.NEURO, epoch::Union{Int64, Vector{Int64}, AbstractRange})::Nothing
    for idx in epoch
        @assert !(idx < 1 || idx > nepochs(obj)) "epoch must be in [1, $(nepochs(obj))]."
    end
    return nothing
end

function _check_segment(obj::NeuroAnalyzer.NEURO, seg::Tuple{Real, Real})::Nothing
    _check_segment(obj, seg[1], seg[2])
    return nothing
end

function _check_segment(obj::NeuroAnalyzer.NEURO, from::Real, to::Real)::Nothing
    @assert to > from "Segment end must be greater than segment start."
    @assert from >= obj.time_pts[1] "Segment start must be ≥ $(obj.time_pts[1])."
    @assert to >= obj.time_pts[1] "Segment end must be ≥ $(obj.time_pts[1])."
    @assert from <= obj.time_pts[end] "Segment start must be ≤ $(obj.time_pts[end])."
    @assert to <= obj.time_pts[end] "Segment end must be ≤ $(obj.time_pts[end])."
    return nothing
end

function _check_segment(signal::AbstractVector, from::Real, to::Real)::Nothing
    @assert from > 0 "Segment start must be > 0."
    @assert to > 0 "Segment end must be > 0."
    @assert to >= from "Segment end must be ≥ $from."
    @assert from <= length(signal) "Segment start must be ≤ $(length(signal))."
    @assert to <= length(signal) "Segment end must be ≤ $(length(signal))."
    return nothing
end

function _check_var(s1::Symbol, s2::Vector{Symbol}, var::String)::Nothing
    if length(s2) > 1
        m = var * " must be "
        for idx in 1:(length(s2) - 2)
            m *= ":" * string(s2[idx]) * ", "
        end
        m *= ":" * string(s2[end - 1]) * " or :" * string(s2[end]) * "."
        @assert s1 in s2 "$m"
    else
        m = var * " must be :" * string(s2[1])
        @assert s1 in s2 "$m"
    end
    return nothing
end

function _check_var(s1::String, s2::Vector{String}, var::String)::Nothing
    if length(s2) > 1
        m = var * " must be "
        for idx in 1:(length(s2) - 2)
            m *= s2[idx] * ", "
        end
        m *= s2[end - 1] * " or " * s2[end] * "."
        @assert s1 in s2 "$m"
    else
        m = var * " must be " * string(s2[1])
        @assert s1 in s2 "$m"
    end
    return nothing
end

function _check_markers(markers::Vector{String}, marker::String)::Nothing
    @assert marker in markers "Marker: $marker not found in markers."
    return nothing
end

function _check_markers(obj::NeuroAnalyzer.NEURO, marker::String)::Nothing
    @assert marker in unique(obj.markers[!, :value]) "Marker: $marker not found in markers."
    return nothing
end

function _check_datatype(obj::NeuroAnalyzer.NEURO, type::Union{String, Vector{String}})::Nothing
    if type isa String
        @assert datatype(obj) == type "This function works only for $(uppercase(string(type))) objects."
    else
        @assert obj.header.recording[:data_type] in type "This function works only for $(replace(uppercase(string(type)), "["=>"", "]"=>"", ":"=>"")) objects."
    end
    return nothing
end

function _check_svec(s::String)::Bool
    s = replace(s, " "=>"")
    s = replace(s, "["=>""; count = 1)
    s = replace(s, "]"=>""; count = 1)
    for idx in eachindex(s)
        string(s[idx]) in vcat(string.(0:9), [","]) || return false
    end
    if !(occursin(",", s) && length(split(s, ",")) > 1 && length(split(s, ",")[end]) > 0)
        return false
    else
        return true
    end
end

function _check_srange(s::String)::Bool
    s = replace(s, " "=>"")
    for idx in eachindex(s)
        string(s[idx]) in vcat(string.(0:9), [":"]) || return false
    end
    if occursin(":", s) &&
        length(split(s, ":")) == 2 &&
        length(s) > 0 &&
        length(split(s, ":")[1]) > 0 &&
        length(split(s, ":")[end]) > 0 &&
        parse(Int64, split(s, ":")[1]) < parse(Int64, split(s, ":")[end])
        return true
    else
        return false
    end
end

function _check_stuplei(s::String)::Bool
    s = replace(s, " "=>"")
    for idx in eachindex(s)
        string(s[idx]) in vcat(string.(0:9), [","], ["("], [")"]) || return false
    end
    if occursin("(", s) && occursin(")", s) && length(split(s, ",")) == 2 && length(s) > 0
        return true
    else
        return false
    end
end

function _check_stuplef(s::String)::Bool
    s = replace(s, " "=>"")
    for idx in eachindex(s)
        string(s[idx]) in vcat(string.(0:9), ["."], [","], ["("], [")"]) || return false
    end
    if occursin("(", s) && occursin(")", s) && length(split(s, ",")) == 2 && length(s) > 0
        return true
    else
        return false
    end
end

function _check_sfloat(s::String)::Bool
    for idx in eachindex(s)
        string(s[idx]) in vcat(string.(0:9), ["."]) || return false
    end
    return true
end

function _check_sint(s::String)::Bool
    for idx in eachindex(s)
        string(s[idx]) in string.(0:9) || return false
    end
    return true
end
