function _check_tuple(t::Tuple{Real, Real}, name::String, range::Union{Nothing, Tuple{Real, Real}}=nothing)
    @assert t == tuple_order(t) "$name must contain two values in ascending order."
    @assert t[1] < t[2] "$name must contain two different values in ascending order."
    if range !== nothing 
        @assert !(t[1] < range[1] || t[2] < range[1] || t[1] > range[2] || t[2] > range[2]) "$name must be in [$(range[1]), $(range[2])]."
    end
    return nothing
end

function _check_channels(s::AbstractArray, ch::Union{Int64, Vector{Int64}, <:AbstractRange})
    for idx in ch
        @assert !(idx < 1 || idx > size(s, 1)) "ch must be in [1, $(size(s, 1))]."
    end
    return nothing
end

function _check_channels(obj::NeuroAnalyzer.NEURO, ch::Union{Int64, Vector{Int64}, <:AbstractRange})
    for idx in ch
        @assert !(idx < 1 || idx > nchannels(obj)) "ch must be in [1, $(nchannels(obj))]."
    end
    return nothing
end

function _check_channels(obj::NeuroAnalyzer.NEURO, ch::Union{Int64, Vector{Int64}, <:AbstractRange}, type::String)
    channels = get_channel_bytype(obj, type=type)
    for idx in ch
        @assert idx in channels "ch $idx does not match type: $(uppercase(string(type))) data channels."
        @assert !(idx < 1 || idx > nchannels(obj)) "ch must be in [1, $(nchannels(obj))]."
    end
    return nothing
end

function _check_channels(channels::Union{Int64, Vector{Int64}, <:AbstractRange}, ch::Union{Int64, Vector{Int64}, <:AbstractRange})
    for idx in ch
        @assert idx in channels "ch must be in $channels."
        @assert !(idx < 1 || idx > sort(channels)[end]) "ch must be in [1, $(nchannels(obj))]."
    end
end

function _check_epochs(obj::NeuroAnalyzer.NEURO, epoch::Union{Int64, Vector{Int64}, <:AbstractRange})
    for idx in epoch
        @assert !(idx < 1 || idx > nepochs(obj)) "epoch must be in [1, $(nepochs(obj))]."
    end
    return nothing
end

function _check_cidx(obj::NeuroAnalyzer.NEURO, c::Symbol, cc::Union{Int64, Vector{Int64}, <:AbstractRange})
    c, _ = _get_component(obj, c)
    if ndims(c) == 1
        @assert cc == 1 "cc must be 1."
    else
        for idx in cc
            @assert !(idx < 1 || idx > size(c, 1)) "cc must be in [1, $(size(c, 1))]."
        end
    end
    return nothing
end

function _check_cidx(c::Union{AbstractVector, AbstractMatrix, AbstractArray}, cc::Union{Int64, Vector{Int64}, <:AbstractRange})
    if ndims(c) == 1
        @assert cc == 1 "cc must be in 1."
    else
        for idx in cc
            @assert !(idx < 1 || idx > size(c, 1)) "cc must be in [1, $(size(c, 1))]."
        end
    end
    return nothing
end

function _check_segment(obj::NeuroAnalyzer.NEURO, seg::Tuple{Real, Real})
    from = seg[1]
    to = seg[2]
    @assert to > from "to must be > from."
    @assert from >= obj.time_pts[1] "from must be ≥ $(obj.time_pts[1])."
    @assert to >= obj.time_pts[1] "to must be ≥ $(obj.time_pts[1])."
    @assert from <= obj.time_pts[end] "from must be ≤ $(obj.time_pts[end])."
    @assert to <= obj.time_pts[end] "to must be ≤ $(obj.time_pts[end])."
    return nothing
end

function _check_segment_topo(obj::NeuroAnalyzer.NEURO, seg::Tuple{Real, Real})
    from = seg[1]
    to = seg[2]
    @assert to >= from "to must be >= from."
    @assert from >= obj.time_pts[1] "from must be ≥ $(obj.time_pts[1])."
    @assert to >= obj.time_pts[1] "to must be ≥ $(obj.time_pts[1])."
    @assert from <= obj.time_pts[end] "from must be ≤ $(obj.time_pts[end])."
    @assert to <= obj.time_pts[end] "to must be ≤ $(obj.time_pts[end])."
    return nothing
end

function _check_segment(obj::NeuroAnalyzer.NEURO, from::Int64, to::Int64)
    @assert from >= 0 "from must be ≥ 0."
    @assert to >= 0 "to must be ≥ 0."
    @assert to >= from "to must be ≥ $(obj.time_pts[vsearch(from / sr(obj), obj.time_pts)])."
    @assert from <= signal_len(obj) "from must be ≤ $(signal_len(obj))."
    @assert to <= signal_len(obj) "to must be ≤ $(signal_len(obj))."
    return nothing
end

function _check_segment(signal::AbstractVector, from::Int64, to::Int64)
    @assert from > 0 "from must be > 0."
    @assert to > 0 "to must be > 0."
    @assert to >= from "to must be ≥ $from."
    @assert from <= length(signal) "from must be ≤ $(length(signal))."
    @assert to <= length(signal) "to must be ≤ $(length(signal))."
    return nothing
end

function _check_var(s1::Symbol, s2::Vector{Symbol}, var::String)
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

function _check_var(s1::String, s2::Vector{String}, var::String)
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

function _check_markers(markers::Vector{String}, marker::String)
    @assert marker in markers "Marker: $marker not found in markers."
    return nothing
end

function _check_markers(obj::NeuroAnalyzer.NEURO, marker::String)
    @assert marker in unique(obj.markers[!, :description]) "Marker: $marker not found in markers."
    return nothing
end

function _check_datatype(obj::NeuroAnalyzer.NEURO, type::Union{String, Vector{String}})
    if type isa String
        @assert obj.header.recording[:data_type] == type "This function works only for $(uppercase(string(type))) objects."
    else
        @assert obj.header.recording[:data_type] in type "This function works only for $(replace(uppercase(string(type)), "["=>"", "]"=>"", ":"=>"")) objects."
    end
end

function _check_svec(s::String)
    s = replace(s, " "=>"")
    s = replace(s, "["=>"", count=1)
    s = replace(s, "]"=>"", count=1)
    for idx in eachindex(s)
        string(s[idx]) in vcat(string.(0:9), [","]) || return false
    end
    if !(occursin(",", s) && length(split(s, ",")) > 1 && length(split(s, ",")[end]) > 0)
        return false
    else
        return true
    end
end

function _check_srange(s::String)
    s = replace(s, " "=>"")
    for idx in eachindex(s)
        string(s[idx]) in vcat(string.(0:9), [":"]) || return false
    end
    if occursin(":", s) && length(split(s, ":")) == 2 && length(s) > 0 && length(split(s, ":")[1]) > 0 && length(split(s, ":")[end]) > 0 && parse(Int64, split(s, ":")[1]) < parse(Int64, split(s, ":")[end])
        return true
    else
        return false
    end
end

function _check_stuplei(s::String)
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

function _check_stuplef(s::String)
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

function _check_sfloat(s::String)
    for idx in eachindex(s)
        string(s[idx]) in vcat(string.(0:9), ["."]) || return false
    end
    return true
end

function _check_sint(s::String)
    for idx in eachindex(s)
        string(s[idx]) in string.(0:9) || return false
    end
    return true
end