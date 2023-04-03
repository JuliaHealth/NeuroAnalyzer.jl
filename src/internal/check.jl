function _check_channels(s::AbstractArray, ch::Union{Int64, Vector{Int64}, <:AbstractRange})
    for idx in ch
        (idx < 1 || idx > size(s, 1)) && throw(ArgumentError("ch must be in [1, $(size(s, 1))]."))
    end
    return nothing
end

function _check_channels(obj::NeuroAnalyzer.NEURO, ch::Union{Int64, Vector{Int64}, <:AbstractRange})
    for idx in ch
        (idx < 1 || idx > channel_n(obj)) && throw(ArgumentError("ch must be in [1, $(channel_n(obj))]."))
    end
    return nothing
end

function _check_channels(obj::NeuroAnalyzer.NEURO, ch::Union{Int64, Vector{Int64}, <:AbstractRange}, type::Symbol)
    channels = get_channel_bytype(obj, type=type)
    for idx in ch
        idx in channels || throw(ArgumentError("ch $idx does not match type: $(uppercase(string(type))) data channels."))
        (idx < 1 || idx > channel_n(obj)) && throw(ArgumentError("ch must be in [1, $(channel_n(obj))]."))
    end
    return nothing
end

function _check_channels(channels::Union{Int64, Vector{Int64}, <:AbstractRange}, ch::Union{Int64, Vector{Int64}, <:AbstractRange})
    for idx in ch
        idx in channels || throw(ArgumentError("ch must be in $channels."))
        (idx < 1 || idx > sort(channels)[end]) && throw(ArgumentError("ch must be in [1, $(channel_n(obj))]."))
    end
end

function _check_epochs(obj::NeuroAnalyzer.NEURO, epoch::Union{Int64, Vector{Int64}, <:AbstractRange})
    for idx in epoch
        (idx < 1 || idx > epoch_n(obj)) && throw(ArgumentError("epoch must be in [1, $(epoch_n(obj))]."))
    end
    return nothing
end

function _check_cidx(obj::NeuroAnalyzer.NEURO, c::Symbol, cc::Union{Int64, Vector{Int64}, <:AbstractRange})
    c, _ = _get_component(obj, c)
    if ndims(c) == 1
        cc != ndims(c) && throw(ArgumentError("cc must be in 1."))
    else
        for idx in cc
            (idx < 1 || idx > size(c, 1)) && throw(ArgumentError("cc must be in [1, $(size(c, 1))]."))
        end
    end
    return nothing
end

function _check_cidx(c::Union{AbstractVector, AbstractMatrix, AbstractArray}, cc::Union{Int64, Vector{Int64}, <:AbstractRange})
    if ndims(c) == 1
        cc != ndims(c) && throw(ArgumentError("cc must be in 1."))
    else
        for idx in cc
            (idx < 1 || idx > size(c, 1)) && throw(ArgumentError("cc must be in [1, $(size(c, 1))]."))
        end
    end
    return nothing
end

function _check_segment(obj::NeuroAnalyzer.NEURO, seg::Tuple{Real, Real})
    from = seg[1]
    to = seg[2]
    to < from && throw(ArgumentError("to must be > from."))
    from < obj.time_pts[1] && throw(ArgumentError("from must be ≥ $(obj.time_pts[1])."))
    to < obj.time_pts[1] && throw(ArgumentError("to must be ≥ $(obj.time_pts[1])."))
    (from > obj.time_pts[end]) && throw(ArgumentError("from must be ≤ $(obj.time_pts[end])."))
    (to > obj.time_pts[end]) && throw(ArgumentError("to must be ≤ $(obj.time_pts[end])."))
    return nothing
end

function _check_segment(obj::NeuroAnalyzer.NEURO, from::Int64, to::Int64)
    from < 1 && throw(ArgumentError("from must be ≥ 0."))
    to < 1 && throw(ArgumentError("to must be ≥ 0."))
    to < from && throw(ArgumentError("to must be ≥ $(obj.time_pts[vsearch(from / sr(obj), obj.time_pts)])."))
    (from > signal_len(obj)) && throw(ArgumentError("from must be ≤ $(obj.time_pts[end])."))
    (to > signal_len(obj)) && throw(ArgumentError("to must be ≤ $(obj.time_pts[end])."))
    return nothing
end

function _check_segment(signal::AbstractVector, from::Int64, to::Int64)
    from < 0 && throw(ArgumentError("from must be > 0."))
    to < 0 && throw(ArgumentError("to must be > 0."))
    to < from && throw(ArgumentError("to must be ≥ $from."))
    from > length(signal) && throw(ArgumentError("from must be ≤ $(length(signal))."))
    to > length(signal) && throw(ArgumentError("to must be ≤ $(length(signal))."))
    return nothing
end

function _check_var(s1::Symbol, s2::Vector{Symbol}, var::String)
    if length(s2) > 1
        m = var * " must be "
        for idx in 1:(length(s2) - 2)
            m *= ":" * string(s2[idx]) * ", "
        end
        m *= ":" * string(s2[end - 1]) * " or :" * string(s2[end]) * "."
        s1 in s2 || throw(ArgumentError(m))
    else
        m = var * " must be :" * string(s2[1])
        s1 in s2 || throw(ArgumentError(m))
    end
    return nothing
end

function _check_var(s1::String, s2::Vector{String}, var::String)
    if length(s2) > 1
        m = var * " must be "
        for idx in 1:(length(s2) - 2)
            m *= ":" * s2[idx] * ", "
        end
        m *= ":" * s2[end - 1] * " or " * s2[end] * "."
        s1 in s2 || throw(ArgumentError(m))
    else
        m = var * " must be " * string(s2[1])
        s1 in s2 || throw(ArgumentError(m))
    end
    return nothing
end

function _check_markers(markers::Vector{String}, marker::String)
    marker in markers || throw(ArgumentError("Marker: $marker not found in markers."))
    return nothing
end

function _check_markers(obj::NeuroAnalyzer.NEURO, marker::String)
    marker in unique(obj.markers[!, :description]) || throw(ArgumentError("Marker: $marker not found in markers."))
    return nothing
end

function _check_datatype(obj::NeuroAnalyzer.NEURO, type::Union{Symbol, Vector{Symbol}})
    if type isa Symbol
        Symbol(obj.header.recording[:data_type]) == type || throw(ArgumentError("This function works only for $(uppercase(string(type))) objects. Think carefully."))
    else
        Symbol(obj.header.recording[:data_type]) in type || throw(ArgumentError("This function works only for $(replace(uppercase(string(type)), "["=>"", "]"=>"", ":"=>"")) objects. Think carefully."))
    end
end
