function _check_channels(record::NeuroAnalyzer.RECORD, channel::Union{Int64, Vector{Int64}, AbstractRange})
    for idx in channel
        (idx < 1 || idx > channel_n(record)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(channel_n(record))."))
    end
end

function _check_channels(eeg::NeuroAnalyzer.RECORD, channel::Union{Int64, Vector{Int64}, AbstractRange}, type::Symbol)
    channels = get_channel_bytype(record, type=type)
    for idx in channel
        idx in channels || throw(ArgumentError("Channel $idx does not match type: $(uppercase(string(type))) data channels."))
        (idx < 1 || idx > channel_n(record)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(channel_n(record))."))
    end
end

function _check_channels(channels::Union{Int64, Vector{Int64}, AbstractRange}, channel::Union{Int64, Vector{Int64}, AbstractRange})
    for idx in channel
        idx in channels || throw(ArgumentError("Channel $idx does not match signal channels."))
        (idx < 1 || idx > length(channels)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(channel_n(record))."))
    end
end

function _check_epochs(record::NeuroAnalyzer.RECORD, epoch::Union{Int64, Vector{Int64}, AbstractRange})
    for idx in epoch
        (idx < 1 || idx > epoch_n(record)) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(epoch_n(record))."))
    end
end

function _check_cidx(record::NeuroAnalyzer.RECORD, c::Symbol, cc::Union{Int64, Vector{Int64}, AbstractRange})
    c, _ = _get_component(record, c)
    for idx in cc
        (idx < 1 || idx > size(c, 1)) && throw(ArgumentError("cc must be ≥ 1 and ≤ $(size(c, 1))."))
    end
end

function _check_cidx(c::Array{Float64, 3}, cc::Union{Int64, Vector{Int64}, AbstractRange})
    for idx in cc
        (idx < 1 || idx > size(c, 1)) && throw(ArgumentError("cc must be ≥ 1 and ≤ $(size(c, 1))."))
    end
end

function _check_segment(record::NeuroAnalyzer.RECORD, from::Int64, to::Int64)
    from < 1 && throw(ArgumentError("from must be > 0."))
    to < 1 && throw(ArgumentError("to must be > 0."))
    to < from && throw(ArgumentError("to must be ≥ $from."))
    (from > signal_len(record)) && throw(ArgumentError("from must be ≤ $(signal_len(record))."))
    (to > signal_len(record)) && throw(ArgumentError("to must be ≤ $(signal_len(record))."))
end

function _check_segment(signal::AbstractVector, from::Int64, to::Int64)
    from < 0 && throw(ArgumentError("from must be > 0."))
    to < 0 && throw(ArgumentError("to must be > 0."))
    to < from && throw(ArgumentError("to must be ≥ $from."))
    from > length(signal) && throw(ArgumentError("from must be ≤ $(length(signal))."))
    to > length(signal) && throw(ArgumentError("to must be ≤ $(length(signal))."))
end

function _check_var(s1::Symbol, s2::Vector{Symbol}, var::String)
    m = var * " must be "
    for idx in 1:(length(s2) - 2)
        m *= ":" * string(s2[idx]) * ", "
    end
    m *= ":" * string(s2[end - 1]) * " or :" * string(s2[end]) * "."
    s1 in s2 || throw(ArgumentError(m))
end

function _check_var(s1::String, s2::Vector{String}, var::String)
    m = var * " must be "
    for idx in 1:(length(s2) - 2)
        m *= ":" * s2[idx] * ", "
    end
    m *= ":" * s2[end - 1] * " or " * s2[end] * "."
    s1 in s2 || throw(ArgumentError(m))
end

function _check_markers(markers::Vector{String}, marker::String)
    marker in markers || throw(ArgumentError("Marker: $marker not found in markers."))
end

function _check_markers(record::NeuroAnalyzer.RECORD, marker::String)
    marker in unique(record.markers[!, :description]) || throw(ArgumentError("Marker: $marker not found in markers."))
end
