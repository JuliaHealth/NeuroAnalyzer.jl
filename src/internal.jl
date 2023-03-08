##################################
#                                #
#  Low-level internal functions  #
#                                #
##################################

_get_range(signal::Union{AbstractVector, AbstractArray}) = round(abs(minimum(signal)) + abs(maximum(signal)), digits=0)

function _gen_clabels(eeg::NeuroAnalyzer.EEG, c::Symbol)
    c, _ = _get_component(eeg, c)
    clabels = Vector{String}()
    for idx in 1:size(c, 1)
        push!(clabels, lpad(string(idx), length(string(size(c, 1))), "0"))
    end
    return clabels
end

function _gen_clabels(c::Array{Float64, 3})
    clabels = Vector{String}()
    for idx in 1:size(c, 1)
        push!(clabels, lpad(string(idx), length(string(size(c, 1))), "0"))
    end
    return clabels
end

function _set_defaults(xl::String, yl::String, tt::String, x::String, y::String, t::String)
    xl == "default" && (xl = x)
    yl == "default" && (yl = y)
    tt == "default" && (tt = t)
    return xl, yl, tt
end

function _len(eeg::NeuroAnalyzer.EEG, len::Int64, def_l::Int64)
    # return default length: one epoch (if epoch_len_seconds < def_l) or def_l seconds
    if len == 0
        if eeg_epoch_len(eeg) > def_l * eeg_sr(eeg)
            len = def_l * eeg_sr(eeg)
        else
            len = eeg_epoch_len(eeg)
        end
    end
    return len
end

function _get_epoch_markers(eeg::NeuroAnalyzer.EEG)
    return round.(s2t.(collect(1:eeg_epoch_len(eeg):eeg_epoch_len(eeg) * eeg_epoch_n(eeg)), eeg_sr(eeg)), digits=2)
end

function _get_t(from::Int64, to::Int64, fs::Int64)
    t = collect((from / fs):(1 / fs):(to / fs))
    t = t[1:(end - 1)]
    t[1] = floor(t[1], digits=2)
    t[2:(end - 1)] = round.(t[2:(end - 1)], digits=3)
    t[end] = ceil(t[end], digits=2)
    return t
end

function _convert_t(t1::Float64, t2::Float64)
    abs(t1) < 1.0 && (ts1 = string(floor(t1 * 1000, digits=2)) * " ms")
    abs(t1) >= 1.0 && (ts1 = string(floor(t1, digits=2)) * " s")
    abs(t2) < 1.0 && (ts2 = string(ceil(t2 * 1000, digits=2)) * " ms")
    abs(t2) >= 1.0 && (ts2 = string(ceil(t2, digits=2)) * " s")
    return t1, ts1, t2, ts2
end

function _tuple_max(t::Tuple{Real, Real})
    abs(t[1]) > abs(t[2]) && (t = (-abs(t[1]), abs(t[1])))
    abs(t[1]) < abs(t[2]) && (t = (-abs(t[2]), abs(t[2])))
    return t
end

function _channel2channel_name(channel::Union{Int64, Vector{Int64}, AbstractRange})
    if typeof(channel) == Int64
        return channel
    else
        if collect(channel[1]:channel[end]) == channel
            channel_name = string(channel[1]) * ":" * string(channel[end])
        else
            channel_name = ""
            for idx in 1:(length(channel) - 1)
                channel_name *= string(channel[idx])
                channel_name *= ", "
            end
            channel_name *= string(channel[end])
        end
    end
    return channel_name
end

function _s2epoch(eeg::NeuroAnalyzer.EEG, from::Int64, to::Int64)
    epoch = floor(Int64, from / eeg_epoch_len(eeg)):ceil(Int64, to / eeg_epoch_len(eeg))
    from / eeg_epoch_len(eeg) > from ÷ eeg_epoch_len(eeg) && (epoch = epoch[1] + 1:epoch[end])
    epoch[1] == 0 && (epoch = 1:epoch[end])
    epoch[1] == epoch[end] && (epoch = epoch[1])
    return epoch
end

function _epoch2s(eeg::NeuroAnalyzer.EEG, epoch::Int64)
    t1 = (epoch - 1) * eeg_epoch_len(eeg) + 1
    t2 = epoch * eeg_epoch_len(eeg)
    return t1, t2
end

function _get_channel_idx(labels::Vector{String}, channel::Union{String, Int64})
    if typeof(channel) == String
        channel_found = nothing
        for idx in eachindex(labels)
            if channel == labels[idx]
                channel_found = idx
            end
        end
        if channel_found === nothing
            throw(ArgumentError("channel name does not match signal labels."))
        end
    else
        channel < 1 || channel > length(labels) && throw(ArgumentError("channel index does not match signal channels."))
        channel_found = channel
    end

    return channel_found
end

function _angle_quadrant(a::Real)
    if a >= 0
        a = mod(a, 360)
        a <= 90 && (q = 1)
        (a > 90 && a <= 180) && (q = 2)
        (a > 180 && a <= 270) && (q = 3)
        (a > 270 && a < 360) && (q = 4)
    else
        a = mod(a, -360)
        a >= -90 && (q = 4)
        (a < -90 && a >= -180) && (q = 3)
        (a < -180 && a >= -270) && (q = 2)
        (a < -270 && a > -360) && (q = 1)
    end
    return q    
end 

function _free_gpumem(threshold::Real=0.95)
    m = CUDA.MemoryInfo()
    usedmem = m.total_bytes - m.free_bytes
    totalmem = m.total_bytes
    if usedmem / totalmem > threshold
        # CUDA.reclaim()
        GC.gc(true)
    end
end 

function _dict2labeled_matrix(d::Dict; rev::Bool=true)
    l = Vector{String}()
    v = Vector{Vector{Float64}}()
    for (kk, vv) in d
        push!(l, kk)
        push!(v, vv)
    end
    rev == true && return reverse!(l), reverse!(c)
    rev == false && return l, c
end

function _labeled_matrix2dict(l::Vector{String}, v::Vector{Vector{Float64}})
    length(l) == length(v) || throw(ArgumentError("Length of labels and values do not match."))
    return Dict(zip(l, v))
end

function _map_channels(channel::Union{Int64, Vector{Int64}, AbstractRange}, channels=Vector{Int64})
    channel_orig = channel
    if typeof(channel) == Int64
        channel = vsearch(channel, channels)
    else
        for idx in eachindex(channel)
            channel[idx] = vsearch(channel[idx], channels)
        end
    end
    return channel, channel_orig
end

_c(n) = collect(1:n)

function _delete_markers(markers::DataFrame, segment::Tuple{Int64, Int64})
    for marker_idx in nrow(markers):-1:1
        markers[marker_idx, :start] in segment[1]:segment[2] && deleteat!(markers, marker_idx)
    end
    return markers
end

function _shift_markers(m::DataFrame, pos::Int64, offset::Int64)
    markers = deepcopy(m)
    for marker_idx in 1:nrow(markers)
        markers[marker_idx, :start] > pos && (markers[marker_idx, :start] -= offset)
    end
    return markers
end

function _s2v(s::Union{<:Number, Vector{<:Number}})
    if typeof(s) <: Number
        return [s]
    else
        return s
    end
end

function _split(df::DataFrame, ratio::Float64=0.8)
    n = nrow(df)
    idx = shuffle(1:n)
    train_idx = view(idx, 1:floor(Int, ratio * n))
    test_idx = view(idx, (floor(Int, ratio * n) + 1):n)
    return df[train_idx, :], df[test_idx, :]
end
