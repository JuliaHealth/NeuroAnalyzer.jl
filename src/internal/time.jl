function _get_t(obj::NeuroAnalyzer.NEURO)
    fs = sr(obj)
    time_pts = round.(collect(0:1/fs:size(obj.data, 2) * size(obj.data, 3) / fs)[1:end-1], digits=3)
    if length(obj.epoch_time) > 0
        epoch_time = round.((collect(0:1/fs:size(obj.data, 2) / fs) .+ obj.epoch_time[1])[1:end-1], digits=3)
    else
        epoch_time = round.((collect(0:1/fs:size(obj.data, 2) / fs))[1:end-1], digits=3)
    end
    return time_pts, epoch_time
end

function _get_t(from::Int64, to::Int64, fs::Int64)
    t = collect((from / fs):(1 / fs):(to / fs))
    t .-= t[1]
    t = round.(t, digits=3)
    #t = t[1:(end - 1)]
    #t[1] = floor(t[1], digits=2)
    #t[2:(end - 1)] = round.(t[2:(end - 1)], digits=3)
    #t[end] = ceil(t[end], digits=2)
end

function _convert_t(t1::Float64, t2::Float64)
    abs(t1) < 1.0 && (ts1 = string(floor(t1, digits=2) * 1000) * " ms")
    abs(t1) >= 1.0 && (ts1 = string(floor(t1, digits=2)) * " s")
    abs(t2) < 1.0 && (ts2 = string(ceil(t2, digits=2) * 1000) * " ms")
    abs(t2) >= 1.0 && (ts2 = string(ceil(t2, digits=2)) * " s")
    return t1, ts1, t2, ts2
end

function _s2epoch(obj::NeuroAnalyzer.NEURO, from::Int64, to::Int64)
    ep = floor(Int64, from / epoch_len(obj)):ceil(Int64, to / epoch_len(obj))
    from / epoch_len(obj) > from ÷ epoch_len(obj) && (ep = ep[1] + 1:ep[end])
    ep[1] == 0 && (ep = 1:ep[end])
    ep[1] == ep[end] && (ep = ep[1])
    return ep
end

function _epoch2s(obj::NeuroAnalyzer.NEURO, ep::Int64)
    t1 = (ep - 1) * epoch_len(obj) + 1
    t2 = ep * epoch_len(obj)
    return t1, t2
end
