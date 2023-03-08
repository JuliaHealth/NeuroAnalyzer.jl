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
