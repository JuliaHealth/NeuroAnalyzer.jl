export reflect
export reflect!
export chop
export chop!

"""
    reflect(obj; n)

Expand signal by adding reflected signal before the signal and after the signal, i.e. a signal 1234 becomes 432112344321. This may reduce edge artifacts, but will also affect amplitude of the filtered signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `n::Int64=sr(obj)`: number of samples to add, default is 1 second

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function reflect(obj::NeuroAnalyzer.NEURO; n::Int64=sr(obj))

    # add up to one epoch
    n > epoch_len(obj) && (n = epoch_len(obj))

    obj_new = deepcopy(obj)
    ch_n = channel_n(obj)
    ep_n = epoch_n(obj)
    s = zeros(ch_n, epoch_len(obj) + 2 * n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s1 = obj_new.data[:, 1:n, ep_idx]
            s2 = obj_new.data[:, end:-1:(end - n + 1), ep_idx]
            @views s[ch_idx, :, ep_idx] = _reflect(obj.data[ch_idx, :, ep_idx], s1[ch_idx, :], s2[ch_idx, :])
        end
    end
    obj_new.data = s

    t = collect(0:(1 / sr(obj)):(size(obj_new.data, 2) / sr(obj)))[1:(end - 1)]
    obj_new.time_pts = t
    obj_new.epoch_time = t .+ obj.epoch_time[1]
    obj_new.header.recording[:duration_samples] = length(t) * ep_n
    obj_new.header.recording[:duration_seconds] = length(t) * ep_n * (1 / sr(obj))
    obj_new.header.recording[:epoch_duration_samples] = size(obj_new.data, 2)
    obj_new.header.recording[:epoch_duration_seconds] = size(obj_new.data, 2) * (1 / sr(obj))

    push!(obj_new.header.history, "reflect(OBJ, n=$n)")

    return obj_new
end

"""
    reflect!(obj; n)

Expand signal by adding reflected signal before the signal and after the signal, i.e. a signal 1234 becomes 432112344321. This may reduce edge artifacts, but will also affect amplitude of the filtered signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `n::Int64=sr(obj)`: number of samples to add, default is 1 second
"""
function reflect!(obj::NeuroAnalyzer.NEURO; n::Int64=sr(obj))

    obj_tmp = reflect(obj, n=n)

    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    obj.time = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time

    reset_components!(obj)

    return nothing
end

"""
    chop(obj; n)

Reduce signal by removing reflected signal before the signal and after the signal, i.e. a signal 432112344321 becomes 1234.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `n::Int64=sr(obj)`: number of samples to remove, default is 1 second

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function chop(obj::NeuroAnalyzer.NEURO; n::Int64=sr(obj))

    # add up to one epoch
    n > epoch_len(obj) && (n = epoch_len(obj))

    obj_new = deepcopy(obj)
    ch_n = channel_n(obj)
    ep_n = epoch_n(obj)
    s = zeros(ch_n, epoch_len(obj) - 2 * n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            @views s[ch_idx, :, ep_idx] = _chop(obj.data[ch_idx, :, ep_idx], n)
        end
    end
    obj_new.data = s

    t = collect(0:(1 / sr(obj)):(size(obj_new.data, 2) / sr(obj)))[1:(end - 1)]
    obj_new.time_pts = t
    obj_new.epoch_time = t .+ obj.epoch_time[1]
    obj_new.header.recording[:duration_samples] = length(t) * ep_n
    obj_new.header.recording[:duration_seconds] = length(t) * ep_n * (1 / sr(obj))
    obj_new.header.recording[:epoch_duration_samples] = size(obj_new.data, 2)
    obj_new.header.recording[:epoch_duration_seconds] = size(obj_new.data, 2) * (1 / sr(obj))

    push!(obj_new.header.history, "chop(OBJ, n=$n)")

    return obj_new
end

"""
    chop!(obj; c, v)

Reduce signal by removing reflected signal before the signal and after the signal, i.e. a signal 432112344321 becomes 1234.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `n::Int64=sr(obj)`: number of samples to remove, default is 1 second
"""
function chop!(obj::NeuroAnalyzer.NEURO; n::Int64=sr(obj))

    obj_tmp = chop(obj, n=n)

    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    obj.time = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time

    reset_components!(obj)

    return nothing
end
