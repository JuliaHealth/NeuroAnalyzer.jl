export reflect
export reflect!
export chop!

"""
    reflect(obj; <keyword arguments>)

Expand signal by adding reflected signal before the signal and after the signal, i.e. a signal 1234 becomes 432112344321. This may reduce edge artifacts, but will also affect amplitude of the filtered signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `n::Int64=sr(obj)`: number of samples to add, default is 1 second

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function reflect(obj::NeuroAnalyzer.NEURO; n::Int64 = sr(obj))::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    reflect!(obj_new; n = n)

    return obj_new
end

"""
    reflect!(obj; <keyword arguments>)

Expand signal by adding reflected signal before the signal and after the signal, i.e. a signal 1234 becomes 432112344321. This may reduce edge artifacts, but will also affect amplitude of the filtered signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `n::Int64=sr(obj)`: number of samples to add, default is 1 second

# Returns

- `Nothing`
"""
function reflect!(obj::NeuroAnalyzer.NEURO; n::Int64 = sr(obj))::Nothing
    # validate
    n >= 1 || throw(ArgumentError("n must be ≥ 1."))
    # add up to one epoch
    n > epoch_len(obj) && (n = epoch_len(obj))

    ch_n = nchannels(obj)
    ep_n = nepochs(obj)
    el = epoch_len(obj)
    s = zeros(ch_n, el + 2 * n, ep_n)

    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        x = @view obj.data[ch_idx, :, ep_idx]
        # mirror the first n samples before and the last n samples after the signal
        # (_reflect reverses the 2nd and 3rd arguments)
        s[ch_idx, :, ep_idx] = _reflect(x, @view(x[1:n]), @view(x[(end - n + 1):end]))
    end

    # epochs begin n samples earlier
    obj.epochs = copy(_epochs_table(obj))
    obj.epochs[!, :onset] .-= n / sr(obj)
    obj.data = s
    obj.time_pts, obj.epoch_time = _get_t(obj)

    # every epoch grows by 2n samples, n of them before the original samples
    obj.markers = _remap_markers(
        obj.markers, sr(obj), s -> ((s - 1) ÷ el) * (el + 2 * n) + n + (s - 1) % el + 1,
    )

    push!(obj.history, "reflect(obj, n=$n)")

    return nothing
end

"""
    chop(obj; <keyword arguments>)

Reduce signal by removing reflected signal before the signal and after the signal, i.e. a signal 432112344321 becomes 1234.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `n::Int64=sr(obj)`: number of samples to remove, default is 1 second

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function chop(obj::NeuroAnalyzer.NEURO; n::Int64 = sr(obj))::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    chop!(obj_new; n = n)

    return obj_new
end

"""
    chop!(obj; <keyword arguments>)

Reduce signal by removing reflected signal before the signal and after the signal, i.e. a signal 432112344321 becomes 1234.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `n::Int64=sr(obj)`: number of samples to remove, default is 1 second

# Returns

- `Nothing`
"""
function chop!(obj::NeuroAnalyzer.NEURO; n::Int64 = sr(obj))::Nothing
    # validate
    n >= 1 || throw(ArgumentError("n must be ≥ 1."))
    el = epoch_len(obj)
    # at least one sample must remain in each epoch
    2 * n < el ||
        throw(ArgumentError("n must be < $(cld(el, 2)) (half of the epoch length)."))

    ch_n = nchannels(obj)
    ep_n = nepochs(obj)
    s = zeros(ch_n, el - 2 * n, ep_n)

    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        s[ch_idx, :, ep_idx] = _chop(@view(obj.data[ch_idx, :, ep_idx]), n)
    end

    # epochs begin n samples later
    obj.epochs = copy(_epochs_table(obj))
    obj.epochs[!, :onset] .+= n / sr(obj)
    obj.data = s
    obj.time_pts, obj.epoch_time = _get_t(obj)

    # every epoch loses n samples at both ends; markers in the removed parts are dropped
    obj.markers = _remap_markers(
        obj.markers, sr(obj),
        s -> begin
            w = (s - 1) % el + 1
            n < w <= el - n ? ((s - 1) ÷ el) * (el - 2 * n) + w - n : nothing
        end,
    )

    push!(obj.history, "chop(obj, n=$n)")

    return nothing
end
