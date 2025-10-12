export filter_poly
export filter_poly!

"""
    filter_poly(s; <keyword arguments>)

Filter using polynomial filter.

# Arguments

- `s::AbstractVector`
- `order::Int64=8`: polynomial order
- `window::Int64=10`: window length

# Returns

- `s_filtered::Vector{Float64}`
"""
function filter_poly(s::AbstractVector; order::Int64=8, window::Int64=10)::Vector{Float64}

    @assert order > 1 "order must be > 1."
    @assert !(window < 1 || window > length(s)) "window must be in [1, $(length(s))]."

    s_filtered = deepcopy(s)

    # TO DO: smooth spikes between windows
    window_n = length(s) ÷ window
    window_last = length(s) - (window_n * window)

    for window_idx in 1:(window_n - 1)
        s_tmp = @views s[((window_idx - 1) * window + 1):window_idx * window]
        t = eachindex(s_tmp)
        p = Polynomials.fit(t, s_tmp, order)
        @inbounds for idx in eachindex(s_tmp)
            s_tmp[idx] = p(t[idx])
        end
        @inbounds s_filtered[((window_idx - 1) * window + 1):window_idx * window] = s_tmp
    end

    s_tmp = s[((window_n - 1) * window + 1):window_n * window + window_last]
    t = eachindex(s_tmp)
    p = Polynomials.fit(t, s_tmp, order)
    @inbounds for idx in eachindex(s_tmp)
        s_tmp[idx] = p(t[idx])
    end
    @inbounds s_filtered[((window_n - 1) * window + 1):window_n * window + window_last] = s_tmp

    # smooth peaks using Loess
    for window_idx in window:window:window * (window_n - 1)
        s_tmp = @views s_filtered[window_idx - window ÷ 4:window_idx + window ÷ 4 - 1]
        t = collect(1.0:1:length(s_tmp))
        model = Loess.loess(t, Vector(s_tmp), span=1.0)
        s_smoothed = Loess.predict(model, t)
        s_smoothed[1] = s_tmp[1]
        s_smoothed[end] = s_tmp[end]
        @inbounds s_filtered[window_idx - window ÷ 4:window_idx + window ÷ 4 - 1] = s_smoothed
    end

    # smooth peaks using median/mean
    # for window_idx in window:window:window * (window_n - 1)
    #     s_filtered[window_idx-2] = mean(s_filtered[window_idx-3:window_idx])
    #     s_filtered[window_idx-1] = mean(s_filtered[window_idx-2:window_idx])
    #     s_filtered[window_idx] = mean(s_filtered[window_idx-2:window_idx+2])
    #     s_filtered[window_idx+1] = mean(s_filtered[window_idx:window_idx+2])
    #     s_filtered[window_idx+2] = mean(s_filtered[window_idx:window_idx+3])
    # end

    return s_filtered

end

"""
    filter_poly(s; <keyword arguments>)

Filter using polynomial filter.

# Arguments

- `s::AbstractArray`
- `order::Int64=8`: polynomial order
- `window::Int64=10`: window length

# Returns

- `s_filtered::Array{Float64, 3}`: convoluted signal
"""
function filter_poly(s::AbstractArray; order::Int64=8, window::Int64=10)::Array{Float64, 3}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_filtered = similar(s)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_filtered[ch_idx, :, ep_idx] = @views filter_poly(s[ch_idx, :, ep_idx], order=order, window=window)
        end
    end

    return s_filtered

end

"""
    filter_poly(obj; <keyword arguments>)

Filter using polynomial filter.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `order::Int64=8`: polynomial order
- `window::Int64=10`: window length

# Returns

- `obj_new::NeuroAnalyzer.NEURO`: convoluted signal
"""
function filter_poly(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, order::Int64=8, window::Int64=10)::NeuroAnalyzer.NEURO

    ch = get_channel(obj, ch=ch)
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = filter_poly(obj.data[ch, :, :], order=order, window=window)
    reset_components!(obj_new)
    push!(obj_new.history, "filter_poly(OBJ, ch=$ch, order=$order, window=$window")

    return obj_new

end

"""
    filter_poly!(obj; <keyword arguments>)

Filter using polynomial filter.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `order::Int64=8`: polynomial order
- `window::Int64=10`: window length

# Returns

Nothing
"""
function filter_poly!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, order::Int64=8, window::Int64=10)::Nothing

    obj_new = filter_poly(obj, ch=ch, order=order, window=window)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
