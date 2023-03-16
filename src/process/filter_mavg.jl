export filter_mavg
export filter_mavg!

"""
    filter_mavg(s; <keyword arguments>)

Filter using moving average filter (with threshold).

# Arguments

- `s::AbstractVector`
- `k::Int64=8`: k-value for (window length = 2 × k + 1)
- `t::Real=0`: threshold (`t = t * std(s) + mean(s)`)
- `window::Union{Nothing, AbstractVector}=nothing`: weighting window

# Returns

- `s_filtered::Vector{Float64}`
"""
function filter_mavg(s::AbstractVector; k::Int64=8, t::Real=0, window::AbstractVector=ones(2 * k + 1))

    # check k
    k < 2 || k > length(s) && throw(ArgumentError("k must be in [2, signal length ($(length(s)))]."))
    mod(k, 2) != 0 && throw(ArgumentError("k must be even."))

    # check window
    length(window) != (2 * k + 1) && throw(ArgumentError("window length must be 2 × k + 1 ($(2 * k + 1))."))

    s_filtered = zeros(length(s))

    @inbounds for idx in (1 + k):(length(s) - k)
        if t > 0
            if s[idx] > t * std(s) + mean(s)
                s_filtered[idx] = @views mean(s[(idx - k):(idx + k)] .* window)
            end
        else
            s_filtered[idx] = @views mean(s[(idx - k):(idx + k)] .* window)
        end
    end

    return s_filtered

end

"""
    filter_mavg(s; k, t, window)

Filter using moving average filter (with threshold).

# Arguments

- `s::AbstractArray`
- `k::Int64=8`: k-value for (window length = 2 × k + 1)
- `t::Real=0`: threshold (`t = t * std(s) + mean(s)`)
- `window::Union{Nothing, AbstractVector}=nothing`: weighting window

# Returns

- `s_filtered::Array{Float64, 3}`: convoluted signal
"""
function filter_mavg(s::AbstractArray; k::Int64=8, t::Real=0, window::AbstractVector=ones(2 * k + 1))

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_filtered = similar(s)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_filtered[ch_idx, :, ep_idx] = @views filter_mavg(s[ch_idx, :, ep_idx], k=k, t=t, window=window)
        end
    end

    return s_filtered

end

"""
    filter_mavg(obj; ch, k, t, window)

Filter using moving average filter (with threshold).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `k::Int64=8`: k-value for (window length = 2 × k + 1)
- `t::Real=0`: threshold (`t = t * std(s) + mean(s)`)
- `window::Union{Nothing, AbstractVector}=nothing`: weighting window

# Returns

- `obj_new::NeuroAnalyzer.NEURO`: convoluted signal
"""
function filter_mavg(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), k::Int64=8, t::Real=0, window::AbstractVector=ones(2 * k + 1))

    _check_channels(obj, ch)

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = filter_mavg(obj.data[ch, :, :], k=k, t=t, window=window)
    reset_components!(obj_new)
    push!(obj_new.history, "filter_mavg(OBJ, ch=$ch, k=$k, t=$t, window=$window")

    return obj_new

end

"""
    filter_mavg!(obj; ch, k, t, window)

Filter using moving average filter (with threshold).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `k::Int64=8`: k-value for (window length = 2 × k + 1)
- `t::Real=0`: threshold (`t = t * std(s) + mean(s)`)
- `window::Union{Nothing, AbstractVector}=nothing`: weighting window

"""
function filter_mavg!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), k::Int64=8, t::Real=0, window::AbstractVector=ones(2 * k + 1))

    obj_new = filter_mavg(obj, ch=ch, k=k, t=t, window=window)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
