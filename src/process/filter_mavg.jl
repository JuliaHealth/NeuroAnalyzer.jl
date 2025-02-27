export filter_mavg
export filter_mavg!

"""
    filter_mavg(s; <keyword arguments>)

Filter using moving average filter (with threshold).

# Arguments

- `s::AbstractVector`
- `k::Int64=8`: window length is `2 × k + 1`; for cut-off frequency f, k is `sqrt(0.196202 + f^2) / f`
- `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`); only samples below/above the threshold are being filtered
- `window::Union{Nothing, AbstractVector}=nothing`: weighting window

# Returns

- `s_filtered::Vector{Float64}`
"""
function filter_mavg(s::AbstractVector; k::Int64=8, t::Real=0, window::AbstractVector=ones(2 * k + 1))::Vector{Float64}

    # check k
    @assert k in 1:length(s) "k must be in [1, signal length ($(length(s)))]."

    # check window
    @assert length(window) == (2 * k + 1) "window length must be `2 × k + 1` ($(2 * k + 1))."

    s_filtered = deepcopy(s)

    @inbounds for idx in (1 + k):(length(s) - k)
        if t > 0
            if s[idx] < mean(s) - t * std(s) || s[idx] > (mean(s) + t * std(s))
                s_filtered[idx] = @views mean(s[(idx - k):(idx + k)] .* window)
            end
        else
            s_filtered[idx] = @views mean(s[(idx - k):(idx + k)] .* window)
        end
    end

    return s_filtered

end

"""
    filter_mavg(s; <keyword arguments>)

Filter using moving average filter (with threshold).

# Arguments

- `s::AbstractArray`
- `k::Int64=8`: window length is `2 × k + 1`; for cut-off frequency f, k is `sqrt(0.196202 + f^2) / f`
- `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`); only samples below/above the threshold are being filtered
- `window::Union{Nothing, AbstractVector}=nothing`: weighting window

# Returns

- `s_filtered::Array{Float64, 3}`
"""
function filter_mavg(s::AbstractArray; k::Int64=8, t::Real=0, window::AbstractVector=ones(2 * k + 1))::Array{Float64, 3}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_filtered = similar(s)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            s_filtered[ch_idx, :, ep_idx] = @views filter_mavg(s[ch_idx, :, ep_idx], k=k, t=t, window=window)
        end
    end

    return s_filtered

end

"""
    filter_mavg(obj; <keyword arguments>)

Filter using moving average filter (with threshold).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `k::Int64=8`: window length is `2 × k + 1`; for cut-off frequency F, k is `sqrt(0.196202 + F^2) / F`, where F is a normalized frequency (`F = f/fs`)
- `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`); only samples below/above the threshold are being filtered
- `window::Union{Nothing, AbstractVector}=nothing`: weighting window

# Returns

- `obj_new::NeuroAnalyzer.NEURO`

# Source

1. https://dsp.stackexchange.com/questions/9966/what-is-the-cut-off-frequency-of-a-moving-average-filter
"""
function filter_mavg(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, k::Int64=8, t::Real=0, window::AbstractVector=ones(2 * k + 1))::NeuroAnalyzer.NEURO

    ch = get_channel(obj, ch=ch)
    _info("Window length: $(2 * k + 1) samples")
    _info("Approximate cut-off frequency: $(round(0.442947 / (sqrt((2 * k + 1)^2 - 1)), digits=2) * sr(obj)) Hz")

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = @views filter_mavg(obj.data[ch, :, :], k=k, t=t, window=window)
    reset_components!(obj_new)
    push!(obj_new.history, "filter_mavg(OBJ, ch=$ch, k=$k, t=$t, window=$window")

    return obj_new

end

"""
    filter_mavg!(obj; <keyword arguments>)

Filter using moving average filter (with threshold).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `k::Int64=8`: window length is `2 × k + 1`; for cut-off frequency f, k is `sqrt(0.196202 + f^2) / f`
- `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`)
- `window::Union{Nothing, AbstractVector}=nothing`: weighting window

# Returns

Nothing
"""
function filter_mavg!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, k::Int64=8, t::Real=0, window::AbstractVector=ones(2 * k + 1))::Nothing

    obj_new = filter_mavg(obj, ch=ch, k=k, t=t, window=window)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
