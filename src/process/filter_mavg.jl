export filter_mavg
export filter_mavg!

"""
    filter_mavg(s; <keyword arguments>)

Filter using moving average filter (with threshold).

# Arguments

- `s::AbstractVector`
- `k::Int64=8`: window length is `2 × k + 1`; for cutoff frequency F, k is `sqrt(0.196202 + F^2) / F`, where F is a normalized frequency (`F = f/fs`)
- `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`); only samples below/above the threshold are being filtered
- `ww::Union{Nothing, AbstractVector}=nothing`: weighting window

# Returns

- `s_filtered::Vector{Float64}`
"""
function filter_mavg(s::AbstractVector; k::Int64=8, t::Real=0, ww::AbstractVector=ones(2 * k + 1))::Vector{Float64}

    # check k
    _in(k, (1, length(s)), "k")

    # check weighting window
    @assert length(ww) == (2 * k + 1) "ww length must be 2 × k + 1 ($(2 * k + 1))."

    s_filtered = zeros(length(s))

    @inbounds for idx in 1:k
        s_tmp = s[1:idx]
        if t > 0
            if s[idx] < mean(s) - t * std(s) || s[idx] > (mean(s) + t * std(s))
                s_filtered[idx] = @views mean(s_tmp .* ww[1:length(s_tmp)])
            end
        else
            s_filtered[idx] = @views mean(s_tmp .* ww[1:length(s_tmp)])
        end
    end

    @inbounds for idx in (1 + k):(length(s) - k)
        if t > 0
            if s[idx] < mean(s) - t * std(s) || s[idx] > (mean(s) + t * std(s))
                s_filtered[idx] = @views mean(s[(idx - k):(idx + k)] .* ww)
            end
        else
            s_filtered[idx] = @views mean(s[(idx - k):(idx + k)] .* ww)
        end
    end

    @inbounds for idx in (length(s) - k + 1):length(s)
        s_tmp = s[idx:length(s)]
        if t > 0
            if s[idx] < mean(s) - t * std(s) || s[idx] > (mean(s) + t * std(s))
                s_filtered[idx] = @views mean(s_tmp .* ww[(end - length(s_tmp)):end])
            end
        else
            s_filtered[idx] = @views mean(s_tmp .* ww[(end - length(s_tmp) + 1):end])
        end
    end

    return s_filtered

end

"""
    filter_mavg(s; <keyword arguments>)

Filter using moving average filter (with threshold).

# Arguments

- `s::AbstractArray`
- `k::Int64=8`: window length is `2 × k + 1`; for cutoff frequency F, k is `sqrt(0.196202 + F^2) / F`, where F is a normalized frequency (`F = f/fs`)
- `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`); only samples below/above the threshold are being filtered
- `ww::Union{Nothing, AbstractVector}=nothing`: weighting window

# Returns

- `s_filtered::Array{Float64, 3}`
"""
function filter_mavg(s::AbstractArray; k::Int64=8, t::Real=0, ww::AbstractVector=ones(2 * k + 1))::Array{Float64, 3}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_filtered = similar(s)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            s_filtered[ch_idx, :, ep_idx] = @views filter_mavg(s[ch_idx, :, ep_idx], k=k, t=t, ww=ww)
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
- `k::Int64=8`: window length is `2 × k + 1`; for cutoff frequency F, k is `sqrt(0.196202 + F^2) / F`, where F is a normalized frequency (`F = f/fs`)
- `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`); only samples below/above the threshold are being filtered
- `ww::Union{Nothing, AbstractVector}=nothing`: weighting window

# Returns

- `obj_new::NeuroAnalyzer.NEURO`

# Source

1. https://dsp.stackexchange.com/questions/9966/what-is-the-cutoff-frequency-of-a-moving-average-filter
"""
function filter_mavg(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, k::Int64=8, t::Real=0, ww::AbstractVector=ones(2 * k + 1))::NeuroAnalyzer.NEURO

    ch = get_channel(obj, ch=ch)
    _info("Window length: $(2 * k + 1) samples")
    _info("Approximate cutoff frequency: $(round(0.442947 / (sqrt((2 * k + 1)^2 - 1)), digits=2) * sr(obj)) Hz")
    _info("1st zero at: $(round(sr(obj) / k, digits=2)) Hz")
    _info("2nd zero at: $(round(2 * sr(obj) / k, digits=2)) Hz")
    _info("3rd zero at: $(round(3 * sr(obj) / k, digits=2)) Hz")
    _info("4th zero at: $(round(4 * sr(obj) / k, digits=2)) Hz")
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = @views filter_mavg(obj.data[ch, :, :], k=k, t=t, ww=ww)
    reset_components!(obj_new)
    push!(obj_new.history, "filter_mavg(OBJ, ch=$ch, k=$k, t=$t, ww=$ww")

    return obj_new

end

"""
    filter_mavg!(obj; <keyword arguments>)

Filter using moving average filter (with threshold).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `k::Int64=8`: window length is `2 × k + 1`; for cutoff frequency F, k is `sqrt(0.196202 + F^2) / F`, where F is a normalized frequency (`F = f/fs`)
- `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`)
- `ww::Union{Nothing, AbstractVector}=nothing`: weighting window

# Returns

Nothing
"""
function filter_mavg!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, k::Int64=8, t::Real=0, ww::AbstractVector=ones(2 * k + 1))::Nothing

    obj_new = filter_mavg(obj, ch=ch, k=k, t=t, ww=ww)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
