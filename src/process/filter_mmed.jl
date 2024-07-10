export filter_mmed
export filter_mmed!

"""
    filter_mmed(s; <keyword arguments>)

Filter using moving median filter (with threshold).

# Arguments

- `s::AbstractVector`
- `k::Int64=8`: window length is `2 × k + 1`
- `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`); only samples below/above the threshold are being filtered
- `window::Union{Nothing, AbstractVector}=nothing`: weighting window

# Returns

- `s_filtered::Vector{Float64}`
"""
function filter_mmed(s::AbstractVector; k::Int64=8, t::Real=0, window::AbstractVector=ones(2 * k + 1))

    # check k
    @assert k in eachindex(s) "k must be in [1, signal length ($(length(s)))]."

    # check window
    @assert length(window) == (2 * k + 1) "window length must be 2 × k + 1 ($(2 * k + 1))."

    s_filtered = deepcopy(s)

    @inbounds for idx in (1 + k):(length(s) - k)
        if t > 0
            if s[idx] < mean(s) - t * std(s) || s[idx] > (mean(s) + t * std(s))
                s_filtered[idx] = @views median(s[(idx - k):(idx + k)] .* window)
            end
        else
            s_filtered[idx] = @views median(s[(idx - k):(idx + k)] .* window)
        end
    end

    return s_filtered

end

"""
    filter_mmed(s; <keyword arguments>)

Filter using moving median filter (with threshold).

# Arguments

- `s::AbstractArray`
- `k::Int64=8`: window length is `2 × k + 1`
- `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`); only samples below/above the threshold are being filtered
- `window::Union{Nothing, AbstractVector}=nothing`: weighting window

# Returns

- `s_filtered::Array{Float64, 3}`: convoluted signal
"""
function filter_mmed(s::AbstractArray; k::Int64=8, t::Real=0, window::AbstractVector=ones(2 * k + 1))

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_filtered = similar(s)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_filtered[ch_idx, :, ep_idx] = @views filter_mmed(s[ch_idx, :, ep_idx], k=k, t=t, window=window)
        end
    end

    return s_filtered

end

"""
    filter_mmed(obj; <keyword arguments>)

Filter using moving median filter (with threshold).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `k::Int64=8`: window length is `2 × k + 1`
- `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`); only samples above the threshold are being filtered
- `window::Union{Nothing, AbstractVector}=nothing`: weighting window

# Returns

- `obj_new::NeuroAnalyzer.NEURO`: convoluted signal
"""
function filter_mmed(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), k::Int64=8, t::Real=0, window::AbstractVector=ones(2 * k + 1))

    _check_channels(obj, ch)
    isa(ch, Int64) && (ch = [ch])

    _info("Window length: $(2 * k + 1) samples")

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = filter_mmed(obj.data[ch, :, :], k=k, t=t, window=window)
    reset_components!(obj_new)
    push!(obj_new.history, "filter_mmed(OBJ, ch=$ch, k=$k, t=$t, window=$window")

    return obj_new

end

"""
    filter_mmed!(obj; <keyword arguments>)

Filter using moving median filter (with threshold).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `k::Int64=8`: window length is `2 × k + 1`
- `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`); only samples above the threshold are being filtered
- `window::Union{Nothing, AbstractVector}=nothing`: weighting window

"""
function filter_mmed!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), k::Int64=8, t::Real=0, window::AbstractVector=ones(2 * k + 1))

    obj_new = filter_mmed(obj, ch=ch, k=k, t=t, window=window)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
