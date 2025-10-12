export filter_mmed
export filter_mmed!

"""
    filter_mmed(s; <keyword arguments>)

Filter using moving median filter (with threshold).

# Arguments

- `s::AbstractVector`
- `k::Int64=8`: window length is `2 × k + 1`
- `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`); only samples below/above the threshold are being filtered
- `ww::Union{Nothing, AbstractVector}=nothing`: weighting window

# Returns

- `s_filtered::Vector{Float64}`
"""
function filter_mmed(s::AbstractVector; k::Int64=8, t::Real=0, ww::AbstractVector=ones(2 * k + 1))::Vector{Float64}

    # check k
    _in(k, (1, length(s)), "k")

    # check weighting window
    @assert length(ww) == (2 * k + 1) "ww length must be 2 × k + 1 ($(2 * k + 1))."

    s_filtered = zeros(length(s))

    @inbounds for idx in 1:k
        s_tmp = s[1:idx]
        if t > 0
            if s[idx] < mean(s) - t * std(s) || s[idx] > (mean(s) + t * std(s))
                s_filtered[idx] = @views median(s_tmp .* ww[1:length(s_tmp)])
            end
        else
            s_filtered[idx] = @views median(s_tmp .* ww[1:length(s_tmp)])
        end
    end

    @inbounds for idx in (1 + k):(length(s) - k)
        if t > 0
            if s[idx] < mean(s) - t * std(s) || s[idx] > (mean(s) + t * std(s))
                s_filtered[idx] = @views median(s[(idx - k):(idx + k)] .* ww)
            end
        else
            s_filtered[idx] = @views median(s[(idx - k):(idx + k)] .* ww)
        end
    end

    @inbounds for idx in (length(s) - k + 1):length(s)
        s_tmp = s[idx:length(s)]
        if t > 0
            if s[idx] < mean(s) - t * std(s) || s[idx] > (mean(s) + t * std(s))
                s_filtered[idx] = @views median(s_tmp .* ww[(end - length(s_tmp)):end])
            end
        else
            s_filtered[idx] = @views median(s_tmp .* ww[(end - length(s_tmp) + 1):end])
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
- `ww::Union{Nothing, AbstractVector}=nothing`: weighting window

# Returns

- `s_filtered::Array{Float64, 3}`
"""
function filter_mmed(s::AbstractArray; k::Int64=8, t::Real=0, ww::AbstractVector=ones(2 * k + 1))::Array{Float64, 3}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_filtered = similar(s)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_filtered[ch_idx, :, ep_idx] = @views filter_mmed(s[ch_idx, :, ep_idx], k=k, t=t, ww=ww)
        end
    end

    return s_filtered

end

"""
    filter_mmed(obj; <keyword arguments>)

Filter using moving median filter (with threshold).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `k::Int64=8`: window length is `2 × k + 1`
- `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`); only samples above the threshold are being filtered
- `ww::Union{Nothing, AbstractVector}=nothing`: weighting window

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function filter_mmed(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, k::Int64=8, t::Real=0, ww::AbstractVector=ones(2 * k + 1))::NeuroAnalyzer.NEURO

    ch = get_channel(obj, ch=ch)
    _info("Window length: $(2 * k + 1) samples")

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = filter_mmed(obj.data[ch, :, :], k=k, t=t, ww=ww)
    reset_components!(obj_new)
    push!(obj_new.history, "filter_mmed(OBJ, ch=$ch, k=$k, t=$t, ww=$ww")

    return obj_new

end

"""
    filter_mmed!(obj; <keyword arguments>)

Filter using moving median filter (with threshold).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `k::Int64=8`: window length is `2 × k + 1`
- `t::Real=0`: threshold (`t = mean(s) - t * std(s):mean(s) + t * std(s)`); only samples above the threshold are being filtered
- `ww::Union{Nothing, AbstractVector}=nothing`: weighting window

# Returns

Nothing
"""
function filter_mmed!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, k::Int64=8, t::Real=0, ww::AbstractVector=ones(2 * k + 1))::Nothing

    obj_new = filter_mmed(obj, ch=ch, k=k, t=t, ww=ww)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
