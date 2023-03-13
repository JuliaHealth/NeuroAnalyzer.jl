export denoise_wien
export denoise_wien!

"""
    denoise_wien(signal)

Perform Wiener deconvolution denoising.

# Arguments

- `signal::AbstractArray`

# Returns

- `signal_new::Vector{Float64}`
"""
function denoise_wien(signal::AbstractArray)

    ch_n, _, ep_n = size(signal)
    signal_new = similar(signal)

    @inbounds @simd for ep_idx in 1:ep_n
        s_m = @views mean(signal[:, :, ep_idx], dims=1)'[:, 1]
        m = mean(s_m)
        noise = rand(Float64, size(s_m)) .* m
        Threads.@threads for ch_idx in 1:ch_n
            signal_new[ch_idx, :, ep_idx] = @views wiener(signal[ch_idx, :, ep_idx], s_m, noise)
        end
    end

    return signal_new
end

"""
    denoise_wien(obj; channel)

Perform Wiener deconvolution denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels

# Returns
- `obj_new::NeuroAnalyzer.NEURO`
"""
function denoise_wien(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)))

    _check_channels(obj, channel)

    obj_new = deepcopy(obj)
    obj_new.data[channel, :, :] = denoise_wien(obj_new.data[channel, :, :])
    reset_components!(obj_new)
    push!(obj_new.header.history, "denoise_wien(OBJ, channel=$channel)")

    return obj_new
end

"""
    denoise_wien!(obj; channel)

Perform Wiener deconvolution denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
"""
function denoise_wien!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)))

    obj_tmp = denoise_wien(obj, channel=channel)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    reset_components!(obj)

    return nothing
end
