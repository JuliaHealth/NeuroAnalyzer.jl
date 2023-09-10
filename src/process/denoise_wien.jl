export denoise_wien
export denoise_wien!

"""
    denoise_wien(s)

Perform Wiener deconvolution denoising.

# Arguments

- `s::AbstractArray`

# Returns

- `s_new::Vector{Float64}`
"""
function denoise_wien(s::AbstractArray)

    ch_n, _, ep_n = size(s)
    s_new = similar(s)

    @inbounds @simd for ep_idx in 1:ep_n
        s_m = @views mean(s[:, :, ep_idx], dims=1)'[:, 1]
        m = mean(s_m)
        noise = rand(Float64, size(s_m)) .* m
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views wiener(s[ch_idx, :, ep_idx], s_m, noise)
        end
    end

    return s_new

end

"""
    denoise_wien(obj; ch)

Perform Wiener deconvolution denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels

# Returns
- `obj_new::NeuroAnalyzer.NEURO`
"""
function denoise_wien(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)))

    _check_channels(obj, ch)

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = denoise_wien(obj_new.data[ch, :, :])
    reset_components!(obj_new)
    push!(obj_new.history, "denoise_wien(OBJ, ch=$ch)")

    return obj_new
end

"""
    denoise_wien!(obj; ch)

Perform Wiener deconvolution denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
"""
function denoise_wien!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)))

    obj_new = denoise_wien(obj, ch=ch)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
