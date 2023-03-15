export denoise_wavelet
export denoise_wavelet!

"""
    denoise_wavelet(s; wt)

Perform wavelet denoising.

# Arguments

- `s::AbstractVector`
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets

# Returns

- `s_new::Vector{Float64}`
"""
function denoise_wavelet(s::AbstractVector; wt::T) where {T<:DiscreteWavelet}

    s_new = denoise(s, wt)

    return s_new

end

"""
    denoise_wavelet(s; wt)

Perform wavelet denoising.

# Arguments

- `s::AbstractArray`
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function denoise_wavelet(s::AbstractArray; wt::T) where {T<:DiscreteWavelet}

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views denoise_wavelet(s[ch_idx, :, ep_idx], wt=wt)
        end
    end

    return s_new

end

"""
    denoise_wavelet(obj; ch, wt)

Perform wavelet denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function denoise_wavelet(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), wt::T) where {T<:DiscreteWavelet}

    _check_channels(obj, ch)

    obj_new = deepcopy(obj)
    obj_new.data = @views denoise_wavelet(obj.data[ch, :, :], wt=wt)
    reset_components!(obj_new)
    push!(obj_new.header.history, "denoise_wavelet(OBJ, wt=$wt)")

    return obj_new

end

"""
    denoise_wavelet!(obj; ch, wt)

Perform wavelet denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
"""
function denoise_wavelet!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), wt::T) where {T<:DiscreteWavelet}

    obj_tmp = denoise_wavelet(obj, ch=ch, wt=wt)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return nothing

end
