export wdenoise
export wdenoise!

"""
    wdenoise(signal; wt)

Perform wavelet denoising.

# Arguments

- `signal::AbstractVector`
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`

# Returns

- `signal_denoised::Vector{Float64}`
"""
function wdenoise(signal::AbstractVector; wt::T) where {T <: DiscreteWavelet}

    # wt in [:db2, :db4, :db8, :db10, :haar, :coif2, :coif4, :coif8] || throw(ArgumentError("wt must be :db2, :db4, :db8, :db10, :haar, :coif2, :coif4, :coif8"))

    # wt === :db2 && (wt = wavelet(WT.db2))
    # wt === :db4 && (wt = wavelet(WT.db4))
    # wt === :db8 && (wt = wavelet(WT.db8))
    # wt === :db10 && (wt = wavelet(WT.db10))
    # wt === :haar && (wt = wavelet(WT.haar))
    # wt === :coif2 && (wt = wavelet(WT.coif2))
    # wt === :coif4 && (wt = wavelet(WT.coif4))
    # wt === :coif8 && (wt = wavelet(WT.coif8))

    return denoise(signal, wt)
end

"""
    wdenoise(obj; channel, wt)

Perform wavelet denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function wdenoise(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), wt::T) where {T <: DiscreteWavelet}

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    obj_new = deepcopy(obj)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            obj_new.data[channel[ch_idx], :, ep_idx] = @views wdenoise(obj_new.data[channel[ch_idx], :, ep_idx], wt=wt)
        end
    end

    reset_components!(obj_new)
    push!(obj_new.header.history, "wdenoise(OBJ, wt=$wt)")

    return obj_new
end

"""
    wdenoise!(obj; channel, wt)

Perform wavelet denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
"""
function wdenoise!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), wt::T) where {T <: DiscreteWavelet}

    obj_tmp = wdenoise(obj, channel=channel, wt=wt)
    obj_tmp = upsample(obj, new_sr=new_sr)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    reset_components!(obj)

    return nothing
end
