export dwtsplit

"""
    dwtsplit(obj; channel, wt, type, n)

Split OBJ channel into bands using discrete wavelet transformation (DWT).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64}`: channel number
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: transformation type: 
    - `:sdwt`: Stationary Wavelet Transforms
    - `:acdwt`: Autocorrelation Wavelet Transforms
- `n::Int64=0`: number of bands, default is maximum number of bands available or total transformation

# Returns
 
- `bands::Array{Float64, 4}`: bands from lowest to highest frequency (by rows)
"""
function dwtsplit(obj::NeuroAnalyzer.NEURO; channel::Int64, wt::T, type::Symbol, n::Int64=0) where {T <: DiscreteWavelet}

    n -= 1
    if n == 0
        n = maxtransformlevels(obj.data[1, :, 1])
        _info("Calculating DWT using maximum level: $n.")
    end
    n < 2 && throw(ArgumentError("n must be ≥ 2."))

    _check_channels(obj, channel)
    ep_n = epoch_n(obj)

    dwt_c = zeros((n + 1), epoch_len(obj), ep_n)
    Threads.@threads for ep_idx in 1:ep_n
        @inbounds dwt_c[:, :, ep_idx] = @views dwt(obj.data[channel, :, ep_idx], wt=wt, type=type, l=n)
    end
    
    bands = similar(dwt_c)
    bands[1, :, :] = dwt_c[1, :, :]
    bands[2:end, :, :] = dwt_c[end:-1:2, :, :]

    return bands
end
