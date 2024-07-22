export dwtsplit

"""
    dwtsplit(obj; <keyword arguments>)

Split into bands using discrete wavelet transformation (DWT).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: transformation type:
    - `:sdwt`: Stationary Wavelet Transforms
    - `:acdwt`: Autocorrelation Wavelet Transforms
- `n::Int64=0`: number of bands, default is maximum number of bands available or total transformation

# Returns

- `b::Array{Float64, 4}`: bands from lowest to highest frequency (by rows)
"""
function dwtsplit(obj::NeuroAnalyzer.NEURO; ch::String, wt::T, type::Symbol, n::Int64=0) where {T<:DiscreteWavelet}

    n -= 1
    if n == -1
        n = maxtransformlevels(obj.data[1, :, 1])
        _info("Calculating DWT using maximum level: $n")
    end
    @assert n >= 2 "n must be ≥ 2."

    ch = _ch_idx(obj, ch)[1]
    ep_n = nepochs(obj)

    dt = zeros((n + 1), epoch_len(obj), ep_n)
    Threads.@threads for ep_idx in 1:ep_n
        @inbounds dt[:, :, ep_idx] = @views dw_trans(obj.data[ch, :, ep_idx], wt=wt, type=type, l=n)
    end

    b = similar(dt)
    b[1, :, :] = dt[1, :, :]
    b[2:end, :, :] = dt[end:-1:2, :, :]

    return b

end
