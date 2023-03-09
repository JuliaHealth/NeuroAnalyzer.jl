export dwt
export idwt

"""
    dwt(signal; wt, type, l)

Perform discrete wavelet transformation (DWT).

# Arguments

- `signal::AbstractVector`
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: transformation type: Stationary Wavelet Transforms (`:sdwt`) or Autocorrelation Wavelet Transforms (`:acdwt`)
- `l::Int64=0`: number of levels, default maximum number of levels available or total transformation

# Returns

- `dwt_c::Array{Float64, 2}`: DWT coefficients cAl, cD1, ..., cDl (by rows)
"""
function dwt(signal::AbstractVector; wt::T, type::Symbol, l::Int64=0) where {T <: DiscreteWavelet}
    _check_var(type, [:sdwt, :acdwt], "type")

    l < 0 && throw(ArgumentError("l must be > 0."))
    l > maxtransformlevels(signal) && throw(ArgumentError("l must be ≤ $(maxtransformlevels(signal))."))

    if l == 0
        l = maxtransformlevels(signal)
        _info("Calculating DWT using maximum level: $l.")
    end

    if type === :sdwt
        dwt_coefs = sdwt(signal, wt, l)
    elseif type === :acdwt
        dwt_coefs = acdwt(signal, wt, l)
    end

    dwt_c = zeros(size(dwt_coefs, 2), size(dwt_coefs, 1))
    dwt_c[1, :] = @view dwt_coefs[:, 1]
    @inbounds @simd for idx in 2:(l + 1)
        dwt_c[idx, :] = @views dwt_coefs[:, (end - idx + 2)]
    end

    return dwt_c
end

"""
    idwt(dwt_coefs; wt, type)

Perform inverse discrete wavelet transformation (iDWT) of the `dwt_coefs`.

# Arguments

- `dwt_coefs::AbstractArray`: DWT coefficients cAl, cD1, ..., cDl (by rows)
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: transformation type: Stationary Wavelet Transforms (`:sdwt`) or Autocorrelation Wavelet Transforms (`:acdwt`)

# Returns

- `signal::Vector{Float64}`: reconstructed signal
"""
function idwt(dwt_coefs::AbstractArray; wt::T, type::Symbol) where {T <: DiscreteWavelet}
    
    _check_var(type, [:sdwt, :acdwt], "type")

    # reconstruct array of DWT coefficients as returned by Wavelets.jl functions
    dwt_c = zeros(size(dwt_coefs, 2), size(dwt_coefs, 1))
    dwt_c[:, 1] = @view dwt_coefs[1, :]
    @inbounds @simd for idx in 2:size(dwt_coefs, 1)
        dwt_c[:, idx] = @views dwt_coefs[(end - idx + 2), :]
    end

    if type === :sdwt
        return isdwt(dwt_c, wt)
    elseif type === :acdwt
        return iacdwt(dwt_c, wt)
    end
end

"""
    dwt(obj; channel, wt, type, l)

Perform discrete wavelet transformation (DWT).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: transformation type: 
    - `:sdwt`: Stationary Wavelet Transforms
    - `:acdwt`: Autocorrelation Wavelet Transforms
- `l::Int64=0`: number of levels, default is maximum number of levels available or total transformation

# Returns
 
- `dwt_c::Array{Float64, 4}`: DWT coefficients cAl, cD1, ..., cDl (by rows)
"""
function dwt(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), wt::T, type::Symbol, l::Int64=0) where {T <: DiscreteWavelet}

    if l == 0
        l = maxtransformlevels(obj.data[1, :, 1])
        _info("Calculating DWT using maximum level: $l.")
    end

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    dwt_c = zeros(ch_n, (l + 1), epoch_len(obj), ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            dwt_c[ch_idx, :, :, ep_idx] = @views NeuroAnalyzer.dwt(obj.data[channel[ch_idx], :, ep_idx], wt=wt, type=type, l=l)
        end
    end

    return dwt_c
end
