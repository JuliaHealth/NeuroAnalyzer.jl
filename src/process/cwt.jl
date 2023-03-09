export cwt
export icwt

"""
    cwt(signal; wt, type, l)

Perform continuous wavelet transformation (CWT).

# Arguments

- `signal::AbstractVector`
- `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

- `cwt_c::Array{Float64, 2}`: CWT coefficients (by rows)
"""
function cwt(signal::AbstractVector; wt::T) where {T <: CWT}
    cwt_coefs = abs.(ContinuousWavelets.cwt(signal, wt))
    cwt_c = zeros(size(cwt_coefs, 2), size(cwt_coefs, 1))
    for idx in 1:size(cwt_coefs, 2)
        cwt_c[idx, :] = @views cwt_coefs[:, idx]
    end
    return cwt_c
end

"""
    icwt(dwt_coefs; wt, type)

Perform inverse continuous wavelet transformation (iCWT) of the `dwt_coefs`.

# Arguments

- `cwt_coefs::AbstractArray`: CWT coefficients (by rows)
- `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `type::Symbol=df`: inverse style type:
    - `:nd`: NaiveDelta
    - `:pd`: PenroseDelta
    - `:df`: DualFrames

# Returns

- `signal::Vector{Float64}`: reconstructed signal
"""
function icwt(cwt_coefs::AbstractArray; wt::T, type::Symbol) where {T <: CWT}

    _check_var(type, [:nd, :pd, :df], "type")

    # reconstruct array of CWT coefficients as returned by ContinuousWavelets.jl functions
    cwt_c = zeros(size(cwt_coefs, 2), size(cwt_coefs, 1))
    for idx in 1:size(cwt_coefs, 2)
        cwt_c[idx, :] = @views cwt_coefs[:, idx]
    end
    type === :nd && return ContinuousWavelets.icwt(cwt_c, wt, NaiveDelta())
    type === :pd && return ContinuousWavelets.icwt(cwt_c, wt, PenroseDelta())
    type === :df && return ContinuousWavelets.icwt(cwt_c, wt, DualFrames())
end

"""
    cwt(obj; channel, wt)

Perform continuous wavelet transformation (CWT).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns
 
- `cwt_c::Array{Float64, 4}`: CWT coefficients (by rows)
"""
function cwt(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), wt::T) where {T <: CWT}

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    l = size(ContinuousWavelets.cwt(obj.data[1, :, 1], wt), 2)
    cwt_c = similar(obj.data)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            cwt_c[ch_idx, :, :, ep_idx] = @views NeuroAnalyzer.cwt(obj.data[channel[ch_idx], :, ep_idx], wt=wt)
        end
    end

    return cwt_c
end
