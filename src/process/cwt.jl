export cw_trans
export icw_trans

"""
    cw_trans(s; wt, type, l)

Perform continuous wavelet transformation (CWT).

# Arguments

- `s::AbstractVector`
- `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

- `ct::Array{Float64, 2}`: CWT coefficients (by rows)
"""
function cw_trans(s::AbstractVector; wt::T) where {T<:CWT}

    cwt_coefs = abs.(ContinuousWavelets.cwt(s, wt))
    ct = zeros(size(cwt_coefs, 2), size(cwt_coefs, 1))

    for idx in 1:size(cwt_coefs, 2)
        ct[idx, :] = @views cwt_coefs[:, idx]
    end

    return ct

end

"""
    icw_trans(ct; wt, type)

Perform inverse continuous wavelet transformation (iCWT).

# Arguments

- `ct::AbstractArray`: CWT coefficients (by rows)
- `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `type::Symbol=df`: inverse style type:
    - `:nd`: NaiveDelta
    - `:pd`: PenroseDelta
    - `:df`: DualFrames

# Returns

- `s::Vector{Float64}`: reconstructed signal
"""
function icw_trans(ct::AbstractArray; wt::T, type::Symbol) where {T<:CWT}

    _check_var(type, [:nd, :pd, :df], "type")

    # reconstruct array of CWT coefficients as returned by ContinuousWavelets.jl functions
    cwt_c = zeros(size(ct, 2), size(ct, 1))
    for idx in 1:size(ct, 2)
        cwt_c[idx, :] = @views ct[:, idx]
    end

    type === :nd && return ContinuousWavelets.icwt(cwt_c, wt, NaiveDelta())
    type === :pd && return ContinuousWavelets.icwt(cwt_c, wt, PenroseDelta())
    type === :df && return ContinuousWavelets.icwt(cwt_c, wt, DualFrames())

end

"""
    cw_trans(s; wt)

Perform continuous wavelet transformation (CWT).

# Arguments

- `s::AbstractArray`
- `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

- `ct::Array{Float64, 4}`: CWT coefficients (by rows)
"""
function cw_trans(s::AbstractArray; wt::T) where {T<:CWT}

    ch_n, ep_len, ep_n = size(s)

    l = size(ContinuousWavelets.cwt(s[1, :, 1], wt), 2)
    ct = zeros(ch_n, l, ep_len, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ct[ch_idx, :, :, ep_idx] = @views cw_trans(s[ch_idx, :, ep_idx], wt=wt)
        end
    end

    return ct

end

"""
    cw_trans(obj; ch, wt)

Perform continuous wavelet transformation (CWT).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

- `ct::Array{Float64, 4}`: CWT coefficients (by rows)
"""
function cw_trans(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), wt::T) where {T<:CWT}

    _check_channels(obj, ch)

    ct = @views cw_trans(obj.data[ch, :, :], wt=wt)

    return ct

end
