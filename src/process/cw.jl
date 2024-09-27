export cw_trans
export icw_trans

"""
    cw_trans(s; <keyword arguments>)

Perform continuous wavelet transformation (CWT).

# Arguments

- `s::AbstractVector`
- `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

- `ct::Matrix{Float64}`: CWT coefficients (by rows)
"""
function cw_trans(s::AbstractVector; wt::T)::Matrix{Float64} where {T<:CWT}

    ct = Matrix(real.(ContinuousWavelets.cwt(s, wt))')

    return ct

end

"""
    icw_trans(ct; <keyword arguments>)

Perform inverse continuous wavelet transformation (iCWT).

# Arguments

- `ct::AbstractArray`: CWT coefficients (by rows)
- `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `type::Symbol=:pd`: inverse style type:
    - `:pd`: PenroseDelta
    - `:nd`: NaiveDelta
    - `:df`: DualFrames

# Returns

- `s::AbstractArray`: reconstructed signal
"""
function icw_trans(ct::AbstractArray; wt::T, type::Symbol=:pd)::AbstractArray where {T<:CWT}

    _check_var(type, [:nd, :pd, :df], "type")

    # reconstruct array of CWT coefficients as returned by ContinuousWavelets.jl functions
    ct = Matrix(ct')

    type === :pd && return ContinuousWavelets.icwt(ct, wt, PenroseDelta())
    type === :nd && return ContinuousWavelets.icwt(ct, wt, NaiveDelta())
    type === :df && return ContinuousWavelets.icwt(ct, wt, DualFrames())

end

"""
    cw_trans(s; <keyword arguments>)

Perform continuous wavelet transformation (CWT).

# Arguments

- `s::AbstractArray`
- `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

- `ct::Array{Float64, 4}`: CWT coefficients (by rows)
"""
function cw_trans(s::AbstractArray; wt::T)::Array{Float64, 4} where {T<:CWT}

    _chk3d(s)
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
    cw_trans(obj; <keyword arguments>)

Perform continuous wavelet transformation (CWT).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channel name or list of channel names
- `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

- `ct::Array{Float64, 4}`: CWT coefficients (by rows)
"""
function cw_trans(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, wt::T)::Array{Float64, 4} where {T<:CWT}

    ch = get_channel(obj, ch=ch)
    ct = @views cw_trans(obj.data[ch, :, :], wt=wt)

    return ct

end
