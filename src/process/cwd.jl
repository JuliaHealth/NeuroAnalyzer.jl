export cwd
export icwd

"""
    cwd(s; <keyword arguments>)

Perform continuous wavelet decomposition (CWD).

# Arguments

- `s::AbstractVector`
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

- `ct::Matrix{Float64}`: CWT coefficients (by rows)
"""
function cwd(s::AbstractVector; wt::T=wavelet(Morlet(2π), β=2))::Matrix{Float64} where {T<:CWT}

    ct = Matrix(real.(ContinuousWavelets.cwt(s, wt))')

    return ct

end

"""
    icwd(ct; <keyword arguments>)

Perform inverse continuous wavelet transformation (iCWT).

# Arguments

- `ct::AbstractArray`: CWT coefficients (by rows)
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `type::Symbol=:pd`: inverse style type:
    - `:pd`: PenroseDelta
    - `:nd`: NaiveDelta
    - `:df`: DualFrames

# Returns

- `s::AbstractArray`: reconstructed signal
"""
function icwd(ct::AbstractArray; wt::T=wavelet(Morlet(2π), β=2), type::Symbol=:pd)::AbstractArray where {T<:CWT}

    _check_var(type, [:nd, :pd, :df], "type")

    # reconstruct array of CWT coefficients as returned by ContinuousWavelets.jl functions
    ct = Matrix(ct')

    if type === :pd
        s = ContinuousWavelets.icwt(ct, wt, PenroseDelta())
    elseif type === :nd
        s = ContinuousWavelets.icwt(ct, wt, NaiveDelta())
    elseif type === :df
        s = ContinuousWavelets.icwt(ct, wt, DualFrames())
    end

    return s

end

"""
    cwd(s; <keyword arguments>)

Perform continuous wavelet decomposition (CWD).

# Arguments

- `s::AbstractArray`
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

- `ct::Array{Float64, 4}`: CWT coefficients (by rows)
"""
function cwd(s::AbstractArray; wt::T=wavelet(Morlet(2π), β=2))::Array{Float64, 4} where {T<:CWT}

    _chk3d(s)
    ch_n, ep_len, ep_n = size(s)

    _log_off()
    l = size(ContinuousWavelets.cwt(s[1, :, 1], wt), 2)
    ct = zeros(ch_n, l, ep_len, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            ct[ch_idx, :, :, ep_idx] = @views cwd(s[ch_idx, :, ep_idx], wt=wt)
        end
    end
    _log_on()

    return ct

end

"""
    cwd(obj; <keyword arguments>)

Perform continuous wavelet decomposition (CWD).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

- `ct::Array{Float64, 4}`: CWT coefficients (by rows)
"""
function cwd(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, wt::T=wavelet(Morlet(2π), β=2))::Array{Float64, 4} where {T<:CWT}

    ch = get_channel(obj, ch=ch)
    ct = @views cwd(obj.data[ch, :, :], wt=wt)

    return ct

end
