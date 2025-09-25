export denoise_cwd
export denoise_cwd!
export denoise_dwd
export denoise_dwd!

"""
    denoise_cwd(s; <keyword arguments>)

Perform denoising using continuous wavelet decomposition (CWD).

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `nf::Real`: noise frequency in Hz
- `w::Int64=5`: width (in Hz) of the area surrounding noise (from `nf - w` to `nf + w`)
- `type::Symbol=:nd`: inverse style type:
    - `:pd`: PenroseDelta
    - `:nd`: NaiveDelta
    - `:df`: DualFrames

# Returns

- `s_new::Vector{Float64}`
"""
function denoise_cwd(s::AbstractVector; fs::Int64, wt::T=wavelet(Morlet(2π), β=2), nf::Real, w::Int64=5, type::Symbol=:nd)::Vector{Float64} where {T <: CWT}

    @assert fs >= 1 "fs must be ≥ 1."
    @assert nf >= 1 "nf must be ≥ 1."
    @assert nf <= fs / 2 "nf must be ≤ $(fs / 2)."

    # perform continuous wavelet transformation
    s_new = cwd(s, wt=wt)
    f = cwtfrq(s, fs=fs, wt=wt)

    f_idx1 = vsearch(nf - w, f)
    f_idx2 = vsearch(nf + w, f)
    s_new[f_idx1:f_idx2, :] .= 0

    # reconstruct
    s_new = vec(icwd(s_new, wt=wt, type=type))

    return s_new

end

"""
    denoise_cwd(s; <keyword arguments>)

Perform denoising using continuous wavelet decomposition (CWD).

# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `nf::Real`: noise frequency in Hz
- `w::Int64=5`: width (in Hz) of the area surrounding noise (from `nf - w` to `nf + w`)
- `type::Symbol=:nd`: inverse style type:
    - `:pd`: PenroseDelta
    - `:nd`: NaiveDelta
    - `:df`: DualFrames

# Returns

- `s_new::Array{Float64, 3}`
"""
function denoise_cwd(s::AbstractArray; fs::Int64, wt::T=wavelet(Morlet(2π), β=2), nf::Real, w::Int64=5, type::Symbol=:nd)::Array{Float64, 3} where {T <: CWT}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)
    s_new = similar(s)

    _log_off()
    f = cwtfrq(s, fs=fs, wt=wt)
    _log_on()
    f_idx = vsearch(nf, f)
    f_idx1 = vsearch(nf - w, f)
    f_idx2 = vsearch(nf + w, f)
    _info("Noise at: $(f[f_idx]) Hz")
    _info("Noise width: $(f[f_idx1]) to $(f[f_idx2]) Hz")

    _log_off()
    # initialize progress bar
    progress_bar && (progbar = Progress(ep_n * ch_n, dt=1, barlen=20, color=:white))
    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views denoise_cwd(s[ch_idx, :, ep_idx], fs=fs, wt=wt, nf=nf, w=w, type=type)
            # update progress bar
            progress_bar && next!(progbar)
        end
    end
    _log_on()

    return s_new

end

"""
    denoise_cwd(obj; <keyword arguments>)

Perform denoising using continuous wavelet decomposition (CWD).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `nf::Real`: noise frequency in Hz
- `w::Int64=5`: width (in Hz) of the area surrounding noise (from `nf - w` to `nf + w`)
- `type::Symbol=:nd`: inverse style type:
    - `:pd`: PenroseDelta
    - `:nd`: NaiveDelta
    - `:df`: DualFrames

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function denoise_cwd(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, wt::T=wavelet(Morlet(2π), β=2), nf::Real, w::Int64=5, type::Symbol=:nd)::NeuroAnalyzer.NEURO where {T <: CWT}

    ch = get_channel(obj, ch=ch)
    obj_new = deepcopy(obj)
    obj_new.data = @views denoise_cwd(obj.data[ch, :, :], fs=sr(obj), wt=wt, nf=nf, w=w, type=type)
    reset_components!(obj_new)
    push!(obj_new.history, "denoise_cwd(OBJ, ch=$ch, wt=$wt, nf=$nf, w=$w, type=$type)")

    return obj_new

end

"""
    denoise_cwd!(obj; <keyword arguments>)

Perform denoising using continuous wavelet decomposition (CWD).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `nf::Real`: noise frequency in Hz
- `w::Int64=5`: width (in Hz) of the area surrounding noise (from `nf - w` to `nf + w`)
- `type::Symbol=:nd`: inverse style type:
    - `:pd`: PenroseDelta
    - `:nd`: NaiveDelta
    - `:df`: DualFrames

# Returns

Nothing
"""
function denoise_cwd!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, wt::T=wavelet(Morlet(2π), β=2), nf::Real, type::Symbol=:nd)::Nothing where {T <: CWT}

    obj_new = denoise_cwd(obj, ch=ch, wt=wt, nf=nf, type=type)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end

"""
    denoise_dwd(s; <keyword arguments>)

Perform denoising using discrete wavelet decomposition (DWD).

# Arguments

- `s::AbstractVector`
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets

# Returns

- `s_new::Vector{Float64}`
"""
function denoise_dwd(s::AbstractVector; wt::T)::Vector{Float64} where {T<:DiscreteWavelet}

    s_new = denoise(s, wt)

    return s_new

end

"""
    denoise_dwd(s; <keyword arguments>)

Perform denoising using discrete wavelet decomposition (DWD).

# Arguments

- `s::AbstractArray`
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets

# Returns

- `s_new::Array{Float64, 3}`
"""
function denoise_dwd(s::AbstractArray; wt::T)::Array{Float64, 3} where {T<:DiscreteWavelet}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views denoise_dwd(s[ch_idx, :, ep_idx], wt=wt)
        end
    end

    return s_new

end

"""
    denoise_dwd(obj; <keyword arguments>)

Perform denoising using discrete wavelet decomposition (DWD).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function denoise_dwd(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, wt::T)::NeuroAnalyzer.NEURO where {T<:DiscreteWavelet}

    ch = get_channel(obj, ch=ch)
    obj_new = deepcopy(obj)
    obj_new.data = @views denoise_dwd(obj.data[ch, :, :], wt=wt)
    reset_components!(obj_new)
    push!(obj_new.history, "denoise_dwd(OBJ, ch=$ch, wt=$wt)")

    return obj_new

end

"""
    denoise_dwd!(obj; <keyword arguments>)

Perform denoising using discrete wavelet decomposition (DWD).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets

# Returns

Nothing
"""
function denoise_dwd!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, wt::T)::Nothing where {T<:DiscreteWavelet}

    obj_new = denoise_dwd(obj, ch=ch, wt=wt)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
