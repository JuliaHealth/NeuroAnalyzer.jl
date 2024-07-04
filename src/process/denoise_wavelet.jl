export denoise_cwt
export denoise_cwt!
export denoise_dwt
export denoise_dwt!

"""
    denoise_cwt(ct; wt, type)

Perform denoising using continuous wavelet transformation (iCWT).

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `wt::T where {T <: CWT}=wavelet(Morlet(π), β=32, Q=128)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `nf::Real`: noise frequency [Hz]
- `type::Symbol=:nd`: inverse style type:
    - `:pd`: PenroseDelta
    - `:nd`: NaiveDelta
    - `:df`: DualFrames

# Returns

- `s_new::Vector{Float64}`: denoised signal
"""
function denoise_cwt(s::AbstractVector; fs::Int64, wt::T=wavelet(Morlet(π), β=32, Q=128), nf::Real, type::Symbol=:nd) where {T <: CWT}

    @assert fs >= 1 "fs must be ≥ 1."
    @assert nf >= 1 "nf must be ≥ 1."
    @assert nf <= fs / 2 "nf must be ≤ $(fs / 2)."

    # perform continuous wavelet transformation
    s_new = cw_trans(s, wt=wt)

    f = round.(ContinuousWavelets.getMeanFreq(length(s), wt, fs), digits=2)

    f_idx1 = vsearch(nf - 10, f)
    f_idx2 = vsearch(nf + 10, f)
    s_new[f_idx1:f_idx2, :] .= 0

    # reconstruct
    s_new = vec(icw_trans(s_new, wt=wt, type=type))

    return s_new

end

"""
    denoise_cwt(s; wt)

Perform denoising using continuous wavelet transformation (CWT).

# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `wt::T where {T <: CWT}=wavelet(Morlet(π), β=32, Q=128)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `nf::Real`: noise frequency [Hz]
- `type::Symbol=:nd`: inverse style type:
    - `:pd`: PenroseDelta
    - `:nd`: NaiveDelta
    - `:df`: DualFrames

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function denoise_cwt(s::AbstractArray; fs::Int64, wt::T=wavelet(Morlet(π), β=32, Q=128), nf::Real, type::Symbol=:nd) where {T <: CWT}

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)

    # initialize progress bar
    progress_bar && (progbar = Progress(ep_n * ch_n, dt=1, barlen=20, color=:white))

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views denoise_cwt(s[ch_idx, :, ep_idx], fs=fs, wt=wt, nf=nf, type=type)
            # update progress bar
            progress_bar && next!(progbar)
        end
    end

    return s_new

end

"""
    denoise_cwt(obj; ch, wt)

Perform denoising using continuous wavelet transformation (CWT).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
- `wt::T where {T <: CWT}=wavelet(Morlet(π), β=32, Q=128)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `nf::Real`: noise frequency [Hz]
- `type::Symbol=:nd`: inverse style type:
    - `:pd`: PenroseDelta
    - `:nd`: NaiveDelta
    - `:df`: DualFrames

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function denoise_cwt(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)), wt::T=wavelet(Morlet(π), β=32, Q=128), nf::Real, type::Symbol=:nd) where {T <: CWT}

    _check_channels(obj, ch)
    isa(ch, Int64) && (ch = [ch])

    obj_new = deepcopy(obj)
    obj_new.data = @views denoise_cwt(obj.data[ch, :, :], fs=sr(obj), wt=wt, nf=nf, type=type)
    reset_components!(obj_new)
    push!(obj_new.history, "denoise_cwt(OBJ, ch=$ch, wt=$wt, nf=$nf, type=$type)")

    return obj_new

end

"""
    denoise_cwt!(obj; ch, wt)

Perform denoising using continuous wavelet transformation (CWT).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
- `wt::T where {T <: CWT}=wavelet(Morlet(π), β=32, Q=128)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `nf::Real`: noise frequency [Hz]
- `type::Symbol=:nd`: inverse style type:
    - `:pd`: PenroseDelta
    - `:nd`: NaiveDelta
    - `:df`: DualFrames
"""
function denoise_cwt!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)), wt::T=wavelet(Morlet(π), β=32, Q=128), nf::Real, type::Symbol=:nd) where {T <: CWT}

    obj_new = denoise_cwt(obj, ch=ch, wt=wt, nf=nf, type=type)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end

"""
    denoise_dwt(s; wt)

Perform wavelet denoising.

# Arguments

- `s::AbstractVector`
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets

# Returns

- `s_new::Vector{Float64}`
"""
function denoise_dwt(s::AbstractVector; wt::T) where {T<:DiscreteWavelet}

    s_new = denoise(s, wt)

    return s_new

end

"""
    denoise_dwt(s; wt)

Perform wavelet denoising.

# Arguments

- `s::AbstractArray`
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function denoise_dwt(s::AbstractArray; wt::T) where {T<:DiscreteWavelet}

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views denoise_dwt(s[ch_idx, :, ep_idx], wt=wt)
        end
    end

    return s_new

end

"""
    denoise_dwt(obj; ch, wt)

Perform wavelet denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function denoise_dwt(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)), wt::T) where {T<:DiscreteWavelet}

    _check_channels(obj, ch)
    isa(ch, Int64) && (ch = [ch])

    obj_new = deepcopy(obj)
    obj_new.data = @views denoise_dwt(obj.data[ch, :, :], wt=wt)
    reset_components!(obj_new)
    push!(obj_new.history, "denoise_dwt(OBJ, ch=$ch, wt=$wt)")

    return obj_new

end

"""
    denoise_dwt!(obj; ch, wt)

Perform wavelet denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
"""
function denoise_dwt!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)), wt::T) where {T<:DiscreteWavelet}

    obj_new = denoise_dwt(obj, ch=ch, wt=wt)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
