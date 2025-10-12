export psd_rel

"""
    psd_rel(s; <keyword arguments>)

Calculate relative power spectrum density. Default method is Welch's periodogram.

# Arguments
- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `db::Bool=false`: normalize do dB
- `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: frequency range to calculate relative power to; if nothing, than calculate relative to total power
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(fs / 2)`
- `gw::Real=5`: Gaussian width in Hz

# Returns

Named tuple containing:
- `p::Vector{Float64}`: powers
- `f::Vector{Float64}`: frequencies
"""
function psd_rel(s::AbstractVector; fs::Int64, db::Bool=false, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing, method::Symbol=:welch, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, gw::Real=5)::@NamedTuple{p::Vector{Float64}, f::Vector{Float64}}

    ref_pw = frq_lim === nothing ? total_power(s, fs=fs, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw) : band_power(s, fs=fs, frq_lim=frq_lim, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw)

    p, f = psd(s, fs=fs, db=db, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw)

    p = p ./ ref_pw

    return (p=p, f=f)

end

"""
    psd_rel(s; <keyword arguments>)

Calculate relative power spectrum density. Default method is Welch's periodogram.

# Arguments
- `s::AbstractMatrix`
- `fs::Int64`: sampling rate
- `db::Bool=false`: normalize do dB
- `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: frequency range to calculate relative power to; if nothing, than calculate relative to total power
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(fs / 2)`
- `gw::Real=5`: Gaussian width in Hz

# Returns

Named tuple containing:
- `p::Matrix{Float64}`: powers
- `f::Vector{Float64}`: frequencies
"""
function psd_rel(s::AbstractMatrix; fs::Int64, db::Bool=false, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing, method::Symbol=:welch, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, gw::Real=5)::@NamedTuple{p::Matrix{Float64}, f::Vector{Float64}}

    ch_n = size(s, 1)

    _, f = psd_rel(s[1, :, 1], fs=fs, db=db, frq_lim=frq_lim, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw)

    p = zeros(ch_n, length(f))

    @inbounds for ch_idx in 1:ch_n
        p[ch_idx, :], _ = psd_rel(s[ch_idx, :], fs=fs, db=db, frq_lim=frq_lim, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw)
    end

    return (p=p, f=f)

end

"""
    psd_rel(s; <keyword arguments>)

Calculate relative power spectrum density. Default method is Welch's periodogram.

# Arguments
- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `db::Bool=false`: normalize do dB
- `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: frequency range to calculate relative power to; if nothing, than calculate relative to total power
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(fs / 2)`
- `gw::Real=5`: Gaussian width in Hz

# Returns

Named tuple containing:
- `p::Array{Float64, 3}`: powers
- `f::Vector{Float64}`: frequencies
"""
function psd_rel(s::AbstractArray; fs::Int64, db::Bool=false, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing, method::Symbol=:welch, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, gw::Real=5)::@NamedTuple{p::Array{Float64, 3}, f::Vector{Float64}}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    _, f = psd_rel(s[1, :, 1], fs=fs, db=db, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw)

    p = zeros(ch_n, length(f), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            p[ch_idx, :, ep_idx], _ = psd_rel(s[ch_idx, :, ep_idx], fs=fs, db=db, frq_lim=frq_lim, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw)
        end
    end

    return (p=p, f=f)

end

"""
    psd_rel(obj; <keyword arguments>)

Calculate relative power spectrum density. Default method is Welch's periodogram.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `db::Bool=false`: normalize do dB
- `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: frequency range to calculate relative power to; if nothing, than calculate relative to total power
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`
- `gw::Real=5`: Gaussian width in Hz

# Returns

Named tuple containing:
- `p::Array{Float64, 3}`: powers
- `f::Vector{Float64}`: frequencies
"""
function psd_rel(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, db::Bool=false, method::Symbol=:welch, nt::Int64=7, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, gw::Real=5)::@NamedTuple{p::Array{Float64, 3}, f::Vector{Float64}}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    p, f = @views psd_rel(obj.data[ch, :, :], fs=sr(obj), frq_lim=frq_lim, db=db, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw)

    return (p=p, f=f)

end
