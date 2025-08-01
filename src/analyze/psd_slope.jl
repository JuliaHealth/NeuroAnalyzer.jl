export psd_slope

"""
    psd_slope(s; <keyword arguments>)

Calculate PSD linear fit and slope. Default method is Welch's periodogram.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `frq_lim::Tuple{Real, Real}=(0, fs / 2)`: calculate slope of the total power (default) or frequency range `frq_lim[1]` to `frq_lim[2]`
- `db::Bool=false`: normalize do dB
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet convolution
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(fs / 2)`
- `gw::Real=5`: Gaussian width in Hz
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

Named tuple containing:
- `lf::Vector{Float64}`: linear fit
- `ls::Float64`: slopes of linear fit
- `pf::Vector{Float64}`: range of frequencies for the linear fit
"""
function psd_slope(s::AbstractVector; fs::Int64, frq_lim::Tuple{Real, Real}=(0, fs / 2), db::Bool=false, method::Symbol=:welch, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, gw::Real=5, wt::T=wavelet(Morlet(2π), β=32, Q=128)) where {T <: CWT}

    _check_tuple(frq_lim, "frq_lim", (0, fs / 2))

    pw, pf = psd(s, fs=fs, db=db, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw, wt=wt)

    f1_idx = vsearch(frq_lim[1], pf)
    f2_idx = vsearch(frq_lim[2], pf)
    lr = NeuroAnalyzer.linreg(pf[f1_idx:f2_idx], pw[f1_idx:f2_idx])
    lf = lr.lf
    ls = lf[2] - lf[1]

    return (lf=lf, ls=ls, pf=pf[f1_idx:f2_idx])

end

"""
    psd_slope(s; <keyword arguments>)

Calculate PSD linear fit and slope. Default method is Welch's periodogram.

# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `frq_lim::Tuple{Real, Real}=(0, fs / 2)`: calculate slope of the total power (default) or frequency range `frq_lim[1]` to `frq_lim[2]`
- `db::Bool=false`: normalize do dB
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet convolution
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(fs / 2)`
- `gw::Real=5`: Gaussian width in Hz
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

Named tuple containing:
- `lf::Matrix{Float64}`: linear fit
- `s::Vector{Float64}`: slope of linear fit
- `pf::Vector{Float64}`: range of frequencies for the linear fit
"""
function psd_slope(s::AbstractArray; fs::Int64, frq_lim::Tuple{Real, Real}=(0, fs / 2), db::Bool=false, method::Symbol=:welch, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, gw::Real=5, wt::T=wavelet(Morlet(2π), β=32, Q=128)) where {T <: CWT}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    lf, ls, pf = psd_slope(s[1, :, 1], fs=fs, frq_lim=frq_lim, db=db, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw, wt=wt)

    lf = zeros(ch_n, length(lf), ep_n)
    ls = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            lf[ch_idx, :, ep_idx], ls[ch_idx, ep_idx], _ = psd_slope(s[ch_idx, :, ep_idx], fs=fs, frq_lim=frq_lim, db=db, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw, wt=wt)
        end
    end

    return (lf=lf, ls=ls, pf=pf)

end

"""
    psd_slope(obj; <keyword arguments>)

Calculate PSD linear fit and slope. Default method is Welch's periodogram.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: calculate slope of the total power (default) or frequency range frq_lim[1] to frq_lim[2]
- `db::Bool=false`: normalize do dB
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet convolution
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`
- `gw::Real=5`: Gaussian width in Hz
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

Named tuple containing:
- `lf::Array{Float64, 3}`: linear fit
- `ls::Matrix{Float64}`: slope of linear fit
- `pf::Vector{Float64}`: range of frequencies for the linear fit
"""
function psd_slope(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), db::Bool=false, method::Symbol=:welch, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, gw::Real=5, wt::T=wavelet(Morlet(2π), β=32, Q=128)) where {T <: CWT}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    _log_off()
    lf, ls, pf = psd_slope(obj.data[ch, :, :], fs=sr(obj), frq_lim=frq_lim, db=db, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw, wt=wt)
    _log_on()

    return (lf=lf, ls=ls, pf=pf)

end
