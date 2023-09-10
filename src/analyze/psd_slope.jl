export psd_slope

"""
    psd_slope(s; fs, f, norm, method, nt, wlen, woverlap, w, frq_lim, frq_n, frq, fs, ncyc)

Calculate PSD linear fit and slope. Default method is Welch periodogram.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `f::Tuple{Real, Real}=(0, fs / 2)`: calculate slope of the total power (default) or frequency range `f[1]` to `f[2]`
- `norm::Bool=false`: normalize do dB
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch periodogram
    - `:fft`: fast-Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=8`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window for Welch and STFT
- `frq_lim::Tuple{Real, Real}=(0, fs / 2)`: frequency bounds
- `frq_n::Int64=length(frq_lim[1]:frq_lim[2])`: number of frequencies
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

# Returns

Named tuple containing:
- `lf::Vector{Float64}`: linear fit
- `ls::Float64`: slopes of linear fit
- `pf::Vector{Float64}`: range of frequencies for the linear fit
"""
function psd_slope(s::AbstractVector; fs::Int64, f::Tuple{Real, Real}=(0, fs / 2), norm::Bool=false, method::Symbol=:welch, nt::Int64=8, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, frq_lim::Tuple{Real, Real}=(0, fs / 2), frq_n::Int64=_tlength(frq_lim), frq::Symbol=:lin, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)


    f = tuple_order(f)
    @assert f[1] >= 0 "Lower frequency bound must be be ≥ 0."
    @assert f[2] <= fs / 2 "Upper frequency bound must be be < $(fs / 2)."

    pw, pf = psd(s, fs=fs, norm=norm, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)

    f1_idx = vsearch(f[1], pf)
    f2_idx = vsearch(f[2], pf)
    _, _, _, _, _, _, lf = linreg(pf[f1_idx:f2_idx], pw[f1_idx:f2_idx])
    ls = lf[2] - lf[1]

    return (lf=lf, ls=ls, pf=pf[f1_idx:f2_idx])

end

"""
    psd_slope(s; fs, f, norm, method, nt, wlen, woverlap, w, frq_lim, frq_n, frq, fs, ncyc)

Calculate PSD linear fit and slope. Default method is Welch periodogram.

# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `f::Tuple{Real, Real}=(0, fs / 2)`: calculate slope of the total power (default) or frequency range `f[1]` to `f[2]`
- `norm::Bool=false`: normalize do dB
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch periodogram
    - `:fft`: fast-Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=8`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window for Welch and STFT
- `frq_lim::Tuple{Real, Real}=(0, fs / 2)`: frequency bounds
- `frq_n::Int64=length(frq_lim[1]:frq_lim[2])`: number of frequencies
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

# Returns

Named tuple containing:
- `lf::Matrix{Float64}`: linear fit
- `s::Vector{Float64}`: slope of linear fit
- `pf::Vector{Float64}`: range of frequencies for the linear fit
"""
function psd_slope(s::AbstractArray; fs::Int64, f::Tuple{Real, Real}=(0, fs / 2), norm::Bool=false, method::Symbol=:welch, nt::Int64=8, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, frq_lim::Tuple{Real, Real}=(0, fs / 2), frq_n::Int64=_tlength(frq_lim), frq::Symbol=:lin, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)


    ch_n = size(s, 1)
    ep_n = size(s, 3)

    lf, ls, pf = psd_slope(s[1, :, 1], fs=fs, f=f, norm=norm, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)

    lf = zeros(ch_n, length(lf), ep_n)
    ls = zeros(ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            lf[ch_idx, :, ep_idx], ls[ch_idx, ep_idx], _ = psd_slope(s[ch_idx, :, ep_idx], fs=fs, f=f, norm=norm, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)
        end
    end

    return (lf=lf, ls=ls, pf=pf)

end

"""
    psd_slope(obj; ch, f, norm, method, nt, wlen, woverlap, w, frq_lim, frq_n, frq, fs, ncyc)

Calculate PSD linear fit and slope. Default method is Welch periodogram.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `f::Tuple{Real, Real}=(0, sr(obj) / 2)`: calculate slope of the total power (default) or frequency range f[1] to f[2]
- `norm::Bool=false`: normalize do dB
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch periodogram
    - `:fft`: fast-Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=8`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window for Welch and STFT
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds
- `frq_n::Int64=length(frq_lim[1]:frq_lim[2])`: number of frequencies
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

# Returns

Named tuple containing:
- `lf::Array{Float64, 3}`: linear fit
- `ls::Array{Float64, 2}`: slope of linear fit
- `pf::Vector{Float64}`: range of frequencies for the linear fit
"""
function psd_slope(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), f::Tuple{Real, Real}=(0, sr(obj) / 2), norm::Bool=false, method::Symbol=:welch, nt::Int64=8, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), frq_n::Int64=_tlength(frq_lim), frq::Symbol=:lin, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)

    _check_channels(obj, ch)

    lf, ls, pf = psd_slope(obj.data[ch, :, :], fs=sr(obj), f=f, norm=norm, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)

    return (lf=lf, ls=ls, pf=pf)

end
