export psd_rel

"""
    psd_rel(s::AbstractVector; fs::Int64, norm::Bool=false, method::Symbol=:welch, nt::Int64=8, f::Union{Tuple{Real, Real}, Nothing}=nothing, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, frq_lim::Tuple{Real, Real}=(0, fs / 2), frq_n::Int64=_tlength(frq_lim), frq::Symbol=:lin, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)

Calculate relative power spectrum density. Default method is Welch periodogram.

# Arguments
- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `norm::Bool=false`: normalize do dB
- `f::Union{Tuple{Real, Real}, Nothing}=nothing`: frequency range to calculate relative power to; if nothing, than calculate relative to total power
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch periodogram
    - `:fft`: fast Fourier transform
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
- `pw::Vector{Float64}`: powers
- `pf::Vector{Float64}`: frequencies
"""
function psd_rel(s::AbstractVector; fs::Int64, norm::Bool=false, f::Union{Tuple{Real, Real}, Nothing}=nothing, method::Symbol=:welch, nt::Int64=8, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, frq_lim::Tuple{Real, Real}=(0, fs / 2), frq_n::Int64=_tlength(frq_lim), frq::Symbol=:lin, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)

    ref_pw = f === nothing ? total_power(s, fs=fs, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w) : band_power(s, fs=fs, f=f, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)

    pw, pf = psd(s, fs=fs, norm=norm, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)

    pw = pw / ref_pw

    return (pw=pw, pf=pf)

end

"""
    psd_rel(s; fs, norm, method, nt, wlen, woverlap, w, frq_lim, frq_n, frq, fs, ncyc)

Calculate relative power spectrum density. Default method is Welch periodogram.

# Arguments
- `s::AbstractMatrix`
- `fs::Int64`: sampling rate
- `norm::Bool=false`: normalize do dB
- `f::Union{Tuple{Real, Real}, Nothing}=nothing`: frequency range to calculate relative power to; if nothing, than calculate relative to total power
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch periodogram
    - `:fft`: fast Fourier transform
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
- `pw::Array{Float64, 3}`: powers
- `pf::Vector{Float64}`: frequencies
"""
function psd_rel(s::AbstractMatrix; fs::Int64, norm::Bool=false, f::Union{Tuple{Real, Real}, Nothing}=nothing, method::Symbol=:welch, nt::Int64=8, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, frq_lim::Tuple{Real, Real}=(0, fs / 2), frq_n::Int64=_tlength(frq_lim), frq::Symbol=:lin, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)

    ch_n = size(s, 1)

    _, pf = psd_rel(s[1, :, 1], fs=fs, norm=norm, f=f, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)

    pw = zeros(ch_n, length(pf))

    @inbounds @simd for ch_idx in 1:ch_n
        pw[ch_idx, :], _ = psd_rel(s[ch_idx, :], fs=fs, norm=norm, f=f, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)
    end
    
    return (pw=pw, pf=pf)

end

"""
    psd_rel(s; fs, norm, method, nt, wlen, woverlap, w, frq_lim, frq_n, frq, fs, ncyc)

Calculate relative power spectrum density. Default method is Welch periodogram.

# Arguments
- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `norm::Bool=false`: normalize do dB
- `f::Union{Tuple{Real, Real}, Nothing}=nothing`: frequency range to calculate relative power to; if nothing, than calculate relative to total power
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch periodogram
    - `:fft`: fast Fourier transform
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
- `pw::Array{Float64, 3}`: powers
- `pf::Vector{Float64}`: frequencies
"""
function psd_rel(s::AbstractArray; fs::Int64, norm::Bool=false, f::Union{Tuple{Real, Real}, Nothing}=nothing, method::Symbol=:welch, nt::Int64=8, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, frq_lim::Tuple{Real, Real}=(0, fs / 2), frq_n::Int64=_tlength(frq_lim), frq::Symbol=:lin, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    _, pf = psd_rel(s[1, :, 1], fs=fs, norm=norm, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)

    pw = zeros(ch_n, length(pf), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            pw[ch_idx, :, ep_idx], _ = psd_rel(s[ch_idx, :, ep_idx], fs=fs, norm=norm, f=f, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)
        end
    end
    
    return (pw=pw, pf=pf)

end

"""
    psd_rel(obj; ch, norm, method, nt, wlen, woverlap, w, frq_lim, frq_n, frq, fs, ncyc)

Calculate relative power spectrum density. Default method is Welch periodogram.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `norm::Bool=false`: normalize do dB
- `f::Union{Tuple{Real, Real}, Nothing}=nothing`: frequency range to calculate relative power to; if nothing, than calculate relative to total power
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch periodogram
    - `:fft`: fast Fourier transform
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
- `pw::Array{Float64, 3}`: powers
- `pf::Array{Float64, 3}`: frequencies
"""
function psd_rel(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), norm::Bool=false, method::Symbol=:welch, nt::Int64=8, f::Union{Tuple{Real, Real}, Nothing}=nothing, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), frq_n::Int64=_tlength(frq_lim), frq::Symbol=:lin, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)

    _check_channels(obj, ch)

    pw, pf = @views psd_rel(obj.data[ch, :, :], fs=sr(obj), f=f, norm=norm, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)

    return (pw=pw, pf=pf)

end
