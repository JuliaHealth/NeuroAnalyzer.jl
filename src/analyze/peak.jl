export peak_frq
export peak_amp
export peak_pow

"""
    peak_frq(s; <keyword arguments>)

Calculate peak frequency in a band.

# Arguments

  - `s::AbstractVector`
  - `fs::Int64`: sampling rate
  - `flim::Tuple{Real, Real}`: lower and upper frequency bounds
  - `method::Symbol=:welch`: method used to calculate PSD:
      + `:welch`: Welch's periodogram
      + `:fft`: fast Fourier transform
      + `:mt`: multi-tapered periodogram
      + `:stft`: short time Fourier transform
      + `:mw`: Morlet wavelet convolution
  - `nt::Int64=16`: number of Slepian tapers
  - `wlen::Int64=fs`: window length (in samples), default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
  - `w::Bool=true`: if true, apply Hanning window
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
  - `demean::Bool=true`: subtract DC before calculating PSD

# Returns

  - `pf::Float64`: peak frequency
"""
function peak_frq(
        s::AbstractVector;
        fs::Int64,
        flim::Tuple{Real, Real},
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = fs,
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        demean::Bool = true,
    )::Float64

    @assert fs >= 1 "fs must be ≥ 1."
    _check_tuple(flim, (0, fs / 2), "flim")

    pw, pf = psd(s,
                fs = fs,
                db = false,
                method = method,
                nt = nt,
                wlen = wlen,
                woverlap = woverlap,
                w = w,
                ncyc = ncyc,
                demean = demean,
            )

    f1_idx = vsearch(flim[1], pf)
    f2_idx = vsearch(flim[2], pf)

    # peak power
    pp = maximum(pw[f1_idx:f2_idx])
    pf = pf[vsearch(pp, pw)]

    return pf

end

"""
    peak_amp(s; <keyword arguments>)

Calculate amplitude at peak frequency in a band.

# Arguments

  - `s::AbstractVector`
  - `fs::Int64`: sampling rate
  - `flim::Tuple{Real, Real}`: lower and upper frequency bounds
  - `method::Symbol=:welch`: method used to calculate PSD:
      + `:welch`: Welch's periodogram
      + `:fft`: fast Fourier transform
      + `:mt`: multi-tapered periodogram
      + `:stft`: short time Fourier transform
      + `:mw`: Morlet wavelet convolution
  - `nt::Int64=16`: number of Slepian tapers
  - `wlen::Int64=fs`: window length (in samples), default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
  - `w::Bool=true`: if true, apply Hanning window
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
  - `demean::Bool=true`: subtract DC before calculating PSD

# Returns

  - `pa::Float64`: amplitude at peak frequency
"""
function peak_amp(
        s::AbstractVector;
        fs::Int64,
        flim::Tuple{Real, Real},
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = fs,
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        demean::Bool = true,
    )::Float64

    @assert fs >= 1 "fs must be ≥ 1."
    _check_tuple(flim, (0, fs / 2), "flim")

    pw, pf = psd(s,
                fs = fs,
                db = false,
                method = method,
                nt = nt,
                wlen = wlen,
                woverlap = woverlap,
                w = w,
                ncyc = ncyc,
                demean = demean,
            )

    f1_idx = vsearch(flim[1], pf)
    f2_idx = vsearch(flim[2], pf)

    # peak amplitude
    pa = sqrt(maximum(pw[f1_idx:f2_idx]))

    return pa

end

"""
    peak_pow(s; <keyword arguments>)

Calculate power at peak frequency in a band.

# Arguments

  - `s::AbstractVector`
  - `fs::Int64`: sampling rate
  - `flim::Tuple{Real, Real}`: lower and upper frequency bounds
  - `method::Symbol=:welch`: method used to calculate PSD:
      + `:welch`: Welch's periodogram
      + `:fft`: fast Fourier transform
      + `:mt`: multi-tapered periodogram
      + `:stft`: short time Fourier transform
      + `:mw`: Morlet wavelet convolution
  - `nt::Int64=16`: number of Slepian tapers
  - `wlen::Int64=fs`: window length (in samples), default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
  - `w::Bool=true`: if true, apply Hanning window
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
  - `demean::Bool=true`: subtract DC before calculating PSD

# Returns

  - `pp::Float64`: power at peak frequency
"""
function peak_pow(
        s::AbstractVector;
        fs::Int64,
        flim::Tuple{Real, Real},
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = fs,
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        demean::Bool = true,
    )::Float64

    @assert fs >= 1 "fs must be ≥ 1."
    _check_tuple(flim, (0, fs / 2), "flim")

    pw, pf = psd(s,
                fs = fs,
                db = false,
                method = method,
                nt = nt,
                wlen = wlen,
                woverlap = woverlap,
                w = w,
                ncyc = ncyc,
                demean = demean,
            )

    f1_idx = vsearch(flim[1], pf)
    f2_idx = vsearch(flim[2], pf)

    # peak power
    pp = maximum(pw[f1_idx:f2_idx])

    return pp

end

"""
    peak_frq(s; <keyword arguments>)

Calculate peak frequency in a band.

# Arguments

  - `s::AbstractArray`
  - `fs::Int64`: sampling rate
  - `flim::Tuple{Real, Real}`: lower and upper frequency bounds
  - `method::Symbol=:welch`: method used to calculate PSD:
      + `:welch`: Welch's periodogram
      + `:fft`: fast Fourier transform
      + `:mt`: multi-tapered periodogram
      + `:stft`: short time Fourier transform
      + `:mw`: Morlet wavelet convolution
  - `nt::Int64=16`: number of Slepian tapers
  - `wlen::Int64=fs`: window length (in samples), default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
  - `w::Bool=true`: if true, apply Hanning window
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
  - `demean::Bool=true`: subtract DC before calculating PSD

# Returns

  - `pf::Matrix{Float64}`: peak frequency
"""
function peak_frq(
        s::AbstractArray;
        fs::Int64,
        flim::Tuple{Real, Real},
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = fs,
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        demean::Bool = true,
    )::Matrix{Float64}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)
    pf = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            pf[ch_idx, ep_idx] = @views peak_frq(
                s[ch_idx, :, ep_idx],
                fs = fs,
                flim = flim,
                method = method,
                nt = nt,
                wlen = wlen,
                woverlap = woverlap,
                w = w,
                ncyc = ncyc,
                demean = demean,
            )
        end
    end

    return pf

end

"""
    peak_amp(s; <keyword arguments>)

Calculate amplitude at peak frequency in a band.

# Arguments

  - `s::AbstractArray`
  - `fs::Int64`: sampling rate
  - `flim::Tuple{Real, Real}`: lower and upper frequency bounds
  - `method::Symbol=:welch`: method used to calculate PSD:
      + `:welch`: Welch's periodogram
      + `:fft`: fast Fourier transform
      + `:mt`: multi-tapered periodogram
      + `:stft`: short time Fourier transform
      + `:mw`: Morlet wavelet convolution
  - `nt::Int64=16`: number of Slepian tapers
  - `wlen::Int64=fs`: window length (in samples), default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
  - `w::Bool=true`: if true, apply Hanning window
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
  - `demean::Bool=true`: subtract DC before calculating PSD

# Returns

  - `pa::Matrix{Float64}`: amplitude at peak frequency
"""
function peak_amp(
        s::AbstractArray;
        fs::Int64,
        flim::Tuple{Real, Real},
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = fs,
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        demean::Bool = true,
    )::Matrix{Float64}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)
    pa = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            pa[ch_idx, ep_idx] = @views peak_amp(
                s[ch_idx, :, ep_idx],
                fs = fs,
                flim = flim,
                method = method,
                nt = nt,
                wlen = wlen,
                woverlap = woverlap,
                w = w,
                ncyc = ncyc,
                demean = demean,
            )
        end
    end

    return pa

end

"""
    peak_pow(s; <keyword arguments>)

Calculate power at peak frequency in a band.

# Arguments

  - `s::AbstractArray`
  - `fs::Int64`: sampling rate
  - `flim::Tuple{Real, Real}`: lower and upper frequency bounds
  - `method::Symbol=:welch`: method used to calculate PSD:
      + `:welch`: Welch's periodogram
      + `:fft`: fast Fourier transform
      + `:mt`: multi-tapered periodogram
      + `:stft`: short time Fourier transform
      + `:mw`: Morlet wavelet convolution
  - `nt::Int64=16`: number of Slepian tapers
  - `wlen::Int64=fs`: window length (in samples), default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
  - `w::Bool=true`: if true, apply Hanning window
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
  - `demean::Bool=true`: subtract DC before calculating PSD

# Returns

  - `pp::Matrix{Float64}`: power at peak frequency
"""
function peak_pow(
        s::AbstractArray;
        fs::Int64,
        flim::Tuple{Real, Real},
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = fs,
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        demean::Bool = true,
    )::Matrix{Float64}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)
    pp = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            pp[ch_idx, ep_idx] = @views peak_pow(
                s[ch_idx, :, ep_idx],
                fs = fs,
                flim = flim,
                method = method,
                nt = nt,
                wlen = wlen,
                woverlap = woverlap,
                w = w,
                ncyc = ncyc,
                demean = demean,
            )
        end
    end

    return pp

end

"""
    peak_frq(obj; <keyword arguments>)

Calculate peak frequency in a band.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `flim::Tuple{Real, Real}`: lower and upper frequency bounds
  - `method::Symbol=:welch`: method used to calculate PSD:
      + `:welch`: Welch's periodogram
      + `:fft`: fast Fourier transform
      + `:mt`: multi-tapered periodogram
      + `:stft`: short time Fourier transform
  - `nt::Int64=16`: number of Slepian tapers
  - `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
  - `w::Bool=true`: if true, apply Hanning window
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(sr(obj) / 2)`
  - `demean::Bool=true`: subtract DC before calculating PSD

# Returns

  - `pf::Matrix{Float64}`: peak frequency
"""
function peak_frq(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex},
        flim::Tuple{Real, Real},
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = sr(obj),
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        demean::Bool = true,
    )::Matrix{Float64}

    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")
    pf = @views peak_frq(
        obj.data[ch, :, :],
        fs = sr(obj),
        flim = flim,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        demean = demean,
    )

    return pf

end

"""
    peak_amp(obj; <keyword arguments>)

Calculate amplitude at peak frequency in a band.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `flim::Tuple{Real, Real}`: lower and upper frequency bounds
  - `method::Symbol=:welch`: method used to calculate PSD:
      + `:welch`: Welch's periodogram
      + `:fft`: fast Fourier transform
      + `:mt`: multi-tapered periodogram
      + `:stft`: short time Fourier transform
  - `nt::Int64=16`: number of Slepian tapers
  - `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
  - `w::Bool=true`: if true, apply Hanning window
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(sr(obj) / 2)`
  - `demean::Bool=true`: subtract DC before calculating PSD

# Returns

  - `pa::Matrix{Float64}`: amplitude at peak frequency
"""
function peak_amp(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex},
        flim::Tuple{Real, Real},
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = sr(obj),
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        demean::Bool = true,
    )::Matrix{Float64}

    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")
    pa = @views peak_amp(
        obj.data[ch, :, :],
        fs = sr(obj),
        flim = flim,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        demean = demean,
    )

    return pa

end

"""
    peak_pow(obj; <keyword arguments>)

Calculate power at peak frequency in a band.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `flim::Tuple{Real, Real}`: lower and upper frequency bounds
  - `method::Symbol=:welch`: method used to calculate PSD:
      + `:welch`: Welch's periodogram
      + `:fft`: fast Fourier transform
      + `:mt`: multi-tapered periodogram
      + `:stft`: short time Fourier transform
  - `nt::Int64=16`: number of Slepian tapers
  - `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
  - `w::Bool=true`: if true, apply Hanning window
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(sr(obj) / 2)`
  - `demean::Bool=true`: subtract DC before calculating PSD

# Returns

  - `pw::Matrix{Float64}`: power at peak frequency
"""
function peak_pow(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex},
        flim::Tuple{Real, Real},
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = sr(obj),
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        demean::Bool = true,
    )::Matrix{Float64}

    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")
    pw = @views peak_amp(
        obj.data[ch, :, :],
        fs = sr(obj),
        flim = flim,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        demean = demean,
    )

    return pw

end
