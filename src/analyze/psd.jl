export psd
export mwpsd
export ghpsd

"""
    psd(s; <keyword arguments>)

Calculate power spectrum density. Default method is Welch's periodogram.

# Arguments

  - `s::Vector{Float64}`
  - `fs::Int64`: sampling rate
  - `db::Bool=false`: normalize do dB
  - `method::Symbol=:welch`: PSD method:
      + `:welch`: Welch's periodogram
      + `:fft`: fast Fourier transform
      + `:mt`: multi-tapered periodogram
      + `:stft`: short-time Fourier transform
      + `:mw`: Morlet wavelet convolution
      + `:gh`: Gaussian and Hilbert transform
  - `nt::Int64=16`: number of Slepian tapers
  - `wlen::Int64=fs`: window length (in samples), default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
  - `w::Bool=true`: if true, apply Hanning window
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
  - `gw::Real=5`: Gaussian width in Hz
  - `demean::Bool=true`: subtract DC before calculating PSD

# Returns

Named tuple containing:

  - `p::Vector{Float64}`: powers
  - `f::Vector{Float64}`: frequencies

# Notes

Setting `demean=true` reduces (but doesn't fully eliminate) DC contamination. Importantly, this mitigates the 0 Hz bin's value, but the bin is still present in the output.
"""
function psd(
        s::AbstractVector;
        fs::Int64,
        db::Bool = false,
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = fs,
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        gw::Real = 5,
        demean::Bool = true,
    )::@NamedTuple{p::Vector{Float64}, f::Vector{Float64}}

    _check_var(method, [:fft, :welch, :mt, :mw, :stft, :gh], "method")
    @assert nt >= 1 "nt must be ≥ 1."
    @assert fs >= 1 "fs must be ≥ 1."
    @assert wlen <= length(s) "wlen must be ≤ $(length(s))."
    @assert wlen >= 2 "wlen must be ≥ 2."
    @assert woverlap < wlen "woverlap must be < $(wlen)."
    @assert woverlap >= 0 "woverlap must be ≥ 0."

    n = length(s)

    if method === :mt
        demean && (s = remove_dc(s))
        w = w ? hanning(n) : ones(n)
        p = mt_pgram(
                s .* w,
                fs = fs,
                nw = ((nt + 1) ÷ 2),
                ntapers = nt,
            )
        pw = power(p)
        f = Vector(freq(p))
        p = pw[1:length(f)]
        db && (p = pow2db.(p))
    elseif method === :stft
        demean && (s = remove_dc(s))
        w = w ? DSP.hanning : nothing
        p = abs.(DSP.stft(
                        s,
                        wlen,
                        woverlap,
                        fs = fs,
                        window = w,
                    )
                )
        # average STFT segments along time
        p = vec(mean(p, dims = 2))
        # create frequencies vector
        f = linspace(0, fs / 2, length(p))
        db && (p = pow2db.(p))
    elseif method === :welch
        demean && (s = remove_dc(s))
        w = w ? DSP.hanning : nothing
        p = DSP.welch_pgram(
                        s,
                        wlen,
                        woverlap,
                        fs = fs,
                        window = w,
                    )
        f = Vector(freq(p))
        p = power(p)
        db && (p = pow2db.(p))
    elseif method === :fft
        w = w ? DSP.hanning : nothing
        p = DSP.periodogram(
                        s,
                        fs = fs,
                        window = w,
                    )
        pw = power(p)
        f = Vector(freq(p))
        p = pw[1:length(f)]
        db && (p = pow2db.(p))
    elseif method === :mw
        p, f = mwpsd(
                    s,
                    db = db,
                    fs = fs,
                    ncyc = ncyc,
                    w = w,
                    demean = demean,
                )
    elseif method === :gh
        p, f = ghpsd(
                    s,
                    fs = fs,
                    db = db,
                    gw = gw,
                    w = w,
                    demean = demean,
                )
    end

    return (p = p, f = f)

end

"""
    psd(s; <keyword arguments>)

Calculate power spectrum density. Default method is Welch's periodogram.

# Arguments

  - `s::AbstractMatrix`
  - `fs::Int64`: sampling rate
  - `db::Bool=false`: normalize do dB
  - `method::Symbol=:welch`: PSD method:
      + `:welch`: Welch's periodogram
      + `:fft`: fast Fourier transform
      + `:mt`: multi-tapered periodogram
      + `:stft`: short-time Fourier transform
      + `:mw`: Morlet wavelet convolution
      + `:gh`: Gaussian and Hilbert transform
  - `nt::Int64=16`: number of Slepian tapers
  - `wlen::Int64=fs`: window length (in samples), default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
  - `w::Bool=true`: if true, apply Hanning window
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
  - `demean::Bool=true`: subtract DC before calculating PSD

# Returns

Named tuple containing:

  - `p::Matrix{Float64}`: powers
  - `f::Vector{Float64}`: frequencies
"""
function psd(
        s::AbstractMatrix;
        fs::Int64,
        db::Bool = false,
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = fs,
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        gw::Real = 5,
        demean::Bool=true,
    )::@NamedTuple{p::Matrix{Float64}, f::Vector{Float64}}

    _, f = @views psd(
        s[1, :],
        fs = fs,
        db = db,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        gw = gw,
        demean = demean,
    )

    p = zeros(size(s, 1), length(f))

    Threads.@threads for ch_idx in axes(s, 1)
        @inbounds p[ch_idx, :], _ = @views psd(
            s[ch_idx, :],
            fs = fs,
            db = db,
            method = method,
            nt = nt,
            wlen = wlen,
            woverlap = woverlap,
            w = w,
            ncyc = ncyc,
            gw = gw,
            demean = demean,
        )
    end

    return (p = p, f = f)

end

"""
    psd(s; <keyword arguments>)

Calculate power spectrum density. Default method is Welch's periodogram.

# Arguments

  - `s::AbstractArray`
  - `fs::Int64`: sampling rate
  - `db::Bool=false`: normalize do dB
  - `method::Symbol=:welch`: PSD method:
      + `:welch`: Welch's periodogram
      + `:fft`: fast Fourier transform
      + `:mt`: multi-tapered periodogram
      + `:stft`: short-time Fourier transform
      + `:mw`: Morlet wavelet convolution
      + `:gh`: Gaussian and Hilbert transform
  - `nt::Int64=16`: number of Slepian tapers
  - `wlen::Int64=fs`: window length (in samples), default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
  - `w::Bool=true`: if true, apply Hanning window
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
  - `gw::Real=5`: Gaussian width in Hz
  - `demean::Bool=true`: subtract DC before calculating PSD

# Returns

Named tuple containing:

  - `p::Array{Float64, 3}`: powers
  - `f::Vector{Float64}`: frequencies
"""
function psd(
        s::AbstractArray;
        fs::Int64,
        db::Bool = false,
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = fs,
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        gw::Real = 5,
        demean::Bool = true,
    )::@NamedTuple{p::Array{Float64, 3}, f::Vector{Float64}}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    _, f = @views psd(
        s[1, :, 1],
        fs = fs,
        db = db,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        gw = gw,
        demean = demean,
    )

    p = zeros(ch_n, length(f), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            p[ch_idx, :, ep_idx], _ = @views psd(
                s[ch_idx, :, ep_idx],
                fs = fs,
                db = db,
                method = method,
                nt = nt,
                wlen = wlen,
                woverlap = woverlap,
                w = w,
                ncyc = ncyc,
                gw = gw,
                demean = demean,
            )
        end
    end

    return (p = p, f = f)

end

"""
    psd(obj; <keyword arguments>)

Calculate power spectrum density. Default method is Welch's periodogram.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name(s)
  - `db::Bool=false`: normalize do dB
  - `method::Symbol=:welch`: PSD method:
      + `:welch`: Welch's periodogram
      + `:fft`: fast Fourier transform
      + `:mt`: multi-tapered periodogram
      + `:stft`: short-time Fourier transform
      + `:mw`: Morlet wavelet convolution
      + `:gh`: Gaussian and Hilbert transform
  - `nt::Int64=16`: number of Slepian tapers
  - `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
  - `w::Bool=true`: if true, apply Hanning window
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
  - `gw::Real=5`: Gaussian width in Hz
  - `flim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds
  - `demean::Bool=true`: subtract DC before calculating PSD

# Returns

Named tuple containing:

  - `p::Array{Float64, 3}`: powers
  - `f::Vector{Float64}`: frequencies
"""
function psd(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex},
        db::Bool = false,
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = sr(obj),
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        gw::Real = 5,
        flim::Tuple{Real, Real} = (0, sr(obj) / 2),
        demean::Bool = true,
    )::@NamedTuple{p::Array{Float64, 3}, f::Vector{Float64}}

    _check_tuple(flim, (0, sr(obj) / 2), "flim")

    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")
    _log_off()
    p, f = @views psd(
        obj.data[ch, :, :],
        fs = sr(obj),
        db = db,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        gw = gw,
        demean = demean,
    )
    _log_on()

    f1 = vsearch(flim[1], f)
    f2 = vsearch(flim[2], f)
    f = f[f1:f2]
    p = p[:, f1:f2, :]

    return (p = p, f = f)

end

"""
    mwpsd(s; <keyword arguments>)

Calculate power spectrum using Morlet wavelet convolution.

# Arguments

  - `s::AbstractVector`
  - `pad::Int64=0`: number of zeros to append
  - `db::Bool=true`: normalize powers to dB
  - `fs::Int64`: sampling rate
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
  - `w::Bool=true`: if true, apply Hanning window
  - `demean::Bool=true`: subtract DC before calculating PSD

# Returns

Named tuple containing:

  - `p::Vector{Float64}`: powers
  - `f::Vector{Float64}`: frequencies
"""
function mwpsd(
        s::AbstractVector;
        pad::Int64 = 0,
        db::Bool = true,
        fs::Int64,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        w::Bool = true,
        demean::Bool = true,
    )::@NamedTuple{p::Vector{Float64}, f::Vector{Float64}}

    @assert fs >= 1 "fs must be ≥ 1."
    @assert pad >= 0 "pad must be ≥ 0."

    demean && (s = remove_dc(s))
    pad > 0 && (s = pad0(s, pad))

    w = w ? hanning(length(s)) : ones(length(s))

    flim = (0, fs / 2)
    nfrq = _tlength(flim)
    f = linspace(flim[1], flim[2], nfrq)

    if ncyc isa Int64
        @assert ncyc >= 1 "ncyc must be ≥ 1"
        ncyc = repeat([ncyc], nfrq)
    else
        @assert ncyc[1] >= 1 "ncyc[1] must be ≥ 1"
        @assert ncyc[2] >= 1 "ncyc[2] must be ≥ 1"
        ncyc = round.(Int64, logspace(ncyc[1], ncyc[2], nfrq))
    end

    p = zeros(length(f))
    @inbounds for frq_idx in 1:nfrq
        kernel = generate_morlet(fs, f[frq_idx], 1, ncyc = ncyc[frq_idx], complex = true)
        # w_conv = tconv(s .* w, kernel=kernel)
        w_conv = fconv(s .* w, kernel = kernel, norm = false)
        p[frq_idx] = median((abs.(w_conv)) .^ 2)
    end

    db && (p = pow2db.(p))

    return (p = p, f = f)

end

"""
    ghpsd(s; <keyword arguments>)

Calculate power spectrum using Gaussian and Hilbert transform.

# Arguments

  - `s::AbstractVector`
  - `fs::Int64`: sampling rate
  - `db::Bool=true`: normalize powers to dB
  - `gw::Real=5`: Gaussian width in Hz
  - `w::Bool=true`: if true, apply Hanning window
  - `demean::Bool=true`: subtract DC before calculating PSD

# Returns

Named tuple containing:

  - `p::Vector{Float64}`: powers
  - `f::Vector{Float64}`: frequencies
"""
function ghpsd(
        s::AbstractVector;
        fs::Int64,
        db::Bool = true,
        gw::Real = 5,
        w::Bool = true,
        demean::Bool = true,
    )::@NamedTuple{p::Vector{Float64}, f::Vector{Float64}}

    @assert fs >= 1 "fs must be ≥ 1."

    demean && (s = remove_dc(s))
    flim = (0, fs / 2)
    nfrq = _tlength(flim)
    f = linspace(flim[1], flim[2], nfrq)

    w = w ? hanning(length(s)) : ones(length(s))

    p = zeros(length(f), length(s))
    ph = zeros(length(f), length(s))

    @inbounds for frq_idx in eachindex(f)
        s_tmp = filter_g(s .* w, fs = fs, f = f[frq_idx], gw = gw)
        p[frq_idx, :] = (abs.(hilbert(s_tmp))) .^ 2
        ph[frq_idx, :] = DSP.angle.(hilbert(s_tmp))
    end

    p[p .== -Inf] .= minimum(p[p .!== -Inf])
    p[p .== +Inf] .= maximum(p[p .!== +Inf])

    p = vec(mean(p, dims = 2))

    db && (p = pow2db.(p))

    return (p = p, f = f)

end
