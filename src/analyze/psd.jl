export psd
export mwpsd

"""
    psd(s; fs, db, method, nt, wlen, woverlap, w, ncyc)

Calculate power spectrum density. Default method is Welch's periodogram.

# Arguments
- `s::Vector{Float64}`
- `fs::Int64`: sampling rate
- `db::Bool=false`: normalize do dB
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(fs / 2)`

# Returns

Named tuple containing:
- `p::Vector{Float64}`: powers
- `f::Vector{Float64}`: frequencies
"""
function psd(s::AbstractVector; fs::Int64, db::Bool=false, method::Symbol=:welch, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32)

    _check_var(method, [:fft, :welch, :mt, :mw, :stft], "method")
    @assert nt >= 1 "nt must be ≥ 1."
    @assert fs >= 1 "fs must be ≥ 1."
    @assert wlen <= length(s) "wlen must be ≤ $(length(s))."
    @assert wlen >= 2 "wlen must be ≥ 2."
    @assert woverlap <= wlen "woverlap must be ≤ $(wlen)."
    @assert woverlap >= 0 "woverlap must be ≥ 0."

    if method === :mt
        w = w ? hanning(length(s)) : ones(length(s))
        p = mt_pgram(s .* w, fs=fs, nw=((nt + 1) ÷ 2), ntapers=nt)
        pw = power(p)
        f = Vector(freq(p))
        p = pw[1:length(f)]
    elseif method === :stft
        w = w ? DSP.hanning : nothing
        p = abs.(DSP.stft(s, wlen, woverlap, fs=fs, window=w))
        # average STFT segments along time
        p = vec(mean(p, dims=2))
        # create frequencies vector
        f = linspace(0, fs / 2, length(pw))
    elseif method === :welch
        w = w ? DSP.hanning : nothing
        p = DSP.welch_pgram(s, wlen, woverlap, fs=fs, window=w)
        pw = power(p)
        f = Vector(freq(p))
        p = pw[1:length(f)]
    elseif method === :fft
        w = w ? DSP.hanning : nothing
        p = DSP.periodogram(s, fs=fs, window=w)
        pw = power(p)
        f = Vector(freq(p))
        p = pw[1:length(f)]
    elseif method === :mw
        p, f = mwpsd(s, db=false, fs=fs, ncyc=ncyc, w=w)
    end

    # replace powers at extreme frequencies
    # pw[1] = pw[2]
    # pw[end] = pw[end - 1]

    db && (p = pow2db.(p))

    return (p=p, f=f)

end

"""
    psd(s; fs, db, method, nt, wlen, woverlap, w, ncyc)

Calculate power spectrum density. Default method is Welch's periodogram.

# Arguments

- `s::AbstractMatrix`
- `fs::Int64`: sampling rate
- `db::Bool=false`: normalize do dB
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(fs / 2)`

# Returns

Named tuple containing:
- `p::Array{Float64, 2}`: powers
- `f::Vector{Float64}`: frequencies
"""
function psd(s::AbstractMatrix; fs::Int64, db::Bool=false, method::Symbol=:welch, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32)

    ch_n = size(s, 1)
    _, f = psd(s[1, :], fs=fs, db=db, method=method, ncyc=ncyc, nt=nt, wlen=wlen, woverlap=woverlap, w=w)

    p = zeros(ch_n, length(f))

    @inbounds for ch_idx in 1:ch_n
        p[ch_idx, :], _ = psd(s[ch_idx, :], fs=fs, db=db, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc)
    end

    return (p=p, f=f)

end

"""
    psd(s; fs, db, method, nt, wlen, woverlap, w, ncyc)

Calculate power spectrum density. Default method is Welch's periodogram.

# Arguments
- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `db::Bool=false`: normalize do dB
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(fs / 2)`

# Returns

Named tuple containing:
- `p::Array{Float64, 3}`: powers
- `f::Vector{Float64}`: frequencies
"""
function psd(s::AbstractArray; fs::Int64, db::Bool=false, method::Symbol=:welch, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    _, f = psd(s[1, :, 1], fs=fs, db=db, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc)

    p = zeros(ch_n, length(f), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            p[ch_idx, :, ep_idx], _ = psd(s[ch_idx, :, ep_idx], fs=fs, db=db, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc)
        end
    end

    return (p=p, f=f)

end

"""
    psd(obj; ch, db, method, nt, wlen, woverlap, w, ncyc)

Calculate power spectrum density. Default method is Welch's periodogram.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `db::Bool=false`: normalize do dB
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`

# Returns

Named tuple containing:
- `p::Array{Float64, 3}`: powers
- `f::Vector{Float64}`: frequencies
"""
function psd(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), db::Bool=false, method::Symbol=:welch, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32)

    _check_channels(obj, ch)
    isa(ch, Int64) && (ch = [ch])

    p, f = psd(obj.data[ch, :, :], fs=sr(obj), db=db, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc)

    return (p=p, f=f)

end

"""
    mwpsd(s; pad, db, fs, ncyc, w)

Calculate power spectrum using Morlet wavelet convolution.

# Arguments

- `s::AbstractVector`
- `pad::Int64=0`: pad with `pad` zeros
- `db::Bool=true`: normalize powers to dB
- `fs::Int64`: sampling rate
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(fs / 2)`
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:
- `p::Matrix{Float64}`: powers
- `f::Vector{Float64}`: frequencies
"""
function mwpsd(s::AbstractVector; pad::Int64=0, db::Bool=true, fs::Int64, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, w::Bool=true)

    @assert fs >= 1 "fs must be ≥ 1."
    @assert pad >= 0 "pad must be ≥ 0."

    pad > 0 && (s = pad0(s, pad))

    w = w ? hanning(length(s)) : ones(length(s))

    frq_lim = (0, fs / 2)
    frq_n = _tlength(frq_lim)
    f = linspace(frq_lim[1], frq_lim[2], frq_n)

    if ncyc isa Int64
        @assert ncyc >= 1 "ncyc must be ≥ 1"
        ncyc = repeat([ncyc], frq_n)
    else
        @assert ncyc[1] >= 1 "ncyc[1] must be ≥ 1"
        @assert ncyc[2] >= 1 "ncyc[2] must be ≥ 1"
        ncyc = round.(Int64, logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n))
    end

    p = zeros(length(f))
    @inbounds for frq_idx in 1:frq_n
        kernel = generate_morlet(fs, f[frq_idx], 1, ncyc=ncyc[frq_idx], complex=true)
        # w_conv = tconv(s .* w, kernel=kernel)
        w_conv = fconv(s .* w, kernel=kernel, norm=true)
        p[frq_idx] = median((abs.(w_conv)).^2)
    end

    db && (p = pow2db.(p))

    return (p=p, f=f)

end

"""
    mwpsd(s; pad, db, fs, ncyc, w)

Calculate power spectrum using Morlet wavelet convolution.

# Arguments

- `s::AbstractMatrix`
- `pad::Int64=0`: pad with `pad` zeros
- `db::Bool=true`: normalize powers to dB
- `fs::Int64`: sampling rate
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(fs / 2)`
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:
- `p::Array{Float64, 2}`: powers
- `f::Vector{Float64}`: frequencies
"""
function mwpsd(s::AbstractMatrix; pad::Int64=0, db::Bool=true, fs::Int64, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, w::Bool=true)

    ch_n = size(s, 1)

    _, f = mwpsd(s[1, :], pad=pad, db=db, fs=fs, ncyc=ncyc, w=w)
    p = zeros(ch_n, length(f))

    @inbounds for ch_idx in 1:ch_n
        p[ch_idx, :], _ = @views mwpsd(s[ch_idx, :], pad=pad, db=db, fs=fs, ncyc=ncyc, w=w)
    end

    return (p=p, f=f)

end

"""
    mwpsd(s; pad, db, fs, ncyc, w)

Calculate power spectrum using Morlet wavelet convolution.

# Arguments

- `s::AbstractArray`
- `pad::Int64=0`: pad with `pad` zeros
- `db::Bool=true`: normalize powers to dB
- `fs::Int64`: sampling rate
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(fs / 2)`
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:
- `p::Array{Float64, 3}`: powers
- `f::Vector{Float64}`: frequencies
"""
function mwpsd(s::AbstractArray; pad::Int64=0, db::Bool=true, fs::Int64, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, w::Bool=true)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    _, f = mwpsd(s[1, :, 1], pad=pad, db=db, fs=fs, ncyc=ncyc, w=w)
    p = zeros(ch_n, length(f), ep_n)

    # initialize progress bar
    progress_bar && (progbar = Progress(ep_n * ch_n, dt=1, barlen=20, color=:white))

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            p[ch_idx, :, ep_idx], _ = @views mwpsd(s[ch_idx, :, ep_idx], pad=pad, db=db, fs=fs, ncyc=ncyc, w=w)

            # update progress bar
            progress_bar && next!(progbar)
        end
    end

    return (p=p, f=f)

end

"""
    mwpsd(obj; ch, pad, db, ncyc, w)

Calculate power spectrum using Morlet wavelet convolution.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all s channels
- `pad::Int64=0`: pad with `pad` zeros
- `db::Bool`=true: normalize powers to dB
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`
- `w::Bool=true`: if true, apply Hanning window

# Returns

- `p::NeuroAnalyzer.POWERSPECTRUM`: power spectrum
"""
function mwpsd(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), pad::Int64=0, db::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, w::Bool=true)

    _check_channels(obj, ch)
    isa(ch, Int64) && (ch = [ch])

    p, f = @views mwpsd(obj.data[ch, :, :], pad=pad, fs=sr(obj), db=db, ncyc=ncyc, w=w)

    return NeuroAnalyzer.POWERSPECTRUM(p, f)

end
