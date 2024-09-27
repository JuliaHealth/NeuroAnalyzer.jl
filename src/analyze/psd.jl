export psd
export mwpsd
export ghpsd
export cwtpsd

"""
    psd(s; <keyword arguments>)

Calculate power spectrum density. Default method is Welch's periodogram.

# Arguments
- `s::Vector{Float64}`
- `fs::Int64`: sampling rate
- `db::Bool=false`: normalize do dB; for CWT power spectrum: normalize to the signal scale so the amplitudes of wavelet coefficients agree with the amplitudes of oscillatory components in a signal
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(fs / 2)`
- `gw::Real=5`: Gaussian width in Hz
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

Named tuple containing:
- `p::Vector{Float64}`: powers
- `f::Vector{Float64}`: frequencies
"""
function psd(s::AbstractVector; fs::Int64, db::Bool=false, method::Symbol=:welch, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, gw::Real=5, wt::T=wavelet(Morlet(2π), β=32, Q=128))::NamedTuple{(:p, :f), Tuple{Vector{Float64}, Vector{Float64}}} where {T <: CWT}

    _check_var(method, [:fft, :welch, :mt, :mw, :stft, :gh, :cwt], "method")
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
        db && (p = pow2db.(p))
    elseif method === :stft
        w = w ? DSP.hanning : nothing
        p = abs.(DSP.stft(s, wlen, woverlap, fs=fs, window=w))
        # average STFT segments along time
        p = vec(mean(p, dims=2))
        # create frequencies vector
        f = linspace(0, fs / 2, length(p))
        db && (p = pow2db.(p))
    elseif method === :welch
        w = w ? DSP.hanning : nothing
        p = DSP.welch_pgram(s, wlen, woverlap, fs=fs, window=w)
        pw = power(p)
        f = Vector(freq(p))
        p = pw[1:length(f)]
        db && (p = pow2db.(p))
    elseif method === :fft
        w = w ? DSP.hanning : nothing
        p = DSP.periodogram(s, fs=fs, window=w)
        pw = power(p)
        f = Vector(freq(p))
        p = pw[1:length(f)]
        db && (p = pow2db.(p))
    elseif method === :mw
        p, f = mwpsd(s, db=db, fs=fs, ncyc=ncyc, w=w)
    elseif method === :gh
        p, f = ghpsd(s, fs=fs, db=db, gw=gw, w=w)
    elseif method === :cwt
        p, f = cwtpsd(s, fs=fs, norm=db, wt=wt)
    end

    # replace powers at extreme frequencies
    # pw[1] = pw[2]
    # pw[end] = pw[end - 1]

    return (p=p, f=f)

end

"""
    psd(s; <keyword arguments>)

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
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(fs / 2)`
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

Named tuple containing:
- `p::Matrix{Float64}`: powers
- `f::Vector{Float64}`: frequencies
"""
function psd(s::AbstractMatrix; fs::Int64, db::Bool=false, method::Symbol=:welch, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, gw::Real=5, wt::T=wavelet(Morlet(2π), β=32, Q=128))::NamedTuple{(:p, :f), Tuple{Matrix{Float64}, Vector{Float64}}} where {T <: CWT}

    ch_n = size(s, 1)
    _, f = psd(s[1, :], fs=fs, db=db, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw, wt=wt)

    p = zeros(ch_n, length(f))

    @inbounds for ch_idx in 1:ch_n
        p[ch_idx, :], _ = psd(s[ch_idx, :], fs=fs, db=db, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw, wt=wt)
    end

    return (p=p, f=f)

end

"""
    psd(s; <keyword arguments>)

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
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(fs / 2)`
- `gw::Real=5`: Gaussian width in Hz
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

Named tuple containing:
- `p::Array{Float64, 3}`: powers
- `f::Vector{Float64}`: frequencies
"""
function psd(s::AbstractArray; fs::Int64, db::Bool=false, method::Symbol=:welch, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, gw::Real=5, wt::T=wavelet(Morlet(2π), β=32, Q=128))::NamedTuple{(:p, :f), Tuple{Array{Float64, 3}, Vector{Float64}}} where {T <: CWT}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    _, f = psd(s[1, :, 1], fs=fs, db=db, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw, wt=wt)

    p = zeros(ch_n, length(f), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            p[ch_idx, :, ep_idx], _ = psd(s[ch_idx, :, ep_idx], fs=fs, db=db, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw, wt=wt)
        end
    end

    return (p=p, f=f)

end

"""
    psd(obj; <keyword arguments>)

Calculate power spectrum density. Default method is Welch's periodogram.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channel name or list of channel names
- `db::Bool=false`: normalize do dB
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`
- `gw::Real=5`: Gaussian width in Hz
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

Named tuple containing:
- `p::Array{Float64, 3}`: powers
- `f::Vector{Float64}`: frequencies
"""
function psd(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, db::Bool=false, method::Symbol=:welch, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, gw::Real=5, wt::T=wavelet(Morlet(2π), β=32, Q=128))::NamedTuple{(:p, :f), Tuple{Array{Float64, 3}, Vector{Float64}}} where {T <: CWT}

    ch = get_channel(obj, ch=ch)
    p, f = psd(obj.data[ch, :, :], fs=sr(obj), db=db, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw, wt=wt)

    return (p=p, f=f)

end

"""
    mwpsd(s; <keyword arguments>)

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
- `p::Vector{Float64}`: powers
- `f::Vector{Float64}`: frequencies
"""
function mwpsd(s::AbstractVector; pad::Int64=0, db::Bool=true, fs::Int64, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, w::Bool=true)::NamedTuple{(:p, :f), Tuple{Vector{Float64}, Vector{Float64}}}

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
        w_conv = fconv(s .* w, kernel=kernel, norm=false)
        p[frq_idx] = median((abs.(w_conv)).^2)
    end

    db && (p = pow2db.(p))

    return (p=p, f=f)

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

# Returns

Named tuple containing:
- `p::Vector{Float64}`: powers
- `f::Vector{Float64}`: frequencies
"""
function ghpsd(s::AbstractVector; fs::Int64, db::Bool=true, gw::Real=5, w::Bool=true)::NamedTuple{(:p, :f), Tuple{Vector{Float64}, Vector{Float64}}}

    @assert fs >= 1 "fs must be ≥ 1."

    frq_lim = (0, fs / 2)
    frq_n = _tlength(frq_lim)
    f = linspace(frq_lim[1], frq_lim[2], frq_n)

    w = w ? hanning(length(s)) : ones(length(s))

    p = zeros(length(f), length(s))
    ph = zeros(length(f), length(s))

    @inbounds for frq_idx in eachindex(f)
        s_tmp = filter_g(s .* w, fs=fs, f=f[frq_idx], gw=gw)
        p[frq_idx, :] = (abs.(hilbert(s_tmp))).^2
        ph[frq_idx, :] = angle.(hilbert(s_tmp))
    end

    p[p .== -Inf] .= minimum(p[p .!== -Inf])
    p[p .== +Inf] .= maximum(p[p .!== +Inf])

    p = vec(mean(p, dims=2))

    db && (p = pow2db.(p))

    return (p=p, f=f)

end

"""
    cwtpsd(s; <keyword arguments>)

Calculate power spectrum using continuous wavelet transformation (CWT).

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `norm::Bool=true`: normalize scaleogram to the signal scale so the amplitudes of wavelet coefficients agree with the amplitudes of oscillatory components in a signal

# Returns

Named tuple containing:
- `p::Vector{Float64}`: powers
- `f::Vector{Float64}`: frequencies
"""
function cwtpsd(s::AbstractVector; fs::Int64, wt::T=wavelet(Morlet(2π), β=32, Q=128), norm::Bool=true)::NamedTuple{(:p, :f), Tuple{Vector{Float64}, Vector{Float64}}} where {T <: CWT}

    @assert fs >= 1 "fs must be ≥ 1."

    p = abs.(ContinuousWavelets.cwt(s, wt)')
    p = vec(mean(p, dims=2))

    # scale
    if norm
        a = amp(s)[1]
        p = normalize_n(p, a)
    end

    f = cwtfrq(s, fs=fs, wt=wt)

    return (p=p, f=f)

end
