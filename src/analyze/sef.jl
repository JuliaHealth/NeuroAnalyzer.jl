export sef

"""
    sef(s; <keyword arguments>)

Calculate spectral edge frequency (SEF) - the frequency below which x percent of the total power of a given signal are located; typically, x is in the range 75 to 95.

# Arguments

- `s::AbstractVector`
- `x::Float64=0.95`: threshold
- `fs::Int64`: sampling rate
- `f::Tuple{Real, Real}=(0, fs / 2)`: lower and upper frequency bounds, default is total power
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

- `sef_frq::Float64`: spectral edge frequency
"""
function sef(s::AbstractVector; x::Float64=0.95, fs::Int64, f::Tuple{Real, Real}=(0, fs / 2), method::Symbol=:welch, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32)::Float64

    @assert fs >= 1 "fs must be ≥ 1."
    _check_tuple(f, "f", (0, fs / 2))

    pw, pf = psd(s, fs=fs, db=false, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc)

    f1_idx = vsearch(f[1], pf)
    f2_idx = vsearch(f[2], pf)

    # dx: frequency resolution
    dx = pf[2] - pf[1]

    pw = pw[f1_idx:f2_idx]
    pf = pf[f1_idx:f2_idx]

    tp = simpson(pw, dx=dx)
    tp_threshold = tp * x

    sef_frq = nothing
    for idx in eachindex(pf)
        if sum(pw[1:idx]) >= tp_threshold
            sef_frq = pf[idx]
            break
        end
    end

    return sef_frq

end

"""
    sef(s; <keyword arguments>)

Calculate spectral edge frequency (SEF) - the frequency below which x percent of the total power of a given signal are located; typically, x is in the range 75 to 95.

# Arguments

- `s::AbstractArray`
- `x::Float64=0.95`: threshold
- `fs::Int64`: sampling rate
- `f::Tuple{Real, Real}=(0, fs / 2)`: lower and upper frequency bounds, default is total power
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

- `sef_frq::Matrix{Float64}`: spectral edge frequency
"""
function sef(s::AbstractArray; x::Float64=0.95, fs::Int64, f::Tuple{Real, Real}=(0, fs / 2), method::Symbol=:welch, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32)::Matrix{Float64}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)
    sef_frq = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            sef_frq[ch_idx, ep_idx] = @views sef(s[ch_idx, :, ep_idx], x=x, fs=fs, f=f, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc)
        end
    end

    return sef_frq

end

"""
    sef(obj; <keyword arguments>)

Calculate spectral edge frequency (SEF) - the frequency below which x percent of the total power of a given signal are located; typically, x is in the range 75 to 95.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `x::Float64=0.95`: threshold
- `f::Tuple{Real, Real}=(0, sr(obj) / 2)`: lower and upper frequency bounds, default is total power
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`

# Returns

- `sef_frq::Matrix{Float64}`: spectral edge frequency
"""
function sef(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, x::Float64=0.95, f::Tuple{Real, Real}=(0, sr(obj) / 2), method::Symbol=:welch, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32)::Matrix{Float64}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    sef_frq = @views sef(obj.data[ch, :, :], x=x, fs=sr(obj), f=f, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc)

    return sef_frq

end
