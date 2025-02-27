export peak_frq

"""
    peak_frq(s; <keyword arguments>)

Calculate peak frequency in a band.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
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

- `pf::Float64`: peak frequency
"""
function peak_frq(s::AbstractVector; fs::Int64, f::Tuple{Real, Real}, method::Symbol=:welch, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32)::Float64

    @assert fs >= 1 "fs must be ≥ 1."
    _check_tuple(f, "f", (0, fs / 2))

    pw, pf = psd(s, fs=fs, db=false, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc)

    f1_idx = vsearch(f[1], pf)
    f2_idx = vsearch(f[2], pf)

    # peak power
    pp = maximum(pw[f1_idx:f2_idx])
    pf = pf[vsearch(pp, pw)]

    return pf

end

"""
    peak_frq(s; <keyword arguments>)

Calculate peak frequency in a band.

# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
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

- `pf::Matrix{Float64}`: peak frequency
"""
function peak_frq(s::AbstractArray; fs::Int64, f::Tuple{Real, Real}, method::Symbol=:welch, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32)::Matrix{Float64}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)
    pf = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            pf[ch_idx, ep_idx] = @views peak_frq(s[ch_idx, :, ep_idx], fs=fs, f=f, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc)
        end
    end

    return pf

end

"""
    peak_frq(obj; <keyword arguments>)

Calculate peak frequency in a band.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
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

- `pf::Matrix{Float64}`: peak frequency
"""
function peak_frq(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, f::Tuple{Real, Real}, method::Symbol=:welch, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32)::Matrix{Float64}

    ch = get_channel(obj, ch=ch)
    pf = @views peak_frq(obj.data[ch, :, :], fs=sr(obj), f=f, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc)

    return pf

end
