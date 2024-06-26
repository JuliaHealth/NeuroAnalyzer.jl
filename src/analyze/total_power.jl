export total_power

"""
    total_power(s; fs, method, nt, wlen, woverlap, w, ncyc)

Calculate total power.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
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

- `tp::Float64`: total power
"""
function total_power(s::AbstractVector; fs::Int64, method::Symbol=:welch, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32)

    pw, pf = psd(s, fs=fs, db=false, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc)

    # dx: frequency resolution
    dx = pf[2] - pf[1]
    tp = simpson(pw, dx=dx)

    return tp

end

"""
    total_power(s; fs, method, nt, wlen, woverlap, w, ncyc)

Calculate total power.

`# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate
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

- `tp::Matrix{Float64}`: total power

"""
function total_power(s::AbstractArray; fs::Int64, method::Symbol=:welch, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    tp = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            tp[ch_idx, ep_idx] = @views total_power(s[ch_idx, :, ep_idx], fs=fs, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc)
        end
    end

    return tp
end

"""
    total_power(obj, ch, method, nt, wlen, woverlap, w, ncyc)

Calculate total power.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(record)`: index of channels, default is all signal channels
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

- `tp::Matrix{Float64}`: total power
"""
function total_power(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), method::Symbol=:welch, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32)

    _check_channels(obj, ch)
    length(ch) == 1 && (ch = [ch])

    tp = @views total_power(obj.data[ch, :, :], fs=sr(obj), method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc)

    return tp

end
