export band_mpower

"""
    band_mpower(s; fs, f, method, nt, wlen, woverlap, w)

Calculate mean and maximum band power and its frequency.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch periodogram
    - `:fft`: fast-Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
- `nt::Int64=8`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window for Welch and STFT

# Returns

Named tuple containing:
- `mbp::Float64`: mean band power
- `maxfrq::Float64`: frequency of maximum band power
- `maxbp::Float64`: power at maximum band frequency
"""
function band_mpower(s::AbstractVector; fs::Int64, f::Tuple{Real, Real}, method::Symbol=:welch, nt::Int64=8, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)

    @assert fs >= 1 "fs must be ≥ 1."
    f = tuple_order(f)
    @assert f[1] >= 0 "Lower frequency bound must be ≥ 0."
    @assert f[2] <= fs / 2 "Upper frequency bound must be ≤ $(fs / 2)."

    pw, pf = psd(s, fs=fs, norm=false, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w)

    f1_idx = vsearch(f[1], pf)
    f2_idx = vsearch(f[2], pf)
    mbp = mean(pw[f1_idx:f2_idx])
    maxfrq = pf[f1_idx:f2_idx][findmax(pw[f1_idx:f2_idx])[2]]
    maxbp = pw[vsearch(maxfrq, pf)]

    return (mbp=mbp, maxfrq=maxfrq, maxbp=maxbp)
end

"""
    band_mpower(s; fs, f, method, nt, wlen, woverlap, w)

Calculate absolute band power between two frequencies.

# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch periodogram
    - `:fft`: fast-Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
- `nt::Int64=8`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window for Welch and STFT

# Returns

Named tuple containing:
- `mbp::Matrix{Float64}`: mean band power per channel per epoch
- `maxfrq::Matrix{Float64}`: frequency of maximum band power per channel per epoch
- `maxbp::Matrix{Float64}`: power at maximum band frequency per channel per epoch
"""
function band_mpower(s::AbstractArray; fs::Int64, f::Tuple{Real, Real}, method::Symbol=:welch, nt::Int64=8, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)

    ch_n = size(s, 1)
    ep_n = size(s, 3)
    mbp = zeros(ch_n, ep_n)
    maxfrq = zeros(ch_n, ep_n)
    maxbp = zeros(ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            mbp[ch_idx, ep_idx], maxfrq[ch_idx, ep_idx], maxbp[ch_idx, ep_idx] = @views band_mpower(s[ch_idx, :, ep_idx], fs=fs, f=f, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
        end
    end

    return (mbp=mbp, maxfrq=maxfrq, maxbp=maxbp)

end

"""
    band_mpower(obj; ch, f, method, nt, wlen, woverlap, w)

Calculate mean and maximum band power and its frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch periodogram
    - `:fft`: fast-Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
- `nt::Int64=8`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window for Welch and STFT

# Returns

Named tuple containing:
- `mbp::Matrix{Float64}`: mean band power per channel per epoch
- `maxfrq::Matrix{Float64}`: frequency of maximum band power per channel per epoch
- `maxbp::Matrix{Float64}`: power at maximum band frequency per channel per epoch
"""
function band_mpower(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), f::Tuple{Real, Real}, method::Symbol=:welch, nt::Int64=8, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)

    _check_channels(obj, ch)
    
    mbp, maxfrq, maxbp = @views band_mpower(obj.data[ch, :, :], fs=sr(obj), f=f, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w)

    return (mbp=mbp, maxfrq=maxfrq, maxbp=maxbp)

end
