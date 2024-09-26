export band_mpower

"""
    band_mpower(s; <keyword arguments>)

Calculate mean and maximum band power and its frequency.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `frq_lim::Tuple{Real, Real}`: lower and upper frequency bounds
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
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(fs / 2)`
- `gw::Real=5`: Gaussian width in Hz
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

Named tuple containing:
- `mbp::Float64`: mean band power
- `maxfrq::Float64`: frequency of maximum band power
- `maxbp::Float64`: power at maximum band frequency
"""
function band_mpower(s::AbstractVector; fs::Int64, frq_lim::Tuple{Real, Real}, method::Symbol=:welch, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, gw::Real=5, wt::T=wavelet(Morlet(2π), β=32, Q=128))::NamedTuple{(:mbp, :maxfrq, :maxbp), Tuple{Float64, Float64, Float64}} where {T <: CWT}

    @assert fs >= 1 "fs must be ≥ 1."
    _check_tuple(frq_lim, "frq_lim", (0, fs / 2))

    pw, pf = psd(s, fs=fs, db=false, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw, wt=wt)

    f1_idx = vsearch(frq_lim[1], pf)
    f2_idx = vsearch(frq_lim[2], pf)
    mbp = mean(pw[f1_idx:f2_idx])
    maxfrq = pf[f1_idx:f2_idx][findmax(pw[f1_idx:f2_idx])[2]]
    maxbp = pw[vsearch(maxfrq, pf)]

    return (mbp=mbp, maxfrq=maxfrq, maxbp=maxbp)
end

"""
    band_mpower(s; <keyword arguments>)

Calculate absolute band power between two frequencies.

# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `frq_lim::Tuple{Real, Real}`: lower and upper frequency bounds
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
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(fs / 2)`
- `gw::Real=5`: Gaussian width in Hz
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

Named tuple containing:
- `mbp::Matrix{Float64}`: mean band power per channel per epoch
- `maxfrq::Matrix{Float64}`: frequency of maximum band power per channel per epoch
- `maxbp::Matrix{Float64}`: power at maximum band frequency per channel per epoch
"""
function band_mpower(s::AbstractArray; fs::Int64, frq_lim::Tuple{Real, Real}, method::Symbol=:welch, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, gw::Real=5, wt::T=wavelet(Morlet(2π), β=32, Q=128))::NamedTuple{(:mbp, :maxfrq, :maxbp), Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}} where {T <: CWT}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)
    mbp = zeros(ch_n, ep_n)
    maxfrq = zeros(ch_n, ep_n)
    maxbp = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            mbp[ch_idx, ep_idx], maxfrq[ch_idx, ep_idx], maxbp[ch_idx, ep_idx] = @views band_mpower(s[ch_idx, :, ep_idx], fs=fs, frq_lim=frq_lim, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw, wt=wt)
        end
    end

    return (mbp=mbp, maxfrq=maxfrq, maxbp=maxbp)

end

"""
    band_mpower(obj; <keyword arguments>)

Calculate mean and maximum band power and its frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channel name or list of channel names
- `frq_lim::Tuple{Real, Real}`: lower and upper frequency bounds
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
- `gw::Real=5`: Gaussian width in Hz
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

Named tuple containing:
- `mbp::Matrix{Float64}`: mean band power per channel per epoch
- `maxfrq::Matrix{Float64}`: frequency of maximum band power per channel per epoch
- `maxbp::Matrix{Float64}`: power at maximum band frequency per channel per epoch
"""
function band_mpower(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, frq_lim::Tuple{Real, Real}, method::Symbol=:welch, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, gw::Real=5, wt::T=wavelet(Morlet(2π), β=32, Q=128))::NamedTuple{(:mbp, :maxfrq, :maxbp), Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}} where {T <: CWT}

    ch = get_channel(obj, ch=ch)

    mbp, maxfrq, maxbp = @views band_mpower(obj.data[ch, :, :], fs=sr(obj), frq_lim=frq_lim, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw, wt=wt)

    return (mbp=mbp, maxfrq=maxfrq, maxbp=maxbp)

end
