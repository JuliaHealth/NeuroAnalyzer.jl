export band_power

"""
    band_power(s; fs, f, mt, st, nt, wlen, woverlap, w)

Calculate absolute band power between two frequencies.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
- `mt::Bool=false`: if true, use multi-tapered periodogram
- `st::Bool=false`: if true, use short time Fourier transform
- `nt::Int64=8`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window for Welch and STFT

# Returns

- `bp::Float64`: band power
"""
function band_power(s::AbstractVector; fs::Int64, f::Tuple{Real, Real}, mt::Bool=false, st::Bool=true, nt::Int64=8, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)

    @assert fs >= 1 "fs must be ≥ 1."
    f = tuple_order(f)
    @assert f[1] >= 0 "Lower frequency bound must be ≥ 0."
    @assert f[2] <= fs / 2 "Upper frequency bound must be ≤ $(fs / 2)."

    pw, pf = psd(s, fs=fs, norm=false, mt=mt, st=st, nt=nt, wlen=wlen, woverlap=woverlap, w=w)

    f1_idx = vsearch(f[1], pf)
    f2_idx = vsearch(f[2], pf)
    frq_idx = [f1_idx, f2_idx]

    # dx: frequency resolution
    dx = pf[2] - pf[1]

    # integrate
    bp = simpson(pw[frq_idx[1]:frq_idx[2]], pf[frq_idx[1]:frq_idx[2]], dx=dx)

    return bp

end

"""
    band_power(s; fs, f, mt, st, nt, wlen, woverlap, w)

Calculate absolute band power between two frequencies.

# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
- `mt::Bool=false`: if true, use multi-tapered periodogram
- `st::Bool=false`: if true, use short time Fourier transform
- `nt::Int64=8`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window for Welch and STFT

# Returns

- `bp::Matrix{Float64}`: band power
"""
function band_power(s::AbstractArray; fs::Int64, f::Tuple{Real, Real}, mt::Bool=false, st::Bool=false, nt::Int64=8, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)

    ch_n = size(s, 1)
    ep_n = size(s, 3)
    bp = zeros(ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            bp[ch_idx, ep_idx] = @views band_power(s[ch_idx, :, ep_idx], fs=fs, f=f, mt=mt, st=st, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
        end
    end

    return bp

end

"""
    band_power(obj; ch, f, mt, st, nt, wlen, woverlap, w)

Calculate absolute band power between two frequencies.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
- `mt::Bool=false`: if true, use multi-tapered periodogram
- `st::Bool=false`: if true, use short time Fourier transform
- `nt::Int64=8`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window for Welch and STFT

# Returns

- `bp::Matrix{Float64}`: band power
"""
function band_power(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), f::Tuple{Real, Real}, mt::Bool=false, st::Bool=false, nt::Int64=8, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)

    _check_channels(obj, ch)

    bp = @views band_power(obj.data[ch, :, :], fs=sr(obj), f=f, mt=mt, st=st, nt=nt, wlen=wlen, woverlap=woverlap, w=w)

    return bp

end
