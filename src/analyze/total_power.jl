export total_power

"""
    total_power(s; fs, mt, st, nt, wlen, woverlap, w)

Calculate total power.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `mt::Bool=false`: if true, use multi-tapered periodogram
- `st::Bool=false`: if true, use short time Fourier transform
- `nt::Int64=8`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window for Welch and STFT

# Returns

- `tp::Float64`: total power
"""
function total_power(s::AbstractVector; fs::Int64, mt::Bool=false, st::Bool=false, nt::Int64=8, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)

    pw, pf = psd(s, fs=fs, norm=false, mt=mt, st=st, nt=nt, wlen=wlen, woverlap=woverlap, w=w)

    # dx: frequency resolution
    dx = pf[2] - pf[1]
    tp = simpson(pw, dx=dx)

    return tp

end

"""
    total_power(s; fs, mt, st, nt, wlen, woverlap, w)

Calculate total power.

`# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `mt::Bool=false`: if true, use multi-tapered periodogram
- `st::Bool=false`: if true, use short time Fourier transform
- `nt::Int64=8`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window for Welch and STFT

# Returns

- `tp::Matrix{Float64}`: total power

"""
function total_power(s::AbstractArray; fs::Int64, mt::Bool=false, st::Bool=false, nt::Int64=8, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    tp = zeros(ch_n, ep_n)
    
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            @views tp[ch_idx, ep_idx] = total_power(s[ch_idx, :, ep_idx], fs=fs, mt=mt, st=st, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
        end
    end

    return tp
end

"""
    total_power(obj, ch, mt, st, nt, wlen, woverlap)

Calculate total power.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(record)`: index of channels, default is all signal channels
- `mt::Bool=false`: if true, use multi-tapered periodogram
- `st::Bool=false`: if true, use short time Fourier transform
- `nt::Int64=8`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window for Welch and STFT

# Returns
 
- `tp::Matrix{Float64}`: total power
"""
function total_power(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), mt::Bool=false, st::Bool=false, nt::Int64=8, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)

    _check_channels(obj, ch)

    tp = @views total_power(obj.data[ch, :, :], fs=sr(obj), mt=mt, st=st, nt=nt, wlen=wlen, woverlap=woverlap, w=w)

    return tp

end
