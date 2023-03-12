export psd

"""
    psd(signal; fs, norm, mt)

Calculate power spectrum density.

# Arguments
- `signal::Vector{Float64}`
- `fs::Int64`: sampling rate
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns

Named tuple containing:
- `psd_pow::Vector{Float64}`
- `psd_frq::Vector{Float64}`
"""
function psd(signal::Vector{Float64}; fs::Int64, norm::Bool=false, mt::Bool=false, nt::Int64=8)

    nt < 1 && throw(ArgumentError("nt must be ≥ 1."))
    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    length(signal) < 4 * fs && (mt = true)

    if mt == true
        p = mt_pgram(signal, fs=fs, nw=(nt÷2+1), ntapers=nt)
    else
        p = welch_pgram(signal, 4*fs, fs=fs)
    end

    psd_pow = power(p)
    psd_frq = Vector(freq(p))
    
    # replace power at 0 Hz
    psd_pow[1] = psd_pow[2]
    
    norm == true && (psd_pow = pow2db.(psd_pow))

    return (psd_pow=psd_pow, psd_frq=psd_frq)
end

"""
    psd(signal; fs, norm, mt)

Calculate power spectrum density.

# Arguments
- `signal::Matrix{Float64}`
- `fs::Int64`: sampling rate
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

Named tuple containing:
- `psd_pow::Matrix{Float64}`
- `psd_frq::Matrix{Float64}`
"""
function psd(signal::Matrix{Float64}; fs::Int64, norm::Bool=false, mt::Bool=false)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    size(signal, 2) < 4 * fs && (mt = true)
    ch_n = size(signal, 1)
    psd_tmp, psd_frq = psd(signal[1, :], fs=fs, norm=norm, mt=mt)
    psd_pow = zeros(ch_n, length(psd_tmp))

    @inbounds @simd for ch_idx in 1:ch_n
        psd_pow[ch_idx, :], _ = psd(signal[ch_idx, :], fs=fs, norm=norm, mt=mt)
    end
    
    return (psd_pow=psd_pow, psd_frq=psd_frq)
end

"""
    psd(signal; fs, norm, mt)

Calculate power spectrum density.

# Arguments
- `signal::AbstractArray`
- `fs::Int64`: sampling rate
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

Named tuple containing:
- `psd_pow::Array{Float64, 3}`
- `psd_frq::Array{Float64, 3}`
"""
function psd(signal::AbstractArray; fs::Int64, norm::Bool=false, mt::Bool=false)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    size(signal, 2) < 4 * fs && (mt = true)
    ch_n = size(signal, 1)
    ep_n = size(signal, 3)
    psd_tmp, psd_frq = psd(signal[1, :, 1], fs=fs, norm=norm, mt=mt)
    psd_pow = zeros(ch_n, length(psd_tmp), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            psd_pow[ch_idx, :, ep_idx], _ = psd(signal[ch_idx, :, ep_idx], fs=fs, norm=norm, mt=mt)
        end
    end
    
    return (psd_pow=psd_pow, psd_frq=psd_frq)
end

"""
    psd(obj; ch, norm, mt)

Calculate power spectrum density.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns

Named tuple containing:
- `psd_pow::Array{Float64, 3}`:powers
- `psd_frq::Vector{Float64}`: frequencies
"""
function psd(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), norm::Bool=false, mt::Bool=false, nt::Int64=8)

    fs = sr(obj)

    _check_channels(obj, ch)
    ch_n = length(ch)
    ep_n = epoch_n(obj)
    
    _, psd_frq = psd(obj.data[1, :, 1], fs=fs, norm=norm, mt=mt, nt=nt)
    psd_pow = zeros(length(ch), length(psd_frq), ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            psd_pow[ch_idx, :, ep_idx], _ = psd(obj.data[ch[ch_idx], :, ep_idx], fs=fs, norm=norm, mt=mt, nt=nt)
        end
    end

    return (psd_pow=psd_pow, psd_frq=psd_frq)
end
