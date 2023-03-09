export psd_rel

"""
    psd_rel(signal; fs, norm, mt, f)

Calculate relative power spectrum density.

# Arguments
- `signal::AbstractVector`
- `fs::Int64`: sampling rate
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `f::Union(Tuple{Real, Real}, Nothing)=nothing`: calculate power relative to frequency range or total power

# Returns

Named tuple containing:
- `psd_pow::Vector{Float64}`
- `psd_frq::Vector{Float64}`
"""
function psd_rel(signal::AbstractVector; fs::Int64, norm::Bool=false, mt::Bool=false, f::Union{Tuple{Real, Real}, Nothing}=nothing)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    if f !== nothing
        f = tuple_order(f)
        f[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0.")) 
        f[2] > fs / 2 && throw(ArgumentError("Lower frequency bound must be ≤ $(fs / 2).")) 
    end
    ref_pow = f === nothing ? total_power(signal, fs=fs, mt=mt) : band_power(signal, fs=fs, mt=mt, f=f)
    if mt == true
        p = mt_pgram(signal, fs=fs)
    else
        p = welch_pgram(signal, 4*fs, fs=fs)
    end

    psd_pow = power(p)
    psd_frq = Vector(freq(p))
    psd_pow[1] = psd_pow[2]
    psd_pow = psd_pow / ref_pow

    norm == true && (psd_pow = pow2db.(psd_pow))

    return (psd_pow=psd_pow, psd_frq=psd_frq)
end

"""
    psd_rel(signal; fs, norm, mt, f)

Calculate relative power spectrum density.

# Arguments
- `signal::Matrix{Float64}`
- `fs::Int64`: sampling rate
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `f::Union(Tuple{Real, Real}, Nothing)=nothing`: calculate power relative to frequency range or total power

# Returns

Named tuple containing:
- `psd_pow::Vector{Float64}`
- `psd_frq::Vector{Float64}`
"""
function psd_rel(signal::Matrix{Float64}; fs::Int64, norm::Bool=false, mt::Bool=false, f::Union{Tuple{Real, Real}, Nothing}=nothing)

    ch_n = size(signal, 1)
    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    f = tuple_order(f)
    f[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0.")) 
    f[2] > fs / 2 && throw(ArgumentError("Lower frequency bound must be ≤ $(fs / 2).")) 

    if mt == true
        psd_tmp = mt_pgram(signal[1, :], fs=fs)
    else
        psd_tmp = welch_pgram(signal[1, :], 4*fs, fs=fs)
    end
    psd_frq = Vector(freq(psd_tmp))
    psd_pow = zeros(ch_n, length(Vector(freq(psd_tmp))))

    Threads.@threads for channel_idx in 1:ch_n
        ref_pow = f === nothing ? total_power(signal[channel_idx, :], fs=fs, mt=mt) : band_power(signal[channel_idx, :], fs=fs, mt=mt, f=f)
        if mt == true
            p = mt_pgram(signal[channel_idx, :], fs=fs)
        else
            p = welch_pgram(signal[channel_idx, :], 4*fs, fs=fs)
        end
        psd_pow[channel_idx, :] = power(p)
        psd_pow[channel_idx, :] = psd_pow[channel_idx, :] / ref_pow
        psd_pow[channel_idx, 1] = psd_pow[channel_idx, 2]
    end

    norm == true && (psd_pow = pow2db.(psd_pow))

    return (psd_pow=psd_pow, psd_frq=psd_frq)
end

"""
    psd_rel(obj; channel, norm, mt, f)

Calculate relative power spectrum density.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `f::Union{Tuple{Real, Real}, Nothing}=nothing`: calculate power relative to frequency range or total power

# Returns

Named tuple containing:
- `psd_pow::Array{Float64, 3}`:powers
- `psd_frq::Array{Float64, 3}`: frequencies
"""
function psd_rel(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), norm::Bool=false, mt::Bool=false, f::Union{Tuple{Real, Real}, Nothing}=nothing)

    fs = sr(obj)
    if f !== nothing
        f = tuple_order(f)
        f[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
        f[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be < $(fs / 2)."))
    end

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    psd_tmp, psd_frq = @views psd_rel(obj.data[1, :, 1], fs=fs, norm=norm, mt=mt, f=f)
    psd_pow = zeros(ch_n, length(psd_tmp), ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        for ch_idx in 1:ch_n
            psd_pow[ch_idx, :, ep_idx], _ = @views psd_rel(obj.data[channel[ch_idx], :, ep_idx], fs=fs, norm=norm, mt=mt, f=f)
        end
    end

    return (psd_pow=psd_pow, psd_frq=psd_frq)
end
