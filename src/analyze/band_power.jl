export band_power
export band_mpower

"""
    band_power(signal; fs, f, mt)

Calculate absolute band power between two frequencies.

# Arguments

- `signal::AbstractVector`
- `fs::Int64`: sampling rate
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

- `sbp::Float64`: signal band power
"""
function band_power(signal::AbstractVector; fs::Int64, f::Tuple{Real, Real}, mt::Bool=false)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    f = tuple_order(f)
    f[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0.")) 
    f[2] > fs / 2 && throw(ArgumentError("Lower frequency bound must be ≤ $(fs / 2).")) 
    if mt == true
        psd = mt_pgram(signal, fs=fs)
    else
        psd = welch_pgram(signal, 4*fs, fs=fs)
    end

    psd_freq = Vector(psd.freq)
    f1_idx = vsearch(f[1], psd_freq)
    f2_idx = vsearch(f[2], psd_freq)
    frq_idx = [f1_idx, f2_idx]

    # dx: frequency resolution
    dx = psd_freq[2] - psd_freq[1]
    return simpson(psd.power[frq_idx[1]:frq_idx[2]], psd_freq[frq_idx[1]:frq_idx[2]], dx=dx)

end

"""
    band_power(signal; fs, f, mt)

Calculate absolute band power between two frequencies.

# Arguments

- `signal::AbstractArray`
- `fs::Int64`: sampling rate
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

- `sbp::Matrix{Float64}`: band power for each channel per epoch
"""
function band_power(signal::AbstractArray; fs::Int64, f::Tuple{Real, Real}, mt::Bool=false)

    ch_n = size(signal, 1)
    ep_n = size(signal, 3)
    sbp = zeros(channel_n, epoch_n)

    @inbounds @simd for epoch_idx in 1:ep_n
        Threads.@threads for channel_idx in 1:ch_n
            @views sbp[channel_idx, epoch_idx] = band_power(signal[channel_idx, :, epoch_idx], fs=fs, f=f, mt=mt)
        end
    end

    return sbp

end

"""
    band_power(obj; channel, f, mt)

Calculate absolute band power between two frequencies.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

- `sbp::Matrix{Float64}`: band power for each channel per epoch
"""
function band_power(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), f::Tuple{Real, Real}, mt::Bool=false)

    _check_channels(obj, channel)

    return @views band_power(obj.data[channel, :, :], fs=sr(obj), f=f, mt=mt)

end

export band_mpower

"""
    band_mpower(signal; fs, f)

Calculate mean and maximum band power and its frequency.

# Arguments

- `signal::AbstractVector`
- `fs::Int64`: sampling rate
- `f::Tuple{Real, Real}`: lower and upper frequency bounds

# Returns

Named tuple containing:
- `mbp::Float64`: mean band power [dB]
- `maxfrq::Float64`: frequency of maximum band power [Hz]
- `maxbp::Float64`: power at maximum band frequency [dB]
"""
function band_mpower(signal::AbstractVector; fs::Int64, f::Tuple{Real, Real}, mt::Bool=false)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    f = tuple_order(f)
    f[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0.")) 
    f[2] > fs / 2 && throw(ArgumentError("Lower frequency bound must be ≤ $(fs / 2).")) 

    if mt == true
        p = mt_pgram(signal, fs=fs)
    else
        p = welch_pgram(signal, 4*fs, fs=fs)
    end

    psd_freq = Vector(p.freq)
    f1_idx = vsearch(f[1], psd_freq)
    f2_idx = vsearch(f[2], psd_freq)
    mbp = mean(p.power[f1_idx:f2_idx])
    maxfrq = psd_freq[f1_idx:f2_idx][findmax(p.power[f1_idx:f2_idx])[2]]
    maxbp = psd.power[vsearch(maxfrq, psd_freq)]

    return (mbp=mbp, maxfrq=maxfrq, maxbp=maxbp)
end

"""
    band_mpower(obj; channel, f, mt)

Calculate mean and maximum band power and its frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

Named tuple containing:
- `mbp::Matrix{Float64}`: mean band power [μV^2/Hz] per channel per epoch
- `maxfrq::Matrix{Float64}`: frequency of maximum band power [Hz] per channel per epoch
- `maxbp::Matrix{Float64}`: power at maximum band frequency [μV^2/Hz] per channel per epoch
"""
function band_mpower(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), f::Tuple{Real, Real}, mt::Bool=false)

    fs = sr(obj)
    f = tuple_order(f)
    f[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
    f[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be < $(fs / 2)."))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    mbp = zeros(ch_n, ep_n)
    maxfrq = zeros(ch_n, ep_n)
    maxbp = zeros(ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            mbp[ch_idx, ep_idx], maxfrq[ch_idx, ep_idx], maxbp[ch_idx, ep_idx] = @views band_mpower(obj.data[channel[ch_idx], :, ep_idx], fs=fs, f=f, mt=mt)
        end
    end

    return (mbp=mbp, maxfrq=maxfrq, maxbp=maxbp)
end
