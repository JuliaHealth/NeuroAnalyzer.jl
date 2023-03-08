export band_power

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

    return 
end

"""
    band_power(obj; channel, f, mt)

Calculate absolute band power between two frequencies.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

- `sbp::Matrix{Float64}`: band power for each channel per epoch
"""
function band_power(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), f::Tuple{Real, Real}, mt::Bool=false)

    _check_channels(obj, channel)

    return @views band_power(obj.data[channel, :, :], fs=sr(obj), f=f, mt=mt)

end