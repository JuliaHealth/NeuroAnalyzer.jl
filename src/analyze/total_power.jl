export total_power

"""
    total_power(signal; fs, mt)

Calculate total power.

# Arguments

- `signal::AbstractVector`
- `fs::Int64`: sampling rate
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

- `stp::Float64`: signal total power
"""
function total_power(signal::AbstractVector; fs::Int64, mt::Bool=false)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    if mt == true
        psd = mt_pgram(signal, fs=fs)
    else
        psd = welch_pgram(signal, 4*fs, fs=fs)
    end
    psd_pow = power(psd)
    psd_pow[1] = psd_pow[2]
    # dx: frequency resolution
    dx = psd.freq[2] - psd.freq[1]
    stp = simpson(psd_pow, dx=dx)

    return stp
end

"""
    total_power(signal; fs, mt)

Calculate total power.

`# Arguments

- `signal::AbstractArray`
- `fs::Int64`: sampling rate
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

- `stp::Float64`: signal total power
"""
function total_power(signal::AbstractArray; fs::Int64, mt::Bool=false)

    ch_n = size(signal, 1)
    ep_n = size(signal, 3)

    stp = zeros(ch_n, ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            @views stp[ch_idx, ep_idx] = total_power(signal[ch_idx, :, ep_idx], fs=fs, mt=mt)
        end
    end

    return stp
end

"""
    total_power(obj, channel, mt)

Calculate total power.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(record)`: index of channels, default is all signal channels
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns
 
- `stp::Matrix{Float64}`: total power for each channel per epoch
"""
function total_power(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), mt::Bool=false)

    fs = sr(obj)
    _check_channels(obj, channel)

    return @views total_power(obj.data[channel, :, :], fs=fs, mt=mt)
end
