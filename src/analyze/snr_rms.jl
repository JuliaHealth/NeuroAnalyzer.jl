export snr
export snr2
export rms
export rmse

"""
    snr(signal)

Calculate SNR.

# Arguments

- `signal::AbstractVector`

# Returns

- `snr::Float64`: SNR

# Source

D. J. Schroeder (1999). Astronomical optics (2nd ed.). Academic Press. ISBN 978-0-12-629810-9, p.278
"""
function snr(signal::AbstractVector)
    return mean(signal) / std(signal)
end

"""
    snr2(signal)

Calculate RMS-based SNR.

# Arguments

- `signal::AbstractVector`

# Returns

- `snr::Float64`: SNR
"""
function snr2(signal::AbstractVector)
    return (maximum(signal) - minimum(signal)) / rms(signal)
end

"""
    rms(signal)

Calculate Root Mean Square.

# Arguments

- `signal::AbstractVector`

# Returns

- rms::Float64`
"""
function rms(signal::AbstractVector)
    # rms = sqrt(mean(signal.^2))    
    return norm(signal) / sqrt(length(signal))
end

"""
    rmse(signal1, signal2)

Calculate RMSE between two signals.

# Arguments

- `signal1::AbstractVector`
- `signal2::AbstractVector`

# Returns

- `r::Float64`
"""
function rmse(signal1::AbstractVector, signal2::AbstractVector)
    length(signal1) == length(signal2) || throw(ArgumentError("Both signals must have the same length."))
    # r = sum(signal1 .* signal2) ./ (sqrt(sum(signal1.^2)) .* sqrt(sum(signal2.^2)))
    return sqrt(mean(signal2 - signal1)^2)
end

"""
    snr(obj; channel, type)

Calculate SNR.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `type::Symbol=:rms`: SNR type:
    - `:mean`: mean-based
    - `:rms`: RMS-based

# Returns

Named tuple containing:
- `snr::Matrix(Float64)`: SNR for each channel over frequencies 1:Nyquist
- `hz::Vector(Float64)`: frequencies
"""
function snr(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), type::Symbol=:rms)

    _check_var(type, [:mean, :rms], "type")
    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    ep_n == 1 && throw(ArgumentError("OBJ must contain ≥ 2 epochs."))

    hz, _ = s_freqs(obj.epoch_time)
    amp = zeros(ch_n, length(hz), ep_n)
    snr = zeros(ch_n, length(hz))

    # create spectrum for each channel
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            _, amp[ch_idx, :, ep_idx], _, _ = @views spectrum(obj.data[channel[ch_idx], :, ep_idx])
        end
    end

    # calculate SNR for each channel spectrum
    @inbounds @simd for hz_idx in 1:length(hz)
        Threads.@threads for ch_idx in 1:ch_n
            if type === :mean
                snr[ch_idx, hz_idx] = @views snr(amp[ch_idx, hz_idx, :])
            else
                snr[ch_idx, hz_idx] = @views snr2(amp[ch_idx, hz_idx, :])
            end
        end
    end

    return (snr=snr, hz=hz)
end