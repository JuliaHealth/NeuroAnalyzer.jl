export tkeo

"""
    tkeo(signal)

Calculate Teager-Kaiser energy-tracking operator: y(t) = x(t)^2 - x(t-1) × x(t+1)

# Arguments

- `signal::AbstractVector`

# Returns

- `out::Vector{Float64}`
"""
function tkeo(signal::AbstractVector)
    out = zeros(length(signal))
    out[1] = signal[1]
    out[end] = signal[end]
    @inbounds @simd for idx in 2:(length(signal) - 1)
        out[idx] = signal[idx]^2 - (signal[idx - 1] * signal[idx + 1])
    end

    return out
end

"""
    tkeo(obj; channel)

Calculate Teager-Kaiser energy-tracking operator: y(t) = x(t)^2 - x(t-1) × x(t+1)

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `tkeo::Array{Float64, 3}`
"""
function tkeo(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    out = zeros(ch_n, epoch_len(obj), ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            out[ch_idx, :, ep_idx] = @views tkeo(obj.data[channel[ch_idx], :, ep_idx])
        end
    end

    return out
end
