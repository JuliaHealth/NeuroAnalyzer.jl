export add_noise
export add_noise!

"""
    add_noise(signal, noise)

Adds noise to signal.

# Arguments

- `signal::AbstractVector`
- `noise::AbstractVector`
1
# Returns

- `s_noisy::AbstractVector`
"""
function add_noise(signal::AbstractVector, noise::AbstractVector)
    length(signal) == length(noise) || throw(ArgumentError("Length of signal and noise must be equal."))
    return signal .+ noise
end

"""
    eeg_add_noise(eeg; channel, noise)

Add random noise to EEG channel(s)

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `noise::AbstractVector`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eeg_add_noise(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), noise::AbstractVector)

    _check_channels(eeg, channel)

    ch_n = length(channel)
    ep_n = eeg_epoch_n(eeg)

    obj_new = deepcopy(eeg)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in eachindex(channel)
            @views obj_new.data[channel[ch_idx], :, ep_idx] = add_noise(obj_new.data[channel[ch_idx], :, ep_idx], noise=noise)
        end
    end

    eeg_reset_components!(obj_new)
    push!(obj_new.header.history, "eeg_add_noise(EEG, channel=$channel)")

    return obj_new
end

"""
    eeg_add_noise!(eeg; channel, noise)

Add random noise to EEG channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg))`: index of channels, default is all channels
- `noise::AbstractVector`
"""
function eeg_add_noise!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(eeg_channel_n(eeg)), noise::AbstractVector)

    obj_tmp = eeg_add_noise(eeg, channel=channel, noise=noise)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    eeg_reset_components!(eeg)

    return nothing
end
