export pli

"""
    pli(signal1, signal2)

Calculate PLI (Phase-Lag Index) between `signal1` and `signal2`.

# Arguments

- `signal1::AbstractVector`
- `signal2::AbstractVector`

# Returns

Named tuple containing:
- `pli_value::Float64`: PLI value
- `signal_diff::Vector{Float64}`: signal difference (signal2 - signal1)
- `phase_diff::Vector{Float64}`: phase difference (signal2 - signal1)
- `s1_phase::Vector{Float64}`: signal 1 phase
- `s2_phase::Vector{Float64}`: signal 2 phase
"""
function pli(signal1::AbstractVector, signal2::AbstractVector)

    length(signal1) == length(signal2) || throw(ArgumentError("Both signals must have the same length."))

    _, _, _, s1_phase = hspectrum(signal1)
    _, _, _, s2_phase = hspectrum(signal2)

    signal_diff = signal2 - signal1
    phase_diff = s2_phase - s1_phase

    pli_value = abs(mean(sign.(imag.(exp.(1im .* phase_diff)))))

    return (pli_value=pli_value, signal_diff=signal_diff, phase_diff=phase_diff, s1_phase=s1_phase, s2_phase=s2_phase)
end

"""
    pli(obj1, obj2; channel1, channel2, epoch1, epoch2)

Calculate PLI (Phase Lag Index) between `obj1` and `obj2`.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type]))`: index of channels, default is all channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type]))`: index of channels, default is all channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `pli_value::Array{Float64, 2}`: PLI value
- `signal_diff::Array{Float64, 3}`: signal difference (signal2 - signal1)
- `phase_diff::Array{Float64, 3}`: phase difference (signal2 - signal1)
- `s1_phase::Array{Float64, 3}`: signal 1 phase
- `s2_phase::Array{Float64, 3}`: signal 2 phase
"""
function pli(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)))

    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("obj1 and obj2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)

    pli_value = zeros(ch_n, ep_n)
    signal_diff = zeros(ch_n, epoch_len(obj1), ep_n)
    phase_diff = zeros(ch_n, epoch_len(obj1), ep_n)
    s1_phase = zeros(ch_n, epoch_len(obj1), ep_n)
    s2_phase = zeros(ch_n, epoch_len(obj1), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            pli_value[ch_idx, ep_idx], signal_diff[ch_idx, :, ep_idx], phase_diff[ch_idx, :, ep_idx], s1_phase[ch_idx, :, ep_idx], s2_phase[ch_idx, :, ep_idx] = @views pli(obj1.data[channel1[ch_idx], :, epoch1[ep_idx]], obj2.data[channel2[ch_idx], :, epoch2[ep_idx]])
        end
    end

    return (pli_value=pli_value, signal_diff=signal_diff, phase_dif=phase_diff, s1_phase=s1_phase, s2_phase=s2_phase)
end

"""
    pli(obj; channel)

Calculate PLIs (Phase Lag Index).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `pli_m::Array{Float64, 3}`: PLI value matrices over epochs
"""
function pli(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    pli_m = zeros(ch_n, ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_idx1
                pli_m[ch_idx1, ch_idx2, ep_idx], _, _, _, _ = @views pli(obj.data[channel[ch_idx1], :, ep_idx], obj.data[channel[ch_idx2], :, ep_idx])
            end
        end
        Threads.@threads for ch_idx1 in 1:(ch_n - 1)
            for ch_idx2 in (ch_idx1 + 1):ch_n
                pli_m[ch_idx1, ch_idx2, ep_idx] = @views pli_m[ch_idx2, ch_idx1, ep_idx]
            end
        end
    end

    return pli_m
end
