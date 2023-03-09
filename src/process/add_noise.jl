export add_noise
export add_noise!

"""
    add_noise(signal, noise)

Adds noise.

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
    add_noise(obj; channel, noise)

Add noise.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `noise::AbstractVector`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function add_noise(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), noise::AbstractVector)

    _check_channels(obj, channel)

    ch_n = length(channel)
    ep_n = epoch_n(obj)

    obj_new = deepcopy(obj)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in eachindex(channel)
            @views obj_new.data[channel[ch_idx], :, ep_idx] = add_noise(obj_new.data[channel[ch_idx], :, ep_idx], noise=noise)
        end
    end

    reset_components!(obj_new)
    push!(obj_new.header.history, "add_noise(OBJ, channel=$channel)")

    return obj_new
end

"""
    add_noise!(obj; channel, noise)

Add noise.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `noise::AbstractVector`
"""
function add_noise!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), noise::AbstractVector)

    obj_tmp = add_noise(obj, channel=channel, noise=noise)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    reset_components!(obj)

    return nothing
end
