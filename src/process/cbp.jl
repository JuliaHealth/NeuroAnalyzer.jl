export cbp

"""
    cbp(signal; pad, frq, fs, demean)

Perform convolution bandpass filtering.

# Arguments

- `signal::AbstractVector`
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `fs::Int64`: sampling rate
- `demean::Bool=true`: demean signal prior to analysis

# Returns

- `signal_new::Vector{Float64}`
"""
function cbp(signal::AbstractVector; pad::Int64=0, frq::Real, fs::Int64, demean::Bool=true)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    frq <= 0 && throw(ArgumentError("frq must be > 0."))
    frq > fs / 2 && throw(ArgumentError("frq must be ≤ $(fs / 2)."))

    pad > 0 && (signal = pad0(signal, pad))

    demean == true && (signal = demean(signal))

    kernel = generate_sine(frq, -1:1/fs:1)

    return tconv(signal, kernel=kernel)
end

"""
    cbp(obj; channel, pad, frq, demean)

Perform convolution bandpass filtering.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `demean::Bool=true`: demean signal prior to analysis

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function cbp(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), pad::Int64=0, frq::Real, demean::Bool=true)

    _check_channels(obj, channel)

    ep_n = epoch_n(obj)
    fs = sr(obj)

    obj_new = deepcopy(obj)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in eachindex(channel)
            obj_new.data[channel[ch_idx], :, ep_idx] = @views s_cbp(obj_new.data[channel[ch_idx], :, ep_idx], pad=pad, frq=frq, fs=fs, demean=demean)
        end
    end

    reset_components!(obj_new)
    push!(obj_new.header.history, "cbp(OBJ, channel=$channel, pad=$pad, frq=$frq, demean=$demean)")

    return obj_new
end

"""
    cbp!(obj; channel, pad, frq, demean)

Perform convolution bandpass filtering.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Tuple{Real, Real}`: filter frequency
- `demean::Bool=true`: demean signal prior to analysis
"""
function cbp!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), pad::Int64=0, frq::Real, demean::Bool=true)

    obj_tmp = cbp(obj, channel=channel, pad=pad, frq=frq, demean=demean)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    reset_components!(obj)

    return nothing
end
