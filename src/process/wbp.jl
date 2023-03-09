export wbp

"""
    wbp(signal; pad, frq, fs, ncyc, demean)

Perform wavelet bandpass filtering.

# Arguments

- `signal::AbstractVector`
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `fs::Int64`: sampling rate
- `ncyc::Int64=6`: number of cycles for Morlet wavelet
- `demean::Bool=true`: demean signal prior to analysis

# Returns

- `signal_new::Vector{Float64}`
"""
function wbp(signal::AbstractVector; pad::Int64=0, frq::Real, fs::Int64, ncyc::Int64=6, demean::Bool=true)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    frq <= 0 && throw(ArgumentError("frq must be > 0."))
    ncyc <= 0 && throw(ArgumentError("ncyc must be > 0."))
    pad < 0 && throw(ArgumentError("pad must be ≥ 0."))
    frq > fs / 2 && throw(ArgumentError("frq must be ≤ $(fs / 2)."))

    pad > 0 && (signal = pad0(signal, pad))

    demean == true && (signal = demean(signal))

    kernel = generate_morlet(fs, frq, 1, ncyc=ncyc, complex=true)

    return real.(fconv(signal, kernel=kernel, norm=true))
end

"""
    wbp(obj; channel, pad, frq, ncyc, demean)

Perform wavelet bandpass filtering.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `ncyc::Int64=6`: number of cycles for Morlet wavelet
- `demean::Bool=true`: demean signal prior to analysis

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function wbp(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), pad::Int64=0, frq::Real, ncyc::Int64=6, demean::Bool=true)

    _check_channels(obj, channel)

    ep_n = epoch_n(obj)
    fs = sr(obj)

    obj_new = deepcopy(obj)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in eachindex(channel)
            obj_new.data[channel[ch_idx], :, ep_idx] = @views wbp(obj_new.data[channel[ch_idx], :, ep_idx], pad=pad, frq=frq, fs=fs, ncyc=ncyc, demean=demean)
        end
    end

    reset_components!(obj_new)
    push!(obj_new.header.history, "wbp(OBJ, channel=$channel, pad=$pad, frq=$frq, ncyc=$ncyc, demean=$demean)")

    return obj_new
end

"""
    wbp!(obj; channel, pad, frq, ncyc, demean)

Perform wavelet bandpass filtering.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `ncyc::Int64=6`: number of cycles for Morlet wavelet
- `demean::Bool=true`: demean signal prior to analysis
"""
function wbp!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), pad::Int64=0, frq::Real, ncyc::Int64=6, demean::Bool=true)

    obj_tmp = wbp(obj, channel=channel, pad=pad, frq=frq, ncyc=ncyc, demean=demean)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    reset_components!(obj)

    return nothing
end

