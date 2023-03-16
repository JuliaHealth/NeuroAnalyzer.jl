export cbp

"""
    cbp(s; pad, frq, fs)

Perform convolution band-pass filtering.

# Arguments

- `s::AbstractVector`
- `pad::Int64`: pad with `pad` zeros
- `frq::Real`: filter frequency
- `fs::Int64`: sampling rate

# Returns

- `cbp::Vector{Float64}`
"""
function cbp(s::AbstractVector; pad::Int64=0, frq::Real, fs::Int64)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    frq <= 0 && throw(ArgumentError("frq must be > 0."))
    frq > fs / 2 && throw(ArgumentError("frq must be ≤ $(fs / 2)."))

    pad > 0 && (s = pad0(s, pad))

    kernel = generate_sine(frq, -1:1/fs:1)

    return tconv(s, kernel=kernel)

end

"""
    cbp(obj; ch, pad, frq)

Perform convolution bandpass filtering.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function cbp(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), pad::Int64=0, frq::Real)

    _check_channels(obj, ch)

    ep_n = epoch_n(obj)
    fs = sr(obj)

    obj_new = deepcopy(obj)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in eachindex(ch)
            obj_new.data[ch[ch_idx], :, ep_idx] = @views cbp(obj_new.data[ch[ch_idx], :, ep_idx], pad=pad, frq=frq, fs=fs)
        end
    end

    reset_components!(obj_new)
    push!(obj_new.history, "cbp(OBJ, ch=$ch, pad=$pad, frq=$frq)")

    return obj_new
    
end

"""
    cbp!(obj; ch, pad, frq)

Perform convolution bandpass filtering.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Tuple{Real, Real}`: filter frequency
"""
function cbp!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), pad::Int64=0, frq::Real)

    obj_new = cbp(obj, ch=ch, pad=pad, frq=frq)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
