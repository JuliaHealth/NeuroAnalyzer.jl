export cbp

"""
    cbp(s; <keyword arguments>)

Perform convolution band-pass filtering.

# Arguments

- `s::AbstractVector`
- `pad::Int64`: pad with `pad` zeros
- `frq::Real`: filter frequency
- `fs::Int64`: sampling rate

# Returns

- `cbp::Vector{Float64}`
"""
function cbp(s::AbstractVector; pad::Int64=0, frq::Real, fs::Int64)::Vector{Float64}

    @assert fs >= 1 "fs must be ≥ 1."
    @assert frq > 0 "frq must be > 0."
    @assert frq <= fs / 2 "frq must be ≤ $(fs / 2)."

    pad > 0 && (s = pad0(s, pad))

    kernel = generate_sine(frq, -1:1/fs:1)

    return tconv(s, kernel=kernel)

end

"""
    cbp(obj; <keyword arguments>)

Perform convolution bandpass filtering.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channel name or list of channel names
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function cbp(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, pad::Int64=0, frq::Real)::NeuroAnalyzer.NEURO

    ch = get_channel(obj, ch=ch)
    ep_n = nepochs(obj)
    fs = sr(obj)

    obj_new = deepcopy(obj)
    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in eachindex(ch)
            obj_new.data[ch[ch_idx], :, ep_idx] = @views cbp(obj_new.data[ch[ch_idx], :, ep_idx], pad=pad, frq=frq, fs=fs)
        end
    end

    reset_components!(obj_new)
    push!(obj_new.history, "cbp(OBJ, ch=$ch, pad=$pad, frq=$frq)")

    return obj_new

end

"""
    cbp!(obj; <keyword arguments>)

Perform convolution bandpass filtering.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channel name or list of channel names
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Tuple{Real, Real}`: filter frequency

# Returns

Nothing
"""
function cbp!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, pad::Int64=0, frq::Real)::Nothing

    obj_new = cbp(obj, ch=ch, pad=pad, frq=frq)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
