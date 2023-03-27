export wbp
export wbp!

"""
    wbp(s; pad, frq, fs, ncyc)

Perform wavelet band-pass filtering.

# Arguments

- `s::AbstractVector`
- `pad::Int64`: pad with `pad` zeros
- `frq::Real`: filter frequency
- `fs::Int64`: sampling rate
- `ncyc::Int64=6`: number of cycles for Morlet wavelet

# Returns

- `s_new::Vector{Float64}`
"""
function wbp(s::AbstractVector; pad::Int64=0, frq::Real, fs::Int64, ncyc::Int64=6)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    frq <= 0 && throw(ArgumentError("frq must be > 0."))
    ncyc <= 0 && throw(ArgumentError("ncyc must be > 0."))
    pad < 0 && throw(ArgumentError("pad must be ≥ 0."))
    frq > fs / 2 && throw(ArgumentError("frq must be ≤ $(fs / 2)."))

    pad > 0 && (s = pad0(s, pad))

    kernel = generate_morlet(fs, frq, 1, ncyc=ncyc, complex=true)

    s_new = real.(fconv(s, kernel=kernel, norm=true))
    
    return s_new

end

"""
    wbp(s; ch, pad, frq, ncyc)

Perform wavelet band-pass filtering.

# Arguments

- `s::AbstractArray`
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `fs::Int64`: sampling rate
- `ncyc::Int64=6`: number of cycles for Morlet wavelet

# Returns

- `s_new::Array{Float64, 3`
"""
function wbp(s::AbstractArray; pad::Int64=0, frq::Real, fs::Int64, ncyc::Int64=6)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views wbp(s[ch_idx, :, ep_idx], pad=pad, frq=frq, fs=fs, ncyc=ncyc)
        end
    end

    return s_new

end

"""
    wbp(obj; ch, pad, frq, ncyc)

Perform wavelet band-pass filtering.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `ncyc::Int64=6`: number of cycles for Morlet wavelet

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function wbp(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), pad::Int64=0, frq::Real, ncyc::Int64=6)

    _check_channels(obj, ch)

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = @views wbp(obj.data[ch, :, :], pad=pad, frq=frq, fs=sr(obj), ncyc=ncyc)
    reset_components!(obj_new)
    push!(obj_new.history, "wbp(OBJ, ch=$ch, pad=$pad, frq=$frq, ncyc=$ncyc)")

    return obj_new

end

"""
    wbp!(obj; ch, pad, frq, ncyc)

Perform wavelet band-pass filtering.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `ncyc::Int64=6`: number of cycles for Morlet wavelet
"""
function wbp!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), pad::Int64=0, frq::Real, ncyc::Int64=6)

    obj_new = wbp(obj, ch=ch, pad=pad, frq=frq, ncyc=ncyc)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end

