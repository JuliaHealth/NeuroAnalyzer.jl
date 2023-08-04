export filter_g

"""
    filter_g(s, fs, pad, f, gw)

Filter using Gaussian in the frequency domain.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `pad::Int64=0`: number of zeros to add
- `f::Real`: filter frequency
- `gw::Real=5`: Gaussian width in Hz

# Returns

Named tuple containing:
- `s_filtered::Vector{Float64}`
"""
function filter_g(s::AbstractVector; fs::Int64, pad::Int64=0, f::Real, gw::Real=5)

    @assert fs >= 1 "fs must be ≥ 1."
    @assert f > 0 "f must be > 0."
    @assert gw > 0 "gw must be > 0."

    # create Gaussian in frequency domain
    gf = linspace(0, fs, length(s))
    gs = (gw * (2 * pi - 1)) / (4 * pi)     # normalized width
    gf .-= f                                # shifted frequencies
    # g = @. exp(-0.5 * (gf / gs)^2)        # Gaussian
    g = @. exp((-gf^2 ) / 2 * gs^2)         # Gaussian
    g ./= abs(maximum(g))                   # gain-normalized

    # filter
    s_filtered = 2 .* ifft0(fft0(s, pad) .* g, pad)
    
    return s_filtered

end

"""
    filter_g(s; pad, f, gw)

Filter using Gaussian in the frequency domain.

# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `pad::Int64=0`: number of zeros to add
- `f::Real`: filter frequency
- `gw::Real=5`: Gaussian width in Hz

# Returns

- `s_filtered::NeuroAnalyzer.NEURO`
"""
function filter_g(s::AbstractArray; fs::Int64, pad::Int64=0, f::Real, gw::Real=5)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_filtered = similar(s)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_filtered[ch_idx, :, ep_idx] = @views filter_g(s[ch_idx, :, ep_idx], fs=fs, pad=pad, f=f, gw=gw)
        end
    end

    return s_filtered

end

"""
    filter_g(obj; ch, pad, f, gw)

Filter using Gaussian in the frequency domain.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `pad::Int64=0`: number of zeros to add
- `f::Real`: filter frequency
- `gw::Real=5`: Gaussian width in Hz

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function filter_g(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), pad::Int64=0, f::Real, gw::Real=5)

    _check_channels(obj, ch)
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = @views filter_g(obj.data[ch, :, :], fs=sr(obj), pad=pad, f=f, gw=gw)
    reset_components!(obj_new)
    push!(obj_new.history, "filter_g(OBJ, ch=$ch, pad=$pad, f=$f)")

    return obj_new

end

"""
    filter_g!(obj; ch, pad, f, gw)

Filter using Gaussian in the frequency domain.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `pad::Int64=0`: number of zeros to add
- `f::Real`: filter frequency
- `gw::Real=5`: Gaussian width in Hz
"""
function filter_g!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), pad::Int64=0, f::Real, gw::Real=5)

    obj_new = filter_g(obj, ch=ch, pad=pad, f=f, gw=gw)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
