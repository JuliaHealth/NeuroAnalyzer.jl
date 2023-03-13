export gfilter

"""
    gfilter(signal, fs, pad, f, gw)

Filter using Gaussian in the frequency domain.

# Arguments

- `signal::AbstractVector`
- `fs::Int64`: sampling rate
- `pad::Int64=0`: number of zeros to add
- `f::Real`: filter frequency
- `gw::Real=5`: Gaussian width in Hz

# Returns

Named tuple containing:
- `s_f::Vector{Float64}`
"""
function gfilter(signal::AbstractVector; fs::Int64, pad::Int64=0, f::Real, gw::Real=5)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    f <= 0 && throw(ArgumentError("f must be > 0."))
    gw <= 0 && throw(ArgumentError("gw must be > 0."))

    # create Gaussian in frequency domain
    gf = linspace(0, fs, length(signal))
    gs = (gw * (2 * pi - 1)) / (4 * pi)     # normalized width
    gf .-= f                                # shifted frequencies
    # g = @. exp(-0.5 * (gf / gs)^2)        # Gaussian
    g = @. exp((-gf^2 ) / 2 * gs^2)         # Gaussian
    g ./= abs(maximum(g))                   # gain-normalized

    # filter
    return 2 .* ifft0(fft0(signal, pad) .* g, pad)
    
end

"""
    gfilter(obj; channel, pad, f, gw)

Filter using Gaussian in the frequency domain.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `pad::Int64=0`: number of zeros to add
- `f::Real`: filter frequency
- `gw::Real=5`: Gaussian width in Hz

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function gfilter(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), pad::Int64=0, f::Real, gw::Real=5)

    _check_channels(obj, channel)

    ep_n = epoch_n(obj)
    fs = sr(obj)

    obj_new = deepcopy(obj)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in eachindex(channel)
            obj_new.data[channel[ch_idx], :, ep_idx] = @views gfilter(obj_new.data[channel[ch_idx], :, ep_idx], fs=fs, pad=pad, f=f, gw=gw)
        end
    end

    reset_components!(obj_new)
    push!(obj_new.header.history, "wbp(OBJ, channel=$channel, pad=$pad, frq=$frq, ncyc=$ncyc, demean=$demean)")

    return obj_new
end

"""
    wbp!(obj; channel, pad, f, gw)

Filter using Gaussian in the frequency domain.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `pad::Int64=0`: number of zeros to add
- `f::Real`: filter frequency
- `gw::Real=5`: Gaussian width in Hz
"""
function wbp!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), pad::Int64=0, f::Real, gw::Real=5)

    obj_tmp = gfilter(obj, channel=channel, pad=pad, f=f, gw=gw)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    reset_components!(obj)

    return nothing
end
