export filter_g

"""
    filter_g(s; <keyword arguments>)

Filter using Gaussian in the frequency domain.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `pad::Int64=0`: number of zeros to add
- `f::Real`: filter frequency
- `gw::Real=5`: Gaussian width in Hz

# Returns

- `s_new::Vector{Float64}`
"""
function filter_g(s::AbstractVector; fs::Int64, pad::Int64=0, f::Real, gw::Real=5)::Vector{Float64}

    @assert fs >= 1 "fs must be ≥ 1."
    @assert f >= 0 "f must be ≥ 0."
    @assert gw > 0 "gw must be > 0."

    # create Gaussian in the frequency domain
    gf = linspace(0, fs, length(s))
    gs = (gw * (2 * pi - 1)) / (4 * pi)     # normalized width
    gf .-= f                                # shifted frequencies
    g = @. exp(-0.5 * (gf / gs)^2)          # Gaussian
    g ./= abs(maximum(g))                   # gain-normalized

    # filter
    s_tmp = fft0(s, pad) .* g
    s_new = abs.(ifft0(s_tmp, pad))

    return s_new

end

"""
    filter_g(s; <keyword arguments>)

Filter using Gaussian in the frequency domain.

# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `pad::Int64=0`: number of zeros to add
- `f::Real`: filter frequency
- `gw::Real=5`: Gaussian width in Hz

# Returns

- `s_new::Array{Float64, 3}`
"""
function filter_g(s::AbstractArray; fs::Int64, pad::Int64=0, f::Real, gw::Real=5)::Array{Float64, 3}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)
    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views filter_g(s[ch_idx, :, ep_idx], fs=fs, pad=pad, f=f, gw=gw)
        end
    end

    return s_new

end

"""
    filter_g(obj; <keyword arguments>)

Filter using Gaussian in the frequency domain.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `pad::Int64=0`: number of zeros to add
- `f::Real`: filter frequency
- `gw::Real=5`: Gaussian width in Hz

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function filter_g(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, pad::Int64=0, f::Real, gw::Real=5)::NeuroAnalyzer.NEURO

    ch = get_channel(obj, ch=ch)
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = @views filter_g(obj.data[ch, :, :], fs=sr(obj), pad=pad, f=f, gw=gw)
    reset_components!(obj_new)
    push!(obj_new.history, "filter_g(OBJ, ch=$ch, pad=$pad, f=$f)")

    return obj_new

end

"""
    filter_g!(obj; <keyword arguments>)

Filter using Gaussian in the frequency domain.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `pad::Int64=0`: number of zeros to add
- `f::Real`: filter frequency
- `gw::Real=5`: Gaussian width in Hz

# Returns

Nothing
"""
function filter_g!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, pad::Int64=0, f::Real, gw::Real=5)::Nothing

    obj_new = filter_g(obj, ch=ch, pad=pad, f=f, gw=gw)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
