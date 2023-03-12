export band_power

"""
    band_power(s; fs, f, mt)

Calculate absolute band power between two frequencies.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

- `bp::Float64`: band power
"""
function band_power(s::AbstractVector; fs::Int64, f::Tuple{Real, Real}, mt::Bool=false)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    f = tuple_order(f)
    f[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0.")) 
    f[2] > fs / 2 && throw(ArgumentError("Lower frequency bound must be ≤ $(fs / 2).")) 

    if mt == true
        p = mt_pgram(s, fs=fs)
    else
        p = welch_pgram(s, 4*fs, fs=fs)
    end

    pf = Vector(p.freq)
    f1_idx = vsearch(f[1], pf)
    f2_idx = vsearch(f[2], pf)
    frq_idx = [f1_idx, f2_idx]

    # dx: frequency resolution
    dx = pf[2] - pf[1]

    # integrate
    return simpson(p.power[frq_idx[1]:frq_idx[2]], pf[frq_idx[1]:frq_idx[2]], dx=dx)

end

"""
    band_power(s; fs, f, mt)

Calculate absolute band power between two frequencies.

# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

- `bp::Matrix{Float64}`: band power for each channel per epoch
"""
function band_power(s::AbstractArray; fs::Int64, f::Tuple{Real, Real}, mt::Bool=false)

    ch_n = size(s, 1)
    ep_n = size(s, 3)
    bp = zeros(ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            bp[ch_idx, ep_idx] = @views band_power(s[ch_idx, :, ep_idx], fs=fs, f=f, mt=mt)
        end
    end

    return bp

end

"""
    band_power(obj; ch, f, mt)

Calculate absolute band power between two frequencies.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

- `bp::Matrix{Float64}`: band power for each ch per epoch
"""
function band_power(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), f::Tuple{Real, Real}, mt::Bool=false)

    _check_channels(obj, ch)

    bp = @views band_power(obj.data[ch, :, :], fs=sr(obj), f=f, mt=mt)

    return bp

end
