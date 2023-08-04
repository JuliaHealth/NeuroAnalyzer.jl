export band_mpower

"""
    band_mpower(s; fs, f)

Calculate mean and maximum band power and its frequency.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `f::Tuple{Real, Real}`: lower and upper frequency bounds

# Returns

Named tuple containing:
- `mbp::Float64`: mean band power
- `maxfrq::Float64`: frequency of maximum band power
- `maxbp::Float64`: power at maximum band frequency
"""
function band_mpower(s::AbstractVector; fs::Int64, f::Tuple{Real, Real}, mt::Bool=false)

    @assert fs >= 1 "fs must be ≥ 1."
    f = tuple_order(f)
    @assert f[1] >= 0 "Lower frequency bound must be ≥ 0."
    @assert f[2] <= fs / 2 "Lower frequency bound must be ≤ $(fs / 2)."

    if mt == true
        p = mt_pgram(s, fs=fs)
    else
        p = welch_pgram(s, 4*fs, fs=fs)
    end

    pf = Vector(p.freq)
    f1_idx = vsearch(f[1], pf)
    f2_idx = vsearch(f[2], pf)
    mbp = mean(p.power[f1_idx:f2_idx])
    maxfrq = pf[f1_idx:f2_idx][findmax(p.power[f1_idx:f2_idx])[2]]
    maxbp = p.power[vsearch(maxfrq, pf)]

    return (mbp=mbp, maxfrq=maxfrq, maxbp=maxbp)
end

"""
    band_mpower(s; fs, f, mt)

Calculate absolute band power between two frequencies.

# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

Named tuple containing:
- `mbp::Matrix{Float64}`: mean band power per channel per epoch
- `maxfrq::Matrix{Float64}`: frequency of maximum band power per channel per epoch
- `maxbp::Matrix{Float64}`: power at maximum band frequency per channel per epoch
"""
function band_mpower(s::AbstractArray; fs::Int64, f::Tuple{Real, Real}, mt::Bool=false)

    ch_n = size(s, 1)
    ep_n = size(s, 3)
    mbp = zeros(ch_n, ep_n)
    maxfrq = zeros(ch_n, ep_n)
    maxbp = zeros(ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            mbp[ch_idx, ep_idx], maxfrq[ch_idx, ep_idx], maxbp[ch_idx, ep_idx] = @views band_mpower(s[ch_idx, :, ep_idx], fs=fs, f=f, mt=mt)
        end
    end

    return (mbp=mbp, maxfrq=maxfrq, maxbp=maxbp)

end

"""
    band_mpower(obj; ch, f, mt)

Calculate mean and maximum band power and its frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `f::Tuple{Real, Real}`: lower and upper frequency bounds
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

Named tuple containing:
- `mbp::Matrix{Float64}`: mean band power per channel per epoch
- `maxfrq::Matrix{Float64}`: frequency of maximum band power per channel per epoch
- `maxbp::Matrix{Float64}`: power at maximum band frequency per channel per epoch
"""
function band_mpower(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), f::Tuple{Real, Real}, mt::Bool=false)

    _check_channels(obj, ch)
    
    mbp, maxfrq, maxbp = @views band_mpower(obj.data[ch, :, :], fs=sr(obj), f=f, mt=mt)

    return (mbp=mbp, maxfrq=maxfrq, maxbp=maxbp)

end
