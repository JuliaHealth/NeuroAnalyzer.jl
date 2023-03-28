export psd_slope

"""
    psd_slope(s; fs, f, norm, mt, nt)

Calculate PSD linear fit and slope.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `f::Tuple{Real, Real}=(0, fs ÷ 2)`: calculate slope of the total power (default) or frequency range `f[1]` to `f[2]`
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns

Named tuple containing:
- `lf::Vector{Float64}`: linear fit
- `ls::Float64`: slopes of linear fit
- `pf::Vector{Float64}`: range of frequencies for the linear fit
"""
function psd_slope(s::AbstractVector; fs::Int64, f::Tuple{Real, Real}=(0, fs ÷ 2), norm::Bool=false, mt::Bool=false, nt::Int64=8)

    f = tuple_order(f)
    f[1] < 0 && throw(ArgumentError("Lower frequency bound must be be ≥ 0."))
    f[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be be < $(fs / 2)."))

    pw, pf = psd(s, fs=fs, norm=norm, mt=mt, nt=nt)

    f1_idx = vsearch(f[1], pf)
    f2_idx = vsearch(f[2], pf)
    _, _, _, _, _, _, lf = linreg(pf[f1_idx:f2_idx], pw[f1_idx:f2_idx])
    ls = lf[2] - lf[1]

    return (lf=lf, ls=ls, pf=pf[f1_idx:f2_idx])

end

"""
    psd_slope(s; fs, f, norm, mt, nt)

Calculate PSD linear fit and slope.

# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `f::Tuple{Real, Real}=(0, fs ÷ 2)`: calculate slope of the total power (default) or frequency range `f[1]` to `f[2]`
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns

Named tuple containing:
- `lf::Matrix{Float64}`: linear fit
- `s::Vector{Float64}`: slope of linear fit
- `pf::Vector{Float64}`: range of frequencies for the linear fit
"""
function psd_slope(s::AbstractArray; fs::Int64, f::Tuple{Real, Real}=(0, fs ÷ 2), norm::Bool=false, mt::Bool=false, nt::Int64=8)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    lf, ls, pf = psd_slope(s[1, :, 1], fs=fs, f=f, norm=norm, mt=mt, nt=nt)

    lf = zeros(ch_n, length(lf), ep_n)
    ls = zeros(ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            lf[ch_idx, :, ep_idx], ls[ch_idx, ep_idx], _ = psd_slope(s[ch_idx, :, ep_idx], fs=fs, f=f, norm=norm, mt=mt, nt=nt)
        end
    end

    return (lf=lf, ls=ls, pf=pf)

end

"""
    psd_slope(obj; ch, f, norm, mt)

Calculate PSD linear fit and slope.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `f::Tuple{Real, Real}=(0, sr(obj) ÷ 2)`: calculate slope of the total power (default) or frequency range f[1] to f[2]
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns

Named tuple containing:
- `lf::Array{Float64, 3}`: linear fit
- `ls::Array{Float64, 2}`: slope of linear fit
- `pf::Vector{Float64}`: range of frequencies for the linear fit
"""
function psd_slope(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), f::Tuple{Real, Real}=(0, sr(obj) ÷ 2), norm::Bool=false, mt::Bool=false, nt::Int64=8)

    _check_channels(obj, ch)

    lf, ls, pf = psd_slope(obj.data[ch, :, :], fs=sr(obj), f=f, norm=norm, mt=mt, nt=nt)

    return (lf=lf, ls=ls, pf=pf)

end
