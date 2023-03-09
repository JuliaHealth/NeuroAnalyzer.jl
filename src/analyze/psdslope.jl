export psdslope

"""
    psdslope(obj; channel, f, norm, mt)

Calculate PSD linear fit and slope.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all channels
- `f::Tuple{Real, Real}=(0, sr(obj)/2)`: calculate slope of the total power (default) or frequency range f[1] to f[2]
- `norm::Bool=false`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns

Named tuple containing:
- `lf::Array{Float64, 3}`: linear fit for each channel and epoch
- `psd_slope::Array{Float64, 2}`: slopes of each linear fit
- `frq::Vector{Float64}`: range of frequencies for the linear fits
"""
function psdslope(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), f::Tuple{Real, Real}=(0, sr(obj)/2), norm::Bool=false, mt::Bool=false, nt::Int64=8)

    fs = sr(obj)
    f = tuple_order(f)
    f[1] < 0 && throw(ArgumentError("Lower frequency bound must be be ≥ 0."))
    f[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be be < $(fs / 2)."))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    _, frq = s_psd(obj.data[1, :, 1], fs=fs, norm=norm, mt=mt, nt=nt)
    f1_idx = vsearch(f[1], frq)
    f2_idx = vsearch(f[2], frq)
    lf = zeros(ch_n, length(frq[f1_idx:f2_idx]), ep_n)
    psd_slope = zeros(ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            pow, _ = s_psd(obj.data[channel[ch_idx], :, ep_idx], fs=fs, norm=norm, mt=mt, nt=nt)
            _, _, _, _, _, _, lf[ch_idx, :, ep_idx] = @views linreg(frq[f1_idx:f2_idx], pow[f1_idx:f2_idx])
            psd_slope[ch_idx, ep_idx] = lf[ch_idx, 2, ep_idx] - lf[ch_idx, 1, ep_idx]
        end
    end

    return (lf=lf, psd_slope=psd_slope, frq=frq[f1_idx:f2_idx])
end
