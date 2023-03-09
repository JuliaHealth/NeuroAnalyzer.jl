export fbsplit

"""
    fbsplit(obj; channel, order, window)

Split signal into frequency bands.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `order::Int64=8`: number of taps for FIR band-pass filter
- `window::Union{Nothing, AbstractVector, Int64}=nothing`: window for `:fir` filter; default is Hamming window, number of taps is calculated using fred harris' rule-of-thumb

# Returns

Named tuple containing:
- `band_names::Vector{Symbol}`
- `band_frq::Vector{Tuple{Real, Real}}`
- `signal_split::Array{Float64, 4}`
"""
function fbsplit(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), order::Int64=8, window::Union{Nothing, AbstractVector, Int64}=nothing)
    
    band = [:delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher]

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    fs = sr(obj)
    signal_split = zeros(length(band), ch_n, epoch_len(obj), ep_n)
    band_frq = Vector{Tuple{Real, Real}}()

    @inbounds for band_idx in eachindex(band)
        band_f = band(obj, band=band[band_idx])
        push!(band_frq, band_f)
        flt = filter_create(fs=fs, fprototype=:fir, ftype=:bp, cutoff=band_f, order=order, window=window, n=epoch_len(obj))
        @inbounds @simd for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
                signal_split[band_idx, ch_idx, :, ep_idx] = @views filter_apply(obj.data[channel[ch_idx], :, ep_idx], flt=flt)
            end
        end
    end

    return (band_names=band, band_frq=band_frq, signal_split=signal_split)
end
