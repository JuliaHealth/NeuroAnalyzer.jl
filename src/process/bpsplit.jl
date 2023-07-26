export bpsplit

"""
    bpsplit(obj; ch, order, window)

Split signal into frequency bands using IIR band-pass filter.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `order::Int64=8`: number of taps for FIR band-pass filter
- `window::Union{Nothing, AbstractVector, Int64}=nothing`: window for `:fir` filter; default is Hamming window, number of taps is calculated using Fred Harris' rule-of-thumb

# Returns

Named tuple containing:
- `s::Array{Float64, 4}`: split signal
- `bn::Vector{Symbol}`: band names
- `bf::Vector{Tuple{Real, Real}}`: band frequencies
"""
function bpsplit(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), order::Int64=8, window::Union{Nothing, AbstractVector, Int64}=nothing)
    
    bn = [:delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher]

    _check_channels(obj, ch)
    ch_n = length(ch)
    ep_n = epoch_n(obj)

    fs = sr(obj)
    s = zeros(length(bn), ch_n, epoch_len(obj), ep_n)
    bf = Vector{Tuple{Real, Real}}()

    @inbounds for band_idx in eachindex(bn)
        band_f = band_frq(obj, band=bn[band_idx])
        push!(bf, band_f)
        flt = filter_create(fs=fs, fprototype=:fir, ftype=:bp, cutoff=band_f, order=order, window=window, n=epoch_len(obj))
        @inbounds @simd for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
                s[band_idx, ch_idx, :, ep_idx] = @views filter_apply(obj.data[ch[ch_idx], :, ep_idx], flt=flt)
            end
        end
    end

    return (s=s, bn=bn, bf=bf)

end
