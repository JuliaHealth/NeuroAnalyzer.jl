export bpsplit

"""
    bpsplit(obj; <keyword arguments>)

Split a signal into standard EEG frequency bands using band-pass filters.

One filter is designed per band (see `filter_create`) and applied to every selected channel and epoch. A band whose upper limit reaches the Nyquist frequency (e.g. `:gamma` at 256 Hz) is extracted with a high-pass filter at its lower limit.

Available bands: `:delta`, `:theta`, `:alpha`, `:alpha_lower`, `:alpha_higher`, `:beta`, `:beta_lower`, `:beta_higher`, `:gamma`, `:gamma_1`, `:gamma_2`, `:gamma_lower`, `:gamma_higher` (frequency limits: see `band_frq`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `band::Vector{Symbol}`: bands to extract; default: all bands listed above
- `fprototype::Symbol=:fir`: filter prototype: `:fir`, `:firls`, `:remez`, `:butterworth`, `:chebyshev1`, `:chebyshev2` or `:elliptic`
- `order::Union{Nothing, Int64}=nothing`: filter order (number of taps for FIR, must be odd); for FIR prototypes, if `nothing` it is calculated from `bw`; required for IIR prototypes
- `bw::Union{Nothing, Real}=nothing`: transition band width in Hz, centred on each cutoff (FIR prototypes only); if `nothing`, it is set per band to `min(2, f1, (f2 - f1) / 2)` Hz, where `(f1, f2)` are the band limits, so that the transition bands neither reach 0 Hz nor overlap
- `rp::Union{Nothing, Real}=nothing`: pass-band ripple in dB (`:chebyshev1`, `:elliptic`)
- `rs::Union{Nothing, Real}=nothing`: stop-band attenuation in dB (`:chebyshev2`, `:elliptic`; target for automatic `:firls`/`:remez` order)
- `w::Union{Nothing, AbstractVector}=nothing`: window vector for `:fir` (its length sets the number of taps) or weight vector for `:firls` (6 elements for band-pass, 4 for high-pass)
- `window::Symbol=:hamming`: window for `:fir` when `w` is not given (`:rect`, `:hann`, `:hamming`, `:blackman`)
- `dir::Symbol=:twopass`: filtering direction:
    - `:twopass`: forward pass followed by reverse pass (zero phase distortion; magnitude response is squared)
    - `:onepass`: single forward pass (introduces phase delay)
    - `:reverse`: single reverse pass

# Returns

Named tuple:

- `s::Array{Float64, 4}`: band-split signal of shape `(n_bands, n_channels, epoch_len, n_epochs)`
- `bn::Vector{Symbol}`: band name symbols (in the order used in `s`)
- `bf::Vector{Tuple{Real, Real}}`: frequency limits `(f_low, f_high)` in Hz for each band, in the same order as `bn`

# Notes

Low bands need long FIR filters (e.g. 1691 taps for `:delta` at 256 Hz with the default `bw = 0.5` Hz). Filter the continuous signal before epoching, or use an IIR prototype for short epochs.
"""
function bpsplit(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    band::Vector{Symbol} = [
        :delta,
        :theta,
        :alpha,
        :alpha_lower,
        :alpha_higher,
        :beta,
        :beta_lower,
        :beta_higher,
        :gamma,
        :gamma_1,
        :gamma_2,
        :gamma_lower,
        :gamma_higher,
    ],
    fprototype::Symbol = :fir,
    order::Union{Nothing, Int64} = nothing,
    bw::Union{Nothing, Real} = nothing,
    rp::Union{Nothing, Real} = nothing,
    rs::Union{Nothing, Real} = nothing,
    w::Union{Nothing, AbstractVector} = nothing,
    window::Symbol = :hamming,
    dir::Symbol = :twopass,
)::@NamedTuple{
    s::Array{Float64, 4},
    bn::Vector{Symbol},
    bf::Vector{Tuple{Real, Real}},
}
    # validate
    _check_var(
        fprototype,
        [:fir, :firls, :remez, :butterworth, :chebyshev1, :chebyshev2, :elliptic],
        "fprototype",
    )
    _check_var(dir, [:twopass, :onepass, :reverse], "dir")
    isempty(band) && throw(ArgumentError("band must not be empty."))
    for b in band
        _check_var(
            b,
            [
                :delta,
                :theta,
                :alpha,
                :alpha_lower,
                :alpha_higher,
                :beta,
                :beta_lower,
                :beta_higher,
                :gamma,
                :gamma_1,
                :gamma_2,
                :gamma_lower,
                :gamma_higher,
            ],
            "band",
        )
    end
    bn = band
    fir = fprototype in (:fir, :firls, :remez)

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)
    # sampling rate
    fs = sr(obj)
    # Nyquist frequency
    nqf = fs / 2
    # epoch lengths
    el = epoch_len(obj)

    # pre-allocate output
    s = zeros(length(bn), ch_n, el, ep_n)
    bf = Vector{Tuple{Real, Real}}(undef, length(bn))

    # design one filter per band, then apply it to all channels and epochs.
    # the outer band loop is sequential (each band uses a different filter)
    # the inner channel loop is parallelized
    @inbounds for band_idx in eachindex(bn)
        band_f = band_frq(obj; band = bn[band_idx])
        bf[band_idx] = band_f
        f1, f2 = band_f
        # default transition width: within (0, f1) and not overlapping the other edge
        bw_band = fir ? something(bw, min(2.0, f1, (f2 - f1) / 2)) : nothing
        # upper limit at (clamped to) Nyquist: high-pass at the lower limit
        hp = f2 >= nqf - 0.1 || (fir && f2 + bw_band / 2 >= nqf)
        flt = filter_create(;
            fs = fs,
            fprototype = fprototype,
            ftype = hp ? :hp : :bp,
            cutoff = hp ? f1 : band_f,
            order = order,
            bw = bw_band,
            rp = rp,
            rs = rs,
            w = w,
            window = window,
        )
        if flt isa Vector{Float64} && el <= 3 * (length(flt) - 1)
            _warn(
                "$(bn[band_idx]): epoch length ($el samples) is short relative to the filter ($(length(flt)) taps): " *
                "edge transients will affect a large part of each epoch.",
            )
        end

        # calculate over channel and epochs
        @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
            ch_idx, ep_idx = idx[1], idx[2]
            s[band_idx, ch_idx, :, ep_idx] = _filter_one(obj, ch[ch_idx], ep_idx, flt, dir)
        end
    end

    return (; s, bn, bf)
end
