export erop

"""
    erop(obj; <keyword arguments>)

Calculate ERO (Event-Related Oscillations) power spectrum.

Wraps psd() to compute an ERO power summary for a single channel:

ERP/ERF objects (two-column output):

- column 1 – power spectrum of the trial-averaged ERP/ERF (epoch 1)
- column 2 – mean power spectrum of the single trials (epochs 2:end)

All other objects (one-column output):

- column 1 – mean power spectrum across all epochs

The two-column layout allows direct comparison between the "evoked" power (column 1, captures phase-locked activity) and the total power (column 2, captures both phase-locked and non-phase-locked activity).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `exclude_bad_segments::Bool=NeuroAnalyzer.exclude_bad_segments`: if `true`, bad segments (markers of type `"BAD"` of all channels or of the selected channels) of a continuous signal are not used: their samples are removed and the remaining parts are joined (default: the `exclude_bad_segments` preference, see `set_exclude_bad_segments()`)
- `exclude_bad_channels::Bool=NeuroAnalyzer.exclude_bad_channels`: if `true`, channels marked as bad are not used (default: the `exclude_bad_channels` preference, see `set_exclude_bad_channels()`)
- `method::Symbol=:welch`: power spectrum method:
    - `:welch`: Welch's periodogram
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:fft`: Fast Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `wlen::Int64=sr(obj)`: window length in samples (default = 1 second); used by `:welch` and `:stft`
- `woverlap::Int64=min(min(round(Int64, wlen * 0.90), wlen - 1), wlen - 1)`: window overlap in samples
- `w::Bool=true`: if `true`, apply Hanning window
- `db::Bool=true`: normalize powers to dB
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of Morlet wavelet cycles; for a tuple the number of cycles is logarithmically (geometrically) spaced from `ncyc[1]` to `ncyc[2]` across frequencies
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`)
- `demean::Bool=true`: subtract DC component before estimating PSD

# Returns

Named tuple:

- `p::Matrix{Float64}`: powers, shape `(frequencies, 1)` or `(frequencies, 2)`
- `f::Vector{Float64}`: frequencies
"""
function erop(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = min(round(Int64, wlen * 0.90), wlen - 1),
    w::Bool = true,
    method::Symbol = :welch,
    db::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool = true,
    exclude_bad_channels::Bool = NeuroAnalyzer.exclude_bad_channels,
    exclude_bad_segments::Bool = NeuroAnalyzer.exclude_bad_segments,
)::@NamedTuple{
    p::Matrix{Float64},
    f::Vector{Float64},
}
    # validate
    isempty(ch) && throw(ArgumentError("No channels selected."))
    # exactly one channel, not a bad channel if excluded
    _get_single_ch(obj, ch, exclude_bad_channels)

    # compute per-epoch power spectra for the selected channel
    psd_data = psd(
        obj;
        ch = ch, exclude_bad_segments = exclude_bad_segments, exclude_bad_channels = exclude_bad_channels,
        db = db,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        gw = gw,
        demean = demean,
    )

    p = psd_data.p[1, :, :]
    f = psd_data.f

    if datatype(obj) in ["erp", "erf"]
        nepochs(obj) > 1 ||
            throw(ArgumentError("ERP/ERF OBJ must contain single trials (epochs 2:end)."))
        p = cat(p[:, 1], mean(p[:, 2:end]; dims = 2); dims = 2)
    else
        p = mean(p; dims = 2)
    end

    return (; p, f)
end
