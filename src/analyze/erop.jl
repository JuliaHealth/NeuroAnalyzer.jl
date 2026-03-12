export erop

"""
    erop(obj; <keyword arguments>)

Calculate ERO (Event-Related Oscillations) power spectrum.

Wraps psd() to compute an ERO power summary for a single channel:

ERP/ERF objects (two-column output):

- column 1 – power spectrum of the trial-averaged ERP/ERF (epoch 1)
- column 2 – mean power spectrum across all epochs (including epoch 1)

All other objects (one-column output):

- column 1 – mean power spectrum across all epochs

The two-column layout allows direct comparison between the "evoked" power (column 1, captures phase-locked activity) and the total power (column 2, captures both phase-locked and non-phase-locked activity).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name
- `method::Symbol=:welch`: power spectrum method:
- `:welch`: Welch's periodogram
- `:stft`: short-time Fourier transform
- `:mt`: multi-tapered periodogram
- `:fft`: Fast Fourier transform
- `:mw`: Morlet wavelet convolution
- `:gh`: Gaussian and Hilbert transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `db::Bool=true`: normalize powers to dB
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `gw::Real=5`: Gaussian width in Hz
- `demean::Bool=true`: subtract DC before calculating PSD

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
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    method::Symbol = :welch,
    db::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool = true,
)::@NamedTuple{p::Matrix{Float64}, f::Vector{Float64}}

    # compute per-epoch power spectra for the selected channel
    _log_off()
    psd_data = psd(
        obj,
        ch = ch,
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
    _log_on()

    p = psd_data.p[1, :, :]
    f = psd_data.f

    if datatype(obj) in ["erp", "erf"]
        p = cat(p[:, 1], mean(p, dims = 2), dims = 2)
    else
        p = mean(p, dims = 2)
    end

    return (p = p, f = f)

end
