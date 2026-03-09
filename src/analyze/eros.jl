export eros

"""
    eros(obj; <keyword arguments>)

Calculate ERO (Event-Related Oscillations) spectrogram.

Wraps NeuroAnalyzer.spectrogram() to compute an ERO time-frequency summary for a single channel:

ERP/ERF objects (two-slice output, along dim 3):

- s[:, :, 1] – spectrogram of the trial-averaged ERP/ERF (epoch 1)
- s[:, :, 2] – mean spectrogram across all epochs

All other objects (one-slice output):

- s[:, :, 1] – mean spectrogram across all epochs

The two-slice layout mirrors erop() and allows comparison between phase-locked (evoked) and total (evoked + induced) time-frequency power.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `method::Symbol=:stft`: spectrogram method:
  - `:stft`: short-time Fourier transform
  - `:mt`: multi-tapered periodogram
  - `:mw`: Morlet wavelet convolution
  - `:gh`: Gaussian and Hilbert transform
  - `:cwt`: continuous wavelet transformation
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length, default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `pad::Int64=0`: number of zeros to append
- `db::Bool=true`: normalize powers to dB
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

Named tuple containing:

- `s::Array{Float64, 3}`: spectrogram(s)
- `f::Vector{Float64}`: frequencies
- `t::Vector{Float64}`: time
"""
function eros(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    pad::Int64 = 0,
    method::Symbol = :stft,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    db::Bool = true,
    gw::Real = 5,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    wt::T = wavelet(Morlet(2π), β = 2),
)::@NamedTuple{s::Array{Float64, 3}, f::Vector{Float64}, t::Vector{Float64}} where {T <: CWT}

    # compute per-epoch power spectra for the selected channel
    _log_off()
    spec_data = NeuroAnalyzer.spectrogram(
        obj,
        ch = ch,
        method = method,
        nt = nt,
        pad = pad,
        db = db,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        gw = gw,
        wt = wt,
    )
    _log_on()
    f = spec_data.f
    t = spec_data.t
    # (freq, time, epochs)
    s = spec_data.p[1, :, :, :]

    if datatype(obj) in ["erp", "erf"]
        s = cat(s[:, :, 1], dropdims(mean(s, dims = 3), dims = 3), dims = 3)
    else
        s = mean(s, dims = 3)
    end

    return (s = s, f = f, t = t)

end
