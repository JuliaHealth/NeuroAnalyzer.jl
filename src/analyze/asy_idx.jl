export asy_idx

"""
    asy_idx(obj; <keyword arguments>)

Computes the log-ratio and normalized difference of mean band power between two sets of channels (typically left vs right hemisphere):

- `asi = log(mean_power(ch1)) - log(mean_power(ch2))`: positive when ch1 has more power than ch2
- `nasi = (mean_power(ch1) - mean_power(ch2)) / (mean_power(ch1) + mean_power(ch2))`: bounded on (-1, 1), analogous to a contrast ratio

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch1::Union{String, Vector{String}, Regex}`: channel name(s), e.g. left frontal channels
- `ch2::Union{String, Vector{String}, Regex}`: channel name(s), e.g. right frontal channels
- `flim::Tuple{Real, Real}`: lower and upper frequency bounds
- `method::Symbol=:welch`: PSD method:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `gw::Real=5`: Gaussian width in Hz
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

Named tuple:

- `asi::Float64`: log band asymmetry
- `nasi::Float64`: normalized band asymmetry
"""
function asy_idx(
    obj::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}, Regex},
    ch2::Union{String, Vector{String}, Regex},
    flim::Tuple{Real, Real},
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool = true,
)::@NamedTuple{asi::Float64, nasi::Float64}

    # resolve channel names to integer indices
    ch1 = get_channel(obj, ch = ch1)
    ch2 = get_channel(obj, ch = ch2)

    # shared keyword arguments for both band_power calls — defined once to
    # avoid duplicating the argument list and to keep the two calls in sync
    bp_kwargs = (
        fs = sr(obj),
        flim = flim,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        gw = gw,
        demean = demean
    )

    _log_off()
    bp1 = band_power(@view(obj.data[ch1, :, :]); bp_kwargs...)
    bp2 = band_power(@view(obj.data[ch2, :, :]); bp_kwargs...)
    _log_on()

    # compute mean band power
    m1 = mean(bp1)
    m2 = mean(bp2)

    # log asymmetry
    asi = log(m1) - log(m2)

    # normalized asymmetry
    nasi = (m1 - m2) / (m1 + m2)

    return (asi = asi, nasi = nasi)

end
