export band_power

"""
    band_power(s; <keyword arguments>)

Calculate the absolute power in a frequency band by:

1. Estimating the PSD with the chosen method (Welch, FFT, MT, STFT, MW, GH)
2. Locating the bin indices that bracket the requested band
3. Integrating the PSD over that band using Simpson's rule

# Arguments

- `s::AbstractVector`: signal vector
- `fs::Int64`: sampling rate
- `flim::Tuple{Real, Real}`: lower and upper frequency bounds
- `method::Symbol=:welch`: PSD method:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length in samples, default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `gw::Real=5`: Gaussian width in Hz
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

- `bp::Float64`: band power
"""
function band_power(
    s::AbstractVector;
    fs::Int64,
    flim::Tuple{Real, Real},
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool = true,
)::Float64

    @assert fs >= 1 "fs must be ≥ 1."
    _check_tuple(flim, (0, fs / 2), "flim")

    # compute the power spectral density over the full frequency range
    psd_data = psd(
        s,
        fs = fs,
        db = false,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        gw = gw,
        demean = demean,
    )
    pow = psd_data.p
    frq = psd_data.f

    # locate the PSD bin indices that bracket the requested frequency band
    f1_idx = vsearch(flim[1], frq)
    f2_idx = vsearch(flim[2], frq)

    # frequency resolution: uniform bin spacing from the PSD
    dx = frq[2] - frq[1]

    # integrate
    bp = simpson(@view(pow[f1_idx:f2_idx]), @view(frq[f1_idx:f2_idx]), dx = dx)

    return bp

end

"""
    band_power(s; <keyword arguments>)

Calculate absolute band power between two frequencies.

# Arguments

- `s::AbstractArray`: signal array (channels, samples, epochs)
- `fs::Int64`: sampling rate
- `flim::Tuple{Real, Real}`: lower and upper frequency bounds
- `method::Symbol=:welch`: PSD method:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length in samples, default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `gw::Real=5`: Gaussian width in Hz
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

- `bp::Matrix{Float64}`: band power, shape `(channels, epochs)`
"""
function band_power(
    s::AbstractArray;
    fs::Int64,
    flim::Tuple{Real, Real},
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean = demean,
)::Matrix{Float64}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    bp = zeros(ch_n, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        bp[ch_idx, ep_idx] = band_power(
            @view(s[ch_idx, :, ep_idx]),
            fs = fs,
            flim = flim,
            method = method,
            nt = nt,
            wlen = wlen,
            woverlap = woverlap,
            w = w,
            ncyc = ncyc,
            gw = gw,
            demean = demean,
        )
    end

    return bp

end

"""
    band_power(obj; <keyword arguments>)

Calculate absolute band power between two frequencies.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `flim::Tuple{Real, Real}`: lower and upper frequency bounds
- `method::Symbol=:welch`: PSD method:
    - `:welch`: Welch's periodogram (default)
    - `:fft`: plain FFT periodogram
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform averaged over segments
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian filter + Hilbert transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `gw::Real=5`: Gaussian width in Hz
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

- `bp::Matrix{Float64}`: band power, shape `(channels, epochs)`
"""
function band_power(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex},
        flim::Tuple{Real, Real},
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = sr(obj),
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        gw::Real = 5,
        demean::Bool = true,
    )::Matrix{Float64}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    _log_off()
    power_data = band_power(
                    @view(obj.data[ch, :, :]),
                    fs = sr(obj),
                    flim = flim,
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

    return power_data

end
