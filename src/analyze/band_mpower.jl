export band_mpower

"""
    band_mpower(s; <keyword arguments>)

Calculate mean and peak band power. For a given frequency band, computes four descriptors:

- `mbp` – mean power across the band
- `maxfrq – frequency of the peak (maximum) power bin within the band
- `maxbp` – power at that peak frequency bin
- `maxba` – amplitude at that peak frequency bin (√maxbp)

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
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `gw::Real=5`: Gaussian width in Hz
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

Named tuple containing:

- `mbp::Float64`: mean band power
- `maxfrq::Float64`: frequency of maximum band power
- `maxbp::Float64`: power at maximum band frequency
- `maxba::Float64`: amplitude at maximum band frequency
"""
function band_mpower(
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
)::@NamedTuple{mbp::Float64, maxfrq::Float64, maxbp::Float64, maxba::Float64}

    @assert fs >= 1 "fs must be ≥ 1."
    _check_tuple(flim, (0, fs / 2), "flim")

    # compute the power spectral density over the full frequency range
    pw, pf = psd(
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

    # find the PSD bin indices that bound the requested frequency band
    f1_idx = vsearch(flim[1], pf)
    f2_idx = vsearch(flim[2], pf)

    # Mean power: @view avoids allocating a copy of the band slice.
    mbp = mean(@view pw[f1_idx:f2_idx])

    # peak power: findmax returns (value, index) within the band view.
    maxbp, peak_offset = findmax(@view pw[f1_idx:f2_idx])
    # absolute index into pw / pf
    peak_idx = f1_idx + peak_offset - 1
    maxfrq = pf[peak_idx]

    # amplitude at peak: square root of power
    maxba = sqrt(maxbp)

    return (mbp = mbp, maxfrq = maxfrq, maxbp = maxbp, maxba = maxba)
end

"""
    band_mpower(s; <keyword arguments>)

Calculate mean and peak band power. For a given frequency band, computes four descriptors:

- `mbp` – mean power across the band
- `maxfrq – frequency of the peak (maximum) power bin within the band
- `maxbp` – power at that peak frequency bin
- `maxba` – amplitude at that peak frequency bin (√maxbp)

# Arguments

- `s::AbstractArray`: signal array (channels × samples × epochs)
- `fs::Int64`: sampling rate
- `flim::Tuple{Real, Real}`: lower and upper frequency bounds
- `method::Symbol=:welch`: PSD method:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `gw::Real=5`: Gaussian width in Hz
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

Named tuple containing:

- `mbp::Matrix{Float64}`: mean band power, shape `(channels, epochs)`
- `maxfrq::Matrix{Float64}`: frequency of maximum band power, shape `(channels, epochs)`
- `maxbp::Matrix{Float64}`: power at maximum band frequency, shape `(channels, epochs)`
- `maxba::Matrix{Float64}`: amplitude at maximum band frequency, shape `(channels, epochs)`
"""
function band_mpower(
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
    demean::Bool = true,
)::@NamedTuple{mbp::Matrix{Float64}, maxfrq::Matrix{Float64}, maxbp::Matrix{Float64}, maxba::Matrix{Float64}}

    # validate that the input is a proper 3-D array (channels × samples × epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    mbp = zeros(ch_n, ep_n)
    maxfrq = zeros(ch_n, ep_n)
    maxbp = zeros(ch_n, ep_n)
    maxba = zeros(ch_n, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        mpower_data = band_mpower(
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
        mbp[ch_idx, ep_idx] = mpower_data.mbp
        maxfrq[ch_idx, ep_idx] = mpower_data.maxfrq
        maxbp[ch_idx, ep_idx] = mpower_data.maxbp
        maxba[ch_idx, ep_idx] = mpower_data.maxba
    end

    return (mbp = mbp, maxfrq = maxfrq, maxbp = maxbp, maxba = maxba)

end

"""
    band_mpower(obj; <keyword arguments>)

Calculate mean and peak band power. For a given frequency band, computes four descriptors:

- `mbp` – mean power across the band
- `maxfrq – frequency of the peak (maximum) power bin within the band
- `maxbp` – power at that peak frequency bin
- `maxba` – amplitude at that peak frequency bin (√maxbp)

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `flim::Tuple{Real, Real}`: lower and upper frequency bounds
- `method::Symbol=:welch`: PSD method:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `gw::Real=5`: Gaussian width in Hz
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

Named tuple containing:

- `mbp::Matrix{Float64}`: mean band power, shape `(channels, epochs)`
- `maxfrq::Matrix{Float64}`: frequency of maximum band power, shape `(channels, epochs)`
- `maxbp::Matrix{Float64}`: power at maximum band frequency, shape `(channels, epochs)`
- `maxba::Matrix{Float64}`: amplitude at maximum band frequency, shape `(channels, epochs)`
"""
function band_mpower(
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
)::@NamedTuple{mbp::Matrix{Float64}, maxfrq::Matrix{Float64}, maxbp::Matrix{Float64}, maxba::Matrix{Float64}}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    _log_off()
    mpower_data = band_mpower(
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

    return mpower_data

end
