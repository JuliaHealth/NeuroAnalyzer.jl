export psd_slope

"""
    psd_slope(s; <keyword arguments>)

Calculate PSD linear fit and slope for a 1-D signal vector.

Default method is Welch's periodogram.

# Arguments

- `s::AbstractVector`: signal vector
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `flim::Tuple{Real, Real}=(0, fs / 2)`: calculate slope of the total power (default) or frequency range `flim[1]` to `flim[2]`
- `db::Bool=false`: if `true`, convert power to dB
- `method::Symbol=:welch`: PSD method:
- `:welch`: Welch's periodogram (default)
- `:fft`: plain FFT periodogram
- `:mt`: multi-tapered periodogram
- `:stft`: short-time Fourier transform averaged over segments
- `:mw`: Morlet wavelet convolution
- `:gh`: Gaussian filter + Hilbert transform
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `wlen::Int64=fs`: window length in samples (default = 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if `true`, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`)
- `demean::Bool=true`: subtract DC component before estimating PSD

# Returns

Named tuple:

- `lf::Vector{Float64}`: linear fit
- `ls::Float64`: slopes of linear fit
- `pf::Vector{Float64}`: range of frequencies for the linear fit
"""
function psd_slope(
    s::AbstractVector;
    fs::Int64,
    flim::Tuple{Real, Real} = (0, fs / 2),
    db::Bool = false,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool = true
)::@NamedTuple{
    lf::Vector{Float64},
    ls::Float64,
    pf::Vector{Float64}
}

    _check_tuple(flim, (0, fs / 2), "flim")

    psd_data = psd(
        s,
        fs = fs,
        db = db,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        gw = gw,
        demean = demean
    )
    pw = psd_data.p
    pf = psd_data.f

    f1_idx = vsearch(flim[1], pf)
    f2_idx = vsearch(flim[2], pf)
    pf = pf[f1_idx:f2_idx]
    pw = pw[f1_idx:f2_idx]

    lr = NeuroAnalyzer.linreg(pf, pw)
    lf = lr.lf
    ls = lf[2] - lf[1]

    return (; lf, ls, pf)

end

"""
    psd_slope(s; <keyword arguments>)

Calculate PSD linear fit and slope for a 3-D signal array.

Default method is Welch's periodogram.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `flim::Tuple{Real, Real}=(0, fs / 2)`: calculate slope of the total power (default) or frequency range `flim[1]` to `flim[2]`
- `db::Bool=false`: if `true`, convert power to dB
- `method::Symbol=:welch`: PSD method:
- `:welch`: Welch's periodogram (default)
- `:fft`: plain FFT periodogram
- `:mt`: multi-tapered periodogram
- `:stft`: short-time Fourier transform averaged over segments
- `:mw`: Morlet wavelet convolution
- `:gh`: Gaussian filter + Hilbert transform
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `wlen::Int64=fs`: window length in samples (default = 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if `true`, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`)
- `demean::Bool=true`: subtract DC component before estimating PSD

# Returns

Named tuple:

- `lf::Array{Float64, 3}`: linear fit
- `ls::Matrix{Float64}`: slope of linear fit
- `pf::Vector{Float64}`: range of frequencies for the linear fit
"""
function psd_slope(
    s::AbstractArray;
    fs::Int64,
    flim::Tuple{Real, Real} = (0, fs / 2),
    db::Bool = false,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool = true
)::@NamedTuple{
    lf::Array{Float64, 3},
    ls::Matrix{Float64},
    pf::Vector{Float64}
}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    lf = psd_slope(
        s[1, :, 1],
        fs = fs,
        flim = flim,
        db = db,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        gw = gw,
        demean = demean
    ).lf

    # pre-allocate outputs
    lf = zeros(ch_n, length(lf), ep_n)
    ls = zeros(ch_n, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        psd_slope_data = psd_slope(
            @view(s[ch_idx, :, ep_idx]),
            fs = fs,
            flim = flim,
            db = db,
            method = method,
            nt = nt,
            wlen = wlen,
            woverlap = woverlap,
            w = w,
            ncyc = ncyc,
            gw = gw,
            demean = demean
        )
        lf[ch_idx, :, ep_idx] = psd_slope_data.lf
        ls[ch_idx, ep_idx] = psd_slope_data.ls
    end

    return (; lf, ls, pf)

end

"""
    psd_slope(obj; <keyword arguments>)

Calculate PSD linear fit and slope for a NEURO object.

Default method is Welch's periodogram.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `flim::Tuple{Real, Real}=(0, sr(obj) / 2)`: calculate slope of the total power (default) or frequency range flim[1] to flim[2]
- `db::Bool=false`: if `true`, convert power to dB
- `method::Symbol=:welch`: PSD method:
- `:welch`: Welch's periodogram (default)
- `:fft`: plain FFT periodogram
- `:mt`: multi-tapered periodogram
- `:stft`: short-time Fourier transform averaged over segments
- `:mw`: Morlet wavelet convolution
- `:gh`: Gaussian filter + Hilbert transform
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if `true`, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`)
- `demean::Bool=true`: subtract DC component before estimating PSD

# Returns

Named tuple:

- `lf::Array{Float64, 3}`: linear fit
- `ls::Matrix{Float64}`: slope of linear fit
- `pf::Vector{Float64}`: range of frequencies for the linear fit
"""
function psd_slope(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    flim::Tuple{Real, Real} = (0, sr(obj) / 2),
    db::Bool = false,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool = true
)::@NamedTuple{
    lf::Array{Float64, 3},
    ls::Matrix{Float64},
    pf::Vector{Float64}
}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    return psd_slope(
        @view(obj.data[ch, :, :]),
        fs = sr(obj),
        flim = flim,
        db = db,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        gw = gw,
        demean = demean
    )

end
