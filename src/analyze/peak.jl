export peak_frq
export peak_amp
export peak_pow

"""
    peak_frq(s; <keyword arguments>)

Calculate peak frequency within a frequency band.

# Arguments

- `s::AbstractVector`: signal vector
- `fs::Int64`: sampling rate
- `flim::Tuple{Real, Real}`: lower and upper frequency bounds in Hz
- `method::Symbol=:welch`: PSD method:
  - `:welch`: Welch's periodogram
  - `:fft`: fast Fourier transform
  - `:mt`: multi-tapered periodogram
  - `:stft`: short-time Fourier transform
  - `:mw`: Morlet wavelet convolution
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length in samples, default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

- `pf::Float64`: peak frequency
"""
function peak_frq(
        s::AbstractVector;
        fs::Int64,
        flim::Tuple{Real, Real},
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = fs,
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        demean::Bool = true,
    )::Float64

    @assert fs >= 1 "fs must be ≥ 1."
    _check_tuple(flim, (0, fs / 2), "flim")

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
        demean = demean
    )
    pw = psd_data.p
    pf = psd_data.f

    f1_idx = vsearch(flim[1], pf)
    f2_idx = vsearch(flim[2], pf)

    peak_offset = argmax(@view(pw[f1_idx:f2_idx]))
    pf = pf[f1_idx + peak_offset - 1]

    return pf


end

"""
    peak_frq(s; <keyword arguments>)

Calculate peak frequency within a frequency band.

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
- `wlen::Int64=fs`: window length in samples, default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

- `pf::Matrix{Float64}`: peak frequency, shape `(channels, epochs)`
"""
function peak_frq(
    s::AbstractArray;
    fs::Int64,
    flim::Tuple{Real, Real},
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    demean::Bool = true,
)::Matrix{Float64}


    # validate that the input is a proper 3-D array (channels × samples × epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate outputs
    pf = zeros(ch_n, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        pf[ch_idx, ep_idx] = peak_frq(
            @view(s[ch_idx, :, ep_idx]),
            fs = fs,
            flim = flim,
            method = method,
            nt = nt,
            wlen = wlen,
            woverlap = woverlap,
            w = w,
            ncyc = ncyc,
            demean = demean,
        )
    end

    return pf

end

"""
    peak_frq(obj; <keyword arguments>)

Calculate peak frequency within a frequency band.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `flim::Tuple{Real, Real}`: lower and upper frequency bounds
- `method::Symbol=:welch`: PSD method:
  - `:welch`: Welch's periodogram
  - `:fft`: fast Fourier transform
  - `:mt`: multi-tapered periodogram
  - `:stft`: short-time Fourier transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

- `pf::Matrix{Float64}`: peak frequency, shape `(channels, epochs)`
"""
function peak_frq(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    flim::Tuple{Real, Real},
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    demean::Bool = true,
)::Matrix{Float64}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    pf = peak_frq(
        @view(obj.data[ch, :, :]),
        fs = sr(obj),
        flim = flim,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        demean = demean,
    )

    return pf

end


"""
    peak_amp(s; <keyword arguments>)

Calculate amplitude at the peak frequency within a frequency band.

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
- `wlen::Int64=fs`: window length in samples, default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

- `pa::Float64`: amplitude at peak frequency
"""
function peak_amp(
    s::AbstractVector;
    fs::Int64,
    flim::Tuple{Real, Real},
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    demean::Bool = true,
)::Float64

    @assert fs >= 1 "fs must be ≥ 1."
    _check_tuple(flim, (0, fs / 2), "flim")

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
        demean = demean
    )
    pw = psd_data.p
    pf = psd_data.f

    f1_idx = vsearch(flim[1], pf)
    f2_idx = vsearch(flim[2], pf)

    # amplitude = √(peak power) in the band
    return sqrt(maximum(@view(pw[f1_idx:f2_idx])))

end

"""
    peak_amp(s; <keyword arguments>)

Calculate amplitude at peak frequency within a frequency band.

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
- `wlen::Int64=fs`: window length in samples, default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

- `pa::Matrix{Float64}`: amplitude at peak frequency, shape `(channels, epochs)`
"""
function peak_amp(
    s::AbstractArray;
    fs::Int64,
    flim::Tuple{Real, Real},
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    demean::Bool = true,
)::Matrix{Float64}

    # validate that the input is a proper 3-D array (channels × samples × epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate outputs
    pa = zeros(ch_n, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        pa[ch_idx, ep_idx] = peak_amp(
            @view(s[ch_idx, :, ep_idx]),
            fs = fs,
            flim = flim,
            method = method,
            nt = nt,
            wlen = wlen,
            woverlap = woverlap,
            w = w,
            ncyc = ncyc,
            demean = demean,
        )
    end

    return pa

end

"""
    peak_amp(obj; <keyword arguments>)

Calculate amplitude at peak frequency within a frequency band.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `flim::Tuple{Real, Real}`: lower and upper frequency bounds
- `method::Symbol=:welch`: PSD method:
  - `:welch`: Welch's periodogram
  - `:fft`: fast Fourier transform
  - `:mt`: multi-tapered periodogram
  - `:stft`: short-time Fourier transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

- `pa::Matrix{Float64}`: amplitude at peak frequency, shape `(channels, epochs)`
"""
function peak_amp(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    flim::Tuple{Real, Real},
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    demean::Bool = true,
)::Matrix{Float64}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    pa = peak_amp(
        @view(obj.data[ch, :, :]),
        fs = sr(obj),
        flim = flim,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        demean = demean,
    )

    return pa

end

"""
    peak_pow(s; <keyword arguments>)

Calculate power at the peak frequency within a frequency band.

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
- `wlen::Int64=fs`: window length in samples, default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

- `pp::Float64`: power at peak frequency
"""
function peak_pow(
        s::AbstractVector;
        fs::Int64,
        flim::Tuple{Real, Real},
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = fs,
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        demean::Bool = true,
    )::Float64

    @assert fs >= 1 "fs must be ≥ 1."
    _check_tuple(flim, (0, fs / 2), "flim")

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
        demean = demean
    )
    pw = psd_data.p
    pf = psd_data.f

    f1_idx = vsearch(flim[1], pf)
    f2_idx = vsearch(flim[2], pf)

    return maximum(@view(pw[f1_idx:f2_idx]))

end


"""
    peak_pow(s; <keyword arguments>)

Calculate power at peak frequency within a frequency band.

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
- `wlen::Int64=fs`: window length in samples, default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

- `pp::Matrix{Float64}`: peak power, shape `(channels, epochs)`
"""
function peak_pow(
    s::AbstractArray;
    fs::Int64,
    flim::Tuple{Real, Real},
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    demean::Bool = true,
)::Matrix{Float64}

    # validate that the input is a proper 3-D array (channels × samples × epochs)
    _chk3d(s)
    @assert size(s, 1) == 1 "s must have 1 channel."

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    pp = zeros(ch_n, ep_n)

    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        pp[ch_idx, ep_idx] = peak_pow(
            @view(s[ch_idx, :, ep_idx]),
            fs = fs,
            flim = flim,
            method = method,
            nt = nt,
            wlen = wlen,
            woverlap = woverlap,
            w = w,
            ncyc = ncyc,
            demean = demean,
        )
    end

    return pp

end

"""
    peak_pow(obj; <keyword arguments>)

Calculate power at peak frequency within a frequency band.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `flim::Tuple{Real, Real}`: lower and upper frequency bounds
- `method::Symbol=:welch`: PSD method:
  - `:welch`: Welch's periodogram
  - `:fft`: fast Fourier transform
  - `:mt`: multi-tapered periodogram
  - `:stft`: short-time Fourier transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

- `pw::Matrix{Float64}`: peak power, shape `(channels, epochs)`
"""
function peak_pow(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    flim::Tuple{Real, Real},
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    demean::Bool = true,
)::Matrix{Float64}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    pw = peak_pow(
        @view(obj.data[ch, :, :]),
        fs = sr(obj),
        flim = flim,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        demean = demean,
    )

    return pw

end
