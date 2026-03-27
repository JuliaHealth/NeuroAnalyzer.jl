export sef

"""
    sef(s; <keyword arguments>)

Calculate spectral edge frequency (SEF) for a 1-D signal vector.

SEF is the frequency below which x percent of the total power of a given signal are located; typically, x is in the range 75 to 95.

# Arguments

- `s::AbstractVector`: signal vector
- `x::Float64=0.95`: threshold
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `f::Tuple{Real, Real}=(0, fs / 2)`: lower and upper frequency bounds, default is total power
- `method::Symbol=:welch`: PSD method:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `wlen::Int64=fs`: window length in samples (default = 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if `true`, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `demean::Bool=true`: subtract DC component before estimating PSD

# Returns

- `Float64`: spectral edge frequency
"""
function sef(
        s::AbstractVector;
        x::Float64 = 0.95,
        fs::Int64,
        f::Tuple{Real, Real} = (0, fs / 2),
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = fs,
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        demean::Bool = true,
    )::Float64

    # validate
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))
    _check_tuple(f, (0, fs / 2), "f")

    pw, pf = psd(
        s;
        fs = fs,
        db = false,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        demean = demean,
    )

    f1_idx = vsearch(f[1], pf)
    f2_idx = vsearch(f[2], pf)

    # dx: frequency resolution
    dx = pf[2] - pf[1]

    pw = pw[f1_idx:f2_idx]
    pf = pf[f1_idx:f2_idx]

    tp = simpson(pw; dx = dx)
    tp_threshold = tp * x

    sef_frq = nothing
    for idx in eachindex(pf)
        if sum(pw[1:idx]) >= tp_threshold
            sef_frq = pf[idx]
            break
        end
    end

    return sef_frq
end

"""
    sef(s; <keyword arguments>)

Calculate spectral edge frequency (SEF) for a 3-D signal array.

SEF is the frequency below which x percent of the total power of a given signal are located; typically, x is in the range 75 to 95.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `x::Float64=0.95`: threshold
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `f::Tuple{Real, Real}=(0, fs / 2)`: lower and upper frequency bounds, default is total power
- `method::Symbol=:welch`: PSD method:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `wlen::Int64=fs`: window length in samples (default = 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if `true`, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `demean::Bool=true`: subtract DC component before estimating PSD

# Returns

- `Matrix{Float64}`: spectral edge frequency
"""
function sef(
        s::AbstractArray;
        x::Float64 = 0.95,
        fs::Int64,
        f::Tuple{Real, Real} = (0, fs / 2),
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = fs,
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        demean::Bool = true,
    )::Matrix{Float64}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    sef_frq = zeros(ch_n, ep_n)

    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        sef_frq[ch_idx, ep_idx] = sef(
            @view(s[ch_idx, :, ep_idx]),
            x = x,
            fs = fs,
            f = f,
            method = method,
            nt = nt,
            wlen = wlen,
            woverlap = woverlap,
            w = w,
            ncyc = ncyc,
            demean = demean,
        )
    end

    return sef_frq
end

"""
    sef(obj; <keyword arguments>)

Calculate spectral edge frequency (SEF) for a NEURO object.

SEF is the frequency below which x percent of the total power of a given signal are located; typically, x is in the range 75 to 95.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `x::Float64=0.95`: threshold
- `f::Tuple{Real, Real}=(0, sr(obj) / 2)`: lower and upper frequency bounds, default is total power
- `method::Symbol=:welch`: PSD method:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if `true`, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `demean::Bool=true`: subtract DC component before estimating PSD

# Returns

- `Matrix{Float64}`: spectral edge frequency
"""
function sef(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex},
        x::Float64 = 0.95,
        f::Tuple{Real, Real} = (0, sr(obj) / 2),
        method::Symbol = :welch,
        nt::Int64 = 7,
        wlen::Int64 = sr(obj),
        woverlap::Int64 = round(Int64, wlen * 0.9),
        w::Bool = true,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        demean::Bool = true,
    )::Matrix{Float64}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ? get_channel(obj; ch = ch, exclude = "bad") :
                       get_channel(obj; ch = ch, exclude = "")

    return sef(
        @view(obj.data[ch, :, :]);
        x = x,
        fs = sr(obj),
        f = f,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        demean = demean,
    )
end
