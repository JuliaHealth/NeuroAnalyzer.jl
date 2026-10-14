export psd
export mwpsd
export ghpsd

"""
    psd(s; <keyword arguments>)

Calculate Power Spectral Density for a 1-D signal vector.

Default method is Welch's periodogram.

# Arguments

- `s::Vector{Float64}`: signal vector
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `db::Bool=false`: normalize powers to dB
- `method::Symbol=:welch`: PSD estimation method:
    - `:welch`: Welch's periodogram (default)
    - `:fft`: plain FFT periodogram
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform; per-segment PSD averaged over segments
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian filter + Hilbert transform
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `wlen::Int64=fs`: window length in samples (default = 1 second); used by `:welch` and `:stft`; used by `:welch` and `:stft`
- `woverlap::Int64=min(min(round(Int64, wlen * 0.90), wlen - 1), wlen - 1)`: window overlap in samples; used by `:welch` and `:stft`
- `w::Bool=true`: if `true`, apply Hanning window (ignored by `:mt`, which uses Slepian tapers)
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of Morlet wavelet cycles; for a tuple the number of cycles is logarithmically (geometrically) spaced from `ncyc[1]` to `ncyc[2]` across frequencies
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`)
- `demean::Bool=true`: subtract DC component before estimating PSD

# Returns

Named tuple:

- `p::Vector{Float64}`: powers
- `f::Vector{Float64}`: frequencies

# Notes

- Setting `demean=true` reduces (but doesn't fully eliminate) DC contamination; the 0 Hz bin is still present in the output.
- `:welch`, `:fft`, `:mt` and `:stft` return a one-sided PSD (units²/Hz). `:mw` and `:gh` return band power (units²) at each frequency and are not directly comparable with the other methods.
"""
function psd(
    s::AbstractVector;
    fs::Int64,
    db::Bool = false,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = min(round(Int64, wlen * 0.9), wlen - 1),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool = true,
)::@NamedTuple{
    p::Vector{Float64},
    f::Vector{Float64},
}
    # validate
    _check_var(method, [:fft, :welch, :mt, :mw, :stft, :gh], "method")
    nt >= 1 || throw(ArgumentError("nt must be ≥ 1."))
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))
    # window parameters only matter for the segment-based methods
    if method in (:welch, :stft)
        wlen <= length(s) || throw(ArgumentError("wlen must be ≤ $(length(s))."))
        wlen >= 2 || throw(ArgumentError("wlen must be ≥ 2."))
        woverlap < wlen || throw(ArgumentError("woverlap must be < $(wlen)."))
        woverlap >= 0 || throw(ArgumentError("woverlap must be ≥ 0."))
    end

    n = length(s)

    if method === :mt
        demean && (s = remove_dc(s))
        # Slepian tapers already window the data; nw so that 2nw - 1 = nt
        pg = DSP.mt_pgram(s; fs = fs, nw = (nt + 1) / 2, ntapers = nt)
        f = Vector(freq(pg))
        p = power(pg)[1:length(f)]

    elseif method === :stft
        demean && (s = remove_dc(s))
        win = w ? DSP.hanning : nothing
        # per-segment PSD (properly scaled), averaged over segments
        sp = DSP.spectrogram(s, wlen, woverlap; fs = fs, window = win)
        f = Vector(freq(sp))
        p = vec(mean(power(sp); dims = 2))

    elseif method === :welch
        demean && (s = remove_dc(s))
        win = w ? DSP.hanning : nothing
        pg = DSP.welch_pgram(s, wlen, woverlap; fs = fs, window = win)
        f = Vector(freq(pg))
        p = power(pg)

    elseif method === :fft
        demean && (s = remove_dc(s))
        win = w ? DSP.hanning : nothing
        pg = DSP.periodogram(s; fs = fs, window = win)
        f = Vector(freq(pg))
        p = power(pg)[1:length(f)]

    elseif method === :mw
        r = mwpsd(s; db = false, fs = fs, ncyc = ncyc, w = w, demean = demean)
        p, f = r.p, r.f

    elseif method === :gh
        r = ghpsd(s; fs = fs, db = false, gw = gw, w = w, demean = demean)
        p, f = r.p, r.f
    end

    f = round.(f, digits = 3)

    db && (p = pow2db.(p))

    return (; p, f)
end

"""
    psd(s; <keyword arguments>)

Calculate Power Spectral Density for each channel of a matrix.

Default method is Welch's periodogram.

# Arguments

- `s::AbstractMatrix`: signal matrix (channels, samples)
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `db::Bool=false`: normalize powers to dB
- `method::Symbol=:welch`: PSD estimation method:
    - `:welch`: Welch's periodogram (default)
    - `:fft`: plain FFT periodogram
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform averaged over segments
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian filter + Hilbert transform
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `wlen::Int64=fs`: window length in samples (default = 1 second); used by `:welch` and `:stft`
- `woverlap::Int64=min(min(round(Int64, wlen * 0.90), wlen - 1), wlen - 1)`: window overlap in samples
- `w::Bool=true`: if `true`, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of Morlet wavelet cycles; for a tuple the number of cycles is logarithmically (geometrically) spaced from `ncyc[1]` to `ncyc[2]` across frequencies
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`)
- `demean::Bool=true`: subtract DC component before estimating PSD

# Returns

Named tuple:

- `p::Matrix{Float64}`: powers, shape (channels, frequencies)
- `f::Vector{Float64}`: frequencies

# Notes

- Setting `demean=true` reduces (but doesn't fully eliminate) DC contamination; the 0 Hz bin is still present in the output.
- `:welch`, `:fft`, `:mt` and `:stft` return a one-sided PSD (units²/Hz). `:mw` and `:gh` return band power (units²) at each frequency and are not directly comparable with the other methods.
"""
function psd(
    s::AbstractMatrix;
    fs::Int64,
    db::Bool = false,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = min(round(Int64, wlen * 0.90), wlen - 1),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool = true,
)::@NamedTuple{
    p::Matrix{Float64},
    f::Vector{Float64},
}
    # validate
    size(s, 1) >= 1 || throw(ArgumentError("s must contain at least one channel."))

    _psd(x) = psd(
        x;
        fs = fs,
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

    # first channel also determines the output data
    psd_data = _psd(@view(s[1, :]))
    f = psd_data.f
    # pre-allocate output
    p = zeros(size(s, 1), length(f))
    p[1, :] = psd_data.p

    Threads.@threads :static for ch_idx = 2:size(s, 1)
        p[ch_idx, :] = _psd(@view(s[ch_idx, :])).p
    end
 
    return (; p, f)
end

"""
    psd(s; <keyword arguments>)

Calculate Power Spectral Density for a 3-D signal array.

Default method is Welch's periodogram.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `db::Bool=false`: normalize powers to dB
- `method::Symbol=:welch`: PSD estimation method:
    - `:welch`: Welch's periodogram (default)
    - `:fft`: plain FFT periodogram
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform averaged over segments
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian filter + Hilbert transform
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `wlen::Int64=fs`: window length in samples (default = 1 second); used by `:welch` and `:stft`
- `woverlap::Int64=min(min(round(Int64, wlen * 0.90), wlen - 1), wlen - 1)`: window overlap in samples
- `w::Bool=true`: if `true`, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of Morlet wavelet cycles; for a tuple the number of cycles is logarithmically (geometrically) spaced from `ncyc[1]` to `ncyc[2]` across frequencies
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`)
- `demean::Bool=true`: subtract DC component before estimating PSD

# Returns

Named tuple:

- `p::Array{Float64, 3}`: powers, shape (channels, frequencies, epochs)
- `f::Vector{Float64}`: frequencies

# Notes

- Setting `demean=true` reduces (but doesn't fully eliminate) DC contamination; the 0 Hz bin is still present in the output.
- `:welch`, `:fft`, `:mt` and `:stft` return a one-sided PSD (units²/Hz). `:mw` and `:gh` return band power (units²) at each frequency and are not directly comparable with the other methods.
"""
function psd(
    s::AbstractArray;
    fs::Int64,
    db::Bool = false,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = min(round(Int64, wlen * 0.90), wlen - 1),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool = true,
)::@NamedTuple{
    p::Array{Float64, 3},
    f::Vector{Float64},
}
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    _psd(x) = psd(
        x;
        fs = fs,
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

    # first channel/epoch also determines the output data
    psd_data = _psd(@view(s[1, :, 1]))
    f = psd_data.f
    # pre-allocate output
    p = zeros(ch_n, length(f), ep_n)
    p[1, :, 1] = psd_data.p

    # calculate over channels and epochs (skipping the one already done)
    idxs = vec(CartesianIndices((ch_n, ep_n)))[2:end]
    Threads.@threads :static for idx in idxs
        ch_idx, ep_idx = idx[1], idx[2]
        p[ch_idx, :, ep_idx] = _psd(@view(s[ch_idx, :, ep_idx])).p
    end
 
    return (; p, f)
end

"""
    psd(obj; <keyword arguments>)

Calculate Power Spectral Density for a NEURO object.

Default method is Welch's periodogram.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `db::Bool=false`: normalize powers to dB
- `method::Symbol=:welch`: PSD estimation method:
    - `:welch`: Welch's periodogram (default)
    - `:fft`: plain FFT periodogram
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform averaged over segments
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian filter + Hilbert transform
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `wlen::Int64=sr(obj)`: window length in samples (default = 1 second); used by `:welch` and `:stft`
- `woverlap::Int64=min(min(round(Int64, wlen * 0.90), wlen - 1), wlen - 1)`: window overlap in samples
- `w::Bool=true`: if `true`, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`)
- `flim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds for output trimming
- `demean::Bool=true`: subtract DC component before estimating PSD

# Returns

Named tuple:

- `p::Array{Float64, 3}`: powers, shape (channels, frequencies, epochs)
- `f::Vector{Float64}`: frequencies (trimmed to `flim`)

# Notes

- Setting `demean=true` reduces (but doesn't fully eliminate) DC contamination; the 0 Hz bin is still present in the output.
- `:welch`, `:fft`, `:mt` and `:stft` return a one-sided PSD (units²/Hz). `:mw` and `:gh` return band power (units²) at each frequency and are not directly comparable with the other methods.
"""
function psd(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    db::Bool = false,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = min(round(Int64, wlen * 0.90), wlen - 1),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    flim::Tuple{Real, Real} = (0, sr(obj) / 2),
    demean::Bool = true,
)::@NamedTuple{
    p::Array{Float64, 3},
    f::Vector{Float64},
}

    # validate
    _check_tuple(flim, (0, sr(obj) / 2), "flim")

    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    isempty(ch) && throw(ArgumentError("No channels selected."))

    psd_data = psd(
        @view(obj.data[ch, :, :]);
        fs = sr(obj),
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

    # trim output to requested frequency band
    f1 = vsearch(flim[1], psd_data.f)
    f2 = vsearch(flim[2], psd_data.f)

    return (p = psd_data.p[:, f1:f2, :], f = psd_data.f[f1:f2])
end

"""
    mwpsd(s; <keyword arguments>)

Calculate power spectrum using Morlet wavelet convolution for a 1-D signal vector.

# Arguments

- `s::AbstractVector`: signal vector
- `pad::Int64=0`: number of zeros to append
- `db::Bool=false`: normalize powers to dB
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of Morlet wavelet cycles; for a tuple the number of cycles is logarithmically (geometrically) spaced from `ncyc[1]` to `ncyc[2]` across frequencies
- `w::Bool=true`: if `true`, apply Hanning window
- `demean::Bool=true`: subtract DC component before estimating PSD

# Returns

Named tuple:

- `p::Vector{Float64}`: powers
- `f::Vector{Float64}`: frequencies

# Notes

- The 0 Hz bin is omitted: a Morlet wavelet is undefined at 0 Hz.
- Each kernel is normalized to unit gain at its centre frequency and spans ±3 SD of its Gaussian (at least ±1 s), so low-frequency kernels are not truncated. Kernels longer than the signal give edge-dominated estimates.
- Power at each frequency is the median over time of the squared wavelet magnitude (robust to transients; ≈ 0.69 × the mean for Gaussian noise). Values are band power, not PSD (units²/Hz); the Hanning window reduces power by ≈ 0.375.
"""
function mwpsd(
    s::AbstractVector;
    pad::Int64 = 0,
    db::Bool = false,
    fs::Int64,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    w::Bool = true,
    demean::Bool = true,
)::@NamedTuple{
    p::Vector{Float64},
    f::Vector{Float64},
}

    # validate
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))
    pad >= 0 || throw(ArgumentError("pad must be ≥ 0."))

    demean && (s = remove_dc(s))
    # window the data first, then pad, so the taper covers only real samples
    s_w = w ? s .* hanning(length(s)) : collect(Float64, s)
    pad > 0 && (s_w = pad0(s_w, pad))

    w = w ? hanning(length(s)) : ones(length(s))

    flim = (0, fs / 2)
    nfrq = _tlength(flim)
    # drop 0 Hz
    f = linspace(flim[1], flim[2], nfrq)[2:end]
    nfrq = length(f)

    if ncyc isa Int64
        ncyc >= 1 || throw(ArgumentError("ncyc must be ≥ 1."))
        ncyc_v = fill(ncyc, nfrq)
    else
        ncyc[1] >= 1 || throw(ArgumentError("ncyc[1] must be ≥ 1."))
        ncyc[2] >= 1 || throw(ArgumentError("ncyc[2] must be ≥ 1."))
        ncyc_v =
            nfrq == 1 ? [ncyc[1]] : round.(Int64, logspace(ncyc[1], ncyc[2], nfrq))
    end

    p = zeros(nfrq)
    for frq_idx = 1:nfrq
        # Gaussian SD of the wavelet; kernel half-length covers ±3 SD
        sd = ncyc_v[frq_idx] / (2 * pi * f[frq_idx])
        t = max(1.0, 3 * sd)
        kernel = generate_morlet(fs, f[frq_idx], t; ncyc = ncyc_v[frq_idx], complex = true)
        w_conv = fconv(s_w; kernel = kernel, norm = true)
        p[frq_idx] = median(abs2.(w_conv))
    end
 
    db && (p = pow2db.(p))
 
    return (; p, f)
end

"""
    ghpsd(s; <keyword arguments>)

Calculate power spectrum using Gaussian filter and Hilbert transform for a 1-D signal vector.

# Arguments

- `s::AbstractVector`: signal vector
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `db::Bool=false`: normalize powers to dB
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`)
- `w::Bool=true`: if `true`, apply Hanning window
- `demean::Bool=true`: subtract DC component before estimating PSD

# Returns

Named tuple:

- `p::Vector{Float64}`: powers
- `f::Vector{Float64}`: frequencies

# Notes

- `filter_g` returns the amplitude envelope of the band-limited (analytic) signal, so power is its square averaged over time.
- Values are band power, not PSD (units²/Hz); the Hanning window reduces power by ≈ 0.375.
"""
function ghpsd(
    s::AbstractVector;
    fs::Int64,
    db::Bool = false,
    gw::Real = 5,
    w::Bool = true,
    demean::Bool = true,
)::@NamedTuple{
    p::Vector{Float64},
    f::Vector{Float64}
}

    # validate
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))

    demean && (s = remove_dc(s))
    flim = (0, fs / 2)
    nfrq = _tlength(flim)
    f = linspace(flim[1], flim[2], nfrq)

    s_w = w ? s .* hanning(length(s)) : collect(Float64, s)

    # mean instantaneous power per frequency, without storing the full (frequencies × samples) matrix
    p = zeros(length(f))
    for frq_idx in eachindex(f)
        env = filter_g(s_w; fs = fs, f = f[frq_idx], gw = gw)
        p[frq_idx] = mean(abs2, env)
    end
 
    db && (p = pow2db.(p))
 
    return (; p, f)
end
