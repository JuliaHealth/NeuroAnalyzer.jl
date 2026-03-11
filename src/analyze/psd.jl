export psd
export mwpsd
export ghpsd

"""
    psd(s; <keyword arguments>)

Calculate Power Spectral Density. Default method is Welch's periodogram.

# Arguments

- `s::Vector{Float64}`: signal vector
- `fs::Int64`: sampling rate
- `db::Bool=false`: normalize powers to dB
- `method::Symbol=:welch`: PSD method:
    - `:welch`: Welch's periodogram (default)
    - `:fft`: plain FFT periodogram
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform averaged over segments
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian filter + Hilbert transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length in samples, default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `gw::Real=5`: Gaussian width in Hz
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

Named tuple:

- `p::Vector{Float64}`: powers
- `f::Vector{Float64}`: frequencies

# Notes

Setting `demean=true` reduces (but doesn't fully eliminate) DC contamination; the 0 Hz bin is still present in the output.
"""
function psd(
    s::AbstractVector;
    fs::Int64,
    db::Bool = false,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool = true,
)::@NamedTuple{p::Vector{Float64}, f::Vector{Float64}}

    _check_var(method, [:fft, :welch, :mt, :mw, :stft, :gh], "method")
    @assert nt >= 1 "nt must be ≥ 1."
    @assert fs >= 1 "fs must be ≥ 1."
    @assert wlen <= length(s) "wlen must be ≤ $(length(s))."
    @assert wlen >= 2 "wlen must be ≥ 2."
    @assert woverlap < wlen "woverlap must be < $(wlen)."
    @assert woverlap >= 0 "woverlap must be ≥ 0."

    n = length(s)

    if method === :mt
        
        demean && (s = remove_dc(s))
        win = w ? hanning(n) : ones(n)
        pg = mt_pgram(s .* win, fs = fs, nw = ((nt + 1) ÷ 2), ntapers = nt)
        f = Vector(freq(pg))
        p = power(pg)[1:length(f)]
        db && (p = pow2db.(p))

    elseif method === :stft

        demean && (s = remove_dc(s))
        win = w ? DSP.hanning : nothing
        sg = abs.(DSP.stft(s, wlen, woverlap, fs = fs, window = win))
        # average STFT segments over time; build a matching frequency vector.
        p = dropdims(mean(sg, dims = 2), dims = 2)
        f = linspace(0, fs / 2, length(p))
        db && (p = pow2db.(p))

    elseif method === :welch

        demean && (s = remove_dc(s))
        win = w ? DSP.hanning : nothing
        pg  = DSP.welch_pgram(s, wlen, woverlap, fs = fs, window = win)
        f   = Vector(freq(pg))
        p   = power(pg)
        db && (p = pow2db.(p))

    elseif method === :fft

        win = w ? DSP.hanning : nothing
        pg  = DSP.periodogram(s, fs = fs, window = win)
        f   = Vector(freq(pg))
        p   = power(pg)[1:length(f)]
        db && (p = pow2db.(p))

    elseif method === :mw

        mwpsd_data = mwpsd(s, db = db, fs = fs, ncyc = ncyc, w = w, demean = demean)
        p = mwpsd_data.p
        f = mwpsd_data.f

    elseif method === :gh

        ghpsd_data = ghpsd(s, fs = fs, db = db, gw = gw, w = w, demean = demean)
        p = ghpsd_data.p
        f = ghpsd_data.f

    end

    return (p = p, f = f)

end

"""
    psd(s; <keyword arguments>)

Calculate Power Spectral Density for each channel of a matrix. Default method is Welch's periodogram.

# Arguments

- `s::AbstractMatrix`: signal matrix (channels, samples)
- `fs::Int64`: sampling rate
- `db::Bool=false`: normalize powers to dB
- `method::Symbol=:welch`: PSD method:
    - `:welch`: Welch's periodogram (default)
    - `:fft`: plain FFT periodogram
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform averaged over segments
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian filter + Hilbert transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length in samples, default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

Named tuple:

- `p::Matrix{Float64}`: powers, shape `(channels, frequencies)`
- `f::Vector{Float64}`: frequencies
"""
function psd(
    s::AbstractMatrix;
    fs::Int64,
    db::Bool = false,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool=true,
)::@NamedTuple{p::Matrix{Float64}, f::Vector{Float64}}

    # pilot call to determine output frequency vector length
    f = psd(
        @view(s[1, :]),
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
    ).f

    p = zeros(size(s, 1), length(f))

    @inbounds Threads.@threads :dynamic for ch_idx in axes(s, 1)
        p[ch_idx, :] = psd(
            @view(s[ch_idx, :]),
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
        ).p
    end

    return (p = p, f = f)

end

"""
    psd(s; <keyword arguments>)

Calculate Power Spectral Density for a 3-D signal array. Default method is Welch's periodogram.

# Arguments

- `s::AbstractArray`: signal array (channels, samples, epochs)
- `fs::Int64`: sampling rate
- `db::Bool=false`: normalize powers to dB
- `method::Symbol=:welch`: PSD method:
    - `:welch`: Welch's periodogram (default)
    - `:fft`: plain FFT periodogram
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform averaged over segments
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian filter + Hilbert transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length in samples, default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `gw::Real=5`: Gaussian width in Hz
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

Named tuple:

- `p::Array{Float64, 3}`: powers
- `f::Vector{Float64}`: frequencies
"""
function psd(
    s::AbstractArray;
    fs::Int64,
    db::Bool = false,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool = true
)::@NamedTuple{p::Array{Float64, 3}, f::Vector{Float64}}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pilot call to determine output frequency vector length
    f = psd(
        @view(s[1, :, 1]),
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
    ).f

    # pre-allocate output
    p = zeros(ch_n, length(f), ep_n)

    # calculate over channels and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        p[ch_idx, :, ep_idx] = psd(
            @view(s[ch_idx, :, ep_idx]),
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
    ).p
    end

    return (p = p, f = f)

end

"""
    psd(obj; <keyword arguments>)

Calculate Power Spectral Density. Default method is Welch's periodogram.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `db::Bool=false`: normalize powers to dB
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
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `gw::Real=5`: Gaussian width in Hz
- `flim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds for output trimming
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

Named tuple:

- `p::Array{Float64, 3}`: powers, shape `(channels, frequencies, epochs)`
- `f::Vector{Float64}`: frequencies (trimmed to `flim`)
"""
function psd(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    db::Bool = false,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    flim::Tuple{Real, Real} = (0, sr(obj) / 2),
    demean::Bool = true,
)::@NamedTuple{p::Array{Float64, 3}, f::Vector{Float64}}

    _check_tuple(flim, (0, sr(obj) / 2), "flim")

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    _log_off()
    psd_data = psd(@view(obj.data[ch, :, :]),
        fs = sr(obj),
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
    _log_on()

    # trim output to requested frequency band
    f1 = vsearch(flim[1], psd_data.f)
    f2 = vsearch(flim[2], psd_data.f)

    return (p = psd_data.p[:, f1:f2, :], f = psd_data.f[f1:f2])

end

"""
    mwpsd(s; <keyword arguments>)

Calculate power spectrum using Morlet wavelet convolution.

# Arguments

- `s::AbstractVector`: signal vector
- `pad::Int64=0`: number of zeros to append
- `db::Bool=true`: normalize powers to dB
- `fs::Int64`: sampling rate
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `w::Bool=true`: if true, apply Hanning window
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

Named tuple:

- `p::Vector{Float64}`: powers
- `f::Vector{Float64}`: frequencies
"""
function mwpsd(
    s::AbstractVector;
    pad::Int64 = 0,
    db::Bool = true,
    fs::Int64,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    w::Bool = true,
    demean::Bool = true,
)::@NamedTuple{p::Vector{Float64}, f::Vector{Float64}}

    @assert fs >= 1 "fs must be ≥ 1."
    @assert pad >= 0 "pad must be ≥ 0."

    demean && (s = remove_dc(s))
    pad > 0 && (s = pad0(s, pad))

    w = w ? hanning(length(s)) : ones(length(s))

    flim = (0, fs / 2)
    nfrq = _tlength(flim)
    f = linspace(flim[1], flim[2], nfrq)

    if ncyc isa Int64
        @assert ncyc >= 1 "ncyc must be >= 1"
        ncyc = repeat([ncyc], nfrq)
    else
        @assert ncyc[1] >= 1 "ncyc[1] must be >= 1"
        @assert ncyc[2] >= 1 "ncyc[2] must be >= 1"
        ncyc = round.(Int64, logspace(ncyc[1], ncyc[2], nfrq))
    end

    p = zeros(length(f))
    @inbounds for frq_idx in 1:nfrq
        kernel = generate_morlet(fs, f[frq_idx], 1, ncyc = ncyc[frq_idx], complex = true)
        w_conv = fconv(s .* win, kernel = kernel, norm = false)
        p[frq_idx] = median(abs2.(w_conv))
    end

    db && (p = pow2db.(p))

    return (p = p, f = f)

end

"""
    ghpsd(s; <keyword arguments>)

Calculate power spectrum using Gaussian filter and Hilbert transform.

# Arguments

- `s::AbstractVector`: signal vector
- `fs::Int64`: sampling rate
- `db::Bool=true`: normalize powers to dB
- `gw::Real=5`: Gaussian width in Hz
- `w::Bool=true`: if true, apply Hanning window
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

Named tuple:

- `p::Vector{Float64}`: powers
- `f::Vector{Float64}`: frequencies
"""
function ghpsd(
    s::AbstractVector;
    fs::Int64,
    db::Bool = true,
    gw::Real = 5,
    w::Bool = true,
    demean::Bool = true,
)::@NamedTuple{p::Vector{Float64}, f::Vector{Float64}}

    @assert fs >= 1 "fs must be ≥ 1."

    demean && (s = remove_dc(s))
    flim = (0, fs / 2)
    nfrq = _tlength(flim)
    f = linspace(flim[1], flim[2], nfrq)

    w = w ? hanning(length(s)) : ones(length(s))

    # pre-allocate output
    p = zeros(length(f), length(s))

    @inbounds for frq_idx in eachindex(f)
        s_filt = filter_g(s .* win, fs = fs, f = f[frq_idx], gw = gw)
        p[frq_idx, :] = abs2.(DSP.hilbert(s_filt))
    end

    p[p .== -Inf] .= minimum(p[p .!= -Inf])
    p[p .== +Inf] .= maximum(p[p .!= +Inf])

    # average instantaneous power across time at each frequency.
    p = dropdims(mean(p, dims = 2), dims = 2)

    db && (p = pow2db.(p))

    return (p = p, f = f)

end
