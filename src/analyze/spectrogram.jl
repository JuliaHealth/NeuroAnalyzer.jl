export spectrogram
export mwspectrogram
export ghtspectrogram
export cwtspectrogram
export hhtspectrogram

"""
    spectrogram(s; <keyword arguments>)

Calculate spectrogram using STFT or multi-tapered method.

# Arguments

- `s::AbstractVector`: signal vector
- `fs::Int64`: sampling frequency
- `db::Bool=true`: normalize powers to dB
- `method::Symbol=:stft`: PSD method:
- `:stft`: short-time Fourier transform
- `:mt`: multi-tapered periodogram
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length in samples
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple:

- `p::Matrix{Float64}`: powers
- `f::Vector{Float64}`: frequencies
- `t::Vector{Float64}`: time points
"""
function spectrogram(
    s::AbstractVector;
    fs::Int64,
    db::Bool = true,
    method::Symbol = :stft,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
)::@NamedTuple{p::Matrix{Float64}, f::Vector{Float64}, t::Vector{Float64}}

    _check_var(method, [:stft, :mt], "method")
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))
    wlen <= length(s) || throw(ArgumentError("wlen must be ≤ $(length(s))."))
    wlen >= 1 || throw(ArgumentError("wlen must be ≥ 1."))
    woverlap < wlen || throw(ArgumentError("woverlap must be < $(wlen)."))
    woverlap >= 0 || throw(ArgumentError("woverlap must be ≥ 0."))

    if method === :stft

        win = w ? hanning : nothing
        if length(s) < fs * 1.1
            wlen = div(wlen, 2)
            woverlap = div(woverlap, 2)
        end
        pg = DSP.spectrogram(s, wlen, woverlap, fs = fs, window = win)

    elseif method === :mt

        win = w ? hanning(length(s)) : ones(length(s))
        pg  = DSP.mt_spectrogram(s .* win, fs = fs, nw = ((nt + 1) ÷ 2), ntapers = nt)

    end

    p = pg.power

    p[p .== -Inf] .= minimum(p[p .!= -Inf])
    p[p .== +Inf] .= maximum(p[p .!= +Inf])
    db && (p = pow2db.(p))

    t = 0:(1 / fs):(length(s) / fs)
    t = linspace(t[1], t[end - 1], size(p, 2))
    f = linspace(0, fs / 2, size(p, 1))

    return (p = p, f = f, t = t)

end

"""
    spectrogram(s; <keyword arguments>)

Calculate spectrogram for each channel of a matrix.

# Arguments

- `s::AbstractMatrix`: signal matrix (channels, samples)
- `fs::Int64`: sampling frequency
- `db::Bool=true`: normalize powers to dB
- `method::Symbol=:stft`: PSD method:
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length, default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple:

- `p::Array{Float64, 3}`: powers, shape `(freq, time, channels)`
- `f::Vector{Float64}`: frequencies
- `t::Vector{Float64}`: time points
"""
function spectrogram(
    s::AbstractMatrix;
    fs::Int64,
    db::Bool = true,
    method::Symbol = :stft,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
)::@NamedTuple{p::Array{Float64, 3}, f::Vector{Float64}, t::Vector{Float64}}

    # pilot call to determine output frequency vector length
    spec_data = NeuroAnalyzer.spectrogram(
        @view(s[1, :]),
        fs = fs,
        db = db,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w
    )
    f = spec_data.f
    t = spec_data.t

    # pre-allocate output
    @inbounds Threads.@threads :dynamic for ch_idx in axes(s, 1)
        p[:, :, ch_idx] = NeuroAnalyzer.spectrogram(
            @view(s[ch_idx, :]),
            fs = fs,
            db = db,
            method = method,
            nt = nt,
            wlen = wlen,
            woverlap = woverlap,
            w = w
        ).p
    end

    return (p = p, f = f, t = t)

end

"""
    spectrogram(obj; <keyword arguments>)

Calculate spectrogram.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pad::Int64=0`: number of zeros to append
- `method::Symbol=:stft`: spectrogram method:
- `:stft`: short-time Fourier transform
- `:mt`: multi-tapered periodogram
- `:mw`: Morlet wavelet convolution
- `:gh`: Gaussian and Hilbert transform
- `:cwt`: continuous wavelet transformation
- `:hht`: Hilbert-Huang transform
- `db::Bool=true`: normalize powers to dB
- `nt::Int64=7`: number of Slepian tapers
- `gw::Real=10`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple:

- `p::Array{Float64, 4}`: powers (magnitudes for `:cwt`), shape `(freq, time, channels, epochs)`
- `f::Vector{Float64}`: frequencies
- `t::Vector{Float64}`: time points
"""
function spectrogram(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    pad::Int64 = 0,
    method::Symbol = :stft,
    db::Bool = true,
    nt::Int64 = 7,
    gw::Real = 10,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    wt::T = wavelet(Morlet(2π), β = 2),
    wlen::Int64 = sr(obj),
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
)::@NamedTuple{
    p::Array{Float64, 4},
    f::Vector{Float64},
    t::Vector{Float64}
} where {T <: CWT}

    _check_var(method, [:stft, :mt, :mw, :gh, :cwt, :hht], "method")

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)
    # sampling rate
    fs = sr(obj)
    # time points
    t = obj.epoch_time

    # pilot call to determine output dimensions
    if method === :stft
        spec_data = NeuroAnalyzer.spectrogram(
            @view(obj.data[1, :, 1]),
            fs = fs,
            db = db,
            method = :stft,
            wlen = wlen,
            woverlap = woverlap,
            w = w
        )
        f = spec_data.f
        p_tmp = spec_data.p
    elseif method === :mt
        spec_data = NeuroAnalyzer.spectrogram(
            @view(obj.data[1, :, 1]),
            fs = fs,
            db = db,
            method = :mt,
            nt = nt,
            wlen = wlen,
            woverlap = woverlap,
            w = w
        )
        f = spec_data.f
        p_tmp = spec_data.p
    elseif method === :mw
        spec_data = NeuroAnalyzer.mwspectrogram(
            @view(obj.data[1, :, 1]),
            pad = pad,
            fs = fs,
            db = db,
            ncyc = ncyc,
            w = w,
        )
        f = spec_data.f
        p_tmp = spec_data.p
    elseif method === :gh
        spec_data = NeuroAnalyzer.ghtspectrogram(
            @view(obj.data[1, :, 1]),
            fs = fs,
            db = db,
            gw = gw,
            w = w
        )
        f = spec_data.f
        p_tmp = spec_data.p
    elseif method === :cwt
        _log_off()
        spec_data = NeuroAnalyzer.cwtspectrogram(@view(obj.data[1, :, 1]), fs = fs, wt = wt)
        _log_on()
        f = spec_data.f
        # cwtspectrogram returns field .m not .p
        p_tmp = spec_data.m
    elseif method === :hht
        spec_data = NeuroAnalyzer.hhtspectrogram(@view(obj.data[1, :, 1]), t, fs = fs, db = db)
        f = spec_data.f
        p_tmp = spec_data.p
    end

    # pre-allocate outputs
    p = zeros(size(p_tmp, 1), size(p_tmp, 2), ch_n, ep_n)

    # initialize progress bar
    progbar = Progress(ep_n * ch_n, dt = 1, barlen = 20, color = :white, enabled = progress_bar)

    # calculate over channels and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_local, ep_idx = idx[1], idx[2]
        ch_idx = ch[ch_local]   # resolve local index to actual channel number

        if method === :stft
            p[:, :, ch_local, ep_idx] = NeuroAnalyzer.spectrogram(
                @view(obj.data[ch_idx, :, ep_idx]),
                fs = fs,
                db = db,
                method = :stft,
                wlen = wlen,
                woverlap = woverlap,
                w = w
            ).p
        elseif method === :mt
            p[:, :, ch_local, ep_idx] = NeuroAnalyzer.spectrogram(
                @view(obj.data[ch_idx, :, ep_idx]),
                fs = fs,
                db = db,
                method = :mt,
                nt = nt,
                wlen = wlen,
                woverlap = woverlap,
                w = w
            ).p
        elseif method === :mw
            p[:, :, ch_local, ep_idx] = NeuroAnalyzer.mwspectrogram(
                @view(obj.data[ch_idx, :, ep_idx]),
                pad = pad,
                fs = fs,
                db = db,
                ncyc = ncyc,
                w = w
            ).p
        elseif method === :gh
            p[:, :, ch_local, ep_idx] = NeuroAnalyzer.ghtspectrogram(
                @view(obj.data[ch_idx, :, ep_idx]),
                fs = fs,
                db = db,
                gw = gw,
                w = w
            ).p
        elseif method === :cwt
            _log_off()
            # cwtspectrogram returns field .m (magnitudes) rather than .p.
            p[:, :, ch_local, ep_idx] = NeuroAnalyzer.cwtspectrogram(
                @view(obj.data[ch_idx, :, ep_idx]),
                fs = fs,
                wt = wt
            ).m
            _log_on()
        elseif method === :hht
            # hhtspectrogram requires time points for EMD
            p[:, :, ch_local, ep_idx] = NeuroAnalyzer.hhtspectrogram(@view(obj.data[ch_idx, :, ep_idx]),
                                                                     t,
                                                                     fs = fs,
                                                                     db = db).p
        end

        progress_bar && next!(progbar)
    end

    f = round.(f, digits = 2)
    t = round.(t, digits = 3)
    t .+= obj.epoch_time[1]

    return (p = p, f = f, t = t)

end

"""
    mwspectrogram(s; <keyword arguments>)

Calculate spectrogram using Morlet wavelet convolution.

# Arguments

- `s::AbstractVector`: signal vector
- `pad::Int64=0`: number of zeros to append
- `db::Bool=true`: normalize powers to dB
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple:

- `cs::Matrix{ComplexF64}`: convolved analytic signal
- `p::Matrix{Float64}`: powers
- `ph::Matrix{Float64}`: phases
- `f::Vector{Float64}`: frequencies
- `t::Vector{Float64}`: time points
"""
function mwspectrogram(
    s::AbstractVector;
    pad::Int64 = 0,
    db::Bool = true,
    fs::Int64,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    w::Bool = true,
)::@NamedTuple{
    cs::Matrix{ComplexF64},
    p::Matrix{Float64},
    ph::Matrix{Float64},
    f::Vector{Float64},
    t::Vector{Float64}
}

    !(fs >= 1) && throw(ArgumentError("fs must be > 1."))

    pad > 0 && (s = pad0(s, pad))

    win = w ? hanning(length(s)) : ones(length(s))

    if ncyc isa Int64
        !(ncyc >= 1) && throw(ArgumentError("ncyc must be >= 1."))
    else
        !(ncyc[1] >= 1) && throw(ArgumentError("ncyc[1] must be >= 1."))
        !(ncyc[2] >= 1) && throw(ArgumentError("ncyc[2] must be >= 1."))
    end

    flim = (0, fs / 2)
    nfrq = _tlength(flim)
    f = linspace(flim[1], flim[2], nfrq)

    cs = zeros(ComplexF64, length(f), length(s))
    p = zeros(length(f), length(s))
    ph = zeros(length(f), length(s))

    if ncyc isa Int64
        ncyc = repeat([ncyc], nfrq)
    else
        ncyc = round.(Int64, logspace(ncyc[1], ncyc[2], nfrq))
    end

    @inbounds for frq_idx in 1:nfrq
        kernel = generate_morlet(fs, f[frq_idx], 1, ncyc = ncyc[frq_idx], complex = true)
        cs[frq_idx, :] = fconv(s .* win, kernel = kernel, norm = true)
        @views p[frq_idx, :]  = abs2.(cs[frq_idx, :])
        @views ph[frq_idx, :] = DSP.angle.(cs[frq_idx, :])
    end

    p[p .== -Inf] .= minimum(p[p .!= -Inf])
    p[p .== +Inf] .= maximum(p[p .!= +Inf])
    db && (p = pow2db.(p))

    t = 0:(1 / fs):(length(s) / fs)
    t = linspace(t[1], t[end - 1], size(p, 2))

    return (cs = cs, p = p, ph = ph, f = f, t = t)

end

"""
    mwspectrogram(s; <keyword arguments>)

Calculate Morlet wavelet spectrogram for each channel of a matrix.

# Arguments

- `s::AbstractMatrix`: signal matrix (channels, samples)
- `pad::Int64=0`: number of zeros to append
- `db::Bool=true`: normalize powers to dB
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple:

- `cs::Array{ComplexF64, 3}`: convolved analytic signals
- `p::Array{Float64, 3}`: powers
- `ph::Array{Float64, 3}`: phases
- `f::Vector{Float64}`: frequencies
- `t::Vector{Float64}`: time points
"""
function mwspectrogram(
        s::AbstractMatrix;
        pad::Int64 = 0,
        db::Bool = true,
        fs::Int64,
        ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
        w::Bool = true,
    )::@NamedTuple{
        cs::Array{ComplexF64, 3}, p::Array{Float64, 3}, ph::Array{Float64, 3}, f::Vector{Float64}, t::Vector{Float64},
    }

    mwspec_data = mwspectrogram(@view(s[1, :]), pad = pad, db = db, fs = fs, ncyc = ncyc, w = w)
    f_tmp = mwspec_data.f
    t_tmp = mwspec_data.t

    cs = zeros(ComplexF64, length(f_tmp), length(t_tmp), size(s, 1))
    p  = zeros(length(f_tmp), length(t_tmp), size(s, 1))
    ph = zeros(length(f_tmp), length(t_tmp), size(s, 1))

    Threads.@threads :dynamic for ch_idx in axes(s, 1)
        @inbounds begin
            mwspec_data = mwspectrogram(
                @view(s[ch_idx, :]),
                pad = pad,
                db = db,
                fs = fs,
                ncyc = ncyc,
                w = w)
            cs[:, :, ch_idx] .= mwspec_data.cs
            p[:, :, ch_idx] .= mwspec_data.p
            ph[:, :, ch_idx] .= mwspec_data.ph
        end
    end

    return (cs = cs, p = p, ph = ph, f = f_tmp, t = t_tmp)

end

"""
    ghtspectrogram(s; <keyword arguments>)

Calculate spectrogram using Gaussian filter and Hilbert transform.

# Arguments

- `s::AbstractVector`: signal vector
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `db::Bool=true`: normalize powers to dB
- `gw::Real=10`: Gaussian width in Hz
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple:

- `p::Matrix{Float64}`: powers
- `ph::Matrix{Float64}`: phases
- `f::Vector{Float64}`: frequencies
- `t::Vector{Float64}`: time points
"""
function ghtspectrogram(
    s::AbstractVector;
    fs::Int64,
    db::Bool = true,
    gw::Real = 10,
    w::Bool = true
)::@NamedTuple{
    p::Matrix{Float64},
    ph::Matrix{Float64},
    f::Vector{Float64},
    t::Vector{Float64}
}

    !(fs >= 1) && throw(ArgumentError("fs must be ≥ 1."))

    flim = (0, fs / 2)
    nfrq = _tlength(flim)
    f = linspace(flim[1], flim[2], nfrq)

    win = w ? hanning(length(s)) : ones(length(s))
    sw = s .* win

    # pre-allocate outputs
    p = zeros(length(f), length(s))
    ph = zeros(length(f), length(s))

    @inbounds for frq_idx in eachindex(f)
        swg    = filter_g(sw, fs = fs, f = f[frq_idx], gw = gw)
        h_data = htransform(swg)
        p[frq_idx, :]  = h_data.p
        ph[frq_idx, :] = h_data.ph
    end

    p[p .== -Inf] .= minimum(p[p .!= -Inf]) - eps()
    p[p .== +Inf] .= maximum(p[p .!= +Inf]) + eps()
    db && (p = pow2db.(p))

    t = 0:(1 / fs):(length(s) / fs)
    t = linspace(t[1], t[end - 1], size(p, 2))

    return (p = p, ph = ph, f = f, t = t)

end

"""
    ghtspectrogram(s; <keyword arguments>)

Calculate spectrogram using Gaussian and Hilbert transform for each channel of a matrix.

# Arguments

- `s::AbstractArray`: signal matrix (channels, samples)
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `db::Bool=true`: normalize powers to dB
- `gw::Real=10`: Gaussian width in Hz
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple:

- `p::Array{Float64, 3}`: powers, shape `(freq, time, channels)`
- `ph::Array{Float64, 3}`: phases, shape `(freq, time, channels)`
- `f::Vector{Float64}`: frequencies
- `t::Vector{Float64}`: time points
"""
function ghtspectrogram(
    s::AbstractMatrix;
    fs::Int64,
    db::Bool = true,
    gw::Real = 10,
    w::Bool = true
)::@NamedTuple{
    p::Array{Float64, 3},
    ph::Array{Float64, 3},
    f::Vector{Float64},
    t::Vector{Float64}
}

    # pilot call to determine output dimensions
    ght_data = ghtspectrogram(
        @view(s[1, :]),
        fs = fs,
        db = db,
        gw = gw,
        w = w
    )
    f_tmp = ght_data.f
    t_tmp = ght_data.t

    # pre-allocate outputs
    p = zeros(length(f_tmp), length(t_tmp), size(s, 1))
    ph = zeros(length(f_tmp), length(t_tmp), size(s, 1))

    # calculate over channels
    @inbounds Threads.@threads :dynamic for ch_idx in axes(s, 1)
        ght_data = ghtspectrogram(
            @view(s[ch_idx, :]),
            fs = fs,
            db = db,
            gw = gw,
            w = w
        )
        p[:, :, ch_idx] = ght_data.p
        ph[:, :, ch_idx] = ght_data.ph
    end

    return (p = p, ph = ph, f = f_tmp, t = t_tmp)

end

"""
    cwtspectrogram(s; <keyword arguments>)

Calculate scaleogram using Continuous Wavelet Transformation (CWT).

# Arguments

- `s::AbstractVector`: signal vector
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

Named tuple:

- `m::Matrix{Float64}`: magnitudes
- `f::Vector{Float64}`: frequencies
- `t::Vector{Float64}`: time points
"""
function cwtspectrogram(
    s::AbstractVector;
    fs::Int64,
    wt::T = wavelet(Morlet(2π), β = 2)
)::@NamedTuple{
    m::Matrix{Float64},
    f::Vector{Float64},
    t::Vector{Float64}
} where {T <: CWT}

    !(fs >= 1) && throw(ArgumentError("fs must be ≥ 1."))

    m = abs.(ContinuousWavelets.cwt(s, wt)')

    # frequencies
    f = cwtfrq(s, fs = fs, wt = wt)

    # sort frequencies in ascending order
    f_idx = sortperm(f)
    f = f[f_idx]
    m = m[f_idx, :]

    t = 0:(1 / fs):(length(s) / fs)
    t = linspace(t[1], t[end - 1], size(m, 2))

    return (m = m, f = f, t = t)

end

"""
    cwtspectrogram(s; <keyword arguments>)

Calculate scaleogram using Continuous Wavelet Transformation (CWT) for each channel of a matrix.

# Arguments

- `s::AbstractMatrix`: signal matrix (channels, samples)
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

Named tuple:

- `m::Array{Float64, 3}`: magnitudes, shape `(time points, frequencies, channels)`
- `f::Vector{Float64}`: frequencies
- `t::Vector{Float64}`: time points
"""
function cwtspectrogram(
    s::AbstractMatrix;
    fs::Int64,
    wt::T = wavelet(Morlet(2π), β = 2)
)::@NamedTuple{
    m::Array{Float64, 3},
    f::Vector{Float64},
    t::Vector{Float64}
} where {T <: CWT}

    # pilot call to determine output dimensions
    cwt_data = cwtspectrogram(@view(s[1, :]), fs = fs, wt = wt)
    f = cwt_data.f
    t = cwt_data.t

    # pre-allocate output
    m = zeros(length(f), length(t), size(s, 1))
    @inbounds Threads.@threads :dynamic for ch_idx in axes(s, 1)
        m[:, :, ch_idx] = cwtspectrogram(@view(s[ch_idx, :]), fs = fs, wt = wt).m
    end

    return (m = m, f = f, t = t)

end

"""
    hhtspectrogram(s, t; <keyword arguments>)

Calculate spectrogram using Hilbert-Huang transform.

1. Empirical Mode Decomposition (EMD) – Decomposes the signal into Intrinsic Mode Functions (IMFs).
2. Hilbert Spectral Analysis (HSA) – Applies the Hilbert transform to each IMF to extract instantaneous frequency and powers, producing a spectrogram-like representation.

# Arguments

- `s::AbstractVector`: signal vector
- `t::AbstractVector`: time points (required for EMD)
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `db::Bool=true`: normalize powers to dB

# Returns

Named tuple:

- `p::Matrix{Float64}`: powers (frequencies × time points)
- `f::Vector{Float64}`: frequencies
- `t::Vector{Float64}`: time points
"""
function hhtspectrogram(
    s::AbstractVector,
    t::AbstractVector;
    fs::Int64,
    db::Bool = true
)::@NamedTuple{
    p::Matrix{Float64},
    f::Vector{Float64},
    t::Vector{Float64}
}

    !(fs >= 1) && throw(ArgumentError("fs must be ≥ 1."))

    # pre-allocate outputs
    imf_p = Vector{Vector{Float64}}()
    imf_fi = Vector{Vector{Float64}}()

    # perform EMD
    imfs = emd(s, t)
    imfs_n = size(imfs, 1)
    # last IMF is a residual, ignore it
    for idx2 in 1:(imfs_n - 1)
        # get instantaneous frequencies of IMF
        push!(imf_fi, frqinst(imfs[idx2, :]) .* fs)
        # perform Hilbert transform and get instantaneous powers of IMF
        push!(imf_p, htransform(imfs[idx2, :]).p)
    end

    # project IMF powers onto the frequency grid using instantaneous frequency

    f = collect(0:0.5:fs ÷ 2) # frequency range of interest (Hz)
    p = zeros(length(f), length(t))

    # interpolate each IMF's frequency and power onto the grid
    @inbounds for (pow, freq) in zip(imf_p, imf_fi)
        for idx_t in eachindex(t)
            # find the closest frequency bin
            f_idx = argmin(abs.(f .- freq[idx_t]))
            if 1 ≤ f_idx ≤ length(f)
                p[f_idx, idx_t] += pow[idx_t] # add power
            end
        end
    end

    db && (p = pow2db.(p))
    p[p .== -Inf] .= -eps()

    return (p = p, f = f, t = t)

end


"""
    hhtspectrogram(s; <keyword arguments>)

Calculate spectrogram using Hilbert-Huang transform.

1. Empirical Mode Decomposition (EMD) – Decomposes the signal into Intrinsic Mode Functions (IMFs).
2. Hilbert Spectral Analysis (HSA) – Applies the Hilbert transform to each IMF to extract instantaneous frequency and powers, producing a spectrogram-like representation.

# Arguments

- `s::AbstractMatrix`: signal matrix
- `t::AbstractVector`: time points (required for EMD)
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `db::Bool=true`: normalize powers to dB

# Returns

Named tuple:

- `p::Matrix{Float64}`: powers (frequencies × time points)
- `f::Vector{Float64}`: frequencies
- `t::Vector{Float64}`: time points
"""
function hhtspectrogram(
    s::AbstractMatrix,
    t::AbstractVector;
    fs::Int64,
    db::Bool = true
)::@NamedTuple{
    p::Matrix{Float64},
    f::Vector{Float64},
    t::Vector{Float64}
}

    !(fs >= 1) && throw(ArgumentError("fs must be ≥ 1."))

    # pre-allocate outputs
    imf_p = Vector{Vector{Float64}}()
    imf_fi = Vector{Vector{Float64}}()

    # perform EMD
    @inbounds for idx1 in axes(s, 1)
        imfs = emd(s[idx1, :], t)
        imfs_n = size(imfs, 1)
        # last IMF is a residual, ignore it
        for idx2 in 1:(imfs_n - 1)
            # get instantaneous frequencies of IMF
            push!(imf_fi, frqinst(imfs[idx2, :]) .* fs)
            # perform Hilbert transform and get instantaneous powers of IMF
            push!(imf_p, htransform(imfs[idx2, :]).p)
        end
    end

    # project IMF powers onto the frequency grid using instantaneous frequency

    f = collect(0:0.5:fs ÷ 2) # frequency range of interest (Hz)
    p = zeros(length(f), length(t))

    # interpolate each IMF's frequency and power onto the grid
    @inbounds for (pow, freq) in zip(imf_p, imf_fi)
        for idx_t in eachindex(t)
            # find the closest frequency bin
            f_idx = argmin(abs.(f .- freq[idx_t]))
            if 1 ≤ f_idx ≤ length(f)
                p[f_idx, idx_t] += pow[idx_t] # add power
            end
        end
    end

    db && (p = pow2db.(p))
    p[p .== -Inf] .= -eps()

    return (p = p, f = f, t = t)

end
