export spectrogram
export mwspectrogram
export ghspectrogram
export cwtspectrogram

"""
    spectrogram(s; <keyword arguments>)

Calculate spectrogram. Default method is short time Fourier transform.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling frequency
- `db::Bool=true`: normalize powers to dB
- `method::Symbol=:stft`: method used to calculate PSD:
    - `:stft`: short time Fourier transform
    - `:mt`: multi-tapered periodogram
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length, default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:
- `p::Matrix{Float64}`: powers
- `f::Vector{Float64}`: frequencies
- `t::Vector{Float64}`: time
"""
function spectrogram(s::AbstractVector; fs::Int64, db::Bool=true, method::Symbol=:stft, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)::@NamedTuple{p::Matrix{Float64}, f::Vector{Float64}, t::Vector{Float64}}

    _check_var(method, [:stft, :mt], "method")
    @assert fs >= 1 "fs must be ≥ 1."
    @assert wlen <= length(s) "wlen must be ≤ $(length(s))."
    @assert wlen >= 1 "wlen must be ≥ 1."
    @assert woverlap <= wlen "woverlap must be ≤ $(wlen)."
    @assert woverlap >= 0 "woverlap must be ≥ 0."

    if method === :stft
        w = w ? hanning : nothing
        if length(s) < fs * 1.1
            wlen = div(wlen, 2)
            woverlap = div(woverlap, 2)
        end
        p = DSP.spectrogram(s, wlen, woverlap, fs=fs, window=w)
    elseif method === :mt
        w = w ? hanning(length(s)) : ones(length(s))
        p = DSP.mt_spectrogram(s .* w, fs=fs, nw=((nt + 1) ÷ 2), ntapers=nt)
    end

    p = p.power
    p[p .== -Inf] .= minimum(p[p .!== -Inf])
    p[p .== +Inf] .= maximum(p[p .!== +Inf])
    db && (p = pow2db.(p))

    t = 0:1/fs:(length(s) / fs)
    t = linspace(t[1], t[end - 1], size(p, 2))
    f = linspace(0, fs/2, size(p, 1))

    return (p=p, f=f, t=t)

end

"""
    mwspectrogram(s; <keyword arguments>)

Calculate spectrogram using wavelet convolution.

# Arguments

- `s::AbstractVector`
- `pad::Int64`: pad with `pad` zeros
- `db::Bool=true`: normalize powers to dB
- `fs::Int64`: sampling rate
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(fs / 2)`
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:
- `cs::Matrix{ComplexF64}`: convoluted signal
- `p::Matrix{Float64}`: powers
- `ph::Matrix{Float64}`: phases
- `f::Vector{Float64}`: frequencies
- `t::Vector{Float64}`: time
"""
function mwspectrogram(s::AbstractVector; pad::Int64=0, db::Bool=true, fs::Int64, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, w::Bool=true)::@NamedTuple{cs::Matrix{ComplexF64}, p::Matrix{Float64}, ph::Matrix{Float64}, f::Vector{Float64}, t::Vector{Float64}}

    @assert fs >= 1 "fs must be > 1."

    pad > 0 && (s = pad0(s, pad))

    w = w ? hanning(length(s)) : ones(length(s))

    if ncyc isa Int64
        @assert ncyc >= 1 "ncyc must be ≥ 1."
    else
        @assert ncyc[1] >= 1 "ncyc[1] must be ≥ 1."
        @assert ncyc[2] >= 1 "ncyc[2] must be ≥ 1."
    end

    # get frequency range
    frq_lim = (0, fs / 2)
    frq_n = _tlength(frq_lim)
    f = linspace(frq_lim[1], frq_lim[2], frq_n)

    cs = zeros(ComplexF64, length(f), length(s))
    p = zeros(length(f), length(s))
    ph = zeros(length(f), length(s))

    if ncyc isa Int64
        ncyc = repeat([ncyc], frq_n)
    else
        ncyc = round.(Int64, log10space(log10(ncyc[1]), log10(ncyc[2]), frq_n))
    end

    @inbounds for frq_idx in 1:frq_n
        kernel = generate_morlet(fs, f[frq_idx], 1, ncyc=ncyc[frq_idx], complex=true)
        # cs[frq_idx, :] = fconv(s .* w, kernel=kernel, db=false)
        cs[frq_idx, :] = fconv(s .* w, kernel=kernel, norm=true)
        # alternative: a[frq_idx, :] = LinearAlgebra.norm.(real.(cs[frq_idx, :]), imag.(cs[frq_idx]))
        p[frq_idx, :] = @views @. abs(cs[frq_idx, :])^2
        ph[frq_idx, :] = @views @. angle(cs[frq_idx, :])
    end

    p[p .== -Inf] .= minimum(p[p .!== -Inf])
    p[p .== +Inf] .= maximum(p[p .!== +Inf])
    db && (p = pow2db.(p))

    t = 0:1/fs:(length(s) / fs)
    t = linspace(t[1], t[end - 1], size(p, 2))

    return (cs=cs, p=p, ph=ph, f=f, t=t)

end

"""
    ghspectrogram(s; <keyword arguments>)

Calculate spectrogram using Gaussian and Hilbert transform.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `db::Bool=true`: normalize powers to dB
- `gw::Real=5`: Gaussian width in Hz
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:
- `p::Matrix{Float64}`: powers
- `ph::Matrix{Float64}`: phases
- `f::Vector{Float64}`: frequencies
- `t::Vector{Float64}`: time
"""
function ghspectrogram(s::AbstractVector; fs::Int64, db::Bool=true, gw::Real=5, w::Bool=true)::@NamedTuple{p::Matrix{Float64}, ph::Matrix{Float64}, f::Vector{Float64}, t::Vector{Float64}}

    @assert fs >= 1 "fs must be ≥ 1."

    frq_lim = (0, fs / 2)
    frq_n = _tlength(frq_lim)
    f = linspace(frq_lim[1], frq_lim[2], frq_n)

    w = w ? hanning(length(s)) : ones(length(s))

    p = zeros(length(f), length(s))
    ph = zeros(length(f), length(s))

    @inbounds for frq_idx in eachindex(f)
        s_tmp = filter_g(s .* w, fs=fs, f=f[frq_idx], gw=gw)
        p[frq_idx, :] = (abs.(hilbert(s_tmp))).^2
        ph[frq_idx, :] = angle.(hilbert(s_tmp))
    end

    p[p .== -Inf] .= minimum(p[p .!== -Inf])
    p[p .== +Inf] .= maximum(p[p .!== +Inf])
    db && (p = pow2db.(p))

    t = 0:1/fs:(length(s) / fs)
    t = linspace(t[1], t[end - 1], size(p, 2))


    return (p=p, ph=ph, f=f, t=t)

end

"""
    cwtspectrogram(s; <keyword arguments>)

Calculate scaleogram using continuous wavelet transformation (CWT).

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `norm::Bool=true`: normalize scaleogram to the signal scale so the amplitudes of wavelet coefficients agree with the amplitudes of oscillatory components in a signal

# Returns

Named tuple containing:
- `p::Matrix{Float64}`: powers
- `f::Vector{Float64}`: frequencies
- `t::Vector{Float64}`: time
"""
function cwtspectrogram(s::AbstractVector; fs::Int64, wt::T=wavelet(Morlet(2π), β=32, Q=128), norm::Bool=true)::@NamedTuple{p::Matrix{Float64}, f::Vector{Float64}, t::Vector{Float64}} where {T <: CWT}

    @assert fs >= 1 "fs must be ≥ 1."

    # w = w ? hanning(length(s)) : ones(length(s))

    p = abs.(ContinuousWavelets.cwt(s, wt)')

    # scale
    if norm
        a = amp(s)[1]
        p = normalize_n(p, a)
    end

    f = cwtfrq(s, fs=fs, wt=wt)

    # reverse order
    f_idx = sortperm(f)
    f = f[f_idx]
    p = p[f_idx, :]

    # p[p .== -Inf] .= minimum(p[p .!== -Inf])
    # p[p .== +Inf] .= maximum(p[p .!== +Inf])
    # db && (p = pow2db.(p))

    t = 0:1/fs:(length(s) / fs)
    t = linspace(t[1], t[end - 1], size(p, 2))

    return (p=p, f=f, t=t)

end

"""
    spectrogram(obj; <keyword arguments>)

Calculate spectrogram. Default method is short time Fourier transform.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `pad::Int64=0`: number of zeros to add
- `method::Symbol=:stft`: method of calculating spectrogram:
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `db::Bool=true`: normalize powers to dB; for CWT scaleogram: normalize to the signal scale so the amplitudes of wavelet coefficients agree with the amplitudes of oscillatory components in a signal
- `nt::Int64=7`: number of Slepian tapers
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:
- `p::Array{Float64, 4}`: powers
- `f::Vector{Float64}`: frequencies
- `t::Vector{Float64}`: time points
"""
function spectrogram(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, pad::Int64=0, method::Symbol=:stft, db::Bool=true, nt::Int64=7, gw::Real=5, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, wt::T=wavelet(Morlet(2π), β=32, Q=128), wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)::@NamedTuple{p::Array{Float64, 4}, f::Vector{Float64}, t::Vector{Float64}} where {T <: CWT}

    _check_var(method, [:stft, :mt, :mw, :gh, :cwt], "method")
    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    ch_n = length(ch)
    ep_n = nepochs(obj)
    fs = sr(obj)

    if method === :stft
        p_tmp, f, _ = @views NeuroAnalyzer.spectrogram(obj.data[1, :, 1], fs=fs, db=db, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
    elseif method === :mt
        p_tmp, f, _ = @views NeuroAnalyzer.spectrogram(obj.data[1, :, 1], fs=fs, db=db, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
    elseif method === :mw
    _, p_tmp, _, f, _ = @views NeuroAnalyzer.mwspectrogram(obj.data[1, :, 1], pad=pad, fs=fs, db=db, ncyc=ncyc, w=w)
    elseif method === :gh
        p_tmp, _, f, _ = @views NeuroAnalyzer.ghspectrogram(obj.data[1, :, 1], fs=fs, db=db, gw=gw, w=w)
    elseif method === :cwt
        _log_off()
        p_tmp, f, _ = @views NeuroAnalyzer.cwtspectrogram(obj.data[1, :, 1], fs=fs, wt=wt, norm=db)
        _log_on()
    end

    t = linspace(0, (epoch_len(obj) / fs), size(p_tmp, 2))
    p = zeros(size(p_tmp, 1), size(p_tmp, 2), ch_n, ep_n)

    # initialize progress bar
    progress_bar && (progbar = Progress(ep_n * ch_n, dt=1, barlen=20, color=:white))

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            if method === :stft
                p[:, :, ch_idx, ep_idx], _, _ = @views NeuroAnalyzer.spectrogram(obj.data[ch[ch_idx], :, ep_idx], fs=fs, db=db, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
            elseif method === :mt
                p[:, :, ch_idx, ep_idx], _, _ = @views NeuroAnalyzer.spectrogram(obj.data[ch[ch_idx], :, ep_idx], fs=fs, db=db, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
            elseif method === :mw
                _, p[:, :, ch_idx, ep_idx], _, _ = @views NeuroAnalyzer.mwspectrogram(obj.data[ch[ch_idx], :, ep_idx], pad=pad, fs=fs, db=db, ncyc=ncyc, w=w)
            elseif method === :gh
                p[:, :, ch_idx, ep_idx], _, _ = @views NeuroAnalyzer.ghspectrogram(obj.data[ch[ch_idx], :, ep_idx], fs=fs, db=db, gw=gw, w=w)
            elseif method === :cwt
                _log_off()
                p[:, :, ch_idx, ep_idx], _ = @views NeuroAnalyzer.cwtspectrogram(obj.data[ch[ch_idx], :, ep_idx], fs=fs, wt=wt, norm=db)
                _log_on()
            end

            # update progress bar
            progress_bar && next!(progbar)
        end
    end

    f = round.(f, digits=2)
    t = round.(t, digits=2)
    t .+= obj.epoch_time[1]

    return (p=p, f=f, t=t)

end
