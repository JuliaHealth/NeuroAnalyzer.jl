export spectrogram
export mwspectrogram
export ghspectrogram
export cwtspectrogram

"""
    spectrogram(s; fs, db, method, nt, wlen, woverlap, w)

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
- `sp::Matrix{Float64}`: powers
- `sf::Vector{Float64}`: frequencies
- `st::Vector{Float64}`: time
"""
function spectrogram(s::AbstractVector; fs::Int64, db::Bool=true, method::Symbol=:stft, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)

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
        sp = DSP.spectrogram(s, wlen, woverlap, fs=fs, window=w)
    elseif method === :mt
        w = w ? hanning(length(s)) : ones(length(s))
        sp = DSP.mt_spectrogram(s .* w, fs=fs, nw=((nt + 1) ÷ 2), ntapers=nt)
    end

    sp = sp.power
    sp[sp .== -Inf] .= minimum(sp[sp .!== -Inf])
    sp[sp .== +Inf] .= maximum(sp[sp .!== +Inf])
    db && (sp = pow2db.(sp))

    t = 0:1/fs:(length(s) / fs)
    st = linspace(t[1], t[end - 1], size(sp, 2))
    sf = linspace(0, fs/2, size(sp, 1))

    return (sp=sp, sf=sf, st=st)

end

"""
    mwspectrogram(s; pad, db, fs, ncyc, w)

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
- `cs::Matrix(ComplexF64}`: convoluted signal
- `sp::Matrix{Float64}`: powers
- `sph::Matrix{Float64}`: phases
- `sf::Vector{Float64}`: frequencies
- `st::Vector{Float64}`: time
"""
function mwspectrogram(s::AbstractVector; pad::Int64=0, db::Bool=true, fs::Int64, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, w::Bool=true)

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
    sf = linspace(frq_lim[1], frq_lim[2], frq_n)

    cs = zeros(ComplexF64, length(sf), length(s))
    sp = zeros(length(sf), length(s))
    sph = zeros(length(sf), length(s))

    if ncyc isa Int64
        ncyc = repeat([ncyc], frq_n)
    else
        ncyc = round.(Int64, linspace(ncyc[1], ncyc[2], frq_n))
    end

    @inbounds for frq_idx in 1:frq_n
        kernel = generate_morlet(fs, sf[frq_idx], 1, ncyc=ncyc[frq_idx], complex=true)
        # cs[frq_idx, :] = fconv(s .* w, kernel=kernel, db=false)
        cs[frq_idx, :] = fconv(s .* w, kernel=kernel, norm=true)
        # alternative: w_amp[frq_idx, :] = LinearAlgebra.db.(real.(cs), imag.(cs), 2)
        sp[frq_idx, :] = @views @. (2 * abs(cs[frq_idx, :]))^2
        sph[frq_idx, :] = @views @. angle(cs[frq_idx, :])
    end

    sp[sp .== -Inf] .= minimum(sp[sp .!== -Inf])
    sp[sp .== +Inf] .= maximum(sp[sp .!== +Inf])
    db && (sp = pow2db.(sp))

    t = 0:1/fs:(length(s) / fs)
    st = linspace(t[1], t[end - 1], size(sp, 2))

    return (cs=cs, sp=sp, sph=sph, sf=sf, st=st)

end

"""
    ghspectrogram(s; fs, db, gw, w)

Calculate spectrogram using Gaussian and Hilbert transform.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `db::Bool=true`: normalize powers to dB
- `gw::Real=5`: Gaussian width in Hz
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:
- `sp::Matrix{Float64}`: powers
- `sph::Matrix{Float64}`: phases
- `sf::Vector{Float64}`: frequencies
- `st::Vector{Float64}`: time
"""
function ghspectrogram(s::AbstractVector; fs::Int64, db::Bool=true, gw::Real=5, w::Bool=true)

    @assert fs >= 1 "fs must be ≥ 1."

    frq_lim = (0, fs / 2)
    frq_n = _tlength(frq_lim)

    w = w ? hanning(length(s)) : ones(length(s))

    sf = linspace(frq_lim[1], frq_lim[2], frq_n)
    sp = zeros(length(sf), length(s))
    sph = zeros(length(sf), length(s))

    @inbounds for frq_idx in eachindex(sf)
        s_tmp = filter_g(s .* w, fs=fs, f=sf[frq_idx], gw=gw)
        sp[frq_idx, :] = (abs.(hilbert(s_tmp))).^2
        sph[frq_idx, :] = angle.(hilbert(s_tmp))
    end

    sp[sp .== -Inf] .= minimum(sp[sp .!== -Inf])
    sp[sp .== +Inf] .= maximum(sp[sp .!== +Inf])
    db && (sp = pow2db.(sp))

    t = 0:1/fs:(length(s) / fs)
    st = linspace(t[1], t[end - 1], size(sp, 2))

    return (sp=sp, sph=sph, sf=sf, st=st)

end

"""
    cwtspectrogram(s; fs, wt, w, db)

Calculate spectrogram using continuous wavelet transformation (CWT).

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `w::Bool=true`: if true, apply Hanning window
- `db::Bool=true`: normalize powers to dB

# Returns

Named tuple containing:
- `sp::Matrix{Float64}`: powers
- `sf::Vector{Float64}`: frequency indices
- `st::Vector{Float64}`: time
"""
function cwtspectrogram(s::AbstractVector; fs::Int64, wt::T=wavelet(Morlet(2π), β=32, Q=128), w::Bool=true, db::Bool=true) where {T <: CWT}

    @assert fs >= 1 "fs must be ≥ 1."

    frq_lim = (0, fs / 2)

    w = w ? hanning(length(s)) : ones(length(s))

    sp = abs.(ContinuousWavelets.cwt(s .* w, wt)') .^ 2
    sf = round.(ContinuousWavelets.getMeanFreq(s .* w, wt, fs), digits=2)
    sf_idx = sortperm(sf)
    sf = sf[sf_idx]
    sp = sp[sf_idx, :]
    sf = sf[vsearch(frq_lim[1], sf):vsearch(frq_lim[2], sf)]
    sp = sp[vsearch(frq_lim[1], sf):vsearch(frq_lim[2], sf), :]
    sp[sp .== -Inf] .= minimum(sp[sp .!== -Inf])
    sp[sp .== +Inf] .= maximum(sp[sp .!== +Inf])
    db && (sp = pow2db.(sp))

    t = 0:1/fs:(length(s) / fs)
    st = linspace(t[1], t[end - 1], size(sp, 2))

    return (sp=sp, sf=sf, st=st)

end

"""
    spectrogram(obj; ch, pad, method, db, nt, gw, ncyc, wt, wlen, woverlap, w, wt, gw)

Calculate spectrogram. Default method is short time Fourier transform.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `pad::Int64=0`: number of zeros to add
- `method::Symbol=:stft`: method of calculating spectrogram:
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `db::Bool=true`: normalize powers to dB
- `nt::Int64=7`: number of Slepian tapers
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:
- `sp::Array{Float64, 3}`: powers
- `sf::Vector{Float64}`: frequencies (frequency indices for continuous wavelet transformation)
- `st::Vector{Float64}`: time points
"""
function spectrogram(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), pad::Int64=0, method::Symbol=:stft, db::Bool=true, nt::Int64=7, gw::Real=5, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, wt::T=wavelet(Morlet(2π), β=32, Q=128), wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true) where {T <: CWT}

    _check_var(method, [:stft, :mt, :mw, :gh, :cwt], "method")
    _check_channels(obj, ch)

    ch_n = length(ch)
    ep_n = nepochs(obj)
    fs = sr(obj)

    if method === :stft
        p_tmp, sf, _ = @views NeuroAnalyzer.spectrogram(obj.data[1, :, 1], fs=fs, db=db, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
    elseif method === :mt
        p_tmp, sf, _ = @views NeuroAnalyzer.spectrogram(obj.data[1, :, 1], fs=fs, db=db, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
    elseif method === :mw
    _, p_tmp, _, sf, _ = @views NeuroAnalyzer.mwspectrogram(obj.data[1, :, 1], pad=pad, fs=fs, db=db, ncyc=ncyc, w=w)
    elseif method === :gh
        p_tmp, _, sf, _ = @views NeuroAnalyzer.ghspectrogram(obj.data[1, :, 1], fs=fs, db=db, gw=gw, w=w)
    elseif method === :cwt
        p_tmp, sf, _ = @views NeuroAnalyzer.cwtspectrogram(obj.data[1, :, 1], fs=fs, db=db, wt=wt, w=w)
    end

    st = linspace(0, (epoch_len(obj) / fs), size(p_tmp, 2))
    sp = zeros(size(p_tmp, 1), size(p_tmp, 2), ch_n, ep_n)

    # initialize progress bar
    progress_bar && (progbar = Progress(ep_n * ch_n, dt=1, barlen=20, color=:white))

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            if method === :stft
                sp[:, :, ch_idx, ep_idx], _, _ = @views NeuroAnalyzer.spectrogram(obj.data[ch[ch_idx], :, ep_idx], fs=fs, db=db, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
            elseif method === :mt
                sp[:, :, ch_idx, ep_idx], _, _ = @views NeuroAnalyzer.spectrogram(obj.data[ch[ch_idx], :, ep_idx], fs=fs, db=db, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
            elseif method === :mw
                _, sp[:, :, ch_idx, ep_idx], _, _ = @views NeuroAnalyzer.mwspectrogram(obj.data[ch[ch_idx], :, ep_idx], pad=pad, fs=fs, db=db, ncyc=ncyc, w=w)
            elseif method === :gh
                sp[:, :, ch_idx, ep_idx], _, _ = @views NeuroAnalyzer.ghspectrogram(obj.data[ch[ch_idx], :, ep_idx], fs=fs, db=db, gw=gw, w=w)
            elseif method === :cwt
                sp[:, :, ch_idx, ep_idx], _ = @views NeuroAnalyzer.cwtspectrogram(obj.data[ch[ch_idx], :, ep_idx], fs=fs, db=db, wt=wt, w=w)
            end

            # update progress bar
            progress_bar && next!(progbar)
        end
    end

    sf = round.(sf, digits=2)
    st = round.(st, digits=2)
    st .+= obj.epoch_time[1]

    return (sp=sp, sf=sf, st=st)

end
