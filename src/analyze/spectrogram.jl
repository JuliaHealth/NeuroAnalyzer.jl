export spectrogram
export wspectrogram
export ghspectrogram
export cwtspectrogram
export specseg

"""
    spectrogram(signal; fs, norm, mt, st, demean)

Calculate spectrogram.

# Arguments

- `signal::AbstractVector`
- `fs::Int64`: sampling frequency
- `norm::Bool=true`: normalize powers to dB
- `mt::Bool=false`: if true use multi-tapered spectrogram
- `st::Bool=false`: if true use short time Fourier transform
- `demean::Bool=true`: demean signal prior to analysis

# Returns

Named tuple containing:
- `s_pow::Matrix{Float64}`: powers
- `s_frq::Vector{Float64}`: frequencies
- `s_t::Vector{Float64}`: time
"""
function spectrogram(signal::AbstractVector; fs::Int64, norm::Bool=true, mt::Bool=false, st::Bool=false, demean::Bool=true)

    (mt == true && st == true) && throw(ArgumentError("Both mt and st must not be true."))
    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))

    demean == true && (signal = demean(signal))

    nfft = length(signal)
    interval = fs
    overlap = round(Int64, fs * 0.75)

    if st == true
        s_pow = abs.(stft(signal, interval, overlap, nfft=nfft, fs=fs, window=hanning))
        norm == true && (s_pow = pow2db.(s_pow))
        t = 0:1/fs:(length(signal) / fs)
        s_t = linspace(t[1], t[end], size(s_pow, 2))
        s_frq = linspace(0, fs/2, size(s_pow, 1))
        return (s_pow=s_pow, s_frq=s_frq, s_t=s_t)
    end

    if mt == true
        spec = mt_spectrogram(signal, fs=fs)
    else
        spec = spectrogram(signal, interval, overlap, nfft=nfft, fs=fs, window=hanning)
    end    
    s_pow = spec.power
    norm == true ? s_pow = pow2db.(spec.power) : s_pow = spec.power

    #s_t = collect(spec.time)
    #s_frq = Vector(spec.freq)
    t = 0:1/fs:(length(signal) / fs)
    s_t = linspace(t[1], t[end], size(s_pow, 2))
    s_frq = linspace(0, fs/2, size(s_pow, 1))
    
    return (s_pow=s_pow, s_frq=s_frq, s_t=s_t)
end

"""
    wspectrogram(signal; pad, norm, frq_lim, frq_n, frq, fs, ncyc, demean)

Calculate spectrogram using wavelet convolution.

# Arguments

- `signal::AbstractVector`
- `pad::Int64`: pad the `signal` with `pad` zeros
- `norm::Bool=true`: normalize powers to dB
- `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
- `frq_n::Int64`: number of frequencies
- `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
- `fs::Int64`: sampling rate
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n) for frq === :log or ncyc = linspace(ncyc[1], ncyc[2], frq_n) for frq === :lin
- `demean::Bool=true`: demean signal prior to analysis

# Returns

Named tuple containing:
- `w_conv::Matrix(ComplexF64}`: convoluted signal
- `w_powers::Matrix{Float64}`
- `w_phases::Matrix{Float64}`
- `frq_list::Vector{Float64}`
"""
function wspectrogram(signal::AbstractVector; pad::Int64=0, norm::Bool=true, frq_lim::Tuple{Real, Real}, frq_n::Int64, frq::Symbol=:lin, fs::Int64, ncyc::Union{Int64, Tuple{Int64, Int64}}=6, demean::Bool=true)

    _check_var(frq, [:log, :lin], "frq")

    pad > 0 && (signal = pad0(signal, pad))
    if typeof(ncyc) == Int64
        ncyc < 1 && throw(ArgumentError("ncyc must be ≥ 1."))
    else
        ncyc[1] < 1 && throw(ArgumentError("ncyc[1] must be ≥ 1."))
        ncyc[2] < 1 && throw(ArgumentError("ncyc[2] must be ≥ 1."))
    end

    # get frequency range
    fs < 1 && throw(ArgumentError("fs must be > 1."))
    frq_lim = tuple_order(frq_lim)
    frq_lim[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
    frq_lim[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be ≤ $(fs / 2)."))
    frq_n < 2 && throw(ArgumentError("frq_n must be ≥ 2."))
    frq_lim[1] == 0 && (frq_lim = (0.1, frq_lim[2]))
    if frq === :log
        # frq_lim = (frq_lim[1], frq_lim[2])
        frq_list = round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), frq_n), digits=1)
    else
        frq_list = linspace(frq_lim[1], frq_lim[2], frq_n)
    end

    demean == true && (signal = demean(signal))
    w_conv = zeros(ComplexF64, length(frq_list), length(signal))
    w_powers = zeros(length(frq_list), length(signal))
    w_amp = zeros(length(frq_list), length(signal))
    w_phases = zeros(length(frq_list), length(signal))

    if typeof(ncyc) != Tuple{Int64, Int64}
        ncyc = repeat([ncyc], frq_n)
    else
        if frq === :log
            ncyc = round.(Int64, logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n))
        else
            ncyc = round.(Int64, linspace(ncyc[1], ncyc[2], frq_n))
        end
    end

    @inbounds @simd for frq_idx in 1:frq_n
        kernel = generate_morlet(fs, frq_list[frq_idx], 1, ncyc=ncyc[frq_idx], complex=true)
        w_conv[frq_idx, :] = fconv(signal, kernel=kernel, norm=false)
        # alternative: w_amp[frq_idx, :] = LinearAlgebra.norm.(real.(w_conv), imag.(w_conv), 2)
        w_powers[frq_idx, :] = @views @. abs(w_conv[frq_idx, :])^2
        w_phases[frq_idx, :] = @views @. angle(w_conv[frq_idx, :])
    end

    # remove reflected part of the signal
    # if reflect == true
    #     w_conv = w_conv[:, (length(signal) ÷ 3 + 1):(2 * length(signal) ÷ 3)]
    #     w_amp = w_amp[:, (length(signal) ÷ 3 + 1):(2 * length(signal) ÷ 3)]
    #     w_powers = w_powers[:, (length(signal) ÷ 3 + 1):(2 * length(signal) ÷ 3)]
    #     w_phases = w_phases[:, (length(signal) ÷ 3 + 1):(2 * length(signal) ÷ 3)]
    # end

    norm == true && (w_powers = pow2db.(w_powers))

    return (w_conv=w_conv, w_powers=w_powers, w_phases=w_phases, frq_list=frq_list)
end

"""
    ghspectrogram(signal; fs, norm, frq_lim, frq_n, frq, fs, demean)

Calculate spectrogram using Gaussian and Hilbert transform.

# Arguments

- `signal::AbstractVector`
- `fs::Int64`: sampling rate
- `norm::Bool=true`: normalize powers to dB
- `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
- `frq_n::Int64`: number of frequencies
- `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
- `gw::Real=5`: Gaussian width in Hz
- `demean::Bool`=true: demean signal prior to analysis

# Returns

Named tuple containing:
- `s_pow::Matrix{Float64}`: powers
- `s_pow::Matrix{Float64}`: phases
- `frq_list::Vector{Float64}`: frequencies
"""
function ghspectrogram(signal::AbstractVector; fs::Int64, norm::Bool=true, frq_lim::Tuple{Real, Real}, frq_n::Int64, frq::Symbol=:lin, gw::Real=5, demean::Bool=true)

    _check_var(frq, [:log, :lin], "frq")

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    frq_lim = tuple_order(frq_lim)
    frq_lim[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
    frq_lim[2] > fs ÷ 2 && throw(ArgumentError("Upper frequency bound must be ≤ $(fs ÷ 2)."))
    frq_n < 2 && throw(ArgumentError("frq_n frequency bound must be ≥ 2."))
    frq_lim[1] == 0 && (frq_lim = (0.1, frq_lim[2]))
    if frq === :log
        frq_lim = (frq_lim[1], frq_lim[2])
        s_frq = round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), frq_n), digits=1)
    else
        s_frq = linspace(frq_lim[1], frq_lim[2], frq_n)
    end

    demean == true && (signal = demean(signal))
    s_pow = zeros(length(s_frq), length(signal))
    s_ph = zeros(length(s_frq), length(signal))
    @inbounds @simd for frq_idx in eachindex(s_frq)
        s = gfilter(signal, fs=fs, f=s_frq[frq_idx], gw=gw)
        s_pow[frq_idx, :] = abs.(hilbert(s)).^2
        s_ph[frq_idx, :] = angle.(hilbert(s))
    end
    norm == true && (s_pow = pow2db.(s_pow))

    return (s_pow=s_pow, s_ph=s_ph, s_frq=s_frq)
end

"""
    cwtspectrogram(signal; wt, pad, norm, frq_lim, fs, demean)

Calculate spectrogram using continuous wavelet transformation (CWT).

# Arguments

- `signal::AbstractVector`
- `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `fs::Int64`: sampling rate
- `norm::Bool=true`: normalize powers to dB
- `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
- `demean::Bool`=true: demean signal prior to analysis

# Returns

Named tuple containing:
- `h_powers::Matrix{Float64}`
- `frq_list::Vector{Float64}`
"""
function cwtspectrogram(signal::AbstractVector; wt::T, fs::Int64, norm::Bool=true, frq_lim::Tuple{Real, Real}, demean::Bool=true) where {T <: CWT}

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    frq_lim = tuple_order(frq_lim)
    frq_lim[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
    frq_lim[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be ≤ $(fs / 2)."))

    demean == true && (signal = demean(signal))

    h_powers = abs.(ContinuousWavelets.cwt(signal, wt)')
    frq_list = ContinuousWavelets.getMeanFreq(ContinuousWavelets.computeWavelets(length(signal), wt)[1])
    frq_list[1] = 0
    frq_list[1] < frq_lim[1] && throw(ArgumentError("Lower frequency bound must be ≥ $(frq_list[1])."))
    frq_list[end] < frq_lim[2] && throw(ArgumentError("Upper frequency bound must be ≤ $(frq_list[end])."))
    frq_list = frq_list[vsearch(frq_lim[1], frq_list):vsearch(frq_lim[2], frq_list)]
    h_powers = h_powers[vsearch(frq_lim[1], frq_list):vsearch(frq_lim[2], frq_list), :]

    norm == true && (h_powers = pow2db.(h_powers))

    return (h_powers=h_powers, frq_list=frq_list)
end

"""
    spectrogram(obj; channel, norm, mt, st, demean)

Calculate spectrogram.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `method::Symbol=:standard`: method of calculating spectrogram:
    - `:standard`: standard
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `pad::Int64=0`: number of zeros to add
- `frq_lim::Tuple{Real, Real}=(0, 0)`: frequency limits
- `frq_n::Int64=0`: number of frequencies, default is length(frq_lim[1]:frq_lim[2])
- `norm::Bool=true`: normalize powers to dB
- `demean::Bool=true`: demean signal prior to analysis
- `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n) for frq === :log or ncyc = linspace(ncyc[1], ncyc[2], frq_n) for frq === :lin
- `wt<:CWT=wavelet(Morlet(π), β=2)`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

Named tuple containing:
- `s_pow::Array{Float64, 3}`
- `s_frq::Vector{Float64}`
- `s_t::Vector{Float64}`
"""
function spectrogram(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), pad::Int64=0, frq_lim::Tuple{Real, Real}=(0, 0), frq_n::Int64=0, method::Symbol=:standard, norm::Bool=true, demean::Bool=true, frq::Symbol=:log, gw::Real=5, ncyc::Union{Int64, Tuple{Int64, Int64}}=6, wt::T=wavelet(Morlet(π), β=2)) where {T <: CWT}

    _check_var(method, [:standard, :stft, :mt, :mw, :gh, :cwt], "method")
    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    # get frequency range
    fs = sr(obj)
    frq_lim == (0, 0) && (frq_lim = (0, div(fs, 2)))
    frq_n == 0 && (frq_n = length(frq_lim[1]:frq_lim[2]))

    if method === :standard
        p_tmp, s_frq, _ = @views spectrogram(obj.data[1, :, 1], fs=fs, norm=norm, mt=false, st=false, demean=demean)
    elseif method === :mt
        p_tmp, s_frq, _ = @views spectrogram(obj.data[1, :, 1], fs=fs, norm=norm, mt=true, st=false, demean=demean)
    elseif method === :stft
        p_tmp, s_frq, _ = @views spectrogram(obj.data[1, :, 1], fs=fs, norm=norm, mt=false, st=true, demean=demean)
    elseif method === :mw
    _, p_tmp, _, s_frq = @views wspectrogram(obj.data[1, :, 1], pad=pad, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc, demean=demean)
    elseif method === :gh
        p_tmp, s_frq = @views ghspectrogram(obj.data[1, :, 1], fs=fs, frq_lim=frq_lim, frq_n=frq_n, norm=norm, frq=frq, demean=demean, gw=gw)
    elseif method === :cwt
        p_tmp, s_frq = @views cwtspectrogram(obj.data[1, :, 1], wt=wt, fs=fs, frq_lim=frq_lim, norm=norm, demean=demean)
    end

    s_t = linspace(0, (epoch_len(obj) / fs), size(p_tmp, 2))
    s_pow = zeros(size(p_tmp, 1), size(p_tmp, 2), ch_n, ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ch_n, 1))

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            if method === :standard
                s_pow[:, :, ch_idx, ep_idx], _, _ = @views spectrogram(obj.data[channel[ch_idx], :, ep_idx], fs=fs, norm=norm, mt=false, st=false, demean=demean)
            elseif method === :mt
                s_pow[:, :, ch_idx, ep_idx], _, _ = @views spectrogram(obj.data[channel[ch_idx], :, ep_idx], fs=fs, norm=norm, mt=true, st=false, demean=demean)
            elseif method === :stft
                s_pow[:, :, ch_idx, ep_idx], _, _ = @views spectrogram(obj.data[channel[ch_idx], :, ep_idx], fs=fs, norm=norm, mt=false, st=true, demean=demean)
            elseif method === :mw
                _, s_pow[:, :, ch_idx, ep_idx], _, _ = @views wspectrogram(obj.data[channel[ch_idx], :, ep_idx], pad=pad, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc, demean=demean)
            elseif method === :gh
                s_pow[:, :, ch_idx, ep_idx], _, _ = @views ghspectrogram(obj.data[channel[ch_idx], :, ep_idx], fs=fs, frq_lim=frq_lim, frq_n=frq_n, norm=norm, frq=frq, demean=demean, gw=gw)
            elseif method === :cwt
                s_pow[:, :, ch_idx, ep_idx], _ = @views cwtspectrogram(obj.data[channel[ch_idx], :, ep_idx], wt=wt, fs=fs, frq_lim=frq_lim, norm=norm, demean=demean)
            end

            # update progress bar
            progress_bar == true && next!(p)
        end
    end

    s_frq = round.(s_frq, digits=2)
    s_t = round.(s_t, digits=2)
    s_t .+= obj.epoch_time[1]

    return (s_pow=s_pow, s_frq=s_frq, s_t=s_t)
end

"""
    specseg(sp, sf, st; t, f)

Return spectrogram segment.

# Arguments

- `sp::Matrix{Float64}`: spectrogram powers
- `st::Vector{Float64}`: spectrogram time
- `sf::Vector{Float64}`: spectrogram frequencies
- `t::Tuple{Real, Real}`: time bounds
- `f::Tuple{Real, Real}`: frequency bounds

# Returns

Named tuple containing:
- `seg_pow::Matrix{Float64}`: powers
- `seg_shape::Shape{Real, Int64}`: shape for plotting
- `t_idx::Tuple{Real, Real}`: time indices
- `f_idx::Tuple{Real, Real}`: frequency indices
"""
function specseg(sp::Matrix{Float64}, st::Vector{Float64}, sf::Vector{Float64}; t::Tuple{Real, Real}, f::Tuple{Real, Real})

    t = tuple_order(t)
    f = tuple_order(f)

    t[1] < st[1] && throw(ArgumentError("t[1] must be ≥ $(st[1])."))
    t[2] > st[end] && throw(ArgumentError("t[2] must be ≤ $(st[end])."))
    f[1] < sf[1] && throw(ArgumentError("f[1] must be ≥ $(sf[1])."))
    f[2] > sf[end] && throw(ArgumentError("f[2] must be ≤ $(sf[end])."))

    f_idx1 = vsearch(f[1], sf)
    f_idx2 = vsearch(f[2], sf)
    t_idx1 = vsearch(t[1], st)
    t_idx2 = vsearch(t[2], st)
    seg_pow = sp[f_idx1:f_idx2, t_idx1:t_idx2]
    seg_shape = Shape([(st[t_idx1], sf[f_idx1]), (st[t_idx2], sf[f_idx1]), (st[t_idx2], sf[f_idx2]), (st[t_idx1], sf[f_idx2])])

    return (seg_pow=seg_pow, seg_shape=seg_shape, t_idx=(t_idx1,t_idx2), f_idx=(f_idx1,f_idx2))
end


"""
    specseg(sp, sf, st; t, f)

Return spectrogram segment.

# Arguments

- `sp::AbstractArray`: spectrogram powers
- `st::AbstractVector`: spectrogram time
- `sf::AbstractVector`: spectrogram frequencies
- `t::Tuple{Real, Real}`: time bounds
- `f::Tuple{Real, Real}`: frequency bounds

# Returns

Named tuple containing:
- `seg_pow::Array{Float64, 3}`: segment of powers
- `seg_shape::Shape{Real, Int64}`: segment coordinates (shape for plotting)
- `t_idx::Tuple{Real, Real}`: time indices
- `f_idx::Tuple{Real, Real}`: frequency indices
"""
function specseg(sp::AbstractArray, st::AbstractVector, sf::AbstractVector; channel::Int64, t::Tuple{Real, Real}, f::Tuple{Real, Real})

    t = tuple_order(t)
    f = tuple_order(f)

    channel < 1 && throw(ArgumentError("channel must be ≥ 1."))
    channel > size(sp, 3) && throw(ArgumentError("channel must be ≤ $(size(sp, 3))."))

    t[1] < st[1] && throw(ArgumentError("t[1] must be ≥ $(st[1])."))
    t[2] > st[end] && throw(ArgumentError("t[2] must be ≤ $(st[end])."))
    f[1] < sf[1] && throw(ArgumentError("f[1] must be ≥ $(sf[1])."))
    f[2] > sf[end] && throw(ArgumentError("f[2] must be ≤ $(sf[end])."))

    f_idx1 = vsearch(f[1], sf)
    f_idx2 = vsearch(f[2], sf)
    t_idx1 = vsearch(t[1], st)
    t_idx2 = vsearch(t[2], st)
    seg_pow = sp[f_idx1:f_idx2, t_idx1:t_idx2, channel, :]
    seg_shape = Shape([(st[t_idx1], sf[f_idx1]), (st[t_idx2], sf[f_idx1]), (st[t_idx2], sf[f_idx2]), (st[t_idx1], sf[f_idx2])])

    return (seg_pow=seg_pow, seg_shape=seg_shape, t_idx=(t_idx1,t_idx2), f_idx=(f_idx1,f_idx2))
end
