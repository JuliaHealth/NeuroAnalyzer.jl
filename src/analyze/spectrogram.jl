export spectrogram
export wspectrogram
export ghspectrogram
export cwtspectrogram

"""
    spectrogram(s; fs, norm, mt, st)

Calculate spectrogram.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling frequency
- `norm::Bool=true`: normalize powers to dB
- `mt::Bool=false`: if true use multi-tapered spectrogram
- `st::Bool=false`: if true use short time Fourier transform

# Returns

Named tuple containing:
- `sp::Matrix{Float64}`: powers
- `sf::Vector{Float64}`: frequencies
- `st::Vector{Float64}`: time
"""
function spectrogram(s::AbstractVector; fs::Int64, norm::Bool=true, mt::Bool=false, st::Bool=false)

    (mt == true && st == true) && throw(ArgumentError("Both mt and st must not be true."))
    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))

    nfft = length(s)
    interval = fs
    overlap = round(Int64, fs * 0.75)

    if st == true
        sp = abs.(stft(s, interval, overlap, nfft=nfft, fs=fs, window=hanning))
        norm == true && (sp = pow2db.(sp))
        t = 0:1/fs:(length(s) / fs)
        st = linspace(t[1], t[end], size(sp, 2))
        sf = linspace(0, fs/2, size(sp, 1))
        return (sp=sp, sf=sf, st=st)
    end

    if mt == true
        sp = DSP.mt_spectrogram(s, fs=fs)
    else
        sp = DSP.spectrogram(s, interval, overlap, nfft=nfft, fs=fs, window=hanning)
    end
    sp = sp.power
    norm == true && (sp = pow2db.(sp))

    #st = collect(spec.time)
    #sf = Vector(spec.freq)
    t = 0:1/fs:(length(s) / fs)
    st = linspace(t[1], t[end], size(sp, 2))
    sf = linspace(0, fs/2, size(sp, 1))
    
    return (sp=sp, sf=sf, st=st)

end

"""
    wspectrogram(s; pad, norm, fs, frq_lim, frq_n, frq, ncyc)

Calculate spectrogram using wavelet convolution.

# Arguments

- `s::AbstractVector`
- `pad::Int64`: pad with `pad` zeros
- `norm::Bool=true`: normalize powers to dB
- `fs::Int64`: sampling rate
- `frq_lim::Tuple{Real, Real}=(0, fs ÷ 2)`: frequency bounds for the spectrogram
- `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

# Returns

Named tuple containing:
- `cs::Matrix(ComplexF64}`: convoluted signal
- `sp::Matrix{Float64}`: powers
- `sph::Matrix{Float64}`: phases
- `sf::Vector{Float64}`: frequencies
"""
function wspectrogram(s::AbstractVector; pad::Int64=0, norm::Bool=true, fs::Int64, frq_lim::Tuple{Real, Real}=(0, fs ÷ 2), frq_n::Int64=_tlength(frq_lim), frq::Symbol=:lin, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)

    _check_var(frq, [:log, :lin], "frq")

    pad > 0 && (s = pad0(s, pad))

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
        sf = round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), frq_n), digits=1)
    else
        sf = linspace(frq_lim[1], frq_lim[2], frq_n)
    end

    cs = zeros(ComplexF64, length(sf), length(s))
    sp = zeros(length(sf), length(s))
    # w_amp = zeros(length(sf), length(s))
    sph = zeros(length(sf), length(s))

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
        kernel = generate_morlet(fs, sf[frq_idx], 1, ncyc=ncyc[frq_idx], complex=true)
        cs[frq_idx, :] = fconv(s, kernel=kernel, norm=false)
        # alternative: w_amp[frq_idx, :] = LinearAlgebra.norm.(real.(cs), imag.(cs), 2)
        sp[frq_idx, :] = @views @. abs(cs[frq_idx, :])^2
        sph[frq_idx, :] = @views @. angle(cs[frq_idx, :])
    end

    # remove reflected part of the signal
    # if reflect == true
    #     cs = cs[:, (length(s) ÷ 3 + 1):(2 * length(s) ÷ 3)]
    #     w_amp = w_amp[:, (length(s) ÷ 3 + 1):(2 * length(s) ÷ 3)]
    #     sp = sp[:, (length(s) ÷ 3 + 1):(2 * length(s) ÷ 3)]
    #     sph = sph[:, (length(s) ÷ 3 + 1):(2 * length(s) ÷ 3)]
    # end

    norm == true && (sp = pow2db.(sp))

    return (cs=cs, sp=sp, sph=sph, sf=sf)
    
end

"""
    ghspectrogram(s; fs, norm, frq_lim, frq_n, frq, fs)

Calculate spectrogram using Gaussian and Hilbert transform.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `norm::Bool=true`: normalize powers to dB
- `frq_lim::Tuple{Real, Real}=(0, fs ÷ 2)`: frequency bounds for the spectrogram
- `frq_n::Int64_tlength(frq_lim)`: number of frequencies
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `gw::Real=5`: Gaussian width in Hz

# Returns

Named tuple containing:
- `sp::Matrix{Float64}`: powers
- `sph::Matrix{Float64}`: phases
- `sf::Vector{Float64}`: frequencies
"""
function ghspectrogram(s::AbstractVector; fs::Int64, norm::Bool=true, frq_lim::Tuple{Real, Real}=(0, fs ÷ 2), frq_n::Int64=_tlength(frq_lim), frq::Symbol=:lin, gw::Real=5)

    _check_var(frq, [:log, :lin], "frq")

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    frq_lim = tuple_order(frq_lim)
    frq_lim[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
    frq_lim[2] > fs ÷ 2 && throw(ArgumentError("Upper frequency bound must be ≤ $(fs ÷ 2)."))
    frq_n < 2 && throw(ArgumentError("frq_n frequency bound must be ≥ 2."))
    frq_lim[1] == 0 && (frq_lim = (0.1, frq_lim[2]))

    if frq === :log
        frq_lim = (frq_lim[1], frq_lim[2])
        sf = round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), frq_n), digits=1)
    else
        sf = linspace(frq_lim[1], frq_lim[2], frq_n)
    end

    sp = zeros(length(sf), length(s))
    sph = zeros(length(sf), length(s))

    @inbounds @simd for frq_idx in 1:length(sf)
        s = gfilter(s, fs=fs, f=sf[frq_idx], gw=gw)
        sp[frq_idx, :] = abs.(hilbert(s)).^2
        sph[frq_idx, :] = angle.(hilbert(s))
    end
    norm == true && (sp = pow2db.(sp))

    return (sp=sp, sph=sph, sf=sf)

end

"""
    cwtspectrogram(s; wt, pad, norm, frq_lim, fs)

Calculate spectrogram using continuous wavelet transformation (CWT).

# Arguments

- `s::AbstractVector`
- `wt<:CWT`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `fs::Int64`: sampling rate
- `norm::Bool=true`: normalize powers to dB
- `frq_lim::Tuple{Real, Real}=(0, fs ÷ 2)`: frequency bounds for the spectrogram

# Returns

Named tuple containing:
- `sp::Matrix{Float64}`: powers
- `sf::Vector{Float64}`: frequencies
"""
function cwtspectrogram(s::AbstractVector; wt::T, fs::Int64, norm::Bool=true, frq_lim::Tuple{Real, Real}=(0, fs ÷ 2)) where {T <: CWT}

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    frq_lim = tuple_order(frq_lim)
    frq_lim[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
    frq_lim[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be ≤ $(fs / 2)."))

    sp = abs.(ContinuousWavelets.cwt(s, wt)')
    sf = ContinuousWavelets.getMeanFreq(ContinuousWavelets.computeWavelets(length(s), wt)[1])
    sf[1] = 0
    sf[1] < frq_lim[1] && throw(ArgumentError("Lower frequency bound must be ≥ $(sf[1])."))
    sf[end] < frq_lim[2] && throw(ArgumentError("Upper frequency bound must be ≤ $(sf[end])."))
    sf = sf[vsearch(frq_lim[1], sf):vsearch(frq_lim[2], sf)]
    sp = sp[vsearch(frq_lim[1], sf):vsearch(frq_lim[2], sf), :]

    norm == true && (sp = pow2db.(sp))

    return (sp=sp, sf=sf)

end

"""
    spectrogram(obj; ch, norm, mt, st)

Calculate spectrogram.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `method::Symbol=:standard`: method of calculating spectrogram:
    - `:standard`: standard
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `pad::Int64=0`: number of zeros to add
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) ÷ 2)`: frequency limits
- `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
- `norm::Bool=true`: normalize powers to dB
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: `ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc = linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`
- `wt<:CWT=wavelet(Morlet(π), β=2)`: continuous wavelet, e.g. `wt = wavelet(Morlet(π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

Named tuple containing:
- `sp::Array{Float64, 3}`
- `sf::Vector{Float64}`
- `st::Vector{Float64}`
"""
function spectrogram(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), pad::Int64=0, frq_lim::Tuple{Real, Real}=(0, sr(obj) ÷ 2), frq_n::Int64=_tlength(frq_lim), method::Symbol=:standard, norm::Bool=true, frq::Symbol=:log, gw::Real=5, ncyc::Union{Int64, Tuple{Int64, Int64}}=6, wt::T=wavelet(Morlet(π), β=2)) where {T <: CWT}

    _check_var(method, [:standard, :stft, :mt, :mw, :gh, :cwt], "method")
    _check_channels(obj, ch)

    ch_n = length(ch)
    ep_n = epoch_n(obj)
    fs = sr(obj)

    if method === :standard
        p_tmp, sf, _ = @views NeuroAnalyzer.spectrogram(obj.data[1, :, 1], fs=fs, norm=norm, mt=false, st=false)
    elseif method === :mt
        p_tmp, sf, _ = @views NeuroAnalyzer.spectrogram(obj.data[1, :, 1], fs=fs, norm=norm, mt=true, st=false)
    elseif method === :stft
        p_tmp, sf, _ = @views NeuroAnalyzer.spectrogram(obj.data[1, :, 1], fs=fs, norm=norm, mt=false, st=true)
    elseif method === :mw
    _, p_tmp, _, sf = @views NeuroAnalyzer.wspectrogram(obj.data[1, :, 1], pad=pad, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)
    elseif method === :gh
        p_tmp, _, sf = @views NeuroAnalyzer.ghspectrogram(obj.data[1, :, 1], fs=fs, frq_lim=frq_lim, frq_n=frq_n, norm=norm, frq=frq, gw=gw)
    elseif method === :cwt
        p_tmp, sf = @views NeuroAnalyzer.cwtspectrogram(obj.data[1, :, 1], wt=wt, fs=fs, frq_lim=frq_lim, norm=norm)
    end

    st = linspace(0, (epoch_len(obj) / fs), size(p_tmp, 2))
    sp = zeros(size(p_tmp, 1), size(p_tmp, 2), ch_n, ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ch_n, 1))

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            if method === :standard
                sp[:, :, ch_idx, ep_idx], _, _ = @views NeuroAnalyzer.spectrogram(obj.data[ch[ch_idx], :, ep_idx], fs=fs, norm=norm, mt=false, st=false)
            elseif method === :mt
                sp[:, :, ch_idx, ep_idx], _, _ = @views NeuroAnalyzer.spectrogram(obj.data[ch[ch_idx], :, ep_idx], fs=fs, norm=norm, mt=true, st=false)
            elseif method === :stft
                sp[:, :, ch_idx, ep_idx], _, _ = @views NeuroAnalyzer.spectrogram(obj.data[ch[ch_idx], :, ep_idx], fs=fs, norm=norm, mt=false, st=true)
            elseif method === :mw
                _, sp[:, :, ch_idx, ep_idx], _, _ = @views NeuroAnalyzer.wspectrogram(obj.data[ch[ch_idx], :, ep_idx], pad=pad, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)
            elseif method === :gh
                sp[:, :, ch_idx, ep_idx], _, _ = @views NeuroAnalyzer.ghspectrogram(obj.data[ch[ch_idx], :, ep_idx], fs=fs, frq_lim=frq_lim, frq_n=frq_n, norm=norm, frq=frq, gw=gw)
            elseif method === :cwt
                sp[:, :, ch_idx, ep_idx], _ = @views NeuroAnalyzer.cwtspectrogram(obj.data[ch[ch_idx], :, ep_idx], wt=wt, fs=fs, frq_lim=frq_lim, norm=norm)
            end

            # update progress bar
            progress_bar == true && next!(p)
        end
    end

    sf = round.(sf, digits=2)
    st = round.(st, digits=2)
    st .+= obj.epoch_time[1]

    return (sp=sp, sf=sf, st=st)

end
