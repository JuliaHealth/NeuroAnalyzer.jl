export env_up
export env_lo
export henv_up
export henv_lo
export tenv
export tenv_mean
export tenv_median
export penv
export penv_mean
export penv_median
export senv
export senv_mean
export senv_median
export henv
export henv_mean
export henv_median
export env_cor

"""
    env_up(s, x; <keyword arguments>)

Calculate upper cubic-spline envelope from local maxima.

# Arguments

- `s::AbstractVector`: signal vector
- `x::AbstractVector`: x-axis points (e.g. time points)
- `d::Int64=32`: minimum distance between peaks in samples; smaller values give a tighter fit

# Returns

- `e::Vector{Float64}`: upper envelope (zeros if fewer than 2 peaks found)
"""
function env_up(s::AbstractVector, x::AbstractVector; d::Int64 = 32)::Vector{Float64}

    !(length(s) == length(x)) && throw(ArgumentError("Lengths of s ($(length(s))) and x ($(length(x))) must be equal."))

    e = zeros(length(s))

    # locate local maxima at least d samples apart
    p_idx = findpeaks(s, d = d)

    if length(p_idx) < 2
        _info("Envelope cannot be interpolated, less than 2 peaks detected")
    else
        # anchor the spline at the first and last sample so the envelope
        # covers the full signal range without extrapolation artifacts
        p_idx[1] != 1 && pushfirst!(p_idx, 1)
        p_idx[end] != length(s) && push!(p_idx, length(s))

        # fit a cubic spline through the peak positions and evaluate on x
        model = Spline1D(x[p_idx], s[p_idx], bc = "extrapolate")
        e = model(x)
    end

    return e

end

"""
    env_lo(s, x; <keyword arguments>)

Calculate lower cubic-spline envelope from local minima.

# Arguments

- `s::AbstractVector`: signal vector
- `x::AbstractVector`: x-axis points (e.g. time points)
- `d::Int64=32`: minimum distance between peaks in samples; smaller values give a tighter fit

# Returns

- `e::Vector{Float64}`: lower envelope (zeros if fewer than 2 troughs found)
"""
function env_lo(s::AbstractVector, x::AbstractVector; d::Int64 = 32)::Vector{Float64}

    !(length(s) == length(x)) && throw(ArgumentError("Lengths of s ($(length(s))) and x ($(length(x))) must be equal."))

    e = zeros(length(s))

    # flip the signal so minima become maxima, then find peaks.
    s_tmp = _flipx(s)
    p_idx = findpeaks(s_tmp, d = d)

    if length(p_idx) < 2
        _info("Envelope cannot be interpolated, less than 2 peaks detected")
    else
        # anchor the spline at the first and last sample so the envelope
        # covers the full signal range without extrapolation artifacts
        p_idx[1] != 1 && pushfirst!(p_idx, 1)
        p_idx[end] != length(s) && push!(p_idx, length(s))

        # fit a cubic spline through the peak positions and evaluate on x
        model = Spline1D(x[p_idx], s[p_idx], bc = "extrapolate")
        e = model(x)
    end

    return e

end

"""
    henv_up(s)

Calculate upper amplitude envelope using the Hilbert transform.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `e::Vector{Float64}`: instantaneous amplitude (upper envelope)

# Notes

The Hilbert transform works best for narrowband signals (energy concentrated around a single frequency).
"""
function henv_up(s::AbstractVector)::Vector{Float64}

    h = htransform(s)
    e = h.a

    return e

end

"""
    henv_lo(s)

Calculate lower amplitude envelope using the Hilbert transform.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `e::Vector{Float64}`: negative instantaneous amplitude (lower envelope)

# Notes

The Hilbert transform works best for narrowband signals (energy concentrated around a single frequency).
"""
function henv_lo(s::AbstractVector)::Vector{Float64}

    h = htransform(-s)
    e = h.a

    return -e

end

"""
    tenv(obj; <keyword arguments>)

Calculate temporal envelope.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `d::Int64=32`: minimum distance between peaks in samples; smaller values give a tighter fit

# Returns

Named tuple:

- `e::Array{Float64, 3}`: temporal envelope, shape `(channels, samples, epochs)`
- `t::Vector{Float64}`: time points
"""
function tenv(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    d::Int64 = 32
)::@NamedTuple{e::Array{Float64, 3}, t::Vector{Float64}}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)
    # epoch time points
    t = obj.epoch_time

    # pre-allocate output
    e = zeros(ch_n, epoch_len(obj), ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        e[ch_idx, :, ep_idx] = env_up(@view(obj.data[ch[ch_idx], :, ep_idx]), t, d = d)
    end

    return (e = e, t = t)

end

"""
    tenv_mean(obj; <keyword arguments>)

Calculate temporal envelope: mean and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `dims::Int64`: mean over channels (`dims=1`), epochs (`dims=2`), or both (`dims=3`)
- `d::Int64=32`: minimum distance between peaks in samples; smaller values give a tighter fit

# Returns

Named tuple:

- `em::Matrix{Float64}`: mean temporal envelope
- `eu::Matrix{Float64}`: 95% CI upper bound
- `el::Matrix{Float64}`: 95% CI lower bound
- `t::Vector{Float64}`: time points
"""
function tenv_mean(
    obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, dims::Int64, d::Int64 = 32
)::@NamedTuple{
    em::Matrix{Float64},
    eu::Matrix{Float64},
    el::Matrix{Float64},
    t::Vector{Float64},
}

    if dims == 1
        !(nchannels(obj) >= 2) && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        !(nepochs(obj) >= 2) && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        !(nchannels(obj) >= 2) && throw(ArgumentError("Number of channels must be ≥ 2."))
        !(nepochs(obj) >= 2) && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    tenv_data = tenv(obj, ch = ch, d = d)
    s_a = tenv_data.e
    t = tenv_data.t

    ch_n = size(s_a, 1)
    ep_n = size(s_a, 3)

    if dims == 1

        # mean over channels at each time point, per epoch

        # pre-allocate outputs
        em = zeros(length(t), ep_n)
        eu = zeros(length(t), ep_n)
        el = zeros(length(t), ep_n)

        @inbounds for ep_idx in 1:ep_n
            em[:, ep_idx] = dropdims(mean(@view(s_a[:, :, ep_idx]), dims = 1), dims = 1)
            ci = 1.96 * std(@view(em[:, ep_idx])) / sqrt(length(t))
            eu[:, ep_idx] = em[:, ep_idx] .+ ci
            el[:, ep_idx] = em[:, ep_idx] .- ci
        end

    elseif dims == 2

        # mean over epochs at each time point, per channel

        # pre-allocate outputs
        em = zeros(length(t), ch_n)
        eu = zeros(length(t), ch_n)
        el = zeros(length(t), ch_n)

        @inbounds for ch_idx in 1:ch_n
            em[:, ch_idx] = dropdims(mean(@view(s_a[ch_idx, :, :]), dims = 2), dims = 2)
            ci = 1.96 * std(@view(em[:, ch_idx])) / sqrt(length(t))
            eu[:, ch_idx] = em[:, ch_idx] .+ ci
            el[:, ch_idx] = em[:, ch_idx] .- ci
        end

    else

        # mean over channels and epochs: first average over channels (dims=1),
        # then average the result over epochs (dims=2 of the intermediate matrix)

        tenv_data = tenv_mean(obj, ch = ch, dims = 1, d = d)

        em = mean(tenv_data.em, dims = 2)
        eu = mean(tenv_data.eu, dims = 2)
        el = mean(tenv_data.el, dims = 2)

    end

    return (em = em, eu = eu, el = el, t = t)

end

"""
    tenv_median(obj; <keyword arguments>)

Calculate temporal envelope: median and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `dims::Int64`: median over channels (`dims=1`), epochs (`dims=2`), or both (`dims=3`)
- `d::Int64=32`: minimum distance between peaks in samples; smaller values give a tighter fit

# Returns

Named tuple:

- `em::Matrix{Float64}`: median temporal envelope
- `eu::Matrix{Float64}`: 95% CI upper bound
- `el::Matrix{Float64}`: 95% CI lower bound
- `t::Vector{Float64}`: time points
"""
function tenv_median(
    obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, dims::Int64, d::Int64 = 32
)::@NamedTuple{
    em::Matrix{Float64},
    eu::Matrix{Float64},
    el::Matrix{Float64},
    t::Vector{Float64},
}

    if dims == 1
        !(nchannels(obj) >= 2) && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        !(nepochs(obj) >= 2) && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        !(nchannels(obj) >= 2) && throw(ArgumentError("Number of channels must be ≥ 2."))
        !(nepochs(obj) >= 2) && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    tenv_data = tenv(obj, ch = ch, d = d)
    s_a = tenv_data.e
    t = tenv_data.t

    ch_n = size(s_a, 1)
    ep_n = size(s_a, 3)

    if dims == 1

        # median over channels at each time point, per epoch

        # pre-allocate outputs
        em = zeros(length(t), ep_n)
        eu = zeros(length(t), ep_n)
        el = zeros(length(t), ep_n)

        @inbounds for ep_idx in 1:ep_n
            em[:, ep_idx] = dropdims(median(@view(s_a[:, :, ep_idx]), dims = 1), dims = 1)
            for m_idx in eachindex(t)
                eu[m_idx, ep_idx], el[m_idx, ep_idx] = cimd(@view(s_a[:, m_idx, ep_idx]))
            end
        end

    elseif dims == 2

        # median over epochs at each time point, per channel

        # pre-allocate outputs
        em = zeros(length(t), ch_n)
        eu = zeros(length(t), ch_n)
        el = zeros(length(t), ch_n)

        @inbounds for ch_idx in 1:ch_n
            em[:, ch_idx] = dropdims(median(@view(s_a[ch_idx, :, :]), dims = 2), dims = 2)
            for m_idx in eachindex(t)
                eu[m_idx, ch_idx], el[m_idx, ch_idx] = cimd(@view(s_a[ch_idx, m_idx, :]))
            end
        end

    else

        # median over channels and epochs: first average over channels (dims=1),
        # then average the result over epochs (dims=2 of the intermediate matrix)

        tenv_data = tenv_median(obj, ch = ch, dims = 1, d = d)
        em = median(tenv_data.em, dims = 2)
        eu = median(tenv_data.eu, dims = 2)
        el = median(tenv_data.el, dims = 2)

    end

    return (em = em, eu = eu, el = el, t = t)

end

"""
    penv(obj; <keyword arguments>)

Calculate power spectrum (in dB) envelope.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `d::Int64=8`: minimum distance between peaks in samples; smaller values give a tighter fit
- `method::Symbol=:welch`: PSD method:
- `:welch`: Welch's periodogram
- `:fft`: fast Fourier transform
- `:mt`: multi-tapered periodogram
- `:stft`: short-time Fourier transform
- `:mw`: Morlet wavelet convolution
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window length in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

Named tuple:

- `e::Array{Float64, 3}`: power spectrum envelope, shape `(channels, frequencies, epochs)`
- `f::Vector{Float64}`: frequencies
"""
function penv(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    d::Int64 = 8,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    demean::Bool = true,
)::@NamedTuple{e::Array{Float64, 3}, f::Vector{Float64}}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)
    # sampling rate
    fs = sr(obj)

    # pilot call to determine the frequency vector length
    _log_off()
    psd_data = psd(
        @view(obj.data[ch[1], :, 1]),
        fs = fs,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        demean = demean,
    )
    f = psd_data.f
    _log_on()

    # pre-allocate output
    e = zeros(ch_n, length(pw), ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        psd_data = psd(
            @view(obj.data[ch[ch_idx], :, ep_idx]),
            fs = fs,
            db = true,
            method = method,
            nt = nt,
            wlen = wlen,
            woverlap = woverlap,
            w = w,
            ncyc = ncyc,
            demean = demean,
        )
        e[ch_idx, :, ep_idx] = env_up(psd_data.p, f, d = d)
    end
    _log_on()

    return (e = e, f = f)

end

"""
    penv_mean(obj; <keyword arguments>)

Calculate power spectrum (in dB) envelope: mean and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `dims::Int64`: mean over channels (`dims=1`), epochs (`dims=2`), or both (`dims=3`)
- `d::Int64=8`: minimum distance between peaks in samples; smaller values give a tighter fit
- `method::Symbol=:welch`: PSD method:
- `:welch`: Welch's periodogram
- `:fft`: fast Fourier transform
- `:mt`: multi-tapered periodogram
- `:stft`: short-time Fourier transform
- `:mw`: Morlet wavelet convolution
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window length in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

Named tuple:

- `em::Matrix{Float64}`: mean power envelope
- `eu::Matrix{Float64}`: 95% CI upper bound
- `el::Matrix{Float64}`: 95% CI lower bound
- `f::Vector{Float64}`: power spectrum envelope (useful for plotting over PSD)
"""
function penv_mean(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    dims::Int64,
    d::Int64 = 8,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    demean::Bool = true,
)::@NamedTuple{
    em::Matrix{Float64},
    eu::Matrix{Float64},
    el::Matrix{Float64},
    f::Vector{Float64},
}

    if dims == 1
        !(nchannels(obj) >= 2) && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        !(nepochs(obj) >= 2) && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        !(nchannels(obj) >= 2) && throw(ArgumentError("Number of channels must be ≥ 2."))
        !(nepochs(obj) >= 2) && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    penv_data = penv(obj, ch = ch, d = d, method = method, nt = nt, wlen = wlen,
                    woverlap = woverlap, w = w, ncyc = ncyc, demean = demean)
    pw = penv_data.e
    f  = penv_data.f

    # number of channels
    ch_n = size(pw, 1)
    # number of epochs
    ep_n = size(pw, 3)

    if dims == 1

        # mean over channels at each time point, per epoch

        # pre-allocate outputs
        em = zeros(length(f), ep_n)
        eu = zeros(length(f), ep_n)
        el = zeros(length(f), ep_n)

        @inbounds for ep_idx in 1:ep_n
            em[:, ep_idx] = dropdims(mean(@view(pw[:, :, ep_idx]), dims = 1), dims = 1)
            ci = 1.96 * std(@view(em[:, ep_idx])) / sqrt(length(f))
            eu[:, ep_idx] = em[:, ep_idx] .+ ci
            el[:, ep_idx] = em[:, ep_idx] .- ci
        end

    elseif dims == 2

        # mean over epochs at each time point, per channel

        # pre-allocate outputs
        em = zeros(length(f), ch_n)
        eu = zeros(length(f), ch_n)
        el = zeros(length(f), ch_n)

        @inbounds for ch_idx in 1:ch_n
            em[:, ch_idx] = dropdims(mean(@view(pw[ch_idx, :, :]), dims = 2), dims = 2)
            ci = 1.96 * std(@view(em[:, ch_idx])) / sqrt(length(f))
            eu[:, ch_idx] = em[:, ch_idx] .+ ci
            el[:, ch_idx] = em[:, ch_idx] .- ci
        end

    else

        # mean over channels and epochs: first average over channels (dims=1),
        # then average the result over epochs (dims=2 of the intermediate matrix)

        penv_data = penv_mean(obj,
            ch = ch,
            dims = 1,
            d = d,
            method = method,
            nt = nt,
            wlen = wlen,
            woverlap = woverlap,
            w = w,
            ncyc = ncyc,
            demean = demean
        )
        em = mean(penv_data.em, dims = 2)
        eu = mean(penv_data.eu, dims = 2)
        el = mean(penv_data.el, dims = 2)

    end

    return (em = em, eu = eu, el = el, f = f)

end

"""
    penv_median(obj; <keyword arguments>)

Calculate power spectrum (in dB) envelope: median and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `dims::Int64`: median over channels (dims = 1) or epochs (dims = 2)
- `d::Int64=8`: minimum distance between peaks in samples; smaller values give a tighter fit
- `method::Symbol=:welch`: PSD method:
- `:welch`: Welch's periodogram
- `:fft`: fast Fourier transform
- `:mt`: multi-tapered periodogram
- `:stft`: short-time Fourier transform
- `:mw`: Morlet wavelet convolution
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window length in samples
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `demean::Bool=true`: subtract DC before calculating PSD

# Returns

Named tuple:

- `em::Matrix{Float64}`: median power envelope
- `eu::Matrix{Float64}`: 95% CI upper bound
- `el::Matrix{Float64}`: power spectrum envelope: 95% CI lower bound
- `f::Vector{Float64}`: 95% CI lower bound
"""
function penv_median(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    dims::Int64,
    d::Int64 = 8,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    demean::Bool = true,
)::@NamedTuple{
    em::Matrix{Float64},
    eu::Matrix{Float64},
    el::Matrix{Float64},
    f::Vector{Float64},
}

    if dims == 1
        !(nchannels(obj) >= 2) && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        !(nepochs(obj) >= 2) && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        !(nchannels(obj) >= 2) && throw(ArgumentError("Number of channels must be ≥ 2."))
        !(nepochs(obj) >= 2) && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    penv_data = penv(obj,
        ch = ch,
        d = d,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        demean = demean,
    )
    pw = penv_data.e
    f  = penv_data.f

    # number of channels
    ch_n = size(pw, 1)
    # number of epochs
    ep_n = size(pw, 3)

    if dims == 1

        # median over channels at each time point, per epoch

        # pre-allocate outputs
        em = zeros(length(f), ep_n)
        eu = zeros(length(f), ep_n)
        el = zeros(length(f), ep_n)

        @inbounds for ep_idx in 1:ep_n
            em[:, ep_idx] = dropdims(median(@view(pw[:, :, ep_idx]), dims = 1), dims = 1)
            for m_idx in eachindex(f)
                eu[m_idx, ep_idx], el[m_idx, ep_idx] = cimd(@view(pw[:, m_idx, ep_idx]))
            end
        end

    elseif dims == 2

        # mean over epochs at each time point, per channel

        # pre-allocate outputs
        em = zeros(length(f), ch_n)
        eu = zeros(length(f), ch_n)
        el = zeros(length(f), ch_n)

        @inbounds for ch_idx in 1:ch_n
            em[:, ch_idx] = dropdims(median(@view(pw[ch_idx, :, :]), dims = 2), dims = 2)
            for m_idx in eachindex(f)
                eu[m_idx, ch_idx], el[m_idx, ch_idx] = cimd(@view(pw[ch_idx, m_idx, :]))
            end
        end

    else

        # mean over channels and epochs: first average over channels (dims=1),
        # then average the result over epochs (dims=2 of the intermediate matrix)

        penv_data = penv_median(obj,
            ch = ch,
            dims = 1,
            d = d,
            method = method,
            nt = nt,
            wlen = wlen,
            woverlap = woverlap,
            w = w,
            ncyc = ncyc,
            demean = demean
        )
        em = median(penv_data.em, dims = 2)
        eu = median(penv_data.eu, dims = 2)
        el = median(penv_data.el, dims = 2)

    end

    return (em = em, eu = eu, el = el, f = f)

end

"""
    senv(obj; <keyword arguments>)

Calculate spectral envelope (dominant frequency over time).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `d::Int64=2`: minimum distance between peaks in samples; smaller values give a tighter fit
- `t::Union{Real, Nothing}=nothing`: spectrogram threshold; powers above `t` are zeroed before finding the dominant frequency
- `method::Symbol=:stft` spectrogram method:
- `:stft`: short-time Fourier transform
- `:mt`: multi-tapered periodogram
- `:mw`: Morlet wavelet convolution
- `:gh`: Gaussian and Hilbert transform
- `:cwt`: continuous wavelet transformation
- `pad::Int64=0`: number of zeros to append
- `db::Bool=true`: normalize powers to dB
- `nt::Int64=7`: number of Slepian tapers
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window length in samples
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple:

- `e::Array{Float64, 3}`: spectral envelope, shape `(channels, samples, epochs)`
- `t::Vector{Float64}`: spectrogram time
"""
function senv(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    d::Int64 = 2,
    t::Union{Real, Nothing} = nothing,
    pad::Int64 = 0,
    method::Symbol = :stft,
    db::Bool = true,
    nt::Int64 = 7,
    gw::Real = 5,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    wt::T = wavelet(Morlet(2π), β = 2),
    wlen::Int64 = sr(obj),
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
)::@NamedTuple{e::Array{Float64, 3}, t::Vector{Float64}} where {T <: CWT}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)
    # sampling rate
    fs = sr(obj)

    # pilot call to determine the spectrogram time vector 
    if method === :stft
        spec_data = NeuroAnalyzer.spectrogram(
            @view(obj.data[ch[1], :, 1]),
            fs = fs,
            db = db,
            method = :stft,
            wlen = wlen,
            woverlap = woverlap,
            w = w,
        )
        sp = spec_data.p
    elseif method === :mt
        spec_data = NeuroAnalyzer.spectrogram(
            @view(obj.data[ch[1], :, 1]),
            fs = fs,
            db = db,
            method = :mt,
            nt = nt,
            wlen = wlen,
            woverlap = woverlap,
            w = w,
        )
        sp = spec_data.p
    elseif method === :mw
        spec_data = NeuroAnalyzer.mwspectrogram(
            @view(obj.data[ch[1], :, 1]),
            pad = pad,
            fs = fs,
            db = db,
            ncyc = ncyc,
            w = w,
        )
        sp = spec_data.p
    elseif method === :gh
        spec_data = NeuroAnalyzer.ghtspectrogram(
            @view(obj.data[ch[1], :, 1]),
            fs = fs,
            db = db,
            gw = gw,
            w = w,
        )
        sp = spec_data.p
    elseif method === :cwt
        _log_off()
        spec_data = NeuroAnalyzer.cwtspectrogram(
            @view(obj.data[ch[1], :, 1]),
            wt = wt,
            fs = fs,
        )
        _log_on()
        sp = spec_data.m
    end

    # build the spectrogram time axis and align with the epoch start
    st = linspace(0, (epoch_len(obj) / fs), size(sp, 2))
    st .+= obj.epoch_time[1]

    # pre-allocate output
    e = zeros(ch_n, length(st), ep_n)

    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        # compute spectrogram for this channel/epoch
        local sp_loc, sf_loc
        if method === :stft
            spec_data = NeuroAnalyzer.spectrogram(
                @view(obj.data[ch[ch_idx], :, ep_idx]), fs = fs, db = db,
                method = :stft, wlen = wlen, woverlap = woverlap, w = w,
            )
            sp_loc = spec_data.p
            sp_loc = spec_data.f
        elseif method === :mt
            spec_data = NeuroAnalyzer.spectrogram(
                @view(obj.data[ch[ch_idx], :, ep_idx]), fs = fs, db = db,
                method = :mt, nt = nt, wlen = wlen, woverlap = woverlap, w = w,
            )
            sp_loc = spec_data.p
            sp_loc = spec_data.f
        elseif method === :mw
            spec_data = NeuroAnalyzer.mwspectrogram(
                @view(obj.data[ch[ch_idx], :, ep_idx]),
                pad = pad, fs = fs, db = db, ncyc = ncyc, w = w,
            )
            sp_loc = spec_data.p
            sp_loc = spec_data.f
        elseif method === :gh
            spec_data = NeuroAnalyzer.ghtspectrogram(
                @view(obj.data[ch[ch_idx], :, ep_idx]), fs = fs, db = db, gw = gw, w = w,
            )
            sp_loc = spec_data.p
            sp_loc = spec_data.f
        elseif method === :cwt
            _log_off()
            spec_data = NeuroAnalyzer.cwtspectrogram(
                @view(obj.data[ch[ch_idx], :, ep_idx]), wt = wt, fs = fs,
            )
            sp_loc = spec_data.m
            sp_loc = spec_data.f
            _log_on()
        end

        # optionally zero out powers above the threshold, then reverse so the
        # highest sub-threshold power becomes the "dominant" frequency
        if t !== nothing
            sp_loc[sp_loc .> t] .= 0
            reverse!(sp_loc)
            reverse!(sf_loc)
        end

        # for each time bin find the frequency with maximum power
        f_idx = zeros(length(st))
        m = vec(maximum(sp_loc, dims = 1))
        for idx2 in eachindex(m)
            f_idx[idx2] = sf_loc[vsearch(m[idx2], @view(sp_loc[:, idx2]))]
        end

        e[ch_idx, :, ep_idx] = env_up(f_idx, st, d = d)
    end

    return (e = e, t = st)

end

"""
    senv_mean(obj; <keyword arguments>)

Calculate spectral envelope: mean and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `dims::Int64`: mean over channels (`dims=1`), epochs (`dims=2`), or both (`dims=3`)
- `d::Int64=2`: minimum distance between peaks in samples; smaller values give a tighter fit
- `t::Union{Real, Nothing}=nothing`: spectrogram threshold; powers above `t` are zeroed before finding the dominant frequency
- `method::Symbol=:stft` spectrogram method:
- `:stft`: short-time Fourier transform
- `:mt`: multi-tapered periodogram
- `:mw`: Morlet wavelet convolution
- `:gh`: Gaussian and Hilbert transform
- `:cwt`: continuous wavelet transformation
- `pad::Int64=0`: number of zeros to append
- `db::Bool=true`: normalize powers to dB
- `nt::Int64=7`: number of Slepian tapers
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window length in samples
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple:

- `em::Matrix{Float64}`: spectral envelope: mean
- `eu::Matrix{Float64}`: spectral envelope: 95% CI upper bound
- `el::Matrix{Float64}`: spectral envelope: 95% CI lower bound
- `t::Vector{Float64}`: spectral envelope (useful for plotting over spectrogram)
"""
function senv_mean(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    dims::Int64,
    d::Int64 = 2,
    t::Union{Real, Nothing} = nothing,
    method::Symbol = :stft,
    pad::Int64 = 0,
    db::Bool = true,
    nt::Int64 = 7,
    gw::Real = 5,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    wt::T = wavelet(Morlet(2π), β = 2),
    wlen::Int64 = sr(obj),
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
)::@NamedTuple{
    em::Matrix{Float64},
    eu::Matrix{Float64},
    el::Matrix{Float64},
    t::Vector{Float64},
} where {T <: CWT}

    if dims == 1
        !(nchannels(obj) >= 2) && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        !(nepochs(obj) >= 2) && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        !(nchannels(obj) >= 2) && throw(ArgumentError("Number of channels must be ≥ 2."))
        !(nepochs(obj) >= 2) && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    env_data = senv(
        obj,
        ch = ch,
        d = d,
        t = t,
        pad = pad,
        method = method,
        db = db,
        nt = nt,
        gw = gw,
        ncyc = ncyc,
        wt = wt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
    )
    sp = env_data.e
    st = env_data.t

    # number of channels
    ch_n = size(pw, 1)
    # number of epochs
    ep_n = size(pw, 3)

    if dims == 1

        # mean over channels at each time point, per epoch

        # pre-allocate outputs
        em = zeros(length(st), ep_n)
        eu = zeros(length(st), ep_n)
        el = zeros(length(st), ep_n)

        @inbounds for ep_idx in 1:ep_n
            em[:, ep_idx] = dropdims(mean(@view(sp[:, :, ep_idx]), dims = 1), dims = 1)
            ci = 1.96 * std(@view(em[:, ep_idx])) / sqrt(length(st))
            eu[:, ep_idx] = em[:, ep_idx] .+ ci
            el[:, ep_idx] = em[:, ep_idx] .- ci
        end

    elseif dims == 2

        # mean over epochs at each time point, per channel

        # pre-allocate outputs
        em = zeros(length(st), ch_n)
        eu = zeros(length(st), ch_n)
        el = zeros(length(st), ch_n)

        @inbounds for ch_idx in 1:ch_n
            em[:, ch_idx] = dropdims(mean(@view(sp[ch_idx, :, :]), dims = 2), dims = 2)
            ci = 1.96 * std(@view(em[:, ch_idx])) / sqrt(length(st))
            eu[:, ch_idx] = em[:, ch_idx] .+ ci
            el[:, ch_idx] = em[:, ch_idx] .- ci
        end

    else

        # mean over channels and epochs: first average over channels (dims=1),
        # then average the result over epochs (dims=2 of the intermediate matrix)

        env_data = senv_mean(
            obj,
            ch = ch,
            dims = 1,
            d = d,
            pad = pad,
            method = method,
            db = db,
            nt = nt,
            gw = gw,
            ncyc = ncyc,
            wt = wt,
            wlen = wlen,
            woverlap = woverlap,
            w = w,
        )
        em = mean(env_data.em, dims = 2)
        eu = mean(env_data.eu, dims = 2)
        el = mean(env_data.el, dims = 2)

    end

    return (em = em, eu = eu, el = el, t = st)

end

"""
    senv_median(obj; <keyword arguments>)

Calculate spectral envelope: median and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `dims::Int64`: median over channels (`dims=1`), epochs (`dims=2`), or both (`dims=3`)
- `d::Int64=2`: minimum distance between peaks in samples; smaller values give a tighter fit
- `t::Union{Real, Nothing}=nothing`: spectrogram threshold; powers above `t` are zeroed before finding the dominant frequency
- `method::Symbol=:stft` spectrogram method:
- `:stft`: short-time Fourier transform
- `:mt`: multi-tapered periodogram
- `:mw`: Morlet wavelet convolution
- `:gh`: Gaussian and Hilbert transform
- `:cwt`: continuous wavelet transformation
- `pad::Int64=0`: number of zeros to append
- `db::Bool=true`: normalize powers to dB
- `nt::Int64=7`: number of Slepian tapers
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles; for a tuple, cycles vary per frequency: `ncyc = linspace(ncyc[1], ncyc[2], nfrq)`
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `wlen::Int64=sr(obj)`: window length in samples (default is 1 second)
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window length in samples
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple:

- `em::Matrix{Float64}`: median spectral envelope
- `eu::Matrix{Float64}`: 95% CI upper bound
- `el::Matrix{Float64}`: 95% CI lower bound
- `t::Vector{Float64}`: spectrogram time points
"""
function senv_median(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    dims::Int64,
    d::Int64 = 2,
    t::Union{Real, Nothing} = nothing,
    method::Symbol = :stft,
    pad::Int64 = 0,
    db::Bool = true,
    nt::Int64 = 7,
    gw::Real = 5,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    wt::T = wavelet(Morlet(2π), β = 2),
    wlen::Int64 = sr(obj),
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
)::@NamedTuple{
    em::Matrix{Float64},
    eu::Matrix{Float64},
    el::Matrix{Float64},
    t::Vector{Float64},
} where {T <: CWT}

    if dims == 1
        !(nchannels(obj) >= 2) && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        !(nepochs(obj) >= 2) && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        !(nchannels(obj) >= 2) && throw(ArgumentError("Number of channels must be ≥ 2."))
        !(nepochs(obj) >= 2) && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    senv_data = senv(obj,
        ch = ch,
        d = d,
        t = t,
        pad = pad,
        method = method,
        nt = nt,
        db = db,
        gw = gw,
        ncyc = ncyc,
        wt = wt,
        wlen = wlen,
        woverlap = woverlap,
        w = w
    )
    sp = senv_data.e
    st = senv_data.t

    # number of channels
    ch_n = size(pw, 1)
    # number of epochs
    ep_n = size(pw, 3)

    if dims == 1

        # median over epochs at each time point, per channel

        # pre-allocate outputs
        em = zeros(length(st), ep_n)
        eu = zeros(length(st), ep_n)
        el = zeros(length(st), ep_n)

        @inbounds for ep_idx in 1:ep_n
            em[:, ep_idx] = @views median(sp[:, :, ep_idx], dims = 1)
            for m_idx in eachindex(st)
                eu[m_idx, ep_idx], el[m_idx, ep_idx] = cimd(sp[:, m_idx, ep_idx])
            end
        end
    elseif dims == 2

        # median over epochs at each time point, per channel

        # pre-allocate outputs
        em = zeros(length(st), ch_n)
        eu = zeros(length(st), ch_n)
        el = zeros(length(st), ch_n)

        @inbounds for ch_idx in 1:ch_n
            em[:, ch_idx] = dropdims(median(@view(sp[ch_idx, :, :]), dims = 2), dims = 2)
            for m_idx in eachindex(st)
                # BUG FIX: was `cimd(sp[ch_idx, :, :])` — passed the entire
                # (time, epochs) matrix instead of the epoch vector at m_idx.
                eu[m_idx, ch_idx], el[m_idx, ch_idx] = cimd(@view(sp[ch_idx, m_idx, :]))
            end
        end

    else

        # median over channels and epochs: first average over channels (dims=1),
        # then average the result over epochs (dims=2 of the intermediate matrix)

        senv_data = senv_median(obj,
            ch = ch,
            dims = 1,
            d = d,
            pad = pad,
            method = method,
            db = db,
            nt = nt,
            gw = gw,
            ncyc = ncyc,
            wt = wt,
            wlen = wlen,
            woverlap = woverlap,
            w = w
        )

        em = median(senv_data.em, dims = 2)
        eu = median(senv_data.eu, dims = 2)
        el = median(senv_data.el, dims = 2)

    end

    return (em = em, eu = eu, el = el, t = st)

end

"""
    henv(obj; <keyword arguments>)

Calculate Hilbert spectrum amplitude envelope.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `d::Int64=32`: minimum distance between peaks in samples; smaller values give a tighter fit

# Returns

Named tuple:

- `e::Array{Float64, 3}`: Hilbert amplitude envelope, shape `(channels, samples, epochs)`
- `t::Vector{Float64}`: time points
"""
function henv(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    d::Int64 = 32
)::@NamedTuple{e::Array{Float64, 3}, t::Vector{Float64}}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")
    _warn("henv() uses Hilbert transform, the signal should be narrowband for best results.")

    henv_data = htransform(@view(obj.data[ch, :, :]))
    a = henv_data.a

    # number of channels
    ch_n = size(hamp, 1)
    # number of epochs
    ep_n = size(hamp, 3)

    # pre-allocate output
    e = similar(hamp, Float64)

    # epoch time points
    t = obj.epoch_time

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        e[ch_idx, :, ep_idx] = env_up(@view(hamp[ch_idx, :, ep_idx]), t, d = d)
    end

    return (e = e, t = t)

end

"""
    henv_mean(obj; <keyword arguments>)

Calculate Hilbert spectrum amplitude envelope: mean and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `dims::Int64`: mean over channels (`dims=1`), epochs (`dims=2`), or both (`dims=3`)
- `d::Int64=32`: minimum distance between peaks in samples; smaller values give a tighter fit

# Returns

Named tuple:

- `em::Matrix{Float64}`: mean Hilbert envelope
- `eu::Matrix{Float64}`: Hilbert spectrum amplitude envelope: 95% CI upper bound
- `el::Matrix{Float64}`: Hilbert spectrum amplitude envelope: 95% CI lower bound
- `t::Vector{Float64}`: time points
"""
function henv_mean(
        obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, dims::Int64, d::Int64 = 32
    )::@NamedTuple{
        em::Matrix{Float64},
        eu::Matrix{Float64},
        el::Matrix{Float64},
        t::Vector{Float64},
    }

    if dims == 1
        !(nchannels(obj) >= 2) && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        !(nepochs(obj) >= 2) && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        !(nchannels(obj) >= 2) && throw(ArgumentError("Number of channels must be ≥ 2."))
        !(nepochs(obj) >= 2) && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    henv_data = henv(obj, ch = ch, d = d)
    s_a = henv_data.a
    t = henv_data.t

    # number of channels
    ch_n = size(s_a, 1)
    # number of epochs
    ep_n = size(s_a, 3)

    if dims == 1

        # mean over channels at each time point, per epoch

        # pre-allocate outputs
        em = zeros(length(t), ep_n)
        eu = zeros(length(t), ep_n)
        el = zeros(length(t), ep_n)

        @inbounds for ep_idx in 1:ep_n
            em[:, ep_idx] = dropdims(mean(@view(s_a[:, :, ep_idx]), dims = 1), dims = 1)
            ci = 1.96 * std(@view(em[:, ep_idx])) / sqrt(length(t))
            eu[:, ep_idx] = em[:, ep_idx] .+ ci
            el[:, ep_idx] = em[:, ep_idx] .- ci
        end

    elseif dims == 2

        # mean over epochs

        em = zeros(length(t), ch_n)
        eu = zeros(length(t), ch_n)
        el = zeros(length(t), ch_n)

        @inbounds for ch_idx in 1:ch_n
            em[:, ch_idx] = dropdims(mean(@view(s_a[ch_idx, :, :]), dims = 2), dims = 2)
            ci = 1.96 * std(@view(em[:, ch_idx])) / sqrt(length(t))
            eu[:, ch_idx] = em[:, ch_idx] .+ ci
            el[:, ch_idx] = em[:, ch_idx] .- ci
        end

    else

        # mean over channels and epochs: first average over channels (dims=1),
        # then average the result over epochs (dims=2 of the intermediate matrix)

        henv_data = henv_mean(obj, ch = ch, dims = 1, d = d)
        em = vec(mean(henv_data.em, dims = 2))
        eu = vec(mean(henv_data.eu, dims = 2))
        el = vec(mean(henv_data.el, dims = 2))

    end

    return (em = em, eu = eu, el = el, t = t)

end

"""
    henv_median(obj; <keyword arguments>)

Calculate Hilbert spectrum amplitude envelope of `obj`: median and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `dims::Int64`: median over channels (`dims=1`), epochs (`dims=2`), or both (`dims=3`)
- `d::Int64=32`: minimum distance between peaks in samples; smaller values give a tighter fit

# Returns

Named tuple:

- `em::Matrix{Float64}`: Hilbert spectrum amplitude envelope: median
- `eu::Matrix{Float64}`: Hilbert spectrum amplitude envelope: 95% CI upper bound
- `el::Matrix{Float64}`: Hilbert spectrum amplitude envelope: 95% CI lower bound
- `t::Vector{Float64}`: time points
"""
function henv_median(
    obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, dims::Int64, d::Int64 = 32
)::@NamedTuple{
    em::Matrix{Float64},
    eu::Matrix{Float64},
    el::Matrix{Float64},
    t::Vector{Float64},
}

    if dims == 1
        !(nchannels(obj) >= 1) && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        !(nepochs(obj) >= 1) && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        !(nchannels(obj) >= 1) && throw(ArgumentError("Number of channels must be ≥ 2."))
        !(nepochs(obj) >= 1) && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_a, t = henv(obj, ch = ch, d = d)

    ch_n = size(s_a, 1)
    ep_n = size(s_a, 3)

    if dims == 1
        # median over channels

        em = zeros(length(t), ep_n)
        eu = zeros(length(t), ep_n)
        el = zeros(length(t), ep_n)

        @inbounds for ep_idx in 1:ep_n
            em[:, ep_idx] = @views median(s_a[:, :, ep_idx], dims = 1)
            for m_idx in eachindex(t)
                eu[m_idx, ep_idx], el[m_idx, ep_idx] = cimd(s_a[:, m_idx, ep_idx])
            end
        end
    elseif dims == 2
        # median over epochs

        em = zeros(length(t), ch_n)
        eu = zeros(length(t), ch_n)
        el = zeros(length(t), ch_n)

        @inbounds for ch_idx in 1:ch_n
            em[:, ch_idx] = median(s_a[ch_idx, :, :], dims = 2)
            for m_idx in eachindex(t)
                eu[m_idx, ch_idx], el[m_idx, ch_idx] = cimd(s_a[ch_idx, m_idx, :])
            end
        end
    else
        # median over channels and epochs

        em, eu, el, _ = henv_median(obj, ch = ch, dims = 1, d = d)
        em = median(em, dims = 2)
        eu = median(eu, dims = 2)
        el = median(el, dims = 2)
        em = reshape(em, size(em, 1))
        eu = reshape(eu, size(eu, 1))
        el = reshape(el, size(el, 1))
    end

    return (em = em, eu = eu, el = el, t = t)
end

"""
    env_cor(env1, env2)

Calculate envelope correlation.

# Arguments

- `env1::Array{Float64, 3}`
- `env2::Array{Float64, 3}`

# Returns

Named tuple:

- `ec::Vector{Float64}`: envelope correlation coefficient
- `p::Vector{Float64}`: p-value
"""
function env_cor(env1::Array{Float64, 3}, env2::Array{Float64, 3})::@NamedTuple{ec::Vector{Float64}, p::Vector{Float64}}

    !(size(env1) == size(env2)) && throw(ArgumentError("Both envelopes must have the same size."))

    ep_n = size(env1, 3)
    ec = zeros(ep_n)
    p = zeros(ep_n)

    # compare envelopes per epochs
    for ep_idx in 1:ep_n
        ctest = @views CorrelationTest(vec(env1[:, :, ep_idx]), vec(env2[:, :, ep_idx]))
        @inbounds ec[ep_idx] = ctest.r
        @inbounds p[ep_idx] = pvalue(ctest)
    end

    return (ec = ec, p = p)

end
