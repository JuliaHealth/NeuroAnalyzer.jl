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

# mean (`med = false`) or median (`med = true`) and 95% CI of a (channels, points, epochs)
# array at each point: over channels for each epoch (`dims = 1`, output (points, epochs)),
# over epochs for each channel (`dims = 2`, output (points, channels)) or over channels and
# epochs (`dims = 3`, output (points, 1)); mean: normal approximation, median: `cimd()`
function _env_ci(a::AbstractArray{<:Real, 3}, dims::Int64, med::Bool)
    dims in (1, 2, 3) || throw(ArgumentError("dims must be 1, 2 or 3."))
    # values (observations, points, groups)
    x = if dims == 1
        a
    elseif dims == 2
        permutedims(a, (3, 2, 1))
    else
        reshape(permutedims(a, (1, 3, 2)), size(a, 1) * size(a, 3), size(a, 2), 1)
    end
    n_obs, n_pts, n_grp = size(x)
    em = zeros(n_pts, n_grp)
    el = zeros(n_pts, n_grp)
    eu = zeros(n_pts, n_grp)
    @inbounds for g = 1:n_grp, pt = 1:n_pts
        v = @view(x[:, pt, g])
        if med
            em[pt, g] = median(v)
            el[pt, g], eu[pt, g] = cimd(v)
        else
            em[pt, g] = mean(v)
            ci = 1.96 * std(v) / sqrt(n_obs)
            el[pt, g] = em[pt, g] - ci
            eu[pt, g] = em[pt, g] + ci
        end
    end
    return em, el, eu
end
_env_mean_ci(a::AbstractArray{<:Real, 3}, dims::Int64) = _env_ci(a, dims, false)
_env_median_ci(a::AbstractArray{<:Real, 3}, dims::Int64) = _env_ci(a, dims, true)

"""
    env_up(s, x; <keyword arguments>)

Calculate upper cubic-spline envelope from local maxima for a 1-D signal vector.

# Arguments

- `s::AbstractVector`: signal vector
- `x::AbstractVector`: x-axis points (e.g. time points)
- `d::Int64=32`: minimum distance between peaks in samples; smaller values give a tighter fit

# Returns

- `Vector{Float64}`: upper envelope (zeros if fewer than 2 peaks found)
"""
function env_up(s::AbstractVector, x::AbstractVector; d::Int64 = 32)::Vector{Float64}
    # validate
    length(s) == length(x) || throw(
        ArgumentError("Lengths of s ($(length(s))) and x ($(length(x))) must be equal."),
    )

    e = zeros(length(s))

    # locate local maxima at least d samples apart
    p_idx = findpeaks(s; d = d)

    if length(p_idx) < 2
        _info("Envelope cannot be interpolated, less than 2 peaks detected")
    else
        # anchor the spline at the first and last sample so the envelope
        # covers the full signal range without extrapolation artifacts
        p_idx[1] != 1 && pushfirst!(p_idx, 1)
        p_idx[end] != length(s) && push!(p_idx, length(s))

        # fit a cubic spline through the peak positions and evaluate on x
        model = Spline1D(x[p_idx], s[p_idx]; bc = "extrapolate")
        e = model(x)
    end

    return e
end

"""
    env_lo(s, x; <keyword arguments>)

Calculate lower cubic-spline envelope from local minima for a 1-D signal vector.

# Arguments

- `s::AbstractVector`: signal vector
- `x::AbstractVector`: x-axis points (e.g. time points)
- `d::Int64=32`: minimum distance between peaks in samples; smaller values give a tighter fit

# Returns

- `Vector{Float64}`: lower envelope (zeros if fewer than 2 troughs found)
"""
function env_lo(s::AbstractVector, x::AbstractVector; d::Int64 = 32)::Vector{Float64}
    # validate
    length(s) == length(x) || throw(
        ArgumentError("Lengths of s ($(length(s))) and x ($(length(x))) must be equal."),
    )

    e = zeros(length(s))

    # flip the signal so minima become maxima, then find peaks.
    s_tmp = _flipx(s)
    p_idx = findpeaks(s_tmp; d = d)

    if length(p_idx) < 2
        _info("Envelope cannot be interpolated, less than 2 peaks detected")
    else
        # anchor the spline at the first and last sample so the envelope
        # covers the full signal range without extrapolation artifacts
        p_idx[1] != 1 && pushfirst!(p_idx, 1)
        p_idx[end] != length(s) && push!(p_idx, length(s))

        # fit a cubic spline through the peak positions and evaluate on x
        model = Spline1D(x[p_idx], s[p_idx]; bc = "extrapolate")
        e = model(x)
    end

    return e
end

"""
    henv_up(s)

Calculate upper amplitude envelope using the Hilbert transform for a 1-D signal vector.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `Vector{Float64}`: instantaneous amplitude (upper envelope)

# Notes

The Hilbert transform works best for narrowband signals (energy concentrated around a single frequency).
"""
function henv_up(s::AbstractVector)::Vector{Float64}
    return htransform(s).a
end

"""
    henv_lo(s)

Calculate lower amplitude envelope using the Hilbert transform for a 1-D signal vector.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `Vector{Float64}`: negative instantaneous amplitude (lower envelope)

# Notes

The Hilbert transform works best for narrowband signals (energy concentrated around a single frequency).
"""
function henv_lo(s::AbstractVector)::Vector{Float64}
    return htransform(-s).a
end

"""
    tenv(obj; <keyword arguments>)

Calculate temporal envelope for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `d::Int64=32`: minimum distance between peaks in samples; smaller values give a tighter fit

# Returns

Named tuple:

- `e::Array{Float64, 3}`: temporal envelope, shape (channels, samples, epochs)
- `t::Vector{Float64}`: time points
"""
function tenv(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    d::Int64 = 32,
)::@NamedTuple{
    e::Array{Float64, 3},
    t::Vector{Float64},
}
    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)
    # epoch time points
    t = obj.epoch_time

    # pre-allocate output
    e = zeros(ch_n, epoch_len(obj), ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        e[ch_idx, :, ep_idx] = env_up(@view(obj.data[ch[ch_idx], :, ep_idx]), t, d = d)
    end

    return (; e, t)
end

"""
    tenv_mean(obj; <keyword arguments>)

Calculate temporal envelope (mean and 95% CI) for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `dims::Int64`: mean over channels (`dims=1`), epochs (`dims=2`), or both (`dims=3`)
- `d::Int64=32`: minimum distance between peaks in samples; smaller values give a tighter fit

# Returns

Named tuple:

- `em::Matrix{Float64}`: mean temporal envelope, shape (samples, epochs) when `dims=1`, (samples, channels) when `dims=2` or (samples, 1) when `dims=3`
- `el::Matrix{Float64}`: 95% CI lower bound
- `eu::Matrix{Float64}`: 95% CI upper bound
- `t::Vector{Float64}`: time points
"""
function tenv_mean(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    dims::Int64,
    d::Int64 = 32,
)::@NamedTuple{
    em::Matrix{Float64},
    el::Matrix{Float64},
    eu::Matrix{Float64},
    t::Vector{Float64},
}
    # validate
    if dims == 1
        nchannels(obj) >= 2 || throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        nepochs(obj) >= 2 || throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        nchannels(obj) >= 2 || throw(ArgumentError("Number of channels must be ≥ 2."))
        nepochs(obj) >= 2 || throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    tenv_data = tenv(obj; ch = ch, d = d)
    a = tenv_data.e
    t = tenv_data.t

    ch_n = size(a, 1)
    ep_n = size(a, 3)

    # mean and 95% CI over channels (dims = 1), epochs (dims = 2) or both (dims = 3)
    em, el, eu = _env_mean_ci(a, dims)

    return (; em, el, eu, t)
end

"""
    tenv_median(obj; <keyword arguments>)

Calculate temporal envelope (median and 95% CI) for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `dims::Int64`: median over channels (`dims=1`), epochs (`dims=2`), or both (`dims=3`)
- `d::Int64=32`: minimum distance between peaks in samples; smaller values give a tighter fit

# Returns

Named tuple:

- `em::Matrix{Float64}`: median temporal envelope, shape (samples, epochs) when `dims=1`, (samples, channels) when `dims=2` or (samples, 1) when `dims=3`
- `el::Matrix{Float64}`: 95% CI lower bound
- `eu::Matrix{Float64}`: 95% CI upper bound
- `t::Vector{Float64}`: time points
"""
function tenv_median(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    dims::Int64,
    d::Int64 = 32,
)::@NamedTuple{
    em::Matrix{Float64},
    el::Matrix{Float64},
    eu::Matrix{Float64},
    t::Vector{Float64},
}
    # validate
    if dims == 1
        nchannels(obj) >= 2 || throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        nepochs(obj) >= 2 || throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        nchannels(obj) >= 2 || throw(ArgumentError("Number of channels must be ≥ 2."))
        nepochs(obj) >= 2 || throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    tenv_data = tenv(obj; ch = ch, d = d)
    a = tenv_data.e
    t = tenv_data.t

    ch_n = size(a, 1)
    ep_n = size(a, 3)

    # median and 95% CI over channels (dims = 1), epochs (dims = 2) or both (dims = 3)
    em, el, eu = _env_median_ci(a, dims)

    return (; em, el, eu, t)
end

"""
    penv(obj; <keyword arguments>)

Calculate power spectrum envelope (linear power) for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `d::Int64=8`: minimum distance between peaks in samples; smaller values give a tighter fit
- `method::Symbol=:welch`: PSD estimation method:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian filter + Hilbert transform
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `wlen::Int64=sr(obj)`: window length in samples (default = 1 second); used by `:welch` and `:stft`
- `woverlap::Int64=min(min(round(Int64, wlen * 0.90), wlen - 1), wlen - 1)`: window length in samples
- `w::Bool=true`: if `true`, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of Morlet wavelet cycles; for a tuple the number of cycles is logarithmically (geometrically) spaced from `ncyc[1]` to `ncyc[2]` across frequencies
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`)
- `demean::Bool=true`: subtract DC component before estimating PSD

# Returns

Named tuple:

- `e::Array{Float64, 3}`: power spectrum envelope, shape (channels, frequencies, epochs)
- `f::Vector{Float64}`: frequencies
"""
function penv(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    d::Int64 = 8,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = min(round(Int64, wlen * 0.90), wlen - 1),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool = true,
)::@NamedTuple{
    e::Array{Float64, 3},
    f::Vector{Float64},
}
    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)
    # sampling rate
    fs = sr(obj)

    psd_1(x) = _psd(
        x;
        fs = fs,
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
    psd_data = psd_1(@view(obj.data[ch[1], :, 1]))
    f = psd_data.f
    # pre-allocate output
    e = zeros(ch_n, length(f), ep_n)
    e[1, :, 1] = env_up(psd_data.p, f, d = d)

    # calculate over channel and epochs (skipping the one already done)
    idxs = vec(CartesianIndices((ch_n, ep_n)))[2:end]
    @inbounds Threads.@threads :static for idx in idxs
        ch_idx, ep_idx = idx[1], idx[2]
        psd_data = psd_1(@view(obj.data[ch[ch_idx], :, ep_idx]))
        e[ch_idx, :, ep_idx] = env_up(psd_data.p, f, d = d)
    end

    return (e = e, f = f)
end

"""
    penv_mean(obj; <keyword arguments>)

Calculate power spectrum envelope (linear power) (mean and 95% CI) for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `dims::Int64`: mean over channels (`dims=1`), epochs (`dims=2`), or both (`dims=3`)
- `d::Int64=8`: minimum distance between peaks in samples; smaller values give a tighter fit
- `method::Symbol=:welch`: PSD estimation method:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian filter + Hilbert transform
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `wlen::Int64=sr(obj)`: window length in samples (default = 1 second); used by `:welch` and `:stft`
- `woverlap::Int64=min(min(round(Int64, wlen * 0.90), wlen - 1), wlen - 1)`: window length in samples
- `w::Bool=true`: if `true`, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of Morlet wavelet cycles; for a tuple the number of cycles is logarithmically (geometrically) spaced from `ncyc[1]` to `ncyc[2]` across frequencies
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`)
- `demean::Bool=true`: subtract DC component before estimating PSD

# Returns

Named tuple:

- `em::Matrix{Float64}`: mean power envelope, shape (samples, epochs) when `dims=1`, (samples, channels) when `dims=2` or (samples, 1) when `dims=3`
- `el::Matrix{Float64}`: 95% CI lower bound
- `eu::Matrix{Float64}`: 95% CI upper bound
- `f::Vector{Float64}`: frequencies
"""
function penv_mean(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    dims::Int64,
    d::Int64 = 8,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = min(round(Int64, wlen * 0.90), wlen - 1),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool = true,
)::@NamedTuple{
    em::Matrix{Float64},
    el::Matrix{Float64},
    eu::Matrix{Float64},
    f::Vector{Float64},
}
    # validate
    if dims == 1
        nchannels(obj) >= 2 || throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        nepochs(obj) >= 2 || throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        nchannels(obj) >= 2 || throw(ArgumentError("Number of channels must be ≥ 2."))
        nepochs(obj) >= 2 || throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    penv_data = penv(
        obj;
        ch = ch,
        d = d,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        gw = gw,
        demean = demean,
    )
    pw = penv_data.e
    f = penv_data.f

    # number of channels
    ch_n = size(pw, 1)
    # number of epochs
    ep_n = size(pw, 3)

    # mean and 95% CI over channels (dims = 1), epochs (dims = 2) or both (dims = 3)
    em, el, eu = _env_mean_ci(pw, dims)

    return (; em, el, eu, f)
end

"""
    penv_median(obj; <keyword arguments>)

Calculate power spectrum envelope (linear power) (median and 95% CI) for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `dims::Int64`: median over channels (dims = 1) or epochs (dims = 2)
- `d::Int64=8`: minimum distance between peaks in samples; smaller values give a tighter fit
- `method::Symbol=:welch`: PSD estimation method:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian filter + Hilbert transform
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `wlen::Int64=sr(obj)`: window length in samples (default = 1 second); used by `:welch` and `:stft`
- `woverlap::Int64=min(min(round(Int64, wlen * 0.90), wlen - 1), wlen - 1)`: window length in samples
- `w::Bool=true`: if `true`, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of Morlet wavelet cycles; for a tuple the number of cycles is logarithmically (geometrically) spaced from `ncyc[1]` to `ncyc[2]` across frequencies
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`)
- `demean::Bool=true`: subtract DC component before estimating PSD

# Returns

Named tuple:

- `em::Matrix{Float64}`: median power envelope, shape (samples, epochs) when `dims=1`, (samples, channels) when `dims=2` or (samples, 1) when `dims=3`
- `el::Matrix{Float64}`: 95% CI lower bound
- `eu::Matrix{Float64}`: 95% CI upper bound
- `f::Vector{Float64}`: frequencies
"""
function penv_median(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    dims::Int64,
    d::Int64 = 8,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = min(round(Int64, wlen * 0.90), wlen - 1),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool = true,
)::@NamedTuple{
    em::Matrix{Float64},
    el::Matrix{Float64},
    eu::Matrix{Float64},
    f::Vector{Float64},
}
    # validate
    if dims == 1
        nchannels(obj) >= 2 || throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        nepochs(obj) >= 2 || throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        nchannels(obj) >= 2 || throw(ArgumentError("Number of channels must be ≥ 2."))
        nepochs(obj) >= 2 || throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    penv_data = penv(
        obj;
        ch = ch,
        d = d,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        gw = gw,
        demean = demean,
    )
    pw = penv_data.e
    f = penv_data.f

    # number of channels
    ch_n = size(pw, 1)
    # number of epochs
    ep_n = size(pw, 3)

    # median and 95% CI over channels (dims = 1), epochs (dims = 2) or both (dims = 3)
    em, el, eu = _env_median_ci(pw, dims)

    return (; em, el, eu, f)
end

"""
    senv(obj; <keyword arguments>)

Calculate spectral envelope (dominant frequency over time) for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `d::Int64=2`: minimum distance between peaks in samples; smaller values give a tighter fit
- `threshold::Union{Real, Nothing}=nothing`: spectrogram threshold; powers above `threshold` are zeroed before finding the dominant frequency
- `method::Symbol=:stft` spectrogram method:
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `pad::Int64=0`: number of zeros to append
- `db::Bool=true`: normalize powers to dB
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`)
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of Morlet wavelet cycles; for a tuple the number of cycles is logarithmically (geometrically) spaced from `ncyc[1]` to `ncyc[2]` across frequencies
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `wlen::Int64=sr(obj)`: window length in samples (default = 1 second); used by `:welch` and `:stft`
- `woverlap::Int64=min(min(round(Int64, wlen * 0.90), wlen - 1), wlen - 1)`: window length in samples
- `w::Bool=true`: if `true`, apply Hanning window

# Returns

Named tuple:

- `e::Array{Float64, 3}`: spectral envelope, shape (channels, samples, epochs)
- `t::Vector{Float64}`: spectrogram time
"""
function senv(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    d::Int64 = 2,
    threshold::Union{Real, Nothing} = nothing,
    pad::Int64 = 0,
    method::Symbol = :stft,
    db::Bool = true,
    nt::Int64 = 7,
    gw::Real = 5,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    wt::T = wavelet(Morlet(2π); β = 2),
    wlen::Int64 = sr(obj),
    woverlap::Int64 = min(round(Int64, wlen * 0.90), wlen - 1),
    w::Bool = true,
)::@NamedTuple{
    e::Array{Float64, 3},
    t::Vector{Float64},
} where {T <: CWT}
    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)
    # sampling rate
    fs = sr(obj)

    # dry run to determine the spectrogram time vector
    if method === :stft
        spec_data = NeuroAnalyzer._spectrogram(
            @view(obj.data[ch[1], :, 1]);
            fs = fs,
            db = db,
            method = :stft,
            wlen = wlen,
            woverlap = woverlap,
            w = w,
        )
        sp = spec_data.p
    elseif method === :mt
        spec_data = NeuroAnalyzer._spectrogram(
            @view(obj.data[ch[1], :, 1]);
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
            @view(obj.data[ch[1], :, 1]);
            pad = pad,
            fs = fs,
            db = db,
            ncyc = ncyc,
        )
        sp = spec_data.p
    elseif method === :gh
        spec_data = NeuroAnalyzer.ghtspectrogram(
            @view(obj.data[ch[1], :, 1]);
            fs = fs,
            db = db,
            gw = gw,
        )
        sp = spec_data.p
    elseif method === :cwt
        _log_off()
        spec_data = NeuroAnalyzer.cwtspectrogram(
            @view(obj.data[ch[1], :, 1]);
            wt = wt,
            fs = fs,
        )
        _log_on()
        sp = spec_data.m
    end

    # spectrogram time axis (segment centers for :stft / :mt), aligned with the epoch start
    t = spec_data.t .+ obj.epoch_time[1]

    # pre-allocate output
    e = zeros(ch_n, length(t), ep_n)

    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        # compute spectrogram for this channel/epoch
        local sp_loc, sf_loc
        if method === :stft
            spec_data = NeuroAnalyzer._spectrogram(
                @view(obj.data[ch[ch_idx], :, ep_idx]),
                fs = fs,
                db = db,
                method = :stft,
                wlen = wlen,
                woverlap = woverlap,
                w = w,
            )
            sp_loc = spec_data.p
            sf_loc = spec_data.f

        elseif method === :mt
            spec_data = NeuroAnalyzer._spectrogram(
                @view(obj.data[ch[ch_idx], :, ep_idx]),
                fs = fs,
                db = db,
                method = :mt,
                nt = nt,
                wlen = wlen,
                woverlap = woverlap,
                w = w,
            )
            sp_loc = spec_data.p
            sf_loc = spec_data.f

        elseif method === :mw
            spec_data = NeuroAnalyzer.mwspectrogram(
                @view(obj.data[ch[ch_idx], :, ep_idx]),
                pad = pad,
                fs = fs,
                db = db,
                ncyc = ncyc,
            )
            sp_loc = spec_data.p
            sf_loc = spec_data.f

        elseif method === :gh
            spec_data = NeuroAnalyzer.ghtspectrogram(
                @view(obj.data[ch[ch_idx], :, ep_idx]),
                fs = fs,
                db = db,
                gw = gw,
            )
            sp_loc = spec_data.p
            sf_loc = spec_data.f

        elseif method === :cwt
            spec_data = NeuroAnalyzer.cwtspectrogram(
                @view(obj.data[ch[ch_idx], :, ep_idx]),
                wt = wt,
                fs = fs,
            )
            sp_loc = spec_data.m
            sf_loc = spec_data.f
        end

        # optionally zero out powers above the threshold, then reverse so the
        # highest sub-threshold power becomes the "dominant" frequency
        # powers above the threshold are excluded
        !isnothing(threshold) && (sp_loc[sp_loc .> threshold] .= -Inf)
        # for each time bin the frequency with maximum power
        f_idx = [sf_loc[argmax(@view(sp_loc[:, idx2]))] for idx2 in axes(sp_loc, 2)]

        e[ch_idx, :, ep_idx] = env_up(f_idx, t, d = d)
    end

    return (; e, t)
end

"""
    senv_mean(obj; <keyword arguments>)

Calculate spectral envelope (mean and 95% CI) for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `dims::Int64`: mean over channels (`dims=1`), epochs (`dims=2`), or both (`dims=3`)
- `d::Int64=2`: minimum distance between peaks in samples; smaller values give a tighter fit
- `threshold::Union{Real, Nothing}=nothing`: spectrogram threshold; powers above `threshold` are zeroed before finding the dominant frequency
- `method::Symbol=:stft` spectrogram method:
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `pad::Int64=0`: number of zeros to append
- `db::Bool=true`: normalize powers to dB
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`)
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of Morlet wavelet cycles; for a tuple the number of cycles is logarithmically (geometrically) spaced from `ncyc[1]` to `ncyc[2]` across frequencies
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `wlen::Int64=sr(obj)`: window length in samples (default = 1 second); used by `:welch` and `:stft`
- `woverlap::Int64=min(min(round(Int64, wlen * 0.90), wlen - 1), wlen - 1)`: window length in samples
- `w::Bool=true`: if `true`, apply Hanning window

# Returns

Named tuple:

- `em::Matrix{Float64}`: mean spectral envelope, shape (samples, epochs) when `dims=1`, (samples, channels) when `dims=2` or (samples, 1) when `dims=3`
- `eu::Matrix{Float64}`: 95% CI upper bound
- `el::Matrix{Float64}`: 95% CI lower bound
- `t::Vector{Float64}`: time points
"""
function senv_mean(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    dims::Int64,
    d::Int64 = 2,
    threshold::Union{Real, Nothing} = nothing,
    method::Symbol = :stft,
    pad::Int64 = 0,
    db::Bool = true,
    nt::Int64 = 7,
    gw::Real = 5,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    wt::T = wavelet(Morlet(2π); β = 2),
    wlen::Int64 = sr(obj),
    woverlap::Int64 = min(round(Int64, wlen * 0.90), wlen - 1),
    w::Bool = true,
)::@NamedTuple{
    em::Matrix{Float64},
    el::Matrix{Float64},
    eu::Matrix{Float64},
    t::Vector{Float64},
} where {T <: CWT}
    # validate
    if dims == 1
        nchannels(obj) >= 2 || throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        nepochs(obj) >= 2 || throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        nchannels(obj) >= 2 || throw(ArgumentError("Number of channels must be ≥ 2."))
        nepochs(obj) >= 2 || throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    env_data = senv(
        obj;
        ch = ch,
        d = d,
        threshold = threshold,
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
    t = env_data.t

    # number of channels
    ch_n = size(sp, 1)
    # number of epochs
    ep_n = size(sp, 3)

    # mean and 95% CI over channels (dims = 1), epochs (dims = 2) or both (dims = 3)
    em, el, eu = _env_mean_ci(sp, dims)

    return (; em, el, eu, t)
end

"""
    senv_median(obj; <keyword arguments>)

Calculate spectral envelope (median and 95% CI) for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `dims::Int64`: median over channels (`dims=1`), epochs (`dims=2`), or both (`dims=3`)
- `d::Int64=2`: minimum distance between peaks in samples; smaller values give a tighter fit
- `threshold::Union{Real, Nothing}=nothing`: spectrogram threshold; powers above `threshold` are zeroed before finding the dominant frequency
- `method::Symbol=:stft` spectrogram method:
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `pad::Int64=0`: number of zeros to append
- `db::Bool=true`: normalize powers to dB
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`)
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of Morlet wavelet cycles; for a tuple the number of cycles is logarithmically (geometrically) spaced from `ncyc[1]` to `ncyc[2]` across frequencies
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `wlen::Int64=sr(obj)`: window length in samples (default = 1 second); used by `:welch` and `:stft`
- `woverlap::Int64=min(min(round(Int64, wlen * 0.90), wlen - 1), wlen - 1)`: window length in samples
- `w::Bool=true`: if `true`, apply Hanning window

# Returns

Named tuple:

- `em::Matrix{Float64}`: median spectral envelope, shape (samples, epochs) when `dims=1`, (samples, channels) when `dims=2` or (samples, 1) when `dims=3`
- `el::Matrix{Float64}`: 95% CI lower bound
- `eu::Matrix{Float64}`: 95% CI upper bound
- `t::Vector{Float64}`: time points
"""
function senv_median(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    dims::Int64,
    d::Int64 = 2,
    threshold::Union{Real, Nothing} = nothing,
    method::Symbol = :stft,
    pad::Int64 = 0,
    db::Bool = true,
    nt::Int64 = 7,
    gw::Real = 5,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    wt::T = wavelet(Morlet(2π); β = 2),
    wlen::Int64 = sr(obj),
    woverlap::Int64 = min(round(Int64, wlen * 0.90), wlen - 1),
    w::Bool = true,
)::@NamedTuple{
    em::Matrix{Float64},
    el::Matrix{Float64},
    eu::Matrix{Float64},
    t::Vector{Float64},
} where {T <: CWT}
    # validate
    if dims == 1
        nchannels(obj) >= 2 || throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        nepochs(obj) >= 2 || throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        nchannels(obj) >= 2 || throw(ArgumentError("Number of channels must be ≥ 2."))
        nepochs(obj) >= 2 || throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    senv_data = senv(
        obj;
        ch = ch,
        d = d,
        threshold = threshold,
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
    sp = senv_data.e
    t = senv_data.t

    # number of channels
    ch_n = size(sp, 1)
    # number of epochs
    ep_n = size(sp, 3)

    # median and 95% CI over channels (dims = 1), epochs (dims = 2) or both (dims = 3)
    em, el, eu = _env_median_ci(sp, dims)

    return (; em, el, eu, t)
end

"""
    henv(obj; <keyword arguments>)

Calculate Hilbert spectrum amplitude envelope for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `d::Int64=32`: minimum distance between peaks in samples; smaller values give a tighter fit

# Returns

Named tuple:

- `e::Array{Float64, 3}`: Hilbert spectrum amplitude envelope, shape (channels, samples, epochs)
- `t::Vector{Float64}`: time points
"""
function henv(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    d::Int64 = 32,
)::@NamedTuple{
    e::Array{Float64, 3},
    t::Vector{Float64},
}
    _warn(
        "henv() uses Hilbert transform, the signal should be narrowband for best results.",
    )

    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    isempty(ch) && throw(ArgumentError("No channels selected."))

    a = htransform(@view(obj.data[ch, :, :])).a

    # number of channels
    ch_n = size(a, 1)
    # number of epochs
    ep_n = size(a, 3)

    # pre-allocate output
    e = similar(a, Float64)

    # epoch time points
    t = obj.epoch_time

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        e[ch_idx, :, ep_idx] = env_up(@view(a[ch_idx, :, ep_idx]), t, d = d)
    end

    return (; e, t)
end

"""
    henv_mean(obj; <keyword arguments>)

Calculate Hilbert spectrum amplitude envelope (mean and 95% CI) for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `dims::Int64`: mean over channels (`dims=1`), epochs (`dims=2`), or both (`dims=3`)
- `d::Int64=32`: minimum distance between peaks in samples; smaller values give a tighter fit

# Returns

Named tuple:

- `em::Matrix{Float64}`: mean Hilbert spectrum amplitude envelope, shape (samples, epochs) when `dims=1`, (samples, channels) when `dims=2` or (samples, 1) when `dims=3`
- `el::Matrix{Float64}`: 95% CI lower bound
- `eu::Matrix{Float64}`: 95% CI upper bound
- `t::Vector{Float64}`: time points
"""
function henv_mean(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    dims::Int64, d::Int64 = 32,
)::@NamedTuple{
    em::Matrix{Float64},
    el::Matrix{Float64},
    eu::Matrix{Float64},
    t::Vector{Float64},
}
    # validate
    if dims == 1
        nchannels(obj) >= 2 || throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        nepochs(obj) >= 2 || throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        nchannels(obj) >= 2 || throw(ArgumentError("Number of channels must be ≥ 2."))
        nepochs(obj) >= 2 || throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    henv_data = henv(obj; ch = ch, d = d)
    a = henv_data.e
    t = henv_data.t

    # number of channels
    ch_n = size(a, 1)
    # number of epochs
    ep_n = size(a, 3)

    # mean and 95% CI over channels (dims = 1), epochs (dims = 2) or both (dims = 3)
    em, el, eu = _env_mean_ci(a, dims)

    return (; em, el, eu, t)
end

"""
    henv_median(obj; <keyword arguments>)

Calculate Hilbert spectrum amplitude envelope (median and 95% CI) for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `dims::Int64`: median over channels (`dims=1`), epochs (`dims=2`), or both (`dims=3`)
- `d::Int64=32`: minimum distance between peaks in samples; smaller values give a tighter fit

# Returns

Named tuple:

- `em::Matrix{Float64}`: median Hilbert spectrum amplitude envelope, shape (samples, epochs) when `dims=1`, (samples, channels) when `dims=2` or (samples, 1) when `dims=3`
- `el::Matrix{Float64}`: 95% CI lower bound
- `eu::Matrix{Float64}`: 95% CI upper bound
- `t::Vector{Float64}`: time points
"""
function henv_median(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    dims::Int64,
    d::Int64 = 32,
)::@NamedTuple{
    em::Matrix{Float64},
    el::Matrix{Float64},
    eu::Matrix{Float64},
    t::Vector{Float64},
}
    # validate
    if dims == 1
        nchannels(obj) >= 2 || throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        nepochs(obj) >= 2 || throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        nchannels(obj) >= 2 || throw(ArgumentError("Number of channels must be ≥ 2."))
        nepochs(obj) >= 2 || throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    henv_data = henv(obj; ch = ch, d = d)
    a = henv_data.e
    t = henv_data.t

    # number of channels
    ch_n = size(a, 1)
    # number of epochs
    ep_n = size(a, 3)

    # median and 95% CI over channels (dims = 1), epochs (dims = 2) or both (dims = 3)
    em, el, eu = _env_median_ci(a, dims)

    return (; em, el, eu, t)
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
function env_cor(
    env1::Array{Float64, 3},
    env2::Array{Float64, 3},
)::@NamedTuple{ec::Vector{Float64}, p::Vector{Float64}}
    # validate
    size(env1) == size(env2) ||
        throw(ArgumentError("Both envelopes must have the same size."))

    # number of epochs
    ep_n = size(env1, 3)

    # pre-allocate outputs
    ec = zeros(ep_n)
    p = zeros(ep_n)

    # compare envelopes per epochs
    for ep_idx = 1:ep_n
        ctest = CorrelationTest(vec(env1[:, :, ep_idx]), vec(env2[:, :, ep_idx]))
        @inbounds ec[ep_idx] = ctest.r
        @inbounds p[ep_idx] = pvalue(ctest)
    end

    return (; ec, p)
end
