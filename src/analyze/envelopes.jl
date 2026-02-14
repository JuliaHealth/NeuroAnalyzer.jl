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

Calculate upper envelope.

# Arguments

  - `s::AbstractVector`
  - `x::AbstractVector`: x-axis points (e.g. time points)
  - `d::Int64=32`: distance between peeks in points, lower values get better envelope fit

# Returns

  - `e::Vector{Float64}`: envelope
"""
function env_up(s::AbstractVector, x::AbstractVector; d::Int64 = 32)::Vector{Float64}

    @assert length(s) == length(x) "Lengths of s ($(length(s))) and x ($(length(x))) must be equal."

    e = zeros(length(s))

    # find peaks
    p_idx = findpeaks(s; d = d)

    if length(p_idx) < 2
        _info("Envelope cannot be not interpolated, less than 2 peaks detected")
    else
        # add first time-point
        p_idx[1] != 1 && pushfirst!(p_idx, 1)

        # add last time-point
        p_idx[end] != length(s) && push!(p_idx, length(s))

        # interpolate peaks using cubic spline
        model = Spline1D(x[p_idx], s[p_idx]; bc = "extrapolate")
        e = model(x)
    end

    return e

end

"""
    env_lo(s, x; <keyword arguments>)

Calculate lower envelope.

# Arguments

  - `s::AbstractVector`
  - `x::AbstractVector`: x-axis points (e.g. time points)
  - `d::Int64=32`: distance between peeks in points, lower values get better envelope fit

# Returns

  - `e::Vector{Float64}`: envelope
"""
function env_lo(s::AbstractVector, x::AbstractVector; d::Int64 = 32)::Vector{Float64}

    @assert length(s) == length(x) "Lengths of s ($(length(s))) and x ($(length(x))) must be equal."

    e = zeros(length(s))

    # flip the signal along the X axis
    s_tmp = _flipx(s)

    # find peaks
    p_idx = findpeaks(s_tmp; d = d)

    if length(p_idx) < 2
        _info("Envelope cannot be not interpolated, less than 2 peaks detected")
    else
        # add first time-point
        p_idx[1] != 1 && pushfirst!(p_idx, 1)

        # add last time-point
        p_idx[end] != length(s) && push!(p_idx, length(s))

        # interpolate peaks using cubic spline
        model = Spline1D(x[p_idx], s[p_idx]; bc = "extrapolate")
        e = model(x)
    end

    return e

end

"""
    henv_up(s)

Calculate upper envelope using Hilbert transform.

# Arguments

  - `s::AbstractVector`

# Returns

  - `e::Vector{Float64}`: envelope

# Notes

Hilbert transform works best for narrowband signals (i.e., signals with all energy centered about a single frequency).
"""
function henv_up(s::AbstractVector)::Vector{Float64}

    _, e, _, _ = htransform(s)

    return e

end

"""
    henv_lo(s)

Calculate lower envelope using Hilbert transform.

# Arguments

  - `s::AbstractVector`

# Returns

  - `e::Vector{Float64}`: envelope

# Notes

Hilbert transform works best for narrowband signals (i.e., signals with all energy centered about a single frequency).
"""
function henv_lo(s::AbstractVector)::Vector{Float64}

    _, e, _, _ = htransform(-s)

    return -e

end

"""
    tenv(obj; <keyword arguments>)

Calculate temporal envelope.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:

  - `e::Array{Float64, 3}`: temporal envelope
  - `t::Vector{Float64}`: time points
"""
function tenv(
    obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, d::Int64 = 32
)::@NamedTuple{e::Array{Float64, 3}, t::Vector{Float64}}

    ch = exclude_bads ? get_channel(obj; ch = ch, exclude = "bad") : get_channel(obj; ch = ch, exclude = "")
    ch_n = length(ch)
    ep_n = nepochs(obj)
    t = obj.epoch_time

    e = zeros(ch_n, epoch_len(obj), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            e[ch_idx, :, ep_idx] = @views env_up(obj.data[ch[ch_idx], :, ep_idx], t, d = d)
        end
    end

    return (e = e, t = t)

end

"""
    tenv_mean(obj; <keyword arguments>)

Calculate temporal envelope (amplitude): mean and 95% CI.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  - `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:

  - `e_m::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: mean
  - `e_u::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI upper bound
  - `e_l::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI lower bound
  - `t::Vector{Float64}`: time points
"""
function tenv_mean(
    obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, dims::Int64, d::Int64 = 32
)::@NamedTuple{
    e_m::Union{Vector{Float64}, Matrix{Float64}},
    e_u::Union{Vector{Float64}, Matrix{Float64}},
    e_l::Union{Vector{Float64}, Matrix{Float64}},
    t::Vector{Float64},
}

    if dims == 1
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
    elseif dims == 2
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    elseif dims == 3
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    end

    s_a, t = tenv(obj; ch = ch, d = d)

    ch_n = size(s_a, 1)
    ep_n = size(s_a, 3)

    if dims == 1
        # mean over channels

        e_m = zeros(length(t), ep_n)
        e_u = zeros(length(t), ep_n)
        e_l = zeros(length(t), ep_n)

        @inbounds for ep_idx in 1:ep_n
            e_m[:, ep_idx] = @views mean(s_a[:, :, ep_idx], dims = 1)
            s = @views 1.96 * std(e_m[:, ep_idx]) / sqrt(length(e_m[:, ep_idx]))
            e_u[:, ep_idx] = @views e_m[:, ep_idx] .+ s
            e_l[:, ep_idx] = @views e_m[:, ep_idx] .- s
        end
    elseif dims == 2
        # mean over epochs

        e_m = zeros(length(t), ch_n)
        e_u = zeros(length(t), ch_n)
        e_l = zeros(length(t), ch_n)

        @inbounds for ch_idx in 1:ch_n
            e_m[:, ch_idx] = @views mean(s_a[ch_idx, :, :], dims = 2)
            s = @views 1.96 * std(e_m[:, ch_idx]) / sqrt(length(e_m[:, ch_idx]))
            e_u[:, ch_idx] = @views e_m[:, ch_idx] .+ s
            e_l[:, ch_idx] = @views e_m[:, ch_idx] .- s
        end
    else
        # mean over channels and epochs

        e_m, e_u, e_l, _ = tenv_mean(obj; ch = ch, dims = 1, d = d)

        e_m = mean(e_m; dims = 2)
        e_u = mean(e_u; dims = 2)
        e_l = mean(e_l; dims = 2)

        e_m = reshape(e_m, size(e_m, 1))
        e_u = reshape(e_u, size(e_u, 1))
        e_l = reshape(e_l, size(e_l, 1))
    end

    return (e_m = e_m, e_u = e_u, e_l = e_l, t = t)

end

"""
    tenv_median(obj; <keyword arguments>)

Calculate temporal envelope (amplitude): median and 95% CI.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `dims::Int64`: median over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  - `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:

  - `e_m::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: median
  - `e_u::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI upper bound
  - `e_l::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI lower bound
  - `t::Vector{Float64}`: time points
"""
function tenv_median(
    obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, dims::Int64, d::Int64 = 32
)::@NamedTuple{
    e_m::Union{Vector{Float64}, Matrix{Float64}},
    e_u::Union{Vector{Float64}, Matrix{Float64}},
    e_l::Union{Vector{Float64}, Matrix{Float64}},
    t::Vector{Float64},
}

    if dims == 1
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
    elseif dims == 2
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    elseif dims == 3
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    end

    s_a, t = tenv(obj; ch = ch, d = d)

    ch_n = size(s_a, 1)
    ep_n = size(s_a, 3)

    if dims == 1
        # median over channels

        e_m = zeros(length(t), ep_n)
        e_u = zeros(length(t), ep_n)
        e_l = zeros(length(t), ep_n)

        @inbounds for ep_idx in 1:ep_n
            e_m[:, ep_idx] = @views median(s_a[:, :, ep_idx], dims = 1)
            for m_idx in eachindex(t)
                e_u[m_idx, ep_idx], e_l[m_idx, ep_idx] = cimd(s_a[:, m_idx, ep_idx])
            end
        end
    elseif dims == 2
        # median over epochs

        e_m = zeros(length(t), ch_n)
        e_u = zeros(length(t), ch_n)
        e_l = zeros(length(t), ch_n)

        @inbounds for ch_idx in 1:ch_n
            e_m[:, ch_idx] = @views median(s_a[ch_idx, :, :], dims = 2)
            for m_idx in eachindex(t)
                e_u[m_idx, ch_idx], e_l[m_idx, ch_idx] = cimd(s_a[ch_idx, m_idx, :])
            end
        end
    else
        # median over channels and epochs

        e_m, e_u, e_l, _ = tenv_median(obj; ch = ch, dims = 1, d = d)
        e_m = median(e_m; dims = 2)
        e_u = median(e_u; dims = 2)
        e_l = median(e_l; dims = 2)
        e_m = reshape(e_m, size(e_m, 1))
        e_u = reshape(e_u, size(e_u, 1))
        e_l = reshape(e_l, size(e_l, 1))
    end

    return (e_m = e_m, e_u = e_u, e_l = e_l, t = t)

end

"""
    penv(obj; <keyword arguments>)

Calculate power spectrum (in dB) envelope.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
  - `method::Symbol=:welch`: method used to calculate PSD:
      + `:welch`: Welch's periodogram
      + `:fft`: fast Fourier transform
      + `:mt`: multi-tapered periodogram
      + `:stft`: short time Fourier transform
      + `:mw`: Morlet wavelet convolution
  - `nt::Int64=7`: number of Slepian tapers
  - `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
  - `w::Bool=true`: if true, apply Hanning window
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(sr(obj) / 2)`

# Returns

Named tuple containing:

  - `e::Array{Float64, 3}`: power spectrum envelope
  - `f::Vector{Float64}`: frequencies for each envelope
"""
function penv(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    d::Int64 = 8,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = round(Int64, wlen * 0.90),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
)::@NamedTuple{e::Array{Float64, 3}, f::Vector{Float64}}

    ch = exclude_bads ? get_channel(obj; ch = ch, exclude = "bad") : get_channel(obj; ch = ch, exclude = "")
    ch_n = length(ch)
    ep_n = nepochs(obj)
    fs = sr(obj)

    _log_off()
    pw, f = psd(
        obj.data[1, :, 1]; fs = fs, method = method, nt = nt, wlen = wlen, woverlap = woverlap, w = w, ncyc = ncyc
    )
    e = zeros(ch_n, length(pw), ep_n)
    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            pw, _ = psd(
                obj.data[ch[ch_idx], :, ep_idx],
                fs = fs,
                db = true,
                method = method,
                nt = nt,
                wlen = wlen,
                woverlap = woverlap,
                w = w,
                ncyc = ncyc,
            )
            e[ch_idx, :, ep_idx] = env_up(pw, f, d = d)
        end
    end
    _log_on()

    return (e = e, f = f)

end

"""
    penv_mean(obj; <keyword arguments>)

Calculate power spectrum (in dB) envelope: mean and 95% CI.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  - `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
  - `method::Symbol=:welch`: method used to calculate PSD:
      + `:welch`: Welch's periodogram
      + `:fft`: fast Fourier transform
      + `:mt`: multi-tapered periodogram
      + `:stft`: short time Fourier transform
      + `:mw`: Morlet wavelet convolution
  - `nt::Int64=7`: number of Slepian tapers
  - `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
  - `w::Bool=true`: if true, apply Hanning window
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(sr(obj) / 2)`

# Returns

Named tuple containing:

  - `e_m::Array{Float64, 3}`: power spectrum envelope: mean
  - `e_u::Array{Float64, 3}`: power spectrum envelope: 95% CI upper bound
  - `e_l::Array{Float64, 3}`: power spectrum envelope: 95% CI lower bound
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
    woverlap::Int64 = round(Int64, wlen * 0.90),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
)::@NamedTuple{
    e_m::Union{Vector{Float64}, Matrix{Float64}},
    e_u::Union{Vector{Float64}, Matrix{Float64}},
    e_l::Union{Vector{Float64}, Matrix{Float64}},
    f::Vector{Float64},
}

    if dims == 1
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
    elseif dims == 2
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    elseif dims == 3
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    end

    pw, f = penv(obj; ch = ch, d = d, method = method, nt = nt, wlen = wlen, woverlap = woverlap, w = w, ncyc = ncyc)

    ch_n = size(pw, 1)
    ep_n = size(pw, 3)

    if dims == 1
        # mean over channels

        e_m = zeros(length(f), ep_n)
        e_u = zeros(length(f), ep_n)
        e_l = zeros(length(f), ep_n)

        @inbounds for ep_idx in 1:ep_n
            e_m[:, ep_idx] = @views mean(pw[:, :, ep_idx], dims = 1)
            s = @views 1.96 * std(e_m[:, ep_idx]) / sqrt(length(e_m[:, ep_idx]))
            e_u[:, ep_idx] = @views e_m[:, ep_idx] .+ s
            e_l[:, ep_idx] = @views e_m[:, ep_idx] .- s
        end
    elseif dims == 2
        # mean over epochs

        e_m = zeros(length(f), ch_n)
        e_u = zeros(length(f), ch_n)
        e_l = zeros(length(f), ch_n)

        @inbounds for ch_idx in 1:ch_n
            e_m[:, ch_idx] = @views mean(pw[ch_idx, :, :], dims = 2)
            s = @views 1.96 * std(e_m[:, ch_idx]) / sqrt(length(e_m[:, ch_idx]))
            e_u[:, ch_idx] = @views e_m[:, ch_idx] .+ s
            e_l[:, ch_idx] = @views e_m[:, ch_idx] .- s
        end
    else
        # mean over channels and epochs

        e_m, e_u, e_l, _ = penv_mean(obj; ch = ch, dims = 1, d = d)
        e_m = mean(e_m; dims = 2)
        e_u = mean(e_u; dims = 2)
        e_l = mean(e_l; dims = 2)
        e_m = reshape(e_m, size(e_m, 1))
        e_u = reshape(e_u, size(e_u, 1))
        e_l = reshape(e_l, size(e_l, 1))
    end

    return (e_m = e_m, e_u = e_u, e_l = e_l, f = f)

end

"""
    penv_median(obj; <keyword arguments>)

Calculate power spectrum (in dB) envelope: median and 95% CI.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `dims::Int64`: median over channels (dims = 1) or epochs (dims = 2)
  - `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
  - `method::Symbol=:welch`: method used to calculate PSD:
      + `:welch`: Welch's periodogram
      + `:fft`: fast Fourier transform
      + `:mt`: multi-tapered periodogram
      + `:stft`: short time Fourier transform
      + `:mw`: Morlet wavelet convolution
  - `nt::Int64=7`: number of Slepian tapers
  - `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
  - `w::Bool=true`: if true, apply Hanning window
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(sr(obj) / 2)`

# Returns

Named tuple containing:

  - `e_m::Array{Float64, 3}`: power spectrum envelope: median
  - `e_u::Array{Float64, 3}`: power spectrum envelope: 95% CI upper bound
  - `e_l::Array{Float64, 3}`: power spectrum envelope: 95% CI lower bound
  - `f::Vector{Float64}`: power spectrum envelope (useful for plotting over PSD)
"""
function penv_median(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    dims::Int64,
    d::Int64 = 8,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = round(Int64, wlen * 0.90),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
)::@NamedTuple{
    e_m::Union{Vector{Float64}, Matrix{Float64}},
    e_u::Union{Vector{Float64}, Matrix{Float64}},
    e_l::Union{Vector{Float64}, Matrix{Float64}},
    f::Vector{Float64},
}

    if dims == 1
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
    elseif dims == 2
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    elseif dims == 3
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    end

    pw, f = penv(obj; ch = ch, d = d, method = method, nt = nt, wlen = wlen, woverlap = woverlap, w = w, ncyc = ncyc)

    ch_n = size(pw, 1)
    ep_n = size(pw, 3)

    if dims == 1
        # median over channels

        e_m = zeros(length(f), ep_n)
        e_u = zeros(length(f), ep_n)
        e_l = zeros(length(f), ep_n)

        @inbounds for ep_idx in 1:ep_n
            e_m[:, ep_idx] = @views median(pw[:, :, ep_idx], dims = 1)
            for m_idx in eachindex(f)
                e_u[m_idx, ep_idx], e_l[m_idx, ep_idx] = cimd(pw[:, m_idx, ep_idx])
            end
        end
    elseif dims == 2
        # median over epochs

        e_m = zeros(length(f), ch_n)
        e_u = zeros(length(f), ch_n)
        e_l = zeros(length(f), ch_n)

        @inbounds for ch_idx in 1:ch_n
            e_m[:, ch_idx] = @views median(pw[ch_idx, :, :], dims = 2)
            for m_idx in eachindex(f)
                e_u[m_idx, ch_idx], e_l[m_idx, ch_idx] = cimd(pw[ch_idx, :, :])
            end
        end
    else
        # median over channels and epochs

        e_m, e_u, e_l, _ = penv_median(obj; ch = ch, dims = 1, d = d)
        e_m = median(e_m; dims = 2)
        e_u = median(e_u; dims = 2)
        e_l = median(e_l; dims = 2)
        e_m = reshape(e_m, size(e_m, 1))
        e_u = reshape(e_u, size(e_u, 1))
        e_l = reshape(e_l, size(e_l, 1))
    end

    return (e_m = e_m, e_u = e_u, e_l = e_l, f = f)

end

"""
    senv(obj; <keyword arguments>)

Calculate spectral envelope.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
  - `t::Union{Real, Nothing}=nothing`: spectrogram threshold (maximize all powers > t)
  - `method::Symbol=:stft`: method of calculating spectrogram:
      + `:stft`: short-time Fourier transform
      + `:mt`: multi-tapered periodogram
      + `:mw`: Morlet wavelet convolution
      + `:gh`: Gaussian and Hilbert transform
      + `:cwt`: continuous wavelet transformation
  - `pad::Int64=0`: number of zeros to add
  - `db::Bool=true`: normalize powers to dB
  - `nt::Int64=7`: number of Slepian tapers
  - `gw::Real=5`: Gaussian width in Hz
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(sr(obj) / 2)`
  - `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
  - `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
  - `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:

  - `e::Array{Float64, 3}`: spectral envelope
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
    woverlap::Int64 = round(Int64, wlen * 0.90),
    w::Bool = true,
)::@NamedTuple{e::Array{Float64, 3}, t::Vector{Float64}} where {T <: CWT}

    ch = exclude_bads ? get_channel(obj; ch = ch, exclude = "bad") : get_channel(obj; ch = ch, exclude = "")
    ch_n = length(ch)
    ep_n = nepochs(obj)
    fs = sr(obj)

    if method === :stft
        sp, _, _ = @views NeuroAnalyzer.spectrogram(
            obj.data[1, :, 1], fs = fs, db = db, method = :stft, wlen = wlen, woverlap = woverlap, w = w
        )
    elseif method === :mt
        sp, _, _ = @views NeuroAnalyzer.spectrogram(
            obj.data[1, :, 1], fs = fs, db = db, method = :mt, nt = nt, wlen = wlen, woverlap = woverlap, w = w
        )
    elseif method === :mw
        _, sp, _, _ = @views NeuroAnalyzer.mwspectrogram(
            obj.data[1, :, 1], pad = pad, fs = fs, db = db, ncyc = ncyc, w = w
        )
    elseif method === :gh
        sp, _, _ = @views NeuroAnalyzer.ghtspectrogram(obj.data[1, :, 1], fs = fs, db = db, gw = gw, w = w)
    elseif method === :cwt
        _log_off()
        sp, _ = @views NeuroAnalyzer.cwtspectrogram(obj.data[1, :, 1], wt = wt, fs = fs)
        _log_on()
    end
    st = linspace(0, (epoch_len(obj) / fs), size(sp, 2))
    st .+= obj.epoch_time[1]

    e = zeros(ch_n, length(st), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            # prepare spectrogram
            if method === :stft
                sp, sf, _ = @views NeuroAnalyzer.spectrogram(
                    obj.data[ch[ch_idx], :, ep_idx],
                    fs = fs,
                    db = db,
                    method = :stft,
                    wlen = wlen,
                    woverlap = woverlap,
                    w = w,
                )
            elseif method === :mt
                sp, sf, _ = @views NeuroAnalyzer.spectrogram(
                    obj.data[ch[ch_idx], :, ep_idx],
                    fs = fs,
                    db = db,
                    method = :mt,
                    nt = nt,
                    wlen = wlen,
                    woverlap = woverlap,
                    w = w,
                )
            elseif method === :mw
                _, sp, _, sf = @views NeuroAnalyzer.mwspectrogram(
                    obj.data[ch[ch_idx], :, ep_idx], pad = pad, fs = fs, db = db, ncyc = ncyc, w = w
                )
            elseif method === :gh
                sp, _, sf = @views NeuroAnalyzer.ghtspectrogram(
                    obj.data[ch[ch_idx], :, ep_idx], fs = fs, db = db, gw = gw, w = w
                )
            elseif method === :cwt
                _log_off()
                sp, sf = @views NeuroAnalyzer.cwtspectrogram(obj.data[ch[ch_idx], :, ep_idx], wt = wt, fs = fs)
                _log_on()
            end

            # maximize all powers above threshold (t)
            if t !== nothing
                sp[sp .> t] .= 0
                reverse!(sp)
                reverse!(sf)
            end

            f_idx = zeros(length(st))
            m = vec(maximum(sp, dims = 1))
            for idx2 in eachindex(m)
                f_idx[idx2] = sf[vsearch(m[idx2], sp[:, idx2])]
            end
            e[ch_idx, :, ep_idx] = env_up(f_idx, st, d = d)
        end
    end

    return (e = e, t = st)

end

"""
    senv_mean(obj; <keyword arguments>)

Calculate spectral envelope: mean and 95% CI.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  - `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
  - `t::Union{Real, Nothing}=nothing`: spectrogram threshold (maximize all powers > t)
  - `method::Symbol=:stft`: method of calculating spectrogram:
      + `:stft`: short-time Fourier transform
      + `:mt`: multi-tapered periodogram
      + `:mw`: Morlet wavelet convolution
      + `:gh`: Gaussian and Hilbert transform
      + `:cwt`: continuous wavelet transformation
  - `pad::Int64=0`: number of zeros to add
  - `db::Bool=true`: normalize powers to dB
  - `nt::Int64=7`: number of Slepian tapers
  - `gw::Real=5`: Gaussian width in Hz
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(sr(obj) / 2)`
  - `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
  - `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
  - `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:

  - `e_m::Array{Float64, 3}`: spectral envelope: mean
  - `e_u::Array{Float64, 3}`: spectral envelope: 95% CI upper bound
  - `e_l::Array{Float64, 3}`: spectral envelope: 95% CI lower bound
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
    woverlap::Int64 = round(Int64, wlen * 0.90),
    w::Bool = true,
)::@NamedTuple{
    e_m::Union{Vector{Float64}, Matrix{Float64}},
    e_u::Union{Vector{Float64}, Matrix{Float64}},
    e_l::Union{Vector{Float64}, Matrix{Float64}},
    t::Vector{Float64},
} where {T <: CWT}

    if dims == 1
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
    elseif dims == 2
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    elseif dims == 3
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    end

    sp, st = senv(
        obj;
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

    ch_n = size(sp, 1)
    ep_n = size(sp, 3)

    if dims == 1
        # mean over channels

        e_m = zeros(length(st), ep_n)
        e_u = zeros(length(st), ep_n)
        e_l = zeros(length(st), ep_n)

        @inbounds for ep_idx in 1:ep_n
            e_m[:, ep_idx] = @views mean(sp[:, :, ep_idx], dims = 1)
            s = @views 1.96 * std(e_m[:, ep_idx]) / sqrt(length(e_m[:, ep_idx]))
            e_u[:, ep_idx] = e_m[:, ep_idx] .+ s
            e_l[:, ep_idx] = e_m[:, ep_idx] .- s
        end
    elseif dims == 2
        # mean over epochs

        e_m = zeros(length(st), ch_n)
        e_u = zeros(length(st), ch_n)
        e_l = zeros(length(st), ch_n)

        @inbounds for ch_idx in 1:ch_n
            e_m[:, ch_idx] = @views mean(sp[ch_idx, :, :], dims = 2)
            s = @views 1.96 * std(e_m[:, ch_idx]) / sqrt(length(e_m[:, ch_idx]))
            e_u[:, ch_idx] = @views e_m[:, ch_idx] .+ s
            e_l[:, ch_idx] = @views e_m[:, ch_idx] .- s
        end
    else
        # mean over channels and epochs

        e_m, e_u, e_l, _ = senv_mean(
            obj;
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
        e_m = mean(e_m; dims = 2)
        e_u = mean(e_u; dims = 2)
        e_l = mean(e_l; dims = 2)
        e_m = reshape(e_m, size(e_m, 1))
        e_u = reshape(e_u, size(e_u, 1))
        e_l = reshape(e_l, size(e_l, 1))
    end

    return (e_m = e_m, e_u = e_u, e_l = e_l, t = st)

end

"""
    senv_median(obj; <keyword arguments>)

Calculate spectral envelope: median and 95% CI.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `dims::Int64`: median over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  - `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
  - `t::Union{Real, Nothing}=nothing`: spectrogram threshold (maximize all powers > t)
  - `method::Symbol=:stft`: method of calculating spectrogram:
      + `:stft`: short-time Fourier transform
      + `:mt`: multi-tapered periodogram
      + `:mw`: Morlet wavelet convolution
      + `:gh`: Gaussian and Hilbert transform
      + `:cwt`: continuous wavelet transformation
  - `pad::Int64=0`: number of zeros to add
  - `db::Bool=true`: normalize powers to dB
  - `nt::Int64=7`: number of Slepian tapers
  - `gw::Real=5`: Gaussian width in Hz
  - `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(sr(obj) / 2)`
  - `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
  - `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
  - `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
  - `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:

  - `e_m::Array{Float64, 3}`: spectral envelope: median
  - `e_u::Array{Float64, 3}`: spectral envelope: 95% CI upper bound
  - `e_l::Array{Float64, 3}`: spectral envelope: 95% CI lower bound
  - `t::Vector{Float64}`: time points
"""
function senv_median(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    dims::Int64,
    d::Int64 = 2,
    t::Union{Real, Nothing} = nothing,
    flim::Tuple{Real, Real} = (0, sr(obj) / 2),
    nfrq::Int64 = _tlength(flim),
    method::Symbol = :stft,
    pad::Int64 = 0,
    db::Bool = true,
    nt::Int64 = 7,
    frq::Symbol = :log,
    gw::Real = 5,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    wt::T = wavelet(Morlet(2π), β = 2),
    wlen::Int64 = sr(obj),
    woverlap::Int64 = round(Int64, wlen * 0.90),
    w::Bool = true,
)::@NamedTuple{
    e_m::Union{Vector{Float64}, Matrix{Float64}},
    e_u::Union{Vector{Float64}, Matrix{Float64}},
    e_l::Union{Vector{Float64}, Matrix{Float64}},
    t::Vector{Float64},
} where {T <: CWT}

    if dims == 1
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
    elseif dims == 2
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    elseif dims == 3
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    end

    sp, st = senv(
        obj;
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
        w = w,
    )

    ch_n = size(sp, 1)
    ep_n = size(sp, 3)

    if dims == 1
        # median over channels

        e_m = zeros(length(st), ep_n)
        e_u = zeros(length(st), ep_n)
        e_l = zeros(length(st), ep_n)

        @inbounds for ep_idx in 1:ep_n
            e_m[:, ep_idx] = @views median(sp[:, :, ep_idx], dims = 1)
            for m_idx in eachindex(st)
                e_u[m_idx, ep_idx], e_l[m_idx, ep_idx] = cimd(sp[:, m_idx, ep_idx])
            end
        end
    elseif dims == 2
        # median over epochs

        e_m = zeros(length(st), ch_n)
        e_u = zeros(length(st), ch_n)
        e_l = zeros(length(st), ch_n)

        @inbounds for ch_idx in 1:ch_n
            e_m[:, ch_idx] = @views median(sp[ch_idx, :, :], dims = 2)
            for m_idx in eachindex(st)
                e_u[m_idx, ch_idx], e_l[m_idx, ch_idx] = cimd(sp[ch_idx, :, :])
            end
        end
    else
        # median over channels and epochs

        e_m, e_u, e_l, _ = senv_median(
            obj;
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
        e_m = median(e_m; dims = 2)
        e_u = median(e_u; dims = 2)
        e_l = median(e_l; dims = 2)
        e_m = reshape(e_m, size(e_m, 1))
        e_u = reshape(e_u, size(e_u, 1))
        e_l = reshape(e_l, size(e_l, 1))
    end

    return (e_m = e_m, e_u = e_u, e_l = e_l, t = st)

end

"""
    henv(obj; <keyword arguments>)

Calculate Hilbert spectrum amplitude envelope.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:

  - `e::Array{Float64, 3}`: Hilbert spectrum amplitude envelope
  - `t::Vector{Float64}`: time points
"""
function henv(
    obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, d::Int64 = 32
)::@NamedTuple{e::Array{Float64, 3}, t::Vector{Float64}}

    ch = exclude_bads ? get_channel(obj; ch = ch, exclude = "bad") : get_channel(obj; ch = ch, exclude = "")
    _warn("henv() uses Hilbert transform, the signal should be narrowband for best results.")

    _, hamp, _, _ = @views htransform(obj.data[ch, :, :])

    ch_n = size(hamp, 1)
    ep_n = size(hamp, 3)
    e = similar(hamp)

    t = obj.epoch_time

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s = @view hamp[ch_idx, :, ep_idx]
            e[ch_idx, :, ep_idx] = env_up(s, t, d = d)
        end
    end

    return (e = e, t = t)
end

"""
    henv_mean(obj; <keyword arguments>)

Calculate Hilbert spectrum amplitude envelope: mean and 95% CI.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  - `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:

  - `e_m::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: mean
  - `e_u::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI upper bound
  - `e_l::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI lower bound
  - `t::Vector{Float64}`: time points
"""
function henv_mean(
    obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, dims::Int64, d::Int64 = 32
)::@NamedTuple{
    e_m::Union{Vector{Float64}, Matrix{Float64}},
    e_u::Union{Vector{Float64}, Matrix{Float64}},
    e_l::Union{Vector{Float64}, Matrix{Float64}},
    t::Vector{Float64},
}

    if dims == 1
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
    elseif dims == 2
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    elseif dims == 3
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    end

    s_a, t = henv(obj; ch = ch, d = d)

    ch_n = size(s_a, 1)
    ep_n = size(s_a, 3)

    if dims == 1
        # mean over channels

        e_m = zeros(length(t), ep_n)
        e_u = zeros(length(t), ep_n)
        e_l = zeros(length(t), ep_n)

        @inbounds for ep_idx in 1:ep_n
            e_m[:, ep_idx] = @views mean(s_a[:, :, ep_idx], dims = 1)
            s = @views 1.96 * std(e_m[:, ep_idx]) / sqrt(length(e_m[:, ep_idx]))
            e_u[:, ep_idx] = e_m[:, ep_idx] .+ s
            e_l[:, ep_idx] = e_m[:, ep_idx] .- s
        end
    elseif dims == 2
        # mean over epochs

        e_m = zeros(length(t), ch_n)
        e_u = zeros(length(t), ch_n)
        e_l = zeros(length(t), ch_n)

        @inbounds for ch_idx in 1:ch_n
            e_m[:, ch_idx] = @views mean(s_a[ch_idx, :, :], dims = 2)
            s = @views 1.96 * std(e_m[:, ch_idx]) / sqrt(length(e_m[:, ch_idx]))
            e_u[:, ch_idx] = @views e_m[:, ch_idx] .+ s
            e_l[:, ch_idx] = @views e_m[:, ch_idx] .- s
        end
    else
        # mean over channels and epochs

        e_m, e_u, e_l, _ = henv_mean(obj; ch = ch, dims = 1, d = d)
        e_m = mean(e_m; dims = 2)
        e_u = mean(e_u; dims = 2)
        e_l = mean(e_l; dims = 2)
        e_m = reshape(e_m, size(e_m, 1))
        e_u = reshape(e_u, size(e_u, 1))
        e_l = reshape(e_l, size(e_l, 1))
    end

    return (e_m = e_m, e_u = e_u, e_l = e_l, t = t)

end

"""
    henv_median(obj; <keyword arguments>)

Calculate Hilbert spectrum amplitude envelope of `obj`: median and 95% CI.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `dims::Int64`: median over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
  - `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:

  - `e_m::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: median
  - `e_u::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI upper bound
  - `e_l::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI lower bound
  - `t::Vector{Float64}`: time points
"""
function henv_median(
    obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, dims::Int64, d::Int64 = 32
)::@NamedTuple{
    e_m::Union{Vector{Float64}, Matrix{Float64}},
    e_u::Union{Vector{Float64}, Matrix{Float64}},
    e_l::Union{Vector{Float64}, Matrix{Float64}},
    t::Vector{Float64},
}

    if dims == 1
        @assert nchannels(obj) >= 1 "Number of channels must be ≥ 2."
    elseif dims == 2
        @assert nepochs(obj) >= 1 "Number of epochs must be ≥ 2."
    elseif dims == 3
        @assert nchannels(obj) >= 1 "Number of channels must be ≥ 2."
        @assert nepochs(obj) >= 1 "Number of epochs must be ≥ 2."
    end

    s_a, t = henv(obj; ch = ch, d = d)

    ch_n = size(s_a, 1)
    ep_n = size(s_a, 3)

    if dims == 1
        # median over channels

        e_m = zeros(length(t), ep_n)
        e_u = zeros(length(t), ep_n)
        e_l = zeros(length(t), ep_n)

        @inbounds for ep_idx in 1:ep_n
            e_m[:, ep_idx] = @views median(s_a[:, :, ep_idx], dims = 1)
            for m_idx in eachindex(t)
                e_u[m_idx, ep_idx], e_l[m_idx, ep_idx] = cimd(s_a[:, m_idx, ep_idx])
            end
        end
    elseif dims == 2
        # median over epochs

        e_m = zeros(length(t), ch_n)
        e_u = zeros(length(t), ch_n)
        e_l = zeros(length(t), ch_n)

        @inbounds for ch_idx in 1:ch_n
            e_m[:, ch_idx] = median(s_a[ch_idx, :, :], dims = 2)
            for m_idx in eachindex(t)
                e_u[m_idx, ch_idx], e_l[m_idx, ch_idx] = cimd(s_a[ch_idx, m_idx, :])
            end
        end
    else
        # median over channels and epochs

        e_m, e_u, e_l, _ = henv_median(obj; ch = ch, dims = 1, d = d)
        e_m = median(e_m; dims = 2)
        e_u = median(e_u; dims = 2)
        e_l = median(e_l; dims = 2)
        e_m = reshape(e_m, size(e_m, 1))
        e_u = reshape(e_u, size(e_u, 1))
        e_l = reshape(e_l, size(e_l, 1))
    end

    return (e_m = e_m, e_u = e_u, e_l = e_l, t = t)
end

"""
    env_cor(env1, env2)

Calculate envelope correlation.

# Arguments

  - `env1::Array{Float64, 3}`
  - `env2::Array{Float64, 3}`

# Returns

Named tuple containing:

  - `ec::Vector{Float64}`: envelope correlation coefficient
  - `p::Vector{Float64}`: p-value
"""
function env_cor(env1::Array{Float64, 3}, env2::Array{Float64, 3})::@NamedTuple{ec::Vector{Float64}, p::Vector{Float64}}

    @assert size(env1) == size(env2) "Both envelopes must have the same size."

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
