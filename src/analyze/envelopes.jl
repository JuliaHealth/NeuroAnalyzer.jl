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
function env_up(s::AbstractVector, x::AbstractVector; d::Int64=32)::Vector{Float64}

    @assert length(s) == length(x) "Lengths of s ($(length(s))) and x ($(length(x))) must be equal."

    e = zeros(length(s))

    # find peaks
    p_idx = findpeaks(s, d=d)

    if length(p_idx) < 2
        _info("Envelope cannot be not interpolated, less than 2 peaks detected")
    else
        # add first time-point
        p_idx[1] != 1 && pushfirst!(p_idx, 1)

        # add last time-point
        p_idx[end] != length(s) && push!(p_idx, length(s))

        # interpolate peaks using cubic spline
        model = Spline1D(x[p_idx], s[p_idx], bc="extrapolate")
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
function env_lo(s::AbstractVector, x::AbstractVector; d::Int64=32)::Vector{Float64}

    @assert length(s) == length(x) "Lengths of s ($(length(s))) and x ($(length(x))) must be equal."

    e = zeros(length(s))

    # flip the signal along the X axis
    s_tmp = _flipx(s)

    # find peaks
    p_idx = findpeaks(s_tmp, d=d)

    if length(p_idx) < 2
        _info("Envelope cannot be not interpolated, less than 2 peaks detected")
    else
        # add first time-point
        p_idx[1] != 1 && pushfirst!(p_idx, 1)

        # add last time-point
        p_idx[end] != length(s) && push!(p_idx, length(s))

        # interpolate peaks using cubic spline
        model = Spline1D(x[p_idx], s[p_idx], bc="extrapolate")
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

Calculate temporal envelope (amplitude).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `t_env::Array{Float64, 3}`: temporal envelope
- `s_t::Vector{Float64}`: signal time
"""
function tenv(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, d::Int64=32)::@NamedTuple{t_env::Array{Float64, 3}, s_t::Vector{Float64}}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    ch_n = length(ch)
    ep_n = nepochs(obj)
    s_t = obj.epoch_time

    t_env = zeros(ch_n, epoch_len(obj), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            t_env[ch_idx, :, ep_idx] = @views env_up(obj.data[ch[ch_idx], :, ep_idx], s_t, d=d)
        end
    end

    return (t_env=t_env, s_t=s_t)

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
- `t_env_m::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: mean
- `t_env_u::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI upper bound
- `t_env_l::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function tenv_mean(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, dims::Int64, d::Int64=32)::@NamedTuple{t_env_m::Union{Vector{Float64}, Matrix{Float64}}, t_env_u::Union{Vector{Float64}, Matrix{Float64}}, t_env_l::Union{Vector{Float64}, Matrix{Float64}}, s_t::Vector{Float64}}

    if dims == 1
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
    elseif dims == 2
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    elseif dims == 3
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    end

    s_a, s_t = tenv(obj, ch=ch, d=d)

    ch_n = size(s_a, 1)
    ep_n = size(s_a, 3)

    if dims == 1
        # mean over channels

        t_env_m = zeros(length(s_t), ep_n)
        t_env_u = zeros(length(s_t), ep_n)
        t_env_l = zeros(length(s_t), ep_n)

        @inbounds for ep_idx in 1:ep_n
            t_env_m[:, ep_idx] = @views mean(s_a[:, :, ep_idx], dims=1)
            s = @views 1.96 * std(t_env_m[:, ep_idx]) / sqrt(length(t_env_m[:, ep_idx]))
            t_env_u[:, ep_idx] = @views t_env_m[:, ep_idx] .+ s
            t_env_l[:, ep_idx] = @views t_env_m[:, ep_idx] .- s
        end
    elseif dims == 2
        # mean over epochs

        t_env_m = zeros(length(s_t), ch_n)
        t_env_u = zeros(length(s_t), ch_n)
        t_env_l = zeros(length(s_t), ch_n)

        @inbounds for ch_idx in 1:ch_n
            t_env_m[:, ch_idx] = @views mean(s_a[ch_idx, :, :], dims=2)
            s = @views 1.96 * std(t_env_m[:, ch_idx]) / sqrt(length(t_env_m[:, ch_idx]))
            t_env_u[:, ch_idx] = @views t_env_m[:, ch_idx] .+ s
            t_env_l[:, ch_idx] = @views t_env_m[:, ch_idx] .- s
        end
    else
        # mean over channels and epochs

        t_env_m, t_env_u, t_env_l, _ = tenv_mean(obj, ch=ch, dims=1, d=d)

        t_env_m = mean(t_env_m, dims=2)
        t_env_u = mean(t_env_u, dims=2)
        t_env_l = mean(t_env_l, dims=2)

        t_env_m = reshape(t_env_m, size(t_env_m, 1))
        t_env_u = reshape(t_env_u, size(t_env_u, 1))
        t_env_l = reshape(t_env_l, size(t_env_l, 1))
    end

    return (t_env_m=t_env_m, t_env_u=t_env_u, t_env_l=t_env_l, s_t=s_t)

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
- `t_env_m::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: median
- `t_env_u::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI upper bound
- `t_env_l::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function tenv_median(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, dims::Int64, d::Int64=32)::@NamedTuple{t_env_m::Union{Vector{Float64}, Matrix{Float64}}, t_env_u::Union{Vector{Float64}, Matrix{Float64}}, t_env_l::Union{Vector{Float64}, Matrix{Float64}}, s_t::Vector{Float64}}

    if dims == 1
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
    elseif dims == 2
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    elseif dims == 3
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    end

    s_a, s_t = tenv(obj, ch=ch, d=d)

    ch_n = size(s_a, 1)
    ep_n = size(s_a, 3)

    if dims == 1
        # median over channels

        t_env_m = zeros(length(s_t), ep_n)
        t_env_u = zeros(length(s_t), ep_n)
        t_env_l = zeros(length(s_t), ep_n)

        @inbounds for ep_idx in 1:ep_n
            t_env_m[:, ep_idx] = @views median(s_a[:, :, ep_idx], dims=1)
            for m_idx in eachindex(s_t)
                t_env_u[m_idx, ep_idx], t_env_l[m_idx, ep_idx] = cimd(s_a[:, m_idx, ep_idx])
            end
        end
    elseif dims == 2
        # median over epochs

        t_env_m = zeros(length(s_t), ch_n)
        t_env_u = zeros(length(s_t), ch_n)
        t_env_l = zeros(length(s_t), ch_n)

        @inbounds for ch_idx in 1:ch_n
            t_env_m[:, ch_idx] = @views median(s_a[ch_idx, :, :], dims=2)
            for m_idx in eachindex(s_t)
                t_env_u[m_idx, ch_idx], t_env_l[m_idx, ch_idx] = cimd(s_a[ch_idx, m_idx, :])
            end
        end
    else
        # median over channels and epochs

        t_env_m, t_env_u, t_env_l, _ = tenv_median(obj, ch=ch, dims=1, d=d)
        t_env_m = median(t_env_m, dims=2)
        t_env_u = median(t_env_u, dims=2)
        t_env_l = median(t_env_l, dims=2)
        t_env_m = reshape(t_env_m, size(t_env_m, 1))
        t_env_u = reshape(t_env_u, size(t_env_u, 1))
        t_env_l = reshape(t_env_l, size(t_env_l, 1))
    end

    return (t_env_m=t_env_m, t_env_u=t_env_u, t_env_l=t_env_l, s_t=s_t)

end

"""
    penv(obj; <keyword arguments>)

Calculate power spectrum (in dB) envelope.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`

# Returns

Named tuple containing:
- `p_env::Array{Float64, 3}`: power spectrum envelope
- `p_env_frq::Vector{Float64}`: frequencies for each envelope
"""
function penv(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, d::Int64=8, method::Symbol=:welch, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32)::@NamedTuple{p_env::Array{Float64, 3}, p_env_frq::Vector{Float64}}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    ch_n = length(ch)
    ep_n = nepochs(obj)
    fs = sr(obj)

    _log_off()
    pw, pf = psd(obj.data[1, :, 1], fs=fs, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc)
    p_env = zeros(ch_n, length(pw), ep_n)
    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            pw, _ = psd(obj.data[ch[ch_idx], :, ep_idx], fs=fs, db=true, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc)
            p_env[ch_idx, :, ep_idx] = env_up(pw, pf, d=d)
        end
    end
    _log_on()

    return (p_env=p_env, p_env_frq=pf)

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
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`

# Returns

Named tuple containing:
- `p_env_m::Array{Float64, 3}`: power spectrum envelope: mean
- `p_env_u::Array{Float64, 3}`: power spectrum envelope: 95% CI upper bound
- `p_env_l::Array{Float64, 3}`: power spectrum envelope: 95% CI lower bound
- `p_env_frq::Vector{Float64}`: power spectrum envelope (useful for plotting over PSD)
"""
function penv_mean(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, dims::Int64, d::Int64=8, method::Symbol=:welch, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32)::@NamedTuple{p_env_m::Union{Vector{Float64}, Matrix{Float64}}, p_env_u::Union{Vector{Float64}, Matrix{Float64}}, p_env_l::Union{Vector{Float64}, Matrix{Float64}}, p_env_frq::Vector{Float64}}

    if dims == 1
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
    elseif dims == 2
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    elseif dims == 3
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    end

    pw, pf = penv(obj, ch=ch, d=d, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc)

    ch_n = size(pw, 1)
    ep_n = size(pw, 3)

    if dims == 1
        # mean over channels

        p_env_m = zeros(length(pf), ep_n)
        p_env_u = zeros(length(pf), ep_n)
        p_env_l = zeros(length(pf), ep_n)

        @inbounds for ep_idx in 1:ep_n
            p_env_m[:, ep_idx] = @views mean(pw[:, :, ep_idx], dims=1)
            s = @views 1.96 * std(p_env_m[:, ep_idx]) / sqrt(length(p_env_m[:, ep_idx]))
            p_env_u[:, ep_idx] = @views p_env_m[:, ep_idx] .+ s
            p_env_l[:, ep_idx] = @views p_env_m[:, ep_idx] .- s
        end
    elseif dims == 2
        # mean over epochs

        p_env_m = zeros(length(pf), ch_n)
        p_env_u = zeros(length(pf), ch_n)
        p_env_l = zeros(length(pf), ch_n)

        @inbounds for ch_idx in 1:ch_n
            p_env_m[:, ch_idx] = @views mean(pw[ch_idx, :, :], dims=2)
            s = @views 1.96 * std(p_env_m[:, ch_idx]) / sqrt(length(p_env_m[:, ch_idx]))
            p_env_u[:, ch_idx] = @views p_env_m[:, ch_idx] .+ s
            p_env_l[:, ch_idx] = @views p_env_m[:, ch_idx] .- s
        end
    else
        # mean over channels and epochs

        p_env_m, p_env_u, p_env_l, _ = penv_mean(obj, ch=ch, dims=1, d=d)
        p_env_m = mean(p_env_m, dims=2)
        p_env_u = mean(p_env_u, dims=2)
        p_env_l = mean(p_env_l, dims=2)
        p_env_m = reshape(p_env_m, size(p_env_m, 1))
        p_env_u = reshape(p_env_u, size(p_env_u, 1))
        p_env_l = reshape(p_env_l, size(p_env_l, 1))
    end

    return (p_env_m=p_env_m, p_env_u=p_env_u, p_env_l=p_env_l, p_env_frq=pf)

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
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`

# Returns

Named tuple containing:
- `p_env_m::Array{Float64, 3}`: power spectrum envelope: median
- `p_env_u::Array{Float64, 3}`: power spectrum envelope: 95% CI upper bound
- `p_env_l::Array{Float64, 3}`: power spectrum envelope: 95% CI lower bound
- `p_env_frq::Vector{Float64}`: power spectrum envelope (useful for plotting over PSD)
"""
function penv_median(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, dims::Int64, d::Int64=8, method::Symbol=:welch, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32)::@NamedTuple{p_env_m::Union{Vector{Float64}, Matrix{Float64}}, p_env_u::Union{Vector{Float64}, Matrix{Float64}}, p_env_l::Union{Vector{Float64}, Matrix{Float64}}, p_env_frq::Vector{Float64}}

    if dims == 1
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
    elseif dims == 2
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    elseif dims == 3
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    end

    pw, pf = penv(obj, ch=ch, d=d, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc)

    ch_n = size(pw, 1)
    ep_n = size(pw, 3)

    if dims == 1
        # median over channels

        p_env_m = zeros(length(pf), ep_n)
        p_env_u = zeros(length(pf), ep_n)
        p_env_l = zeros(length(pf), ep_n)

        @inbounds for ep_idx in 1:ep_n
            p_env_m[:, ep_idx] = @views median(pw[:, :, ep_idx], dims=1)
            for m_idx in eachindex(pf)
                p_env_u[m_idx, ep_idx], p_env_l[m_idx, ep_idx] = cimd(pw[:, m_idx, ep_idx])
            end
        end
    elseif dims == 2
        # median over epochs

        p_env_m = zeros(length(pf), ch_n)
        p_env_u = zeros(length(pf), ch_n)
        p_env_l = zeros(length(pf), ch_n)

        @inbounds for ch_idx in 1:ch_n
            p_env_m[:, ch_idx] = @views median(pw[ch_idx, :, :], dims=2)
            for m_idx in eachindex(pf)
                p_env_u[m_idx, ch_idx], p_env_l[m_idx, ch_idx] = cimd(pw[ch_idx, :, :])
            end
        end
    else
        # median over channels and epochs

        p_env_m, p_env_u, p_env_l, _ = penv_median(obj, ch=ch, dims=1, d=d)
        p_env_m = median(p_env_m, dims=2)
        p_env_u = median(p_env_u, dims=2)
        p_env_l = median(p_env_l, dims=2)
        p_env_m = reshape(p_env_m, size(p_env_m, 1))
        p_env_u = reshape(p_env_u, size(p_env_u, 1))
        p_env_l = reshape(p_env_l, size(p_env_l, 1))
    end

    return (p_env_m=p_env_m, p_env_u=p_env_u, p_env_l=p_env_l, p_env_frq=pf)

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
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `pad::Int64=0`: number of zeros to add
- `db::Bool=true`: normalize powers to dB
- `nt::Int64=7`: number of Slepian tapers
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:
- `s_env::Array{Float64, 3}`: spectral envelope
- `s_env_t::Vector{Float64}`: spectrogram time
"""
function senv(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, d::Int64=2, t::Union{Real, Nothing}=nothing, pad::Int64=0, method::Symbol=:stft, db::Bool=true, nt::Int64=7, gw::Real=5, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, wt::T=wavelet(Morlet(2π), β=2), wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)::@NamedTuple{s_env::Array{Float64, 3}, s_env_t::Vector{Float64}} where {T <: CWT}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    ch_n = length(ch)
    ep_n = nepochs(obj)
    fs = sr(obj)

    if method === :stft
        sp, _, _ = @views NeuroAnalyzer.spectrogram(obj.data[1, :, 1], fs=fs, db=db, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
    elseif method === :mt
        sp, _, _ = @views NeuroAnalyzer.spectrogram(obj.data[1, :, 1], fs=fs, db=db, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
    elseif method === :mw
    _, sp, _, _ = @views NeuroAnalyzer.mwspectrogram(obj.data[1, :, 1], pad=pad, fs=fs, db=db, ncyc=ncyc, w=w)
    elseif method === :gh
        sp, _, _ = @views NeuroAnalyzer.ghtspectrogram(obj.data[1, :, 1], fs=fs, db=db, gw=gw, w=w)
    elseif method === :cwt
        _log_off()
        sp, _ = @views NeuroAnalyzer.cwtspectrogram(obj.data[1, :, 1], wt=wt, fs=fs)
        _log_on()
    end
    st = linspace(0, (epoch_len(obj) / fs), size(sp, 2))
    st .+= obj.epoch_time[1]

    s_env = zeros(ch_n, length(st), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            # prepare spectrogram
            if method === :stft
                sp, sf, _ = @views NeuroAnalyzer.spectrogram(obj.data[ch[ch_idx], :, ep_idx], fs=fs, db=db, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
            elseif method === :mt
                sp, sf, _ = @views NeuroAnalyzer.spectrogram(obj.data[ch[ch_idx], :, ep_idx], fs=fs, db=db, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
            elseif method === :mw
                _, sp, _, sf = @views NeuroAnalyzer.mwspectrogram(obj.data[ch[ch_idx], :, ep_idx], pad=pad, fs=fs, db=db, ncyc=ncyc, w=w)
            elseif method === :gh
                sp, _, sf = @views NeuroAnalyzer.ghtspectrogram(obj.data[ch[ch_idx], :, ep_idx], fs=fs, db=db, gw=gw, w=w)
            elseif method === :cwt
                _log_off()
                sp, sf = @views NeuroAnalyzer.cwtspectrogram(obj.data[ch[ch_idx], :, ep_idx], wt=wt, fs=fs)
                _log_on()
            end

            # maximize all powers above threshold (t)
            if t !== nothing
                sp[sp .> t] .= 0
                reverse!(sp)
                reverse!(sf)
            end

            f_idx = zeros(length(st))
            m = vec(maximum(sp, dims=1))
            for idx2 in eachindex(m)
                f_idx[idx2] = sf[vsearch(m[idx2], sp[:, idx2])]
            end
            s_env[ch_idx, :, ep_idx] = env_up(f_idx, st, d=d)
        end
    end

    return (s_env=s_env, s_env_t=st)

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
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `pad::Int64=0`: number of zeros to add
- `db::Bool=true`: normalize powers to dB
- `nt::Int64=7`: number of Slepian tapers
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:
- `s_env_m::Array{Float64, 3}`: spectral envelope: mean
- `s_env_u::Array{Float64, 3}`: spectral envelope: 95% CI upper bound
- `s_env_l::Array{Float64, 3}`: spectral envelope: 95% CI lower bound
- `s_env_t::Vector{Float64}`: spectral envelope (useful for plotting over spectrogram)
"""
function senv_mean(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, dims::Int64, d::Int64=2, t::Union{Real, Nothing}=nothing, method::Symbol=:stft, pad::Int64=0, db::Bool=true, nt::Int64=7, gw::Real=5, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, wt::T=wavelet(Morlet(2π), β=2), wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)::@NamedTuple{s_env_m::Union{Vector{Float64}, Matrix{Float64}}, s_env_u::Union{Vector{Float64}, Matrix{Float64}}, s_env_l::Union{Vector{Float64}, Matrix{Float64}}, s_env_t::Vector{Float64}} where {T <: CWT}

    if dims == 1
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
    elseif dims == 2
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    elseif dims == 3
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    end

    sp, st = senv(obj, ch=ch, d=d, t=t, pad=pad, method=method, db=db, nt=nt, gw=gw, ncyc=ncyc, wt=wt, wlen=wlen, woverlap=woverlap, w=w)

    ch_n = size(sp, 1)
    ep_n = size(sp, 3)

    if dims == 1
        # mean over channels

        s_env_m = zeros(length(st), ep_n)
        s_env_u = zeros(length(st), ep_n)
        s_env_l = zeros(length(st), ep_n)

        @inbounds for ep_idx in 1:ep_n
            s_env_m[:, ep_idx] = @views mean(sp[:, :, ep_idx], dims=1)
            s = @views 1.96 * std(s_env_m[:, ep_idx]) / sqrt(length(s_env_m[:, ep_idx]))
            s_env_u[:, ep_idx] = s_env_m[:, ep_idx] .+ s
            s_env_l[:, ep_idx] = s_env_m[:, ep_idx] .- s
        end
    elseif dims == 2
        # mean over epochs

        s_env_m = zeros(length(st), ch_n)
        s_env_u = zeros(length(st), ch_n)
        s_env_l = zeros(length(st), ch_n)

        @inbounds for ch_idx in 1:ch_n
            s_env_m[:, ch_idx] = @views mean(sp[ch_idx, :, :], dims=2)
            s = @views 1.96 * std(s_env_m[:, ch_idx]) / sqrt(length(s_env_m[:, ch_idx]))
            s_env_u[:, ch_idx] = @views s_env_m[:, ch_idx] .+ s
            s_env_l[:, ch_idx] = @views s_env_m[:, ch_idx] .- s
        end
    else
        # mean over channels and epochs

        s_env_m, s_env_u, s_env_l, _ = senv_mean(obj, ch=ch, dims=1, d=d, pad=pad, method=method, db=db, nt=nt, gw=gw, ncyc=ncyc, wt=wt, wlen=wlen, woverlap=woverlap, w=w)
        s_env_m = mean(s_env_m, dims=2)
        s_env_u = mean(s_env_u, dims=2)
        s_env_l = mean(s_env_l, dims=2)
        s_env_m = reshape(s_env_m, size(s_env_m, 1))
        s_env_u = reshape(s_env_u, size(s_env_u, 1))
        s_env_l = reshape(s_env_l, size(s_env_l, 1))
    end

    return (s_env_m=s_env_m, s_env_u=s_env_u, s_env_l=s_env_l, s_env_t=st)

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
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `pad::Int64=0`: number of zeros to add
- `db::Bool=true`: normalize powers to dB
- `nt::Int64=7`: number of Slepian tapers
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:
- `s_env_m::Array{Float64, 3}`: spectral envelope: median
- `s_env_u::Array{Float64, 3}`: spectral envelope: 95% CI upper bound
- `s_env_l::Array{Float64, 3}`: spectral envelope: 95% CI lower bound
- `s_env_t::Vector{Float64}`: spectral envelope (useful for plotting over spectrogram)
"""
function senv_median(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, dims::Int64, d::Int64=2, t::Union{Real, Nothing}=nothing, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), frq_n::Int64=_tlength(frq_lim), method::Symbol=:stft, pad::Int64=0, db::Bool=true, nt::Int64=7, frq::Symbol=:log, gw::Real=5, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, wt::T=wavelet(Morlet(2π), β=2), wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)::@NamedTuple{s_env_m::Union{Vector{Float64}, Matrix{Float64}}, s_env_u::Union{Vector{Float64}, Matrix{Float64}}, s_env_l::Union{Vector{Float64}, Matrix{Float64}}, s_env_t::Vector{Float64}} where {T <: CWT}

    if dims == 1
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
    elseif dims == 2
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    elseif dims == 3
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    end

    sp, st = senv(obj, ch=ch, d=d, t=t, pad=pad, method=method, nt=nt, db=db, gw=gw, ncyc=ncyc, wt=wt, wlen=wlen, woverlap=woverlap, w=w)

    ch_n = size(sp, 1)
    ep_n = size(sp, 3)

    if dims == 1
        # median over channels

        s_env_m = zeros(length(st), ep_n)
        s_env_u = zeros(length(st), ep_n)
        s_env_l = zeros(length(st), ep_n)

        @inbounds for ep_idx in 1:ep_n
            s_env_m[:, ep_idx] = @views median(sp[:, :, ep_idx], dims=1)
            for m_idx in eachindex(st)
                s_env_u[m_idx, ep_idx], s_env_l[m_idx, ep_idx] = cimd(sp[:, m_idx, ep_idx])
            end
        end
    elseif dims == 2
        # median over epochs

        s_env_m = zeros(length(st), ch_n)
        s_env_u = zeros(length(st), ch_n)
        s_env_l = zeros(length(st), ch_n)

        @inbounds for ch_idx in 1:ch_n
            s_env_m[:, ch_idx] = @views median(sp[ch_idx, :, :], dims=2)
            for m_idx in eachindex(st)
                s_env_u[m_idx, ch_idx], s_env_l[m_idx, ch_idx] = cimd(sp[ch_idx, :, :])
            end
        end
    else
        # median over channels and epochs

        s_env_m, s_env_u, s_env_l, _ = senv_median(obj, ch=ch, dims=1, d=d, pad=pad, method=method, db=db, nt=nt, gw=gw, ncyc=ncyc, wt=wt, wlen=wlen, woverlap=woverlap, w=w)
        s_env_m = median(s_env_m, dims=2)
        s_env_u = median(s_env_u, dims=2)
        s_env_l = median(s_env_l, dims=2)
        s_env_m = reshape(s_env_m, size(s_env_m, 1))
        s_env_u = reshape(s_env_u, size(s_env_u, 1))
        s_env_l = reshape(s_env_l, size(s_env_l, 1))
    end

    return (s_env_m=s_env_m, s_env_u=s_env_u, s_env_l=s_env_l, s_env_t=st)

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
- `h_env::Array{Float64, 3}`: Hilbert spectrum amplitude envelope
- `s_t::Vector{Float64}`: signal time
"""
function henv(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, d::Int64=32)::@NamedTuple{h_env::Array{Float64, 3}, s_t::Vector{Float64}}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    _warn("henv() uses Hilbert transform, the signal should be narrowband for best results.")

    _, hamp, _, _ = @views htransform(obj.data[ch, :, :])

    ch_n = size(hamp, 1)
    ep_n = size(hamp, 3)
    h_env = similar(hamp)

    s_t = obj.epoch_time

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s = @view hamp[ch_idx, :, ep_idx]
            h_env[ch_idx, :, ep_idx] = env_up(s, s_t, d=d)
        end
    end

    return (h_env=h_env, s_t=s_t)
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
- `h_env_m::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: mean
- `h_env_u::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI upper bound
- `h_env_l::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function henv_mean(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, dims::Int64, d::Int64=32)::@NamedTuple{h_env_m::Union{Vector{Float64}, Matrix{Float64}}, h_env_u::Union{Vector{Float64}, Matrix{Float64}}, h_env_l::Union{Vector{Float64}, Matrix{Float64}}, s_t::Vector{Float64}}

    if dims == 1
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
    elseif dims == 2
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    elseif dims == 3
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    end

    s_a, s_t = henv(obj, ch=ch, d=d)

    ch_n = size(s_a, 1)
    ep_n = size(s_a, 3)

    if dims == 1
        # mean over channels

        h_env_m = zeros(length(s_t), ep_n)
        h_env_u = zeros(length(s_t), ep_n)
        h_env_l = zeros(length(s_t), ep_n)

        @inbounds for ep_idx in 1:ep_n
            h_env_m[:, ep_idx] = @views mean(s_a[:, :, ep_idx], dims=1)
            s = @views 1.96 * std(h_env_m[:, ep_idx]) / sqrt(length(h_env_m[:, ep_idx]))
            h_env_u[:, ep_idx] = h_env_m[:, ep_idx] .+ s
            h_env_l[:, ep_idx] = h_env_m[:, ep_idx] .- s
        end
    elseif dims == 2
        # mean over epochs

        h_env_m = zeros(length(s_t), ch_n)
        h_env_u = zeros(length(s_t), ch_n)
        h_env_l = zeros(length(s_t), ch_n)

        @inbounds for ch_idx in 1:ch_n
            h_env_m[:, ch_idx] = @views mean(s_a[ch_idx, :, :], dims=2)
            s = @views 1.96 * std(h_env_m[:, ch_idx]) / sqrt(length(h_env_m[:, ch_idx]))
            h_env_u[:, ch_idx] = @views h_env_m[:, ch_idx] .+ s
            h_env_l[:, ch_idx] = @views h_env_m[:, ch_idx] .- s
        end
    else
        # mean over channels and epochs

        h_env_m, h_env_u, h_env_l, _ = henv_mean(obj, ch=ch, dims=1, d=d)
        h_env_m = mean(h_env_m, dims=2)
        h_env_u = mean(h_env_u, dims=2)
        h_env_l = mean(h_env_l, dims=2)
        h_env_m = reshape(h_env_m, size(h_env_m, 1))
        h_env_u = reshape(h_env_u, size(h_env_u, 1))
        h_env_l = reshape(h_env_l, size(h_env_l, 1))
    end

    return (h_env_m=h_env_m, h_env_u=h_env_u, h_env_l=h_env_l, s_t=s_t)

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
- `h_env_m::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: median
- `h_env_u::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI upper bound
- `h_env_l::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function henv_median(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, dims::Int64, d::Int64=32)::@NamedTuple{h_env_m::Union{Vector{Float64}, Matrix{Float64}}, h_env_u::Union{Vector{Float64}, Matrix{Float64}}, h_env_l::Union{Vector{Float64}, Matrix{Float64}}, s_t::Vector{Float64}}

    if dims == 1
        @assert nchannels(obj) >= 1 "Number of channels must be ≥ 2."
    elseif dims == 2
        @assert nepochs(obj) >= 1 "Number of epochs must be ≥ 2."
    elseif dims == 3
        @assert nchannels(obj) >= 1 "Number of channels must be ≥ 2."
        @assert nepochs(obj) >= 1 "Number of epochs must be ≥ 2."
    end

    s_a, s_t = henv(obj, ch=ch, d=d)

    ch_n = size(s_a, 1)
    ep_n = size(s_a, 3)

    if dims == 1
        # median over channels

        h_env_m = zeros(length(s_t), ep_n)
        h_env_u = zeros(length(s_t), ep_n)
        h_env_l = zeros(length(s_t), ep_n)

        @inbounds for ep_idx in 1:ep_n
            h_env_m[:, ep_idx] = @views median(s_a[:, :, ep_idx], dims=1)
            for m_idx in eachindex(s_t)
                h_env_u[m_idx, ep_idx], h_env_l[m_idx, ep_idx] = cimd(s_a[:, m_idx, ep_idx])
            end
        end
    elseif dims == 2
        # median over epochs

        h_env_m = zeros(length(s_t), ch_n)
        h_env_u = zeros(length(s_t), ch_n)
        h_env_l = zeros(length(s_t), ch_n)

        @inbounds for ch_idx in 1:ch_n
            h_env_m[:, ch_idx] = median(s_a[ch_idx, :, :], dims=2)
            for m_idx in eachindex(s_t)
                h_env_u[m_idx, ch_idx], h_env_l[m_idx, ch_idx] = cimd(s_a[ch_idx, m_idx, :])
            end
        end
    else
        # median over channels and epochs

        h_env_m, h_env_u, h_env_l, _ = henv_median(obj, ch=ch, dims=1, d=d)
        h_env_m = median(h_env_m, dims=2)
        h_env_u = median(h_env_u, dims=2)
        h_env_l = median(h_env_l, dims=2)
        h_env_m = reshape(h_env_m, size(h_env_m, 1))
        h_env_u = reshape(h_env_u, size(h_env_u, 1))
        h_env_l = reshape(h_env_l, size(h_env_l, 1))
    end

    return (h_env_m=h_env_m, h_env_u=h_env_u, h_env_l=h_env_l, s_t=s_t)
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

    return (ec=ec, p=p)

end
