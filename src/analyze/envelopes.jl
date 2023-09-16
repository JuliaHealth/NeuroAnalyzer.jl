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
    tenv(obj; ch, d)

Calculate temporal envelope (amplitude).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `t_env::Array{Float64, 3}`: temporal envelope
- `s_t::Vector{Float64}`: signal time
"""
function tenv(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), d::Int64=32)
    
    _check_channels(obj, ch)

    ch_n = length(ch)
    ep_n = nepochs(obj)
    s_t = obj.epoch_time

    t_env = zeros(ch_n, epoch_len(obj), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s = @view obj.data[ch[ch_idx], :, ep_idx]
            # find peaks
            p_idx = findpeaks(s, d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(s))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(s_t[p_idx], s[p_idx])
                try
                    t_env[ch_idx, :, ep_idx] = model(s_t)
                catch
                    @warn "CubicSpline error, using Loess."
                    model = Loess.loess(s_t[p_idx], s[p_idx], span=0.5)
                    t_env[ch_idx, :, ep_idx] = Loess.predict(model, s_t)
                end
            else
                _info("Less than 5 peaks detected, using Loess.")
                model = Loess.loess(s_t[p_idx], s[p_idx], span=0.5)
                t_env[ch_idx, :, ep_idx] = Loess.predict(model, s_t)
            end
            t_env[ch_idx, 1, ep_idx] = t_env[ch_idx, 2, ep_idx]
        end
    end
    
    return (t_env=t_env, s_t=s_t)

end

"""
    tenv_mean(obj; ch, dims, d)

Calculate temporal envelope (amplitude): mean and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `t_env_m::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: mean
- `t_env_u::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI upper bound
- `t_env_l::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function tenv_mean(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=32)
    
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

        @inbounds @simd for ep_idx in 1:ep_n
            t_env_m[:, ep_idx] = mean(s_a[:, :, ep_idx], dims=1)
            s = std(t_env_m[:, ep_idx]) / sqrt(length(t_env_m[:, ep_idx]))
            t_env_u[:, ep_idx] = @. t_env_m[:, ep_idx] + 1.96 * s
            t_env_l[:, ep_idx] = @. t_env_m[:, ep_idx] - 1.96 * s
        end
    elseif dims == 2
        # mean over epochs

        t_env_m = zeros(length(s_t), ch_n)
        t_env_u = zeros(length(s_t), ch_n)
        t_env_l = zeros(length(s_t), ch_n)

        @inbounds @simd for ch_idx in 1:ch_n
            t_env_m[:, ch_idx] = mean(s_a[ch_idx, :, :], dims=2)
            s = std(t_env_m[:, ch_idx]) / sqrt(length(t_env_m[:, ch_idx]))
            t_env_u[:, ch_idx] = @views @. t_env_m[:, ch_idx] + 1.96 * s
            t_env_l[:, ch_idx] = @views @. t_env_m[:, ch_idx] - 1.96 * s
        end
    else
        # mean over channels and epochs

        t_env_m, t_env_u, t_env_l, _ = tenv_mean(obj, dims=1, d=d)

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
    tenv_median(obj; ch, dims, d)

Calculate temporal envelope (amplitude): median and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `dims::Int64`: median over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `t_env_m::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: median
- `t_env_u::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI upper bound
- `t_env_l::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function tenv_median(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=32)
    
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

        @inbounds @simd for ep_idx in 1:ep_n
            t_env_m[:, ep_idx] = median(s_a[:, :, ep_idx], dims=1)
            t_idx = findpeaks(t_env_m[:, ep_idx], d=d)
            pushfirst!(t_idx, 1)
            push!(t_idx, length(t_env_m[:, ep_idx]))
            if length(t_idx) > 4
                model = CubicSpline(s_t[t_idx], t_env_m[t_idx])
                try
                    t_env_m[:, ep_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(t_env_m[:, ep_idx]) / sqrt(length(t_env_m[:, ep_idx]))
            t_env_u[:, ep_idx] = @. t_env_m[:, ep_idx] + 1.96 * s
            t_env_l[:, ep_idx] = @. t_env_m[:, ep_idx] - 1.96 * s
        end
    elseif dims == 2
        # median over epochs

        t_env_m = zeros(length(s_t), ch_n)
        t_env_u = zeros(length(s_t), ch_n)
        t_env_l = zeros(length(s_t), ch_n)

        @inbounds @simd for ch_idx in 1:ch_n
            t_env_m[:, idx] = median(s_a[ch_idx, :, :], dims=2)
            t_idx = findpeaks(t_env_m[:, ch_idx], d=d)
            pushfirst!(t_idx, 1)
            push!(t_idx, length(t_env_m[:, ch_idx]))
            if length(t_idx) > 4
                model = CubicSpline(s_t[t_idx], t_env_m[t_idx])
                try
                    t_env_m[:, ch_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(t_env_m[:, ch_idx]) / sqrt(length(t_env_m[:, ch_idx]))
            t_env_u[:, ch_idx] = @views @. t_env_m[:, ch_idx] + 1.96 * s
            t_env_l[:, ch_idx] = @views @. t_env_m[:, ch_idx] - 1.96 * s
        end
    else
        # median over channels and epochs

        t_env_m, t_env_u, t_env_l, _ = tenv_median(obj, dims=1, d=d)

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
    penv(obj; ch, d, method, nt, wlen, woverlap, w)

Calculate power (in dB) envelope.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=8`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds
- `frq_n::Int64=length(frq_lim[1]:frq_lim[2])`: number of frequencies
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

# Returns

Named tuple containing:
- `p_env::Array{Float64, 3}`: power spectrum envelope
- `p_env_frq::Vector{Float64}`: frequencies for each envelope
"""
function penv(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), d::Int64=8, method::Symbol=:welch, nt::Int64=8, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), frq_n::Int64=_tlength(frq_lim), frq::Symbol=:lin, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)
    
    _check_channels(obj, ch)

    ch_n = length(ch)
    ep_n = nepochs(obj)
    fs = sr(obj)

    pw, pf = psd(obj.data[1, :, 1], fs=fs, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)

    p_env = zeros(ch_n, length(pw), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            pw, _ = psd(obj.data[ch[ch_idx], :, ep_idx], fs=fs, norm=true, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)
            # find peaks
            p_idx = findpeaks(pw, d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(pw))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(pf[p_idx], pw[p_idx])
                try
                    p_env[ch_idx, :, ep_idx] = model(pf)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            else
                p_env[ch_idx, :, ep_idx] = pw
            end
            p_env[ch_idx, 1, ep_idx] = p_env[ch_idx, 2, ep_idx]
        end
    end
    
    return (p_env=p_env, p_env_frq=pf)

end

"""
    penv_mean(obj; ch, dims, d, method, nt, wlen, woverlap, w)

Calculate power (in dB) envelope: mean and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=8`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds
- `frq_n::Int64=length(frq_lim[1]:frq_lim[2])`: number of frequencies
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

# Returns

Named tuple containing:
- `p_env_m::Array{Float64, 3}`: power spectrum envelope: mean
- `p_env_u::Array{Float64, 3}`: power spectrum envelope: 95% CI upper bound
- `p_env_l::Array{Float64, 3}`: power spectrum envelope: 95% CI lower bound
- `p_env_frq::Vector{Float64}`: power spectrum envelope (useful for plotting over PSD)
"""
function penv_mean(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=8, method::Symbol=:welch, nt::Int64=8, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), frq_n::Int64=_tlength(frq_lim), frq::Symbol=:lin, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)

    if dims == 1
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
    elseif dims == 2
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    elseif dims == 3
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    end

    pw, pf = psd(obj, ch=ch, norm=true, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)

    ch_n = size(pw, 1)
    ep_n = size(pw, 3)

    if dims == 1
        # mean over channels

        p_env_m = zeros(length(pf), ep_n)
        p_env_u = zeros(length(pf), ep_n)
        p_env_l = zeros(length(pf), ep_n)

        @inbounds @simd for ep_idx in 1:ep_n
            p_env_m[:, ep_idx] = mean(pw[:, :, ep_idx], dims=1)
            # find peaks
            p_idx = findpeaks(p_env_m[:, ep_idx], d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(p_env_m[:, ep_idx]))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(pf[p_idx], p_env_m[p_idx])
                try
                    p_env_m[:, ep_idx] = model(pf)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = std(p_env_m[:, ep_idx]) / sqrt(length(p_env_m[:, ep_idx]))
            p_env_u[:, ep_idx] = @. p_env_m[:, ep_idx] + 1.96 * s
            p_env_l[:, ep_idx] = @. p_env_m[:, ep_idx] - 1.96 * s
        end
    elseif dims == 2
        # mean over epochs

        p_env_m = zeros(length(pf), ch_n)
        p_env_u = zeros(length(pf), ch_n)
        p_env_l = zeros(length(pf), ch_n)

        @inbounds @simd for ch_idx in 1:ch_n
            p_env_m[:, ch_idx] = mean(pw[ch_idx, :, :], dims=2)
            # find peaks
            p_idx = findpeaks(p_env_m[:, ch_idx], d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(p_env_m[:, ch_idx]))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(pf[p_idx], p_env_m[p_idx])
                try
                    p_env_m[:, ch_idx] = model(pf)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = std(p_env_m[:, ch_idx]) / sqrt(length(p_env_m[:, ch_idx]))
            p_env_u[:, ch_idx] = @views @. p_env_m[:, ch_idx] + 1.96 * s
            p_env_l[:, ch_idx] = @views @. p_env_m[:, ch_idx] - 1.96 * s
        end
    else
        # mean over channels and epochs

        p_env_m, p_env_u, p_env_l, _ = penv_mean(obj, dims=1, d=d)
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
    penv_median(obj; ch, dims, d, method, nt, wlen, woverlap, w)

Calculate power (in dB) envelope: median and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `dims::Int64`: median over channels (dims = 1) or epochs (dims = 2)
- `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=8`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds
- `frq_n::Int64=length(frq_lim[1]:frq_lim[2])`: number of frequencies
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

# Returns

Named tuple containing:
- `p_env_m::Array{Float64, 3}`: power spectrum envelope: median
- `p_env_u::Array{Float64, 3}`: power spectrum envelope: 95% CI upper bound
- `p_env_l::Array{Float64, 3}`: power spectrum envelope: 95% CI lower bound
- `p_env_frq::Vector{Float64}`: power spectrum envelope (useful for plotting over PSD)
"""
function penv_median(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=8, method::Symbol=:welch, nt::Int64=8, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), frq_n::Int64=_tlength(frq_lim), frq::Symbol=:lin, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)

    if dims == 1
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
    elseif dims == 2
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    elseif dims == 3
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    end

    pw, pf = psd(obj, ch=ch, norm=true, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)

    ch_n = size(pw, 1)
    ep_n = size(pw, 3)

    if dims == 1
        # median over channels

        p_env_m = zeros(length(pf), ep_n)
        p_env_u = zeros(length(pf), ep_n)
        p_env_l = zeros(length(pf), ep_n)

        @inbounds @simd for ep_idx in 1:ep_n
            p_env_m[:, ep_idx] = median(pw[:, :, ep_idx], dims=1)
            # find peaks
            p_idx = findpeaks(p_env_m[:, ep_idx], d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(p_env_m[:, ep_idx]))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(pf[p_idx], p_env_m[p_idx])
                try
                    p_env_m[:, ep_idx] = model(pf)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(p_env_m[:, ep_idx]) / sqrt(length(p_env_m[:, ep_idx]))
            p_env_u[:, ep_idx] = @. p_env_m[:, ep_idx] + 1.96 * s
            p_env_l[:, ep_idx] = @. p_env_m[:, ep_idx] - 1.96 * s
        end
    elseif dims == 2
        # median over epochs

        p_env_m = zeros(length(pf), ch_n)
        p_env_u = zeros(length(pf), ch_n)
        p_env_l = zeros(length(pf), ch_n)

        @inbounds @simd for ch_idx in 1:ch_n
            p_env_m[:, ch_idx] = median(pw[ch_idx, :, :], dims=2)
            # find peaks
            p_idx = findpeaks(p_env_m[:, ch_idx], d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(p_env_m[:, ch_idx]))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(pf[p_idx], p_env_m[p_idx])
                try
                    p_env_m[:, ch_idx] = model(pf)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(p_env_m[:, ch_idx]) / sqrt(length(p_env_m[:, ch_idx]))
            p_env_u[:, ch_idx] = @views @. p_env_m[:, ch_idx] + 1.96 * s
            p_env_l[:, ch_idx] = @views @. p_env_m[:, ch_idx] - 1.96 * s
        end
    else
        # median over channels and epochs
        
        p_env_m, p_env_u, p_env_l, _ = penv_median(obj, dims=1, d=d)
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
    senv(obj; ch, d, t, pad, method, frq_lim, frq_n, norm, nt, frq, gw, ncyc, wt, wlen, woverlap, w, wt, gw)

Calculate spectral envelope.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
- `t::Union{Real, Nothing}=nothing`: spectrogram threshold (maximize all powers > t)
- `method::Symbol=:stft`: method of calculating spectrogram:
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `pad::Int64=0`: number of zeros to add
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency limits
- `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
- `norm::Bool=true`: normalize powers to dB
- `nt::Int64=8`: number of Slepian tapers
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: `ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc = linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, e.g. `wt = wavelet(Morlet(2π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:
- `s_env::Array{Float64, 3}`: spectral envelope
- `s_env_t::Vector{Float64}`: spectrogram time
"""
function senv(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), d::Int64=2, t::Union{Real, Nothing}=nothing, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), frq_n::Int64=_tlength(frq_lim), pad::Int64=0, method::Symbol=:stft, norm::Bool=true, nt::Int64=8, frq::Symbol=:log, gw::Real=5, ncyc::Union{Int64, Tuple{Int64, Int64}}=6, wt::T=wavelet(Morlet(2π), β=2), wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true) where {T <: CWT}
    
    _check_channels(obj, ch)

    ch_n = length(ch)
    ep_n = nepochs(obj)
    fs = sr(obj)

    if method === :stft
        sp, _, _ = @views NeuroAnalyzer.spectrogram(obj.data[1, :, 1], fs=fs, norm=norm, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
    elseif method === :mt
        sp, _, _ = @views NeuroAnalyzer.spectrogram(obj.data[1, :, 1], fs=fs, norm=norm, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
    elseif method === :mw
    _, sp, _, _ = @views NeuroAnalyzer.mwspectrogram(obj.data[1, :, 1], pad=pad, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)
    elseif method === :gh
        sp, _, _ = @views NeuroAnalyzer.ghspectrogram(obj.data[1, :, 1], fs=fs, frq_lim=frq_lim, frq_n=frq_n, norm=norm, frq=frq, gw=gw)
    elseif method === :cwt
        sp, _ = @views NeuroAnalyzer.cwtspectrogram(obj.data[1, :, 1], wt=wt, fs=fs, frq_lim=frq_lim, norm=norm)
    end
    st = linspace(0, (epoch_len(obj) / fs), size(sp, 2))
    st .+= obj.epoch_time[1]

    s_env = zeros(ch_n, length(st), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            # prepare spectrogram
            if method === :stft
                sp, sf, _ = @views NeuroAnalyzer.spectrogram(obj.data[ch[ch_idx], :, ep_idx], fs=fs, norm=norm, method=:stft, wlen=wlen, woverlap=woverlap, w=w)
            elseif method === :mt
                sp, sf, _ = @views NeuroAnalyzer.spectrogram(obj.data[ch[ch_idx], :, ep_idx], fs=fs, norm=norm, method=:mt, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
            elseif method === :mw
            _, sp, _, sf = @views NeuroAnalyzer.mwspectrogram(obj.data[ch[ch_idx], :, ep_idx], pad=pad, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)
            elseif method === :gh
                sp, _, sf = @views NeuroAnalyzer.ghspectrogram(obj.data[ch[ch_idx], :, ep_idx], fs=fs, frq_lim=frq_lim, frq_n=frq_n, norm=norm, frq=frq, gw=gw)
            elseif method === :cwt
                sp, sf = @views NeuroAnalyzer.cwtspectrogram(obj.data[ch[ch_idx], :, ep_idx], wt=wt, fs=fs, frq_lim=frq_lim, norm=norm)
            end

            # maximize all powers above threshold (t)
            if t !== nothing
                sp[sp .> t] .= 0
                reverse!(sp)
                reverse!(sf)
            end
            
            f_idx = zeros(length(st))
            m = maximum(sp, dims=1)
            for idx2 in eachindex(m)
                f_idx[idx2] = sf[vsearch(m[idx2], sp[:, idx2])]
            end
            # find peaks
            p_idx = findpeaks(f_idx, d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(st))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(st[p_idx], f_idx[p_idx])
                try
                    s_env[ch_idx, :, ep_idx] = model(st)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            else
                s_env[ch_idx, :, ep_idx] = f_idx
            end
            s_env[ch_idx, 1, ep_idx] = s_env[ch_idx, 2, ep_idx]
        end
    end
    
    return (s_env=s_env, senv_t=st)

end

"""
    senv_mean(obj; ch, dims, d, t, pad, method, frq_lim, frq_n, norm, nt, frq, gw, ncyc, wt, wlen, woverlap, w, wt, gw)

Calculate spectral envelope: mean and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
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
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency limits
- `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
- `norm::Bool=true`: normalize powers to dB
- `nt::Int64=8`: number of Slepian tapers
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: `ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc = linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, e.g. `wt = wavelet(Morlet(2π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
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
function senv_mean(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=2, t::Union{Real, Nothing}=nothing, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), frq_n::Int64=_tlength(frq_lim), method::Symbol=:stft, pad::Int64=0, norm::Bool=true, nt::Int64=8, frq::Symbol=:log, gw::Real=5, ncyc::Union{Int64, Tuple{Int64, Int64}}=6, wt::T=wavelet(Morlet(2π), β=2), wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true) where {T <: CWT}

    if dims == 1
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
    elseif dims == 2
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    elseif dims == 3
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    end

    sp, st = senv(obj, ch=ch, d=d, t=t, pad=pad, method=method, frq_lim=frq_lim, frq_n=frq_n, norm=norm, nt=nt, frq=frq, gw=gw, ncyc=ncyc, wt=wt, wlen=wlen, woverlap=woverlap, w=w)

    ch_n = size(sp, 1)
    ep_n = size(sp, 3)

    if dims == 1
        # mean over channels

        s_env_m = zeros(length(st), ep_n)
        s_env_u = zeros(length(st), ep_n)
        s_env_l = zeros(length(st), ep_n)

        @inbounds @simd for ep_idx in 1:ep_n
            s_env_m[:, ep_idx] = mean(sp[:, :, ep_idx], dims=1)
            # find peaks
            s_idx = findpeaks(s_env_m[:, ep_idx], d=d)
            # add first time-point
            pushfirst!(s_idx, 1)
            # interpolate peaks using cubic spline or loess
            push!(s_idx, length(s_env_m[:, ep_idx]))
            if length(s_idx) > 4
                model = CubicSpline(st[s_idx], s_env_m[s_idx])
                try
                    s_env_m[:, ep_idx] = model(st)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = std(s_env_m[:, ep_idx]) / sqrt(length(s_env_m[:, ep_idx]))
            s_env_u[:, ep_idx] = @. s_env_m[:, ep_idx] + 1.96 * s
            s_env_l[:, ep_idx] = @. s_env_m[:, ep_idx] - 1.96 * s
        end
    elseif dims == 2
        # mean over epochs

        s_env_m = zeros(length(st), ch_n)
        s_env_u = zeros(length(st), ch_n)
        s_env_l = zeros(length(st), ch_n)

        @inbounds @simd for ch_idx in 1:ch_n
            s_env_m[:, ch_idx] = mean(sp[ch_idx, :, :], dims=2)
            # find peaks
            s_idx = findpeaks(s_env_m[:, ch_idx], d=d)
            # add first time-point
            pushfirst!(s_idx, 1)
            # add last time-point
            push!(s_idx, length(s_env_m[:, ch_idx]))
            # interpolate peaks using cubic spline or loess
            if length(s_idx) > 4
                model = CubicSpline(st[s_idx], s_env_m[s_idx])
                try
                    s_env_m[:, ch_idx] = model(st)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = std(s_env_m[:, ch_idx]) / sqrt(length(s_env_m[:, ch_idx]))
            s_env_u[:, ch_idx] = @views @. s_env_m[:, ch_idx] + 1.96 * s
            s_env_l[:, ch_idx] = @views @. s_env_m[:, ch_idx] - 1.96 * s
        end
    else
        # mean over channels and epochs

        s_env_m, s_env_u, s_env_l, _ = senv_mean(obj, dims=1, d=d, pad=pad, method=method, frq_lim=frq_lim, frq_n=frq_n, norm=norm, nt=nt, frq=frq, gw=gw, ncyc=ncyc, wt=wt, wlen=wlen, woverlap=woverlap, w=w)
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
    senv_median(obj; ch, dims, d, t, pad, method, frq_lim, frq_n, norm, nt, frq, gw, ncyc, wt, wlen, woverlap, w, wt, gw)

Calculate spectral envelope: median and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
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
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency limits
- `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
- `norm::Bool=true`: normalize powers to dB
- `nt::Int64=8`: number of Slepian tapers
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: `ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc = linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, e.g. `wt = wavelet(Morlet(2π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
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
function senv_median(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=2, t::Union{Real, Nothing}=nothing, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), frq_n::Int64=_tlength(frq_lim), method::Symbol=:stft, pad::Int64=0, norm::Bool=true, nt::Int64=8, frq::Symbol=:log, gw::Real=5, ncyc::Union{Int64, Tuple{Int64, Int64}}=6, wt::T=wavelet(Morlet(2π), β=2), wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true) where {T <: CWT}

    if dims == 1
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
    elseif dims == 2
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    elseif dims == 3
        @assert nchannels(obj) >= 2 "Number of channels must be ≥ 2."
        @assert nepochs(obj) >= 2 "Number of epochs must be ≥ 2."
    end

    sp, st = senv(obj, ch=ch, d=d, t=t, pad=pad, method=method, nt=nt, frq_lim=frq_lim, frq_n=frq_n, norm=norm, frq=frq, gw=gw, ncyc=ncyc, wt=wt, wlen=wlen, woverlap=woverlap, w=w)

    ch_n = size(sp, 1)
    ep_n = size(sp, 3)

    if dims == 1
        # median over channels

        s_env_m = zeros(length(st), ep_n)
        s_env_u = zeros(length(st), ep_n)
        s_env_l = zeros(length(st), ep_n)

        @inbounds @simd for ep_idx in 1:ep_n
            s_env_m[:, ep_idx] = median(sp[:, :, ep_idx], dims=1)
            # find peaks
            s_idx = findpeaks(s_env_m[:, ep_idx], d=d)
            # add first time-point
            pushfirst!(s_idx, 1)
            # add last time-point
            push!(s_idx, length(s_env_m[:, ep_idx]))
            # interpolate peaks using cubic spline or loess
            if length(s_idx) > 4
                model = CubicSpline(st[s_idx], s_env_m[s_idx])
                try
                    s_env_m[:, ep_idx] = model(st)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(s_env_m[:, ep_idx]) / sqrt(length(s_env_m[:, ep_idx]))
            s_env_u[:, ep_idx] = @. s_env_m[:, ep_idx] + 1.96 * s
            s_env_l[:, ep_idx] = @. s_env_m[:, ep_idx] - 1.96 * s
        end
    elseif dims == 2
        # median over epochs

        s_env_m = zeros(length(st), ch_n)
        s_env_u = zeros(length(st), ch_n)
        s_env_l = zeros(length(st), ch_n)

        @inbounds @simd for ch_idx in 1:ch_n
            s_env_m[:, ch_idx] = median(sp[ch_idx, :, :], dims=2)
            # find peaks
            s_idx = findpeaks(s_env_m[:, ch_idx], d=d)
            # add first time-point
            pushfirst!(s_idx, 1)
            # add last time-point
            push!(s_idx, length(s_env_m[:, ch_idx]))
            # interpolate peaks using cubic spline or loess
            if length(s_idx) > 4
                model = CubicSpline(st[s_idx], s_env_m[s_idx])
                try
                    s_env_m[:, ch_idx] = model(st)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(s_env_m[:, ch_idx]) / sqrt(length(s_env_m[:, ch_idx]))
            s_env_u[:, ch_idx] = @views @. s_env_m[:, ch_idx] + 1.96 * s
            s_env_l[:, ch_idx] = @views @. s_env_m[:, ch_idx] - 1.96 * s
        end
    else
        # median over channels and epochs

        s_env_m, s_env_u, s_env_l, _ = senv_median(obj, dims=1, d=d, pad=pad, method=method, frq_lim=frq_lim, frq_n=frq_n, norm=norm, nt=nt, frq=frq, gw=gw, ncyc=ncyc, wt=wt, wlen=wlen, woverlap=woverlap, w=w)
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
    henv(obj; ch, d)

Calculate Hilbert spectrum amplitude envelope.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `h_env::Array{Float64, 3}`: Hilbert spectrum amplitude envelope
- `s_t::Vector{Float64}`: signal time
"""
function henv(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), d::Int64=32)

    _check_channels(obj, ch)

    _, hamp, _, _ = @views hspectrum(obj.data[ch, :, :])

    ch_n = size(hamp, 1)
    ep_n = size(hamp, 3)
    h_env = similar(hamp)

    s_t = obj.epoch_time

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s = @view hamp[ch_idx, :, ep_idx]
            # find peaks
            p_idx = findpeaks(s, d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(s))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(s_t[p_idx], s[p_idx])
                try
                    h_env[ch_idx, :, ep_idx] = model(s_t)
                catch
                    @warn "CubicSpline error, using Loess."
                    model = Loess.loess(s_t[p_idx], s[p_idx], span=0.5)
                    h_env[ch_idx, :, ep_idx] = Loess.predict(model, s_t)
                end
            else
                _info("Less than 5 peaks detected, using Loess.")
                model = Loess.loess(s_t[p_idx], s[p_idx], span=0.5)
                h_env[ch_idx, :, ep_idx] = Loess.predict(model, s_t)
            end
            h_env[ch_idx, 1, ep_idx] = h_env[ch_idx, 2, ep_idx]
        end
    end
    
    return (h_env=h_env, s_t=s_t)
end

"""
    henv_mean(obj; ch, dims, d)

Calculate Hilbert spectrum amplitude envelope: mean and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `h_env_m::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: mean
- `h_env_u::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI upper bound
- `h_env_l::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function henv_mean(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=32)
    
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

        @inbounds @simd for ep_idx in 1:ep_n
            h_env_m[:, ep_idx] = mean(s_a[:, :, ep_idx], dims=1)
            s = std(h_env_m[:, ep_idx]) / sqrt(length(h_env_m[:, ep_idx]))
            h_env_u[:, ep_idx] = @. h_env_m[:, ep_idx] + 1.96 * s
            h_env_l[:, ep_idx] = @. h_env_m[:, ep_idx] - 1.96 * s
        end
    elseif dims == 2
        # mean over epochs

        h_env_m = zeros(length(s_t), ch_n)
        h_env_u = zeros(length(s_t), ch_n)
        h_env_l = zeros(length(s_t), ch_n)

        @inbounds @simd for ch_idx in 1:ch_n
            h_env_m[:, ch_idx] = mean(s_a[ch_idx, :, :], dims=2)
            s = std(h_env_m[:, ch_idx]) / sqrt(length(h_env_m[:, ch_idx]))
            h_env_u[:, ch_idx] = @views @. h_env_m[:, ch_idx] + 1.96 * s
            h_env_l[:, ch_idx] = @views @. h_env_m[:, ch_idx] - 1.96 * s
        end
    else
        # mean over channels and epochs

        h_env_m, h_env_u, h_env_l, _ = henv_mean(obj, dims=1, d=d)
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
    henv_median(obj; ch, dims, d)

Calculate Hilbert spectrum amplitude envelope of `obj`: median and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `dims::Int64`: median over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `h_env_m::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: median
- `h_env_u::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI upper bound
- `h_env_l::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function henv_median(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=32)
    
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

        @inbounds @simd for ep_idx in 1:ep_n
            h_env_m[:, ep_idx] = median(s_a[:, :, ep_idx], dims=1)
            t_idx = findpeaks(h_env_m[:, ep_idx], d=d)
            pushfirst!(t_idx, 1)
            push!(t_idx, length(h_env_m[:, ep_idx]))
            if length(t_idx) > 4
                model = CubicSpline(s_t[t_idx], h_env_m[t_idx])
                try
                    h_env_m[:, ep_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(h_env_m[:, ep_idx]) / sqrt(length(h_env_m[:, ep_idx]))
            h_env_u[:, ep_idx] = @views @. h_env_m[:, ep_idx] + 1.96 * s
            h_env_l[:, ep_idx] = @views @. h_env_m[:, ep_idx] - 1.96 * s
        end
    elseif dims == 2
        # median over epochs

        h_env_m = zeros(length(s_t), ch_n)
        h_env_u = zeros(length(s_t), ch_n)
        h_env_l = zeros(length(s_t), ch_n)

        @inbounds @simd for ch_idx in 1:ch_n
            h_env_m[:, ch_idx] = median(s_a[ch_idx, :, :], dims=2)
            t_idx = findpeaks(h_env_m[:, ch_idx], d=d)
            pushfirst!(t_idx, 1)
            push!(t_idx, length(h_env_m[:, ch_idx]))
            if length(t_idx) > 4
                model = CubicSpline(s_t[t_idx], h_env_m[t_idx])
                try
                    h_env_m[:, ch_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(h_env_m[:, ch_idx]) / sqrt(length(h_env_m[:, ch_idx]))
            h_env_u[:, ch_idx] = @views @. h_env_m[:, ch_idx] + 1.96 * s
            h_env_l[:, ch_idx] = @views @. h_env_m[:, ch_idx] - 1.96 * s
        end
    else
        # median over channels and epochs

        h_env_m, h_env_u, h_env_l, _ = henv_median(obj, dims=1, d=d)
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
- `ec::Vector{Float64}`: power correlation value
- `p::Vector{Float64}`: p-value
"""
function env_cor(env1::Array{Float64, 3}, env2::Array{Float64, 3})

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
