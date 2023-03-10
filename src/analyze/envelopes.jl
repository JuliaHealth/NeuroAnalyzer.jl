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
    tenv(obj; channel, d)

Calculate temporal envelope (amplitude).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `t_env::Array{Float64, 3}`: temporal envelope
- `s_t::Vector{Float64}`: signal time
"""
function tenv(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), d::Int64=32)
    
    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    t_env = zeros(ch_n, epoch_len(obj), ep_n)
    s_t = obj.epoch_time

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s = @view obj.data[channel[ch_idx], :, ep_idx]
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
                    @error "CubicSpline error, using Loess."
                    model = loess(s_t[p_idx], s[p_idx], span=0.5)
                    t_env[ch_idx, :, ep_idx] = Loess.predict(model, s_t)
                end
            else
                _info("Less than 5 peaks detected, using Loess.")
                model = loess(s_t[p_idx], s[p_idx], span=0.5)
                t_env[ch_idx, :, ep_idx] = Loess.predict(model, s_t)
            end
            t_env[ch_idx, 1, ep_idx] = t_env[ch_idx, 2, ep_idx]
        end
    end
    
    return (t_env=t_env, s_t=s_t)
end

"""
    tenv_mean(obj; channel, dims, d)

Calculate temporal envelope (amplitude): mean and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `t_env_m::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: mean
- `t_env_u::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI upper bound
- `t_env_l::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function tenv_mean(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=32)
    
    if dims == 1
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_a, s_t = tenv(obj, channel=channel, d=d)
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
            t_env_u[:, ch_idx] = @. t_env_m[:, ch_idx] + 1.96 * s
            t_env_l[:, ch_idx] = @. t_env_m[:, ch_idx] - 1.96 * s
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
    tenv_median(obj; channel, dims, d)

Calculate temporal envelope (amplitude): median and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `dims::Int64`: median over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `t_env_m::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: median
- `t_env_u::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI upper bound
- `t_env_l::Union{Vector{Float64}, Matrix{Float64}}`: temporal envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function tenv_median(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=32)
    
    if dims == 1
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_a, s_t = tenv(obj, channel=channel, d=d)
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
            t_env_u[:, ch_idx] = @. t_env_m[:, ch_idx] + 1.96 * s
            t_env_l[:, ch_idx] = @. t_env_m[:, ch_idx] - 1.96 * s
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
    penv(obj; channel, d)

Calculate power (in dB) envelope.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
- `mt::Bool=false`: if true use multi-tapered periodogram
- `nt::Int64=8`: number of Slepian tapers

# Returns

Named tuple containing:
- `p_env::Array{Float64, 3}`: power spectrum envelope
- `p_env_frq::Vector{Float64}`: frequencies for each envelope
"""
function penv(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), d::Int64=8, mt::Bool=false, nt::Int64=8)
    
    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    fs = sr(obj)

    psd_tmp, frq = psd(obj.data[1, :, 1], fs=fs, mt=mt, nt=nt)
    p_env = zeros(ch_n, length(psd_tmp), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            psd_pow, _ = psd(obj.data[channel[ch_idx], :, ep_idx], fs=fs, mt=mt, norm=true, nt=nt)
            # find peaks
            p_idx = findpeaks(psd_pow, d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(psd_pow))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(frq[p_idx], psd_pow[p_idx])
                try
                    p_env[ch_idx, :, ep_idx] = model(frq)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            else
                p_env[ch_idx, :, ep_idx] = psd_pow
            end
            p_env[ch_idx, 1, ep_idx] = p_env[ch_idx, 2, ep_idx]
        end
    end
    
    return (p_env=p_env, p_env_frq=frq)
end

"""
    penv_mean(obj; channel, dims, d)

Calculate power (in dB) envelope: mean and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

Named tuple containing:
- `p_env_m::Array{Float64, 3}`: power spectrum envelope: mean
- `p_env_u::Array{Float64, 3}`: power spectrum envelope: 95% CI upper bound
- `p_env_l::Array{Float64, 3}`: power spectrum envelope: 95% CI lower bound
- `p_env_frq::Vector{Float64}`: power spectrum envelope (useful for plotting over PSD)
"""
function penv_mean(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=8, mt::Bool=false)
    
    if dims == 1
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_p, s_f = psd(obj, channel=channel, norm=true, mt=mt)
    ch_n = size(s_p, 1)
    ep_n = size(s_p, 3)

    if dims == 1
        # mean over channels

        p_env_m = zeros(length(s_f), ep_n)
        p_env_u = zeros(length(s_f), ep_n)
        p_env_l = zeros(length(s_f), ep_n)

        @inbounds @simd for ep_idx in 1:ep_n
            p_env_m[:, ep_idx] = mean(s_p[:, :, ep_idx], dims=1)
            # find peaks
            p_idx = findpeaks(p_env_m[:, ep_idx], d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(p_env_m[:, ep_idx]))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(s_f[p_idx], p_env_m[p_idx])
                try
                    p_env_m[:, ep_idx] = model(s_f)
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

        p_env_m = zeros(length(s_f), ch_n)
        p_env_u = zeros(length(s_f), ch_n)
        p_env_l = zeros(length(s_f), ch_n)

        @inbounds @simd for ch_idx in 1:ch_n
            p_env_m[:, ch_idx] = mean(s_p[ch_idx, :, :], dims=2)
            # find peaks
            p_idx = findpeaks(p_env_m[:, ch_idx], d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(p_env_m[:, ch_idx]))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(s_f[p_idx], p_env_m[p_idx])
                try
                    p_env_m[:, ch_idx] = model(s_f)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = std(p_env_m[:, ch_idx]) / sqrt(length(p_env_m[:, ch_idx]))
            p_env_u[:, ch_idx] = @. p_env_m[:, ch_idx] + 1.96 * s
            p_env_l[:, ch_idx] = @. p_env_m[:, ch_idx] - 1.96 * s
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
    
    return (p_env_m=p_env_m, p_env_u=p_env_u, p_env_l=p_env_l, p_env_frq=s_f)
end

"""
    penv_median(obj; channel, dims, d)

Calculate power (in dB) envelope: median and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `dims::Int64`: median over channels (dims = 1) or epochs (dims = 2)
- `d::Int64=8`: distance between peeks in samples, lower values get better envelope fit
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

Named tuple containing:
- `p_env_m::Array{Float64, 3}`: power spectrum envelope: median
- `p_env_u::Array{Float64, 3}`: power spectrum envelope: 95% CI upper bound
- `p_env_l::Array{Float64, 3}`: power spectrum envelope: 95% CI lower bound
- `p_env_frq::Vector{Float64}`: power spectrum envelope (useful for plotting over PSD)
"""
function penv_median(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=8, mt::Bool=false)
    
    if dims == 1
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_p, s_f = psd(obj, channel=channel, norm=true, mt=mt)
    ch_n = size(s_p, 1)
    ep_n = size(s_p, 3)

    if dims == 1
        # median over channels

        p_env_m = zeros(length(s_f), ep_n)
        p_env_u = zeros(length(s_f), ep_n)
        p_env_l = zeros(length(s_f), ep_n)

        @inbounds @simd for ep_idx in 1:ep_n
            p_env_m[:, ep_idx] = median(s_p[:, :, ep_idx], dims=1)
            # find peaks
            p_idx = findpeaks(p_env_m[:, ep_idx], d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(p_env_m[:, ep_idx]))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(s_f[p_idx], p_env_m[p_idx])
                try
                    p_env_m[:, ep_idx] = model(s_f)
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

        p_env_m = zeros(length(s_f), ch_n)
        p_env_u = zeros(length(s_f), ch_n)
        p_env_l = zeros(length(s_f), ch_n)

        @inbounds @simd for ch_idx in 1:ch_n
            p_env_m[:, ch_idx] = median(s_p[ch_idx, :, :], dims=2)
            # find peaks
            p_idx = findpeaks(p_env_m[:, ch_idx], d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(p_env_m[:, ch_idx]))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(s_f[p_idx], p_env_m[p_idx])
                try
                    p_env_m[:, ch_idx] = model(s_f)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(p_env_m[:, ch_idx]) / sqrt(length(p_env_m[:, ch_idx]))
            p_env_u[:, ch_idx] = @. p_env_m[:, ch_idx] + 1.96 * s
            p_env_l[:, ch_idx] = @. p_env_m[:, ch_idx] - 1.96 * s
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
    
    return (p_env_m=p_env_m, p_env_u=p_env_u, p_env_l=p_env_l, p_env_frq=s_f)
end

"""
    senv(obj; channel, d, mt, t)

Calculate spectral envelope.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
- `mt::Bool=false`: if true use multi-tapered spectrogram
- `t::Union{Real, Nothing}=nothing`: spectrogram threshold (maximize all powers > t)

# Returns

Named tuple containing:
- `s_env::Array{Float64, 3}`: spectral envelope
- `s_env_t::Vector{Float64}`: spectrogram time
"""
function senv(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), d::Int64=2, mt::Bool=false, t::Union{Real, Nothing}=nothing)
    
    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    fs = sr(obj)
    s_tmp = @view obj.data[1, :, 1]
    interval = fs
    overlap = round(Int64, fs * 0.75)
    # for short signals always use multi-taper
    length(s_tmp) < 4 * fs && (mt = true)
    if mt == true
        spec_tmp = mt_spectrogram(s_tmp, fs=fs)
    else
        spec_tmp = DSP.spectrogram(s_tmp, interval, overlap, nfft=length(s_tmp), fs=fs, window=hanning)
    end
    sp_t = collect(spec_tmp.time)
    sp_t .+= obj.epoch_time[1]

    s_env = zeros(ch_n, length(sp_t), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            # prepare spectrogram
            if mt == true
                spec = @views mt_spectrogram(obj.data[channel[ch_idx], :, ep_idx], fs=fs)
            else
                spec = @views DSP.spectrogram(obj.data[channel[ch_idx], :, ep_idx], interval, overlap, nfft=length(s_tmp), fs=fs, window=hanning)
            end
            s_frq = Vector(spec.freq)
            s_p = pow2db.(spec.power)

            # maximize all powers above threshold (t)
            if t !== nothing
                s_p[s_p .> t] .= 0
                reverse!(s_p)
                reverse!(s_frq)
            end
            
            f_idx = zeros(length(spec.time))
            m = maximum(s_p, dims=1)
            for idx2 in eachindex(m)
                f_idx[idx2] = s_frq[vsearch(m[idx2], s_p[:, idx2])]
            end
            # find peaks
            p_idx = findpeaks(f_idx, d=d)
            # add first time-point
            pushfirst!(p_idx, 1)
            # add last time-point
            push!(p_idx, length(spec.time))
            # interpolate peaks using cubic spline or loess
            if length(p_idx) >= 5
                model = CubicSpline(sp_t[p_idx], f_idx[p_idx])
                try
                    s_env[ch_idx, :, ep_idx] = model(sp_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            else
                s_env[ch_idx, :, ep_idx] = f_idx
            end
            s_env[ch_idx, 1, ep_idx] = s_env[ch_idx, 2, ep_idx]
        end
    end
    
    return (s_env=s_env, senv_t=sp_t)
end

"""
    senv_mean(obj; channel, dims, d, mt, t)

Calculate spectral envelope: mean and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
- `mt::Bool=false`: if true use multi-tapered spectrogram
- `t::Union{Real, Nothing}=nothing`: spectrogram threshold (maximize all powers > t)

# Returns

Named tuple containing:
- `s_env_m::Array{Float64, 3}`: spectral envelope: mean
- `s_env_u::Array{Float64, 3}`: spectral envelope: 95% CI upper bound
- `s_env_l::Array{Float64, 3}`: spectral envelope: 95% CI lower bound
- `s_env_t::Vector{Float64}`: spectral envelope (useful for plotting over spectrogram)
"""
function senv_mean(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=2, mt::Bool=false, t::Union{Real, Nothing}=nothing)

    if dims == 1
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_p, s_t = senv(obj, channel=channel, d=d, mt=mt, t=t)
    ch_n = size(s_p, 1)
    ep_n = size(s_p, 3)

    if dims == 1
        # mean over channels

        s_env_m = zeros(length(s_t), ep_n)
        s_env_u = zeros(length(s_t), ep_n)
        s_env_l = zeros(length(s_t), ep_n)

        @inbounds @simd for ep_idx in 1:ep_n
            s_env_m[:, ep_idx] = mean(s_p[:, :, ep_idx], dims=1)
            # find peaks
            s_idx = findpeaks(s_env_m[:, ep_idx], d=d)
            # add first time-point
            pushfirst!(s_idx, 1)
            # interpolate peaks using cubic spline or loess
            push!(s_idx, length(s_env_m[:, ep_idx]))
            if length(s_idx) > 4
                model = CubicSpline(s_t[s_idx], s_env_m[s_idx])
                try
                    s_env_m[:, ep_idx] = model(s_t)
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

        s_env_m = zeros(length(s_t), ch_n)
        s_env_u = zeros(length(s_t), ch_n)
        s_env_l = zeros(length(s_t), ch_n)

        @inbounds @simd for ch_idx in 1:ch_n
            s_env_m[:, ch_idx] = mean(s_p[ch_idx, :, :], dims=2)
            # find peaks
            s_idx = findpeaks(s_env_m[:, ch_idx], d=d)
            # add first time-point
            pushfirst!(s_idx, 1)
            # add last time-point
            push!(s_idx, length(s_env_m[:, ch_idx]))
            # interpolate peaks using cubic spline or loess
            if length(s_idx) > 4
                model = CubicSpline(s_t[s_idx], s_env_m[s_idx])
                try
                    s_env_m[:, ch_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = std(s_env_m[:, ch_idx]) / sqrt(length(s_env_m[:, ch_idx]))
            s_env_u[:, ch_idx] = @. s_env_m[:, ch_idx] + 1.96 * s
            s_env_l[:, ch_idx] = @. s_env_m[:, ch_idx] - 1.96 * s
        end
    else
        # mean over channels and epochs

        s_env_m, s_env_u, s_env_l, _ = senv_mean(obj, dims=1, d=d, mt=mt)
        s_env_m = mean(s_env_m, dims=2)
        s_env_u = mean(s_env_u, dims=2)
        s_env_l = mean(s_env_l, dims=2)
        s_env_m = reshape(s_env_m, size(s_env_m, 1))
        s_env_u = reshape(s_env_u, size(s_env_u, 1))
        s_env_l = reshape(s_env_l, size(s_env_l, 1))
    end
    
    return (s_env_m=s_env_m, s_env_u=s_env_u, s_env_l=s_env_l, s_env_t=s_t)
end

"""
    senv_median(obj; channel, dims, d, mt)

Calculate spectral envelope: median and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `dims::Int64`: median over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=2`: distance between peeks in samples, lower values get better envelope fit
- `mt::Bool=false`: if true use multi-tapered spectrogram
- `t::Union{Real, Nothing}=nothing`: spectrogram threshold (maximize all powers > t)

# Returns

Named tuple containing:
- `s_env_m::Array{Float64, 3}`: spectral envelope: median
- `s_env_u::Array{Float64, 3}`: spectral envelope: 95% CI upper bound
- `s_env_l::Array{Float64, 3}`: spectral envelope: 95% CI lower bound
- `s_env_t::Vector{Float64}`: spectral envelope (useful for plotting over spectrogram)
"""
function senv_median(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=2, mt::Bool=false, t::Union{Real, Nothing}=nothing)
    
    if dims == 1
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_p, s_t = senv(obj, channel=channel, d=d, mt=mt, t=t)
    ch_n = size(s_p, 1)
    ep_n = size(s_p, 3)

    if dims == 1
        # median over channels

        s_env_m = zeros(length(s_t), ep_n)
        s_env_u = zeros(length(s_t), ep_n)
        s_env_l = zeros(length(s_t), ep_n)

        @inbounds @simd for ep_idx in 1:ep_n
            s_env_m[:, ep_idx] = median(s_p[:, :, ep_idx], dims=1)
            # find peaks
            s_idx = findpeaks(s_env_m[:, ep_idx], d=d)
            # add first time-point
            pushfirst!(s_idx, 1)
            # add last time-point
            push!(s_idx, length(s_env_m[:, ep_idx]))
            # interpolate peaks using cubic spline or loess
            if length(s_idx) > 4
                model = CubicSpline(s_t[s_idx], s_env_m[s_idx])
                try
                    s_env_m[:, ep_idx] = model(s_t)
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

        s_env_m = zeros(length(s_t), ch_n)
        s_env_u = zeros(length(s_t), ch_n)
        s_env_l = zeros(length(s_t), ch_n)

        @inbounds @simd for ch_idx in 1:ch_n
            s_env_m[:, ch_idx] = median(s_p[ch_idx, :, :], dims=2)
            # find peaks
            s_idx = findpeaks(s_env_m[:, ch_idx], d=d)
            # add first time-point
            pushfirst!(s_idx, 1)
            # add last time-point
            push!(s_idx, length(s_env_m[:, ch_idx]))
            # interpolate peaks using cubic spline or loess
            if length(s_idx) > 4
                model = CubicSpline(s_t[s_idx], s_env_m[s_idx])
                try
                    s_env_m[:, ch_idx] = model(s_t)
                catch
                    _info("CubicSpline could not be calculated, using non-smoothed variant instead.")
                end
            end
            s = iqr(s_env_m[:, ch_idx]) / sqrt(length(s_env_m[:, ch_idx]))
            s_env_u[:, ch_idx] = @. s_env_m[:, ch_idx] + 1.96 * s
            s_env_l[:, ch_idx] = @. s_env_m[:, ch_idx] - 1.96 * s
        end
    else
        # median over channels and epochs

        s_env_m, s_env_u, s_env_l, _ = senv_median(obj, dims=1, d=d, mt=mt)
        s_env_m = median(s_env_m, dims=2)
        s_env_u = median(s_env_u, dims=2)
        s_env_l = median(s_env_l, dims=2)
        s_env_m = reshape(s_env_m, size(s_env_m, 1))
        s_env_u = reshape(s_env_u, size(s_env_u, 1))
        s_env_l = reshape(s_env_l, size(s_env_l, 1))
    end
    
    return (s_env_m=s_env_m, s_env_u=s_env_u, s_env_l=s_env_l, s_env_t=s_t)
end

"""
    henv(obj; channel, d)

Calculate Hilbert spectrum amplitude envelope.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `h_env::Array{Float64, 3}`: Hilbert spectrum amplitude envelope
- `s_t::Vector{Float64}`: signal time
"""
function henv(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), d::Int64=32)

    _check_channels(obj, channel)
    _, signal, _, _ = @views spectrum(keep_channel(obj, channel=channel), h=true)

    ch_n = size(signal, 1)
    ep_n = size(signal, 3)
    h_env = similar(signal)
    s_t = obj.epoch_time

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s = @view signal[ch_idx, :, ep_idx]
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
                    @error "CubicSpline error, using Loess."
                    model = loess(s_t[p_idx], s[p_idx], span=0.5)
                    h_env[ch_idx, :, ep_idx] = Loess.predict(model, s_t)
                end
            else
                _info("Less than 5 peaks detected, using Loess.")
                model = loess(s_t[p_idx], s[p_idx], span=0.5)
                h_env[ch_idx, :, ep_idx] = Loess.predict(model, s_t)
            end
            h_env[ch_idx, 1, ep_idx] = h_env[ch_idx, 2, ep_idx]
        end
    end
    
    return (h_env=h_env, s_t=s_t)
end

"""
    henv_mean(obj; channel, dims, d)

Calculate Hilbert spectrum amplitude envelope: mean and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `dims::Int64`: mean over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `h_env_m::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: mean
- `h_env_u::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI upper bound
- `h_env_l::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function henv_mean(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=32)
    
    if dims == 1
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_a, s_t = henv(obj, channel=channel, d=d)
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
            h_env_u[:, ch_idx] = @. h_env_m[:, ch_idx] + 1.96 * s
            h_env_l[:, ch_idx] = @. h_env_m[:, ch_idx] - 1.96 * s
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
    henv_median(obj; channel, dims, d)

Calculate Hilbert spectrum amplitude envelope of `obj`: median and 95% CI.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `dims::Int64`: median over channels (dims = 1), epochs (dims = 2) or channels and epochs (dims = 3)
- `d::Int64=32`: distance between peeks in samples, lower values get better envelope fit

# Returns

Named tuple containing:
- `h_env_m::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: median
- `h_env_u::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI upper bound
- `h_env_l::Union{Vector{Float64}, Matrix{Float64}}`: Hilbert spectrum amplitude envelope: 95% CI lower bound
- `s_t::Vector{Float64}`: signal time
"""
function henv_median(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), dims::Int64, d::Int64=32)
    
    if dims == 1
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
    elseif dims == 2
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    elseif dims == 3
        channel_n(obj) == 1 && throw(ArgumentError("Number of channels must be ≥ 2."))
        epoch_n(obj) == 1 && throw(ArgumentError("Number of epochs must be ≥ 2."))
    end

    s_a, s_t = henv(obj, channel=channel, d=d)
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
            h_env_u[:, ep_idx] = @. h_env_m[:, ep_idx] + 1.96 * s
            h_env_l[:, ep_idx] = @. h_env_m[:, ep_idx] - 1.96 * s
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
            h_env_u[:, ch_idx] = @. h_env_m[:, ch_idx] + 1.96 * s
            h_env_l[:, ch_idx] = @. h_env_m[:, ch_idx] - 1.96 * s
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
    env_cor(obj1, obj2; type, channel1, channel2, epoch1, epoch2)

Calculate envelope correlation.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `type::Symbol=:amp`: envelope type:
    - `:amp`: amplitude
    - `:pow`: power
    - `:spec`: spectrogram
    - `:hamp`: Hilbert spectrum amplitude
- `channel1::Int64`
- `channel2::Int64`
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `ec::Vector{Float64}`: power correlation value
- `ec_p::Vector{Float64}`: power correlation p-value
"""
function env_cor(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; type::Symbol=:amp, channel1::Int64, channel2::Int64, epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)))

    _check_var(type, [:amp, :pow, :spec, :hamp], "type")

    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))

    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("OBJ1 and OBJ2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    
    ec_r = zeros(ep_n)
    ec_p = zeros(ep_n)

    # calculate envelopes
    if type === :amp
        s1, _ = tenv(obj1, channel=channel1)
        s2, _ = tenv(obj2, channel=channel2)
    elseif type === :pow
        s1, _ = penv(obj1, channel=channel1)
        s2, _ = penv(obj2, channel=channel2)
    elseif type === :spec
        s1, _ = senv(obj1, channel=channel1)
        s2, _ = senv(obj2, channel=channel2)
    elseif type === :hamp
        s1, _ = henv(obj1, channel=channel1)
        s2, _ = henv(obj2, channel=channel2)
    end
    s1 = s1[:, :, epoch1]
    s2 = s2[:, :, epoch2]
    
    # compare envelopes per epochs
    Threads.@threads for ep_idx in 1:ep_n
        ec = CorrelationTest(vec(s1[:, :, ep_idx]), vec(s2[:, :, ep_idx]))
        @inbounds ec_r[ep_idx] = ec.r
        @inbounds ec_p[ep_idx] = pvalue(ec)
    end

    return (ec=ec_r, ec_p=ec_p)
end
