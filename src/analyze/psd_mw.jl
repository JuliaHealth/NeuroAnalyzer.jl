export psd_mw

"""
    psd_mw(signal; pad, norm, frq_lim, frq_n, frq, fs, ncyc)

Calculate power spectrum using Morlet wavelet convolution.

# Arguments

- `signal::AbstractVector`
- `pad::Int64`: pad the `signal` with `pad` zeros
- `norm::Bool=true`: normalize powers to dB
- `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
- `frq_n::Int64=0`: number of frequencies
- `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
- `fs::Int64`: sampling rate
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n) for frq === :log or ncyc = linspace(ncyc[1], ncyc[2], frq_n) for frq === :lin

# Returns

Named tuple containing:
- `w_powers::Matrix{Float64}`
- `frq_list::Vector{Float64}`
"""
function psd_mw(signal::AbstractVector; pad::Int64=0, norm::Bool=true, frq_lim::Tuple{Real, Real}=(0, 0), frq_n::Int64=0, frq::Symbol=:lin, fs::Int64, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)

    _check_var(frq, [:log, :lin], "frq")
    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    pad < 0 && throw(ArgumentError("pad must be ≥ 0."))
    frq_lim == (0, 0) && (frq_lim = (0, fs / 2))
    frq_lim = tuple_order(frq_lim)
    frq_lim[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
    frq_lim[2] > fs ÷ 2 && throw(ArgumentError("Upper frequency bound must be ≤ $(fs ÷ 2)."))
    frq_n == 0 && (frq_n = length(frq_lim[1]:frq_lim[2]))
    frq_n < 2 && throw(ArgumentError("frq_n must be ≥ 2."))
    frq_lim[1] == 0 && (frq_lim = (0.1, frq_lim[2]))
    if frq === :log
        frq_lim = (frq_lim[1], frq_lim[2])
        frq_list = round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), frq_n), digits=1)
    else
        frq_list = linspace(frq_lim[1], frq_lim[2], frq_n)
    end

    if typeof(ncyc) != Tuple{Int64, Int64}
        ncyc = repeat([ncyc], frq_n)
    else
        if frq === :log
            ncyc = round.(Int64, logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n))
        else
            ncyc = round.(Int64, linspace(ncyc[1], ncyc[2], frq_n))
        end
    end

    ch_n = size(signal, 1)
    w_powers = zeros(ch_n, )

    pad > 0 && (signal = pad0(signal, pad))
    @inbounds @simd for frq_idx in 1:frq_n
        kernel = generate_morlet(fs, frq_list[frq_idx], 1, ncyc=ncyc[frq_idx], complex=true)
        w_conv = fconv(signal, kernel=kernel, norm=true)
        w_powers[frq_idx] = mean(@. abs(w_conv)^2)
    end

    w_powers = w_powers[1:length(frq_list)]
    norm == true && (w_powers = pow2db.(w_powers))

    return (w_powers=w_powers, frq_list=frq_list)
end


"""
    psd_mw(signal; pad, norm, frq_lim, frq_n, frq, fs, ncyc)

Calculate power spectrum using Morlet wavelet convolution.

# Arguments

- `signal::Array{Float64, 2}`
- `pad::Int64`: pad the `signal` with `pad` zeros
- `norm::Bool=true`: normalize powers to dB
- `frq_lim::Tuple{Real, Real}=(0, 0)`: frequency bounds for the spectrogram
- `frq_n::Int64=10`: number of frequencies
- `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
- `fs::Int64`: sampling rate
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n) for frq === :log or ncyc = linspace(ncyc[1], ncyc[2], frq_n) for frq === :lin

# Returns

Named tuple containing:
- `w_powers::Matrix{Float64}`
- `frq_list::Vector{Float64}`
"""
function psd_mw(signal::Matrix{Float64}; pad::Int64=0, norm::Bool=true, frq_lim::Tuple{Real, Real}=(0, 0), frq_n::Int64=10, frq::Symbol=:lin, fs::Int64, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)

    ch_n = size(signal, 1)
    w_powers, frq_list = psd_mw(signal[1, :], pad=pad, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, fs=fs, ncyc=ncyc)
    w_powers = zeros(ch_n, length(frq_list))
    frq_list = zeros(length(frq_list))

    Threads.@threads for channel_idx in 1:ch_n
        @inbounds w_powers[channel_idx, :], frq_list = @views psd_mw(signal[channel_idx, :], pad=pad, norm=false, frq_lim=frq_lim, frq_n=frq_n, frq=frq, fs=fs, ncyc=ncyc)
    end

    norm == true && (w_powers = pow2db.(w_powers))

    return (w_powers=w_powers, frq_list=frq_list)
end

"""
    psd_mw(obj; channel, pad, norm, frq_lim, frq_n, frq, ncyc)

Calculate power spectrum using Morlet wavelet convolution.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `pad::Int64`: pad the `signal` with `pad` zeros
- `norm::Bool`=true: normalize powers to dB
- `frq_lim::Tuple{Real, Real}=(0, 0)`: frequency bounds for the spectrogram
- `frq_n::Int64=10`: number of frequencies
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n) for frq === :log or ncyc = linspace(ncyc[1], ncyc[2], frq_n) for frq === :lin

# Returns

Named tuple containing:
- `w_pow::Array{Float64, 4}`
- `w_frq::Matrix{Float64}`
"""
function psd_mw(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), pad::Int64=0, norm::Bool=true, frq_lim::Tuple{Real, Real}=(0, 0), frq_n::Int64=0, frq::Symbol=:lin, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    fs = sr(obj)
    p_tmp, w_frq = @views mwpsd(obj.data[1, :, 1], fs=fs, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)
    w_pow = zeros(ch_n, length(p_tmp), ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ch_n, 1))

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            w_pow[ch_idx, :, ep_idx], _ = @views psd_mw(obj.data[channel[ch_idx], :, ep_idx], pad=pad, fs=fs, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)

            # update progress bar
            progress_bar == true && next!(p)
        end
    end

    return (w_pow=w_pow, w_frq=w_frq)
end

