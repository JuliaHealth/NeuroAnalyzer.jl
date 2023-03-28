export psd_mw

"""
    psd_mw(s; pad, norm, frq_lim, frq_n, frq, fs, ncyc)

Calculate power spectrum using Morlet wavelet convolution.

# Arguments

- `s::AbstractVector`
- `pad::Int64=0`: pad with `pad` zeros
- `norm::Bool=true`: normalize powers to dB
- `fs::Int64`: sampling rate
- `frq_lim::Tuple{Real, Real}=(0, fs ÷ 2)`: frequency bounds for the spectrogram
- `frq_n::Int64=length(frq_lim[1]:frq_lim[2])`: number of frequencies
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

# Returns

Named tuple containing:
- `pw::Matrix{Float64}`: powers
- `pf::Vector{Float64}`: frequencies
"""
function psd_mw(s::AbstractVector; pad::Int64=0, norm::Bool=true, fs::Int64, frq_lim::Tuple{Real, Real}=(0, fs ÷ 2), frq_n::Int64=_tlength(frq_lim), frq::Symbol=:lin, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)

    _check_var(frq, [:log, :lin], "frq")
    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    pad < 0 && throw(ArgumentError("pad must be ≥ 0."))
    frq_lim = tuple_order(frq_lim)
    frq_lim[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
    frq_lim[2] > fs ÷ 2 && throw(ArgumentError("Upper frequency bound must be ≤ $(fs ÷ 2)."))
    frq_n < 2 && throw(ArgumentError("frq_n must be ≥ 2."))
    frq_lim[1] == 0 && (frq_lim = (0.1, frq_lim[2]))

    if frq === :log
        frq_lim = (frq_lim[1], frq_lim[2])
        pf = round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), frq_n), digits=1)
    else
        pf = linspace(frq_lim[1], frq_lim[2], frq_n)
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

    ch_n = size(s, 1)
    pw = zeros(ch_n, length(pf))

    pad > 0 && (s = pad0(s, pad))
    @inbounds @simd for frq_idx in 1:frq_n
        kernel = generate_morlet(fs, pf[frq_idx], 1, ncyc=ncyc[frq_idx], complex=true)
        w_conv = fconv(s, kernel=kernel, norm=true)
        pw[frq_idx] = mean(@. abs(w_conv)^2)
    end

    pw = pw[1:length(pf)]
    norm == true && (pw = pow2db.(pw))

    return (pw=pw, pf=pf)

end

"""
    psd_mw(s; pad, norm, frq_lim, frq_n, frq, fs, ncyc)

Calculate power spectrum using Morlet wavelet convolution.

# Arguments

- `s::AbstractMatrix`
- `pad::Int64=0`: pad with `pad` zeros
- `norm::Bool=true`: normalize powers to dB
- `fs::Int64`: sampling rate
- `frq_lim::Tuple{Real, Real}=(0, fs ÷ 2)`: frequency bounds for the spectrogram
- `frq_n::Int64=10`: number of frequencies
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc = linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

# Returns

Named tuple containing:
- `pw::Array{Float64, 2}`: powers
- `pf::Vector{Float64}`: frequencies
"""
function psd_mw(s::AbstractMatrix; pad::Int64=0, norm::Bool=true, fs::Int64, frq_lim::Tuple{Real, Real}=(0, fs ÷ 2), frq_n::Int64=_tlength(frq_lim), frq::Symbol=:lin, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)

    ch_n = size(s, 1)

    _, pf = psd_mw(s[1, :], pad=pad, norm=norm, fs=fs, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)
    pw = zeros(ch_n, length(pf))

    @inbounds @simd for ch_idx in 1:ch_n
        pw[ch_idx, :], _ = @views psd_mw(s[ch_idx, :], pad=pad, norm=norm, fs=fs, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)
    end

    return (pw=pw, pf=pf)

end

"""
    psd_mw(s; pad, norm, frq_lim, frq_n, frq, fs, ncyc)

Calculate power spectrum using Morlet wavelet convolution.

# Arguments

- `s::AbstractArray`
- `pad::Int64=0`: pad with `pad` zeros
- `norm::Bool=true`: normalize powers to dB
- `fs::Int64`: sampling rate
- `frq_lim::Tuple{Real, Real}=(0, fs ÷ 2)`: frequency bounds for the spectrogram
- `frq_n::Int64=10`: number of frequencies
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc = linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

# Returns

Named tuple containing:
- `pw::Array{Float64, 3}`: powers
- `pf::Vector{Float64}`: frequencies
"""
function psd_mw(s::AbstractArray; pad::Int64=0, norm::Bool=true, fs::Int64, frq_lim::Tuple{Real, Real}=(0, fs ÷ 2), frq_n::Int64=_tlength(frq_lim), frq::Symbol=:lin, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    _, pf = psd_mw(s[1, :, 1], pad=pad, norm=norm, fs=fs, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)
    pw = zeros(ch_n, length(pf), ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ch_n, 1))

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            pw[ch_idx, :, ep_idx], _ = @views psd_mw(s[ch_idx, :, ep_idx], pad=pad, norm=norm, fs=fs, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)

            # update progress bar
            progress_bar == true && next!(p)
        end
    end

    return (pw=pw, pf=pf)

end

"""
    psd_mw(obj; ch, pad, norm, frq_lim, frq_n, frq, ncyc)

Calculate power spectrum using Morlet wavelet convolution.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all s channels
- `pad::Int64=0`: pad with `pad` zeros
- `norm::Bool`=true: normalize powers to dB
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) ÷ 2)`: frequency bounds for the spectrogram
- `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc=linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`
# Returns

Named tuple containing:
- `pw::Array{Float64, 3}`: powers
- `pf::Vector{Float64}`: frequencies
"""
function psd_mw(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), pad::Int64=0, norm::Bool=true, frq_lim::Tuple{Real, Real}=(0, sr(obj) ÷ 2), frq_n::Int64=_tlength(frq_lim), frq::Symbol=:lin, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)

    _check_channels(obj, ch)

    pw, pf = @views psd_mw(obj.data[ch, :, :], pad=pad, fs=sr(obj), norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)

    return (pw=pw, pf=pf)

end
