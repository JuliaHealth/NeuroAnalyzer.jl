export coherence

"""
    coherence(s1, s2; <keyword arguments>)

Calculate coherence and MSC (magnitude-squared coherence).

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`
- `method::Symbol=:mt`: method used to calculate CPSD:
    - `:mt`: multi-tapered cross-power spectra
    - `:fft`: fast Fourier transformation
- `fs::Int64`: sampling rate
- `frq_lim::Tuple{Real, Real}=(0, fs / 2)`: frequency bounds
- `demean::Bool=false`: if true, the channel-wise mean will be subtracted from the input signals before the cross spectral powers are computed
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window

# Returns

- `coh::Vector{Float64}`: coherence
- `mscoh::Vector{Float64}`: magnitude-squared coherence
- `p::Vector{Float64}`: frequencies
"""
function coherence(s1::AbstractVector, s2::AbstractVector; method::Symbol=:mt, fs::Int64, frq_lim::Tuple{Real, Real}=(0, fs / 2), demean::Bool=false, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)

    _check_var(method, [:mt, :fft], "method")
    s1, s2 = _veqlen(s1, s2)
    @assert nt >= 1 "nt must be ≥ 1."
    @assert fs >= 1 "fs must be ≥ 1."
    @assert wlen <= length(s1) "wlen must be ≤ $(length(s1))."
    @assert wlen >= 2 "wlen must be ≥ 2."
    @assert woverlap <= wlen "woverlap must be ≤ $(wlen)."
    @assert woverlap >= 0 "woverlap must be ≥ 0."
    _check_tuple(frq_lim, "frq_lim", (0, fs / 2))

    if method === :mt
        w = w ? hanning(length(s1)) : nothing
        # multitaper
        s = hcat(s1, s2)'
        n_samples = length(s1)
        c = mt_coherence(s, fs=fs, demean=demean, nfft=nextpow(2, n_samples), window=nothing, nw=((nt + 1) ÷ 2), ntapers=nt)
        f = DSP.freq(c)
        coh = DSP.coherence(c)
        f1_idx = vsearch(frq_lim[1], f)
        f2_idx = vsearch(frq_lim[2], f)
        f = f[f1_idx:f2_idx]
        coh = @views coh[1, 2, f1_idx:f2_idx]
        mscoh = coh.^2
    elseif method === :fft
        # fft
        s1s1, f = cpsd(s1, s1, fs=fs, wlen=wlen, woverlap=woverlap, w=w, demean=demean)
        s1s2, f = cpsd(s1, s2, fs=fs, wlen=wlen, woverlap=woverlap, w=w, demean=demean)
        s2s2, f = cpsd(s2, s2, fs=fs, wlen=wlen, woverlap=woverlap, w=w, demean=demean)
        f1_idx = vsearch(frq_lim[1], f)
        f2_idx = vsearch(frq_lim[2], f)
        f = f[f1_idx:f2_idx]
        s1s1 = @views s1s1[f1_idx:f2_idx]
        s2s2 = @views s2s2[f1_idx:f2_idx]
        s1s2 = @views s1s2[f1_idx:f2_idx]
        coh = @. abs(s1s2) / sqrt(s1s1 * s2s2)
        mscoh = @. (abs(s1s2))^2 / (s1s1 * s2s2)
    end

    return (coh=coh, mscoh=mscoh, f=f)

end

"""
    coherence(s1, s2; <keyword arguments>)

Calculate coherence and MSC (magnitude-squared coherence).

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`
- `method::Symbol=:mt`: method used to calculate CPSD:
    - `:mt`: multi-tapered cross-power spectra
    - `:fft`: fast Fourier transformation
- `fs::Int64`: sampling rate
- `frq_lim::Tuple{Real, Real}=(0, fs / 2)`: frequency bounds
- `demean::Bool=false`: if true, the channel-wise mean will be subtracted from the input signals before the cross spectral powers are computed
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window

# Returns

- `coh::Array{Float64, 3}`: coherence
- `mscoh::Array{Float64, 3}`: magnitude-squared coherence
- `p::Vector{Float64}`: frequencies
"""
function coherence(s1::AbstractArray, s2::AbstractArray; method::Symbol=:mt, fs::Int64, frq_lim::Tuple{Real, Real}=(0, fs / 2), demean::Bool=false, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."

    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    _, _, f = coherence(s1[1, :, 1], s2[1, :, 1]; method=method, fs=fs, frq_lim=frq_lim, demean=demean, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
    if frq_lim !== nothing
        _check_tuple(frq_lim, "frq_lim", (0, fs / 2))
        idx1 = vsearch(frq_lim[1], f)
        idx2 = vsearch(frq_lim[2], f)
        f = f[idx1:idx2]
    end

    coh = zeros(ch_n, length(f), ep_n)
    mscoh = zeros(ch_n, length(f), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            coh[ch_idx, :, ep_idx], mscoh[ch_idx, :, ep_idx], _ = @views coherence(s1[ch_idx, :, ep_idx], s2[ch_idx, :, ep_idx], method=method, fs=fs, frq_lim=frq_lim, demean=demean, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
        end
    end

    return (coh=coh, mscoh=mscoh, f=f)
end

"""
    coherence(obj1, obj2; <keyword arguments>)

Calculate coherence and MSC (magnitude-squared coherence).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs
- `method::Symbol=:mt`: method used to calculate CPSD:
    - `:mt`: multi-tapered cross-power spectra
    - `:fft`: fast Fourier transformation
- `fs::Int64`: sampling rate
- `frq_lim::Tuple{Real, Real}=(0, fs / 2)`: frequency bounds
- `demean::Bool=false`: if true, the channel-wise mean will be subtracted from the input signals before the cross spectral powers are computed
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window

# Returns

- `coh::Array{Float64, 3}`: coherence
- `mscoh::Array{Float64, 3}`: magnitude-squared coherence
- `p::Vector{Float64}`: frequencies
"""
function coherence(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1), ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2), ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2)), method::Symbol=:mt, frq_lim::Tuple{Real, Real}=(0, sr(obj1) / 2), demean::Bool=false, nt::Int64=7, wlen::Int64=sr(obj1), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)

    _check_channels(obj1, ch1)
    _check_channels(obj2, ch2)
    @assert length(ch1) == length(ch2) "ch1 and ch2 must have the same length."

    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    @assert length(ep1) == length(ep2) "ep1 and ep2 must have the same length."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    @assert sr(obj1) == sr(obj2) "OBJ1 and OBJ2 must have the same sampling rate."

    isa(ch1, Int64) && (ch1 = [ch1])
    isa(ch2, Int64) && (ch2 = [ch2])
    length(ep1) == 1 && (ep1 = [ep1])
    length(ep2) == 1 && (ep2 = [ep2])

    coh, mscoh, f = @views coherence(obj1.data[ch1, :, ep1], obj2.data[ch2, :, ep2], method=method, fs=sr(obj1), frq_lim=frq_lim, demean=demean, nt=nt, wlen=wlen, woverlap=woverlap, w=w)

    return (coh=coh, mscoh=mscoh, f=f)

end
