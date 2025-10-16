export cpsd

"""
    cpsd(s1, s2; <keyword arguments>)

Calculate cross power spectral density (CPSD).

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`
- `method::Symbol=:mt`: method used to calculate CPSD:
    - `:mt`: multi-tapered cross-power spectra
    - `:fft`: fast Fourier transformation
    - `:stft`: short time Fourier transformation
- `fs::Int64`: sampling rate
- `frq_lim::Tuple{Real, Real}=(0, fs / 2)`: frequency bounds
- `demean::Bool=false`: if true, the channel-wise mean will be subtracted from the input signals before the cross spectral powers are computed
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:
- `pxy::Vector{ComplexF64}`: cross-power spectrum
- `f::Vector{Float64}`: frequencies
"""
function cpsd(s1::AbstractVector, s2::AbstractVector; method::Symbol=:mt, fs::Int64, frq_lim::Tuple{Real, Real}=(0, fs / 2), demean::Bool=false, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)::@NamedTuple{pxy::Vector{ComplexF64}, f::Vector{Float64}}

    _check_var(method, [:mt, :fft, :stft], "method")
    s1, s2 = _veqlen(s1, s2)
    @assert nt >= 1 "nt must be ≥ 1."
    @assert fs >= 1 "fs must be ≥ 1."
    @assert wlen <= length(s1) "wlen must be ≤ $(length(s1))."
    @assert wlen >= 2 "wlen must be ≥ 2."
    @assert woverlap < wlen "woverlap must be < $(wlen)."
    @assert woverlap >= 0 "woverlap must be ≥ 0."
    _check_tuple(frq_lim, "frq_lim", (0, fs / 2))

    n_samples = length(s1)

    if method === :mt
        w = w ? DSP.hanning(n_samples) : ones(n_samples)
        s = hcat(s1 .* w, s2 .* w)'
        pxy = DSP.mt_cross_power_spectra(s, fs=fs, demean=demean, nfft=nextpow(2, n_samples), nw=((nt + 1) ÷ 2), ntapers=nt)
        f = DSP.freq(pxy)
        pxy = DSP.power(pxy)
        f1_idx = vsearch(frq_lim[1], f)
        f2_idx = vsearch(frq_lim[2], f)
        f = f[f1_idx:f2_idx]
        pxy = pxy[1, 2, f1_idx:f2_idx]
    elseif method === :stft
        chunks_idx = NeuroAnalyzer._fchunks(length(s1), wlen=wlen, woverlap=woverlap)
        pxy = zeros(ComplexF64, wlen)
        w = w ? hanning(wlen) : ones(wlen)
        for idx in axes(chunks_idx, 1)
            if demean
                s1_tmp = @views detrend(s1[chunks_idx[idx, 1]:chunks_idx[idx, 2]], type=:mean) .* w
                s2_tmp = @views detrend(s2[chunks_idx[idx, 1]:chunks_idx[idx, 2]], type=:mean) .* w
            else
                s1_tmp = @views s1[chunks_idx[idx, 1]:chunks_idx[idx, 2]] .* w
                s2_tmp = @views s2[chunks_idx[idx, 1]:chunks_idx[idx, 2]] .* w
            end
            ss1 = fft(s1_tmp) / length(s1_tmp)
            ss2 = fft(s2_tmp) / length(s2_tmp)
            pxy .+= (conj.(ss1) .* ss2)
        end
        # get mean over segments
        pxy /= size(chunks_idx, 1)
        f, _ = freqs(wlen, fs)
        f1_idx = vsearch(frq_lim[1], f)
        f2_idx = vsearch(frq_lim[2], f)
        f = f[f1_idx:f2_idx]
        pxy = pxy[f1_idx:f2_idx]
    elseif method === :fft
        w = w ? hanning(n_samples) : ones(n_samples)
        if demean
            s1_tmp = detrend(s1, type=:mean) .* w
            s2_tmp = detrend(s2, type=:mean) .* w
        else
            s1_tmp = s1 .* w
            s2_tmp = s2 .* w
        end
        # fft
        ss1 = fft(s1_tmp) / n_samples
        ss2 = fft(s2_tmp) / n_samples
        f, _ = freqs(s1, fs)
        pxy = conj.(ss1) .* ss2
        f1_idx = vsearch(frq_lim[1], f)
        f2_idx = vsearch(frq_lim[2], f)
        f = f[f1_idx:f2_idx]
        pxy = pxy[f1_idx:f2_idx]
    end

    return (pxy=pxy, f=f)

end

"""
    cpsd(s1, s2; <keyword arguments>)

Calculate cross power spectral density (CPSD).

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`
- `method::Symbol=:mt`: method used to calculate CPSD:
    - `:mt`: multi-tapered cross-power spectra
    - `:fft`: fast Fourier transformation
    - `:stft`: short time Fourier transformation
- `fs::Int64`: sampling rate
- `frq_lim::Tuple{Real, Real}=(0, fs / 2)`: frequency bounds
- `demean::Bool=false`: if true, the channel-wise mean will be subtracted from the input signals before the cross spectral powers are computed
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:
- `pxy::Array{ComplexF64, 3}`: cross-power spectrum
- `f::Vector{Float64}`: frequencies
"""
function cpsd(s1::AbstractArray, s2::AbstractArray; method::Symbol=:mt, fs::Int64, frq_lim::Tuple{Real, Real}=(0, fs / 2), demean::Bool=false, nt::Int64=7, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)::@NamedTuple{pxy::Array{ComplexF64, 3}, f::Vector{Float64}}

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."
    _chk3d(s1)
    _chk3d(s2)

    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    _, f = cpsd(s1[1, :, 1], s2[1, :, 1]; method=method, fs=fs, frq_lim=frq_lim, demean=demean, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
    if frq_lim !== nothing
        _check_tuple(frq_lim, "frq_lim", (0, fs / 2))
        idx1 = vsearch(frq_lim[1], f)
        idx2 = vsearch(frq_lim[2], f)
        f = f[idx1:idx2]
    end

    pxy = zeros(ComplexF64, ch_n, length(f), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            pxy[ch_idx, :, ep_idx], _ = @views cpsd(s1[ch_idx, :, ep_idx], s2[ch_idx, :, ep_idx], method=method, fs=fs, frq_lim=frq_lim, demean=demean, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
        end
    end

    return (pxy=pxy, f=f)
end

"""
    cpsd(obj1, obj2; <keyword arguments>)

Calculate cross power spectral density (CPSD).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}}`: list of channels
- `ch2::Union{String, Vector{String}}`: list of channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs
- `method::Symbol=:mt`: method used to calculate CPSD:
    - `:mt`: multi-tapered cross-power spectra
    - `:fft`: fast Fourier transformation
    - `:stft`: short time Fourier transformation
- `frq_lim::Tuple{Real, Real}=(0, sr(obj1) / 2)`: frequency bounds
- `demean::Bool=false`: if true, the channel-wise mean will be subtracted from the input signals before the cross spectral powers are computed
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj1)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:
- `pxy::Array{ComplexF64, 3}`: cross-power spectrum
- `f::Vector{Float64}`: frequencies
"""
function cpsd(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{String, Vector{String}}, ch2::Union{String, Vector{String}}, ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2)), method::Symbol=:mt, frq_lim::Tuple{Real, Real}=(0, sr(obj1) / 2), demean::Bool=false, nt::Int64=7, wlen::Int64=sr(obj1), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true)::@NamedTuple{pxy::Array{ComplexF64, 3}, f::Vector{Float64}}

    ch1 = exclude_bads ? get_channel(obj1, ch=ch1, exclude="bad") : get_channel(obj1, ch=ch1, exclude="")
    ch2 = exclude_bads ? get_channel(obj2, ch=ch2, exclude="bad") : get_channel(obj2, ch=ch2, exclude="")
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])

    @assert sr(obj1) == sr(obj2) "OBJ1 and OBJ2 must have the same sampling rate."
    @assert length(ch1) == length(ch2) "Lengths of ch1 ($(length(ch1))) and ch2 ($(length(ch2))) must be equal."
    @assert length(ep1) == length(ep2) "Lengths of ep1 ($(length(ep1))) and ep2 ($(length(ep2))) must be equal."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    pxy, f = @views cpsd(obj1.data[ch1, :, ep1], obj2.data[ch2, :, ep2], method=method, fs=sr(obj1), frq_lim=frq_lim, demean=demean, nt=nt, wlen=wlen, woverlap=woverlap, w=w)

    return (pxy=pxy, f=f)

end
