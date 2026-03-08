export coherence

"""
    coherence(s1, s2; <keyword arguments>)

Calculate coherence, imaginary part of coherence and magnitude-squared coherence (MSC).

For two signals `s1`, `s2` and their cross-power spectra:
- coh = S12 / √(S11 · S22) (complex coherence)
- imcoh = Im(coh) (imaginary part — insensitive to zero-lag volume conduction)
- msc = |coh|² (magnitude-squared coherence ∈ [0,1])

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector
- `method::Symbol=:mt`: method used to calculate CPSD:
  - `:mt`: multi-tapered cross-power spectra
  - `:fft`: fast Fourier transformation
  - `:stft`: short time Fourier transformation
- `fs::Int64`: sampling rate
- `flim::Tuple{Real, Real}=(0, fs / 2)`: frequency bounds
- `demean::Bool=false`: if true, the channel-wise mean will be subtracted from the input signals before the cross spectral powers are computed
- `nt::Int64=16`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:

- `coh::Vector{ComplexF64}`: coherence
- `imcoh::Vector{Float64}`: imaginary part of coherence
- `msc::Vector{Float64}`: magnitude-squared coherence
- `f::Vector{Float64}`: frequencies
"""
function coherence(
    s1::AbstractVector,
    s2::AbstractVector;
    method::Symbol = :mt,
    fs::Int64,
    flim::Tuple{Real, Real} = (0, fs / 2),
    demean::Bool = false,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
)::@NamedTuple{coh::Vector{ComplexF64}, imcoh::Vector{Float64}, msc::Vector{Float64}, f::Vector{Float64}}

    # check parameters
    _check_var(method, [:mt, :fft, :stft], "method")
    s1, s2 = _veqlen(s1, s2)
    @assert nt >= 1 "nt must be ≥ 1."
    @assert fs >= 1 "fs must be ≥ 1."
    @assert wlen <= length(s1) "wlen must be ≤ $(length(s1))."
    @assert wlen >= 2 "wlen must be ≥ 2."
    @assert woverlap < wlen "woverlap must be < $(wlen)."
    @assert woverlap >= 0 "woverlap must be ≥ 0."
    _check_tuple(flim, (0, fs / 2), "flim")

    # shared kwargs for all three cpsd calls — defined once to keep them in sync
    cpsd_kwargs = (
        fs = fs,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        demean = demean,
        method = method,
        flim = flim,
    )

    # compute the three cross-power spectra needed for coherence
    # S11 = auto-spectrum of s1, S12 = cross-spectrum, S22 = auto-spectrum of s2.
    s1s1, f = cpsd(s1, s1; cpsd_kwargs...)
    # f is the same for all three calls
    s1s2, _ = cpsd(s1, s2; cpsd_kwargs...)
    s2s2, _ = cpsd(s2, s2; cpsd_kwargs...)

    # complex coherence: S12 / √(S11 · S22).
    coh = @. s1s2 / sqrt(s1s1 * s2s2)

    # imaginary coherence: insensitive to instantaneous zero-lag coupling
    # (e.g. volume conduction), captures only time-lagged interactions
    imcoh = imag.(coh)

    # magnitude-squared coherence: real-valued, bounded on [0, 1].
    msc = abs2.(coh)

    return (coh = coh, imcoh = imcoh, msc = msc, f = f)

end

"""
    coherence(s1, s2; <keyword arguments>)

Calculate coherence, imaginary part of coherence and magnitude-squared coherence (MSC).

For two signals `s1`, `s2` and their cross-power spectra:
- coh = S12 / √(S11 · S22) (complex coherence)
- imcoh = Im(coh) (imaginary part — insensitive to zero-lag volume conduction)
- msc = |coh|² (magnitude-squared coherence ∈ [0,1])

# Arguments

- `s1::AbstractArray`: signal array
- `s2::AbstractArray`: signal array
- `method::Symbol=:mt`: method used to calculate CPSD:
  - `:mt`: multi-tapered cross-power spectra
  - `:fft`: fast Fourier transformation
  - `:stft`: short time Fourier transformation
- `fs::Int64`: sampling rate
- `flim::Tuple{Real, Real}=(0, fs / 2)`: frequency bounds
- `demean::Bool=false`: if true, the channel-wise mean will be subtracted from the input signals before the cross spectral powers are computed
- `nt::Int64=16`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window

# Returns

- `coh::Array{ComplexF64, 3}`: coherence of shape `(channels, frequencies, epochs)`
- `imcoh::Array{Float64, 3}`: imaginary part of coherence of shape `(channels, frequencies, epochs)`
- `msc::Array{Float64, 3}`: magnitude-squared coherence of shape `(channels, frequencies, epochs)`
- `f::Vector{Float64}`: frequencies in Hz
"""
function coherence(
    s1::AbstractArray,
    s2::AbstractArray;
    method::Symbol = :mt,
    fs::Int64,
    flim::Tuple{Real, Real} = (0, fs / 2),
    demean::Bool = false,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
)::@NamedTuple{coh::Array{ComplexF64, 3}, imcoh::Array{Float64, 3}, msc::Array{Float64, 3}, f::Vector{Float64}}

    # validate shape
    @assert size(s1) == size(s2) "s1 and s2 must have the same size."
    _chk3d(s1)
    _chk3d(s2)

    # number of channels
    ch_n = size(s1, 1)
    # number of epochs
    ep_n = size(s1, 3)

    # pre-compute the frequency vector with a single pilot call on the first channel/epoch pair
    result = NeuroAnalyzer.coherence(
        @view(s1[1, :, 1]),
        @view(s2[1, :, 1]),
        method = method,
        fs = fs,
        flim = flim,
        demean = demean,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
    )
    f = result.f

    # pre-allocate outputs
    coh = zeros(ComplexF64, ch_n, length(f), ep_n)
    imcoh = zeros(ch_n, length(f), ep_n)
    msc = zeros(ch_n, length(f), ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        result = coherence(
            @view(s1[ch_idx, :, ep_idx]),
            @view(s2[ch_idx, :, ep_idx]),
            method = method,
            fs = fs,
            flim = flim,
            demean = demean,
            nt = nt,
            wlen = wlen,
            woverlap = woverlap,
            w = w,
        )
        coh[ch_idx, :, ep_idx] = result.coh
        imcoh[ch_idx, :, ep_idx] = result.imcoh
        msc[ch_idx, :, ep_idx] = result.msc
    end

    return (coh = coh, imcoh = imcoh, msc = msc, f = f)
end

"""
    coherence(obj1, obj2; <keyword arguments>)

Calculate coherence, imaginary part of coherence and magnitude-squared coherence (MSC).

For two signals `s1`, `s2` and their cross-power spectra:
- coh = S12 / √(S11 · S22) (complex coherence)
- imcoh = Im(coh) (imaginary part — insensitive to zero-lag volume conduction)
- msc = |coh|² (magnitude-squared coherence ∈ [0,1])

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}}`: channel name(s)
- `ch2::Union{String, Vector{String}}`: channel name(s)
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))` epoch number(s)
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))` epoch number(s)
- `method::Symbol=:mt`: method used to calculate CPSD:
  - `:mt`: multi-tapered cross-power spectra
  - `:fft`: fast Fourier transformation
  - `:stft`: short time Fourier transformation
- `fs::Int64`: sampling rate
- `flim::Tuple{Real, Real}=(0, fs / 2)`: frequency bounds
- `demean::Bool=false`: if true, the channel-wise mean will be subtracted from the input signals before the cross spectral powers are computed
- `nt::Int64=16`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window

# Returns

- `coh::Array{ComplexF64, 3}`: coherence of shape `(channels, frequencies, epochs)`
- `imcoh::Array{Float64, 3}`: imaginary part of coherence of shape `(channels, frequencies, epochs)`
- `msc::Array{Float64, 3}`: magnitude-squared coherence of shape `(channels, frequencies, epochs)`
- `f::Vector{Float64}`: frequencies
"""
function coherence(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}},
    ch2::Union{String, Vector{String}},
    ep1::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj1)),
    ep2::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj2)),
    method::Symbol = :mt,
    flim::Tuple{Real, Real} = (0, sr(obj1) / 2),
    demean::Bool = false,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj1),
    woverlap::Int64 = round(Int64, wlen * 0.9),
    w::Bool = true,
)::@NamedTuple{coh::Array{ComplexF64, 3}, imcoh::Array{Float64, 3}, msc::Array{Float64, 3}, f::Vector{Float64}}

    # validate objects
    @assert sr(obj1) == sr(obj2) "OBJ1 and OBJ2 must have the same sampling rate."
    @assert length(ch1) == length(ch2) "Lengths of ch1 ($(length(ch1))) and ch2 ($(length(ch2))) must be equal."
    @assert length(ep1) == length(ep2) "Lengths of ep1 ($(length(ep1))) and ep2 ($(length(ep2))) must be equal."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    # resolve channel names to integer indices, optionally skipping bad channels
    ch1 = exclude_bads ? get_channel(obj1, ch = ch1, exclude = "bad") : get_channel(obj1, ch = ch1, exclude = "")
    ch2 = exclude_bads ? get_channel(obj2, ch = ch2, exclude = "bad") : get_channel(obj2, ch = ch2, exclude = "")

    # validate epoch numbers
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    # normalize scalar epoch arguments to vectors so indexing is uniform
    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])

    return coherence(
        @view(obj1.data[ch1, :, ep1]),
        @view(obj2.data[ch2, :, ep2]),
        method = method,
        fs = sr(obj1),
        flim = flim,
        demean = demean,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
    )

end
