export cpsd

"""
    cpsd(s1, s2; <keyword arguments>)

Calculate the complex cross power spectral density (CPSD) between two signals via one of three estimators:
- `:mt` – multi-taper (DSP.mt_cross_power_spectra)
- `:fft` – single-window FFT
- `:stft` – segmented FFT averaged over overlapping windows

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`
- `method::Symbol=:mt`: method used to calculate CPSD:
  - `:mt`: multi-tapered cross-power spectra
  - `:fft`: fast Fourier transformation
  - `:stft`: short-time Fourier transformation
- `fs::Int64`: sampling rate
- `flim::Tuple{Real, Real}=(0, fs / 2)`: frequency bounds
- `demean::Bool=false`: if true, the channel-wise mean will be subtracted from the input signals before the cross spectral powers are computed
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:

- `pxy::Vector{ComplexF64}`: cross-power spectrum
- `f::Vector{Float64}`: frequencies
"""
function cpsd(
    s1::AbstractVector,
    s2::AbstractVector;
    method::Symbol = :mt,
    fs::Int64,
    flim::Tuple{Real, Real} = (0, fs / 2),
    demean::Bool = false,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = round(Int64, wlen * 0.90),
    w::Bool = true,
)::@NamedTuple{pxy::Vector{ComplexF64}, f::Vector{Float64}}

    # check arguments
    _check_var(method, [:mt, :fft, :stft], "method")
    s1, s2 = _veqlen(s1, s2)
    @assert nt >= 1 "nt must be ≥ 1."
    @assert fs >= 1 "fs must be ≥ 1."
    @assert wlen <= length(s1) "wlen must be ≤ $(length(s1))."
    @assert wlen >= 2 "wlen must be ≥ 2."
    @assert woverlap < wlen "woverlap must be < $(wlen)."
    @assert woverlap >= 0 "woverlap must be ≥ 0."
    _check_tuple(flim, (0, fs / 2), "flim")

    n_samples = length(s1)

    if method === :mt

        # apply Hanning window (or unit window)
        win = w ? DSP.hanning(n_samples) : ones(n_samples)

        # stack signals as rows
        s = hcat(s1 .* win, s2 .* win)'

        # compute the full cross-power spectra matrix.
        pxy_mt = DSP.mt_cross_power_spectra(
            s,
            fs = fs,
            demean = demean,
            nfft = nextfastfft(n_samples),
            nw = ((nt + 1) ÷ 2),
            ntapers = nt,
        )
        f = DSP.freq(pxy_mt)
        pxy = DSP.power(pxy_mt)

        # trim to the requested frequency band.
        f1_idx = vsearch(flim[1], f)
        f2_idx = vsearch(flim[2], f)
        f = f[f1_idx:f2_idx]

        # extract the s1→s2 cross-spectrum (off-diagonal element [1,2])
        pxy = pxy[1, 2, f1_idx:f2_idx]

    elseif method === :stft

        # segment the signals into overlapping windows
        chunks_idx = _fchunks(length(s1), wlen = wlen, woverlap = woverlap)
        pxy = zeros(ComplexF64, nextpow(2, wlen + 1))

        # apply Hanning window (or unit window)
        win = w ? hanning(wlen) : ones(wlen)

        for idx in axes(chunks_idx, 1)
            if demean
                s1_tmp = remove_dc(@view(s1[chunks_idx[idx, 1]:chunks_idx[idx, 2]]))
                s2_tmp = remove_dc(@view(s2[chunks_idx[idx, 1]:chunks_idx[idx, 2]]))
            else
                s1_tmp = @view s1[chunks_idx[idx, 1]:chunks_idx[idx, 2]]
                s2_tmp = @view s2[chunks_idx[idx, 1]:chunks_idx[idx, 2]]
            end
            # FFT each segment
            # zero-pad to next power of 2, normalize by window length
            ss1 = fft0(s1_tmp .* win, nextfastfft(wlen) - wlen) / length(s1_tmp)
            ss2 = fft0(s2_tmp .* win, nextfastfft(wlen) - wlen) / length(s2_tmp)
            # accumulate: CPSD = conj(S1) * S2 for each segment
            pxy .+= conj.(ss1) .* ss2
        end

        # average cross-power over all segments
        pxy /= size(chunks_idx, 1)

        f = freqs(nextfastfft(wlen), fs)[1]

        # trim to the requested frequency band
        f1_idx = vsearch(flim[1], f)
        f2_idx = vsearch(flim[2], f)
        f = f[f1_idx:f2_idx]
        pxy = pxy[f1_idx:f2_idx]

    elseif method === :fft

        # apply Hanning window (or unit window)
        win = w ? hanning(n_samples) : ones(n_samples)
        if demean
            s1_tmp = remove_dc(s1) .* win
            s2_tmp = remove_dc(s2) .* win
        else
            s1_tmp = s1 .* win
            s2_tmp = s2 .* win
        end
        # FFT each signal
        ss1 = fft0(s1_tmp, nextfastfft(n_samples) - n_samples) / n_samples
        ss2 = fft0(s2_tmp, nextfastfft(n_samples) - n_samples) / n_samples
        f, _ = freqs(nextfastfft(n_samples), fs)
        # CPSD = conj(S1) * S2
        pxy = conj.(ss1) .* ss2
        
        # trim to the requested frequency band
        f1_idx = vsearch(flim[1], f)
        f2_idx = vsearch(flim[2], f)
        f = f[f1_idx:f2_idx]
        pxy = pxy[f1_idx:f2_idx]

    end

    return (pxy = pxy, f = f)

end

"""
    cpsd(s1, s2; <keyword arguments>)

Calculate the complex cross power spectral density (CPSD) between two signals via one of three estimators:
- `:mt` – multi-taper (DSP.mt_cross_power_spectra)
- `:fft` – single-window FFT
- `:stft` – segmented FFT averaged over overlapping windows

# Arguments

- `s1::AbstractArray`: signal array (channels × samples × epochs)
- `s2::AbstractArray`: signal array (channels × samples × epochs)
- `method::Symbol=:mt`: method used to calculate CPSD:
  - `:mt`: multi-tapered cross-power spectra
  - `:fft`: fast Fourier transformation
  - `:stft`: short-time Fourier transformation
- `fs::Int64`: sampling rate
- `flim::Tuple{Real, Real}=(0, fs / 2)`: frequency bounds
- `demean::Bool=false`: if true, the channel-wise mean will be subtracted from the input signals before the cross spectral powers are computed
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=fs`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:

- `pxy::Array{ComplexF64, 3}`: cross-power spectrum, shape `(channels, frequencies, epochs)`
- `f::Vector{Float64}`: frequencies
"""
function cpsd(
    s1::AbstractArray,
    s2::AbstractArray;
    method::Symbol = :mt,
    fs::Int64,
    flim::Tuple{Real, Real} = (0, fs / 2),
    demean::Bool = false,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = round(Int64, wlen * 0.90),
    w::Bool = true,
)::@NamedTuple{pxy::Array{ComplexF64, 3}, f::Vector{Float64}}

    # validate that the input is a proper 3-D array (channels × samples × epochs)
    @assert size(s1) == size(s2) "s1 and s2 must have the same size."
    _chk3d(s1)
    _chk3d(s2)

    # number of channels
    ch_n = size(s1, 1)
    # number of epochs
    ep_n = size(s1, 3)

    # pre-compute the frequency vector with a single pilot call on the first channel/epoch pair
    cpsd_data = cpsd(
        @view(s1[1, :, 1]),
        @view(s2[1, :, 1]);
        method = method,
        fs = fs,
        flim = flim,
        demean = demean,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
    )
    f = cpsd_data.f
    
    # pre-allocate output
    pxy = zeros(ComplexF64, ch_n, length(f), ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        pxy[ch_idx, :, ep_idx], _ = cpsd(
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
    end

    return (pxy = pxy, f = f)
end

"""
    cpsd(obj1, obj2; <keyword arguments>)

Calculate the complex cross power spectral density (CPSD) between paired channels of two objects via one of three estimators:
- `:mt` – multi-taper (DSP.mt_cross_power_spectra)
- `:fft` – single-window FFT
- `:stft` – segmented FFT averaged over overlapping windows

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}, Regex}`: channel name(s)
- `ch2::Union{String, Vector{String}, Regex}`: channel name(s)
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: epoch number(s)
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: epoch number(s)
- `method::Symbol=:mt`: method used to calculate CPSD:
  - `:mt`: multi-tapered cross-power spectra
  - `:fft`: fast Fourier transformation
  - `:stft`: short-time Fourier transformation
- `flim::Tuple{Real, Real}=(0, sr(obj1) / 2)`: frequency bounds
- `demean::Bool=false`: if true, the channel-wise mean will be subtracted from the input signals before the cross spectral powers are computed
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj1)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.90)`: window overlap in samples
- `w::Bool=true`: if true, apply Hanning window

# Returns

Named tuple containing:

- `pxy::Array{ComplexF64, 3}`: cross-power spectrum, shape `(channels, frequencies, epochs)`
- `f::Vector{Float64}`: frequencies
"""
function cpsd(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}, Regex},
    ch2::Union{String, Vector{String}, Regex},
    ep1::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj1)),
    ep2::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj2)),
    method::Symbol = :mt,
    flim::Tuple{Real, Real} = (0, sr(obj1) / 2),
    demean::Bool = false,
    nt::Int64 = 7,
    wlen::Int64 = sr(obj1),
    woverlap::Int64 = round(Int64, wlen * 0.90),
    w::Bool = true,
)::@NamedTuple{pxy::Array{ComplexF64, 3}, f::Vector{Float64}}

    # validate objects
    @assert sr(obj1) == sr(obj2) "OBJ1 and OBJ2 must have the same sampling rate."

    # resolve channel names to integer indices, optionally skipping bad channels
    ch1 = exclude_bads ? get_channel(obj1, ch = ch1, exclude = "bad") : get_channel(obj1, ch = ch1, exclude = "")
    ch2 = exclude_bads ? get_channel(obj2, ch = ch2, exclude = "bad") : get_channel(obj2, ch = ch2, exclude = "")
    @assert length(ch1) == length(ch2) "Lengths of ch1 ($(length(ch1))) and ch2 ($(length(ch2))) must be equal."

    # validate epoch indices and ensure both objects have matching epoch structure
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    # normalize scalar epoch arguments to vectors so indexing is uniform
    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])
    @assert length(ep1) == length(ep2) "Lengths of ep1 ($(length(ep1))) and ep2 ($(length(ep2))) must be equal."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    cpsd_data = cpsd(
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

    return cpsd_data

end
