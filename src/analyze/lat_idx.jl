export lat_idx

"""
    lat_idx(obj; <keyword arguments>)

Calculate lateralization index (log(A / B), where A is average power at given frequency (default is 10 Hz, α) for the right hemisphere and B is average power at that frequency for the left hemisphere).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `frq::Union{Real, Tuple{<:Real, <:Real}}=10`: frequency at which the index is calculated; if range is provided, than averaged index across the range is calculated
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet convolution
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`
- `gw::Real=5`: Gaussian width in Hz
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

- `lidx::Float64`: lateralization index
"""
function lat_idx(obj::NeuroAnalyzer.NEURO; frq::Union{Real, Tuple{<:Real, <:Real}}=10, method::Symbol=:welch, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, gw::Real=5, wt::T=wavelet(Morlet(2π), β=32, Q=128))::Float64 where {T <: CWT}

    _check_datatype(obj, ["meg", "eeg", "erp", "erf"])

    @assert length(channel_pick(obj, p=:l)) > 0 "Could not detect left hemisphere channels, check OBJ labels."
    @assert length(channel_pick(obj, p=:r)) > 0 "Could not detect right hemisphere channels, check OBJ labels."

    ch_l = channel_pick(obj, p=:l)
    ch_r = channel_pick(obj, p=:r)

    _log_off()

    # left PSDs
    if datatype(obj) in ["erp", "erf"]
        p_left, f = psd(obj.data[ch_l, :, 1], fs=sr(obj), db=false, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw, wt=wt)
    else
        p_left, f = psd(obj.data[ch_l, :, :], fs=sr(obj), db=false, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw, wt=wt)
    end

    # right PSDs
    if datatype(obj) in ["erp", "erf"]
        p_right, _ = psd(obj.data[ch_r, :, 1], fs=sr(obj), db=false, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw, wt=wt)
    else
        p_right, _ = psd(obj.data[ch_r, :, :], fs=sr(obj), db=false, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw, wt=wt)
    end

    _log_on()

    # average across epochs
    size(p_left, 3) > 1 && (p_left = mean(p_left, dims=3))
    size(p_right, 3) > 1 && (p_right = mean(p_right, dims=3))
    # average across channels
    p_left = mean(p_left, dims=1)
    p_right = mean(p_right, dims=1)

    if length(frq) == 1
        frq_idx = vsearch(frq, f)
        p_left = p_left[1, frq_idx]
        p_right = p_right[1, frq_idx]
    else
        _check_tuple(frq, "frq")
        @assert frq[1] >= 0 "Lower frequency bound must be ≥ 0 Hz."
        @assert frq[2] <= sr(obj) / 2 "Upper frequency bound must be ≤ $(sr(obj) / 2) Hz."
        frq_idx1 = vsearch(frq[1], f)
        frq_idx2 = vsearch(frq[2], f)
        p_left = mean(p_left[1, frq_idx1:frq_idx2])
        p_right = mean(p_right[1, frq_idx1:frq_idx2])
    end

    lidx = log(p_right / p_left)

    return lidx

end
