export psd_rel

"""
    _psd_rel(s; <keyword arguments>)

Calculate relative power spectral density (PSD) for a 1-D signal vector.

Power at each frequency bin is expressed as a fraction of either the total broadband power (when `flim` is `nothing`) or the power within the specified frequency band.

Default method is Welch's periodogram.

# Arguments

- `s::AbstractVector`: signal vector
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `db::Bool=false`: if `true`, convert power to dB
- `flim::Union{Tuple{Real, Real}, Nothing}=nothing`: frequency band `(f_low, f_high)` used as the reference power; `nothing` uses total broadband power
- `method::Symbol=:welch`: PSD estimation method:
    - `:welch`: Welch's periodogram (default)
    - `:fft`: plain FFT periodogram
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform averaged over segments
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian filter + Hilbert transform
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `wlen::Int64=fs`: window length in samples (default = 1 second); used by `:welch` and `:stft`
- `woverlap::Int64=min(min(round(Int64, wlen * 0.90), wlen - 1), wlen - 1)`: window overlap in samples
- `w::Bool=true`: if `true`, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`)
- `demean::Bool=true`: subtract DC component before estimating PSD

# Returns

Named tuple:

- `p::Vector{Float64}`: relative powers (one per frequency bin)
- `f::Vector{Float64}`: corresponding frequencies in Hz
"""
function _psd_rel(
    s::AbstractVector;
    fs::Int64,
    db::Bool = false,
    flim::Union{Tuple{Real, Real}, Nothing} = nothing,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = min(round(Int64, wlen * 0.90), wlen - 1),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool = true,
)::@NamedTuple{
    p::Vector{Float64},
    f::Vector{Float64},
}

    # shared keyword arguments forwarded to every internal PSD/band-power call
    psd_kwargs = (
        fs = fs,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        gw = gw,
        demean = demean,
    )

    # compute the reference power: either within the specified band or broadband
    ref_pw = if flim === nothing
        _total_power(s; psd_kwargs...)
    else
        _band_power(s; flim = flim, psd_kwargs...)
    end

    # relative power on the linear scale, then optionally in dB
    psd_data = _psd(s; db = false, psd_kwargs...)
    p = psd_data.p ./ ref_pw
    db && (p = pow2db.(p))
    f = psd_data.f

    return (; p, f)
end

"""
    _psd_rel(s; <keyword arguments>)

Calculate relative power spectral density for a 2-D signal matrix.

Default method is Welch's periodogram.

# Arguments

- `s::AbstractMatrix`: signal matrix, shape (channels, samples)
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `db::Bool=false`: if `true`, convert power to dB
- `flim::Union{Tuple{Real, Real}, Nothing}=nothing`: frequency range to calculate relative power to; if nothing, than calculate relative to total power
- `method::Symbol=:welch`: PSD estimation method:
    - `:welch`: Welch's periodogram (default)
    - `:fft`: plain FFT periodogram
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform averaged over segments
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian filter + Hilbert transform
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `wlen::Int64=fs`: window length in samples (default = 1 second); used by `:welch` and `:stft`
- `woverlap::Int64=min(min(round(Int64, wlen * 0.90), wlen - 1), wlen - 1)`: window overlap in samples
- `w::Bool=true`: if `true`, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`)
- `demean::Bool=true`: subtract DC component before estimating PSD

# Returns

Named tuple:

- `p::Matrix{Float64}`: relative powers, shape (channels, frequency_bins)
- `f::Vector{Float64}`: corresponding frequencies in Hz
"""
function _psd_rel(
    s::AbstractMatrix;
    fs::Int64,
    db::Bool = false,
    flim::Union{Tuple{Real, Real}, Nothing} = nothing,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = min(round(Int64, wlen * 0.90), wlen - 1),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool = true,
)::@NamedTuple{
    p::Matrix{Float64},
    f::Vector{Float64},
}
    # number of channels
    ch_n = size(s, 1)

    psd_rel_1(x) = _psd_rel(
        x;
        fs = fs,
        db = db,
        flim = flim,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        gw = gw,
        demean = demean,
    )

    # first channel also determines the frequency vector
    r1 = psd_rel_1(@view(s[1, :]))
    f = r1.f
    # pre-allocate output
    p = zeros(ch_n, length(f))
    p[1, :] = r1.p

    @inbounds for ch_idx = 2:ch_n
        p[ch_idx, :] = psd_rel_1(@view(s[ch_idx, :])).p
    end

    return (; p, f)
end

"""
    _psd_rel(s; <keyword arguments>)

Calculate relative power spectral density for a 3-D signal array.

Default method is Welch's periodogram.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `db::Bool=false`: if `true`, convert power to dB
- `flim::Union{Tuple{Real, Real}, Nothing}=nothing`: frequency range to calculate relative power to; if nothing, than calculate relative to total power
- `method::Symbol=:welch`: PSD estimation method:
    - `:welch`: Welch's periodogram (default)
    - `:fft`: plain FFT periodogram
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform averaged over segments
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian filter + Hilbert transform
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `wlen::Int64=fs`: window length in samples (default = 1 second); used by `:welch` and `:stft`
- `woverlap::Int64=min(min(round(Int64, wlen * 0.90), wlen - 1), wlen - 1)`: window overlap in samples
- `w::Bool=true`: if `true`, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: Morlet wavelet cycles, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], nfrq)`, where `nfrq` is the length of `0:(fs / 2)`
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`)
- `demean::Bool=true`: subtract DC component before estimating PSD

# Returns

Named tuple:

- `p::Array{Float64, 3}`: relative powers, shape (channels, frequency_bins, epochs)
- `f::Vector{Float64}`: corresponding frequencies in Hz
"""
function _psd_rel(
    s::AbstractArray;
    fs::Int64,
    db::Bool = false,
    flim::Union{Tuple{Real, Real}, Nothing} = nothing,
    method::Symbol = :welch,
    nt::Int64 = 7,
    wlen::Int64 = fs,
    woverlap::Int64 = min(round(Int64, wlen * 0.90), wlen - 1),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool = true,
)::@NamedTuple{
    p::Array{Float64, 3},
    f::Vector{Float64},
}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    psd_rel_1(x) = _psd_rel(
        x;
        fs = fs,
        db = db,
        flim = flim,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        gw = gw,
        demean = demean,
    )

    # first channel/epoch also determines the frequency vector
    psd_data = psd_rel_1(@view(s[1, :, 1]))
    f = psd_data.f
    # pre-allocate output
    p = zeros(ch_n, length(f), ep_n)
    p[1, :, 1] = psd_data.p

    # calculate over channels and epochs (skipping the one already done)
    idxs = vec(CartesianIndices((ch_n, ep_n)))[2:end]
    @inbounds Threads.@threads :static for idx in idxs
        ch_idx, ep_idx = idx[1], idx[2]
        p[ch_idx, :, ep_idx] = psd_rel_1(@view(s[ch_idx, :, ep_idx])).p
    end

    return (; p, f)
end

"""
    psd_rel(obj; <keyword arguments>)

Calculate relative power spectral density for a NEURO object.

Default method is Welch's periodogram.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `db::Bool=false`: if `true`, convert power to dB
- `flim::Union{Tuple{Real, Real}, Nothing}=nothing`: frequency range to calculate relative power to; if nothing, than calculate relative to total power
- `method::Symbol=:welch`: PSD estimation method:
    - `:welch`: Welch's periodogram (default)
    - `:fft`: plain FFT periodogram
    - `:mt`: multi-tapered periodogram
    - `:stft`: short-time Fourier transform averaged over segments
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian filter + Hilbert transform
- `nt::Int64=7`: number of Slepian tapers (used by `:mt`)
- `wlen::Int64=sr(obj)`: window length in samples (default = 1 second); used by `:welch` and `:stft`
- `woverlap::Int64=min(min(round(Int64, wlen * 0.90), wlen - 1), wlen - 1)`: window overlap in samples
- `w::Bool=true`: if `true`, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of Morlet wavelet cycles; for a tuple the number of cycles is logarithmically (geometrically) spaced from `ncyc[1]` to `ncyc[2]` across frequencies
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`)
- `demean::Bool=true`: subtract DC component before estimating PSD

# Returns

Named tuple:

- `p::Array{Float64, 3}`: relative powers, shape (channels, frequency_bins, epochs)
- `f::Vector{Float64}`: corresponding frequencies in Hz
"""
function psd_rel(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    db::Bool = false,
    method::Symbol = :welch,
    nt::Int64 = 7,
    flim::Union{Tuple{Real, Real}, Nothing} = nothing,
    wlen::Int64 = sr(obj),
    woverlap::Int64 = min(round(Int64, wlen * 0.90), wlen - 1),
    w::Bool = true,
    ncyc::Union{Int64, Tuple{Int64, Int64}} = 32,
    gw::Real = 5,
    demean::Bool = true,
)::@NamedTuple{
    p::Array{Float64, 3},
    f::Vector{Float64},
}
    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    isempty(ch) && throw(ArgumentError("No channels selected."))

    return _psd_rel(
        @view(obj.data[ch, :, :]);
        fs = sr(obj),
        flim = flim,
        db = db,
        method = method,
        nt = nt,
        wlen = wlen,
        woverlap = woverlap,
        w = w,
        ncyc = ncyc,
        gw = gw,
        demean = demean,
    )

end
