module PhaseSlopeIndex

using Einsum
using StatsKit
using FFTW

# Exports
#---
export data2psi

"""
    int(x) = trunc(Int, x)
"""
int(x) = trunc(Int, x)

"""
    dropmean(X, d) = dropdims(mean(X, dims=d), dims=d)
"""
function dropmean(X, d)
    if ndims(X) == 1
        mean(X; dims=d)
    else
        dropdims(mean(X; dims=d); dims=d)
    end
end

"""
    squeeze(X::AbstractArray)

removing singleton dimensions
"""
function squeeze(X::AbstractArray)
    keepd = Tuple(d for d in size(X) if d != 1)
    return reshape(X, keepd)
end

"""
    detrend!(data, n)

(in place) Linear detrend of signals along first axis removing the n-th order polynomial.
This detrend function is limited to linear orders (0th and 1st order).

### Arguments

  - `data::AbstractArray`: N-dim array where signal is in column-major order
  - `n::Integer`: `n = 0` subtracts mean from data, `n = 1` removes linear trend

**Note**: shape of data must be (signal length, ...)
"""
function detrend!(data::AbstractArray, n::Integer)
    original_shape = size(data)
    nsamp = size(data, 1)  # number of samples

    A = Array{Float64}([ones(nsamp) Array(1:nsamp)])

    data = reshape(data, (nsamp, :))  # reshaping data
    if n == 0
        data .-= mean(data; dims=1)
    elseif n == 1
        data .-= A * (A \ data)
    end
    return reshape(data, original_shape)
end

"""
    window = hanning_fun(N)

Hanning window similar to MATLAB `hanning` implementation
"""
function hanning_fun(N::Integer)
    x = [range(0.0, 1.0; length=N + 2);]
    window = 0.5 .* (1 .- cospi.(2 .* x))
    window = (window + window[end:-1:1]) ./ 2  # forcing symmetry
    return window[2:(end - 1)]  # excluding the zero values
end

"""
    data2para(data, seglen, segshift, eplen, freqlist, method)

Extracts and builds a named tuple of parameters.

### Arguments

  - `data::AbstractArray`: NxM array for N data points in M channels.
  - `seglen::Integer`: segment length (determines the frequency resolution).
  - `segshift::Integer`: number of bins by which neighboring segments are shifted
    (e.g. `segshift = seglen / 2` makes overlapping segments).
  - `eplen::Integer`: length of epochs
  - `freqlist::AbstractArray`: a UnitRange or 2D-Array where each column is a frequency band
  - `method::String`: standard deviation estimation method
  - `subave::Bool`: if `true`, subtract average from CS segments
    (for continuous data, `subave = false`).
  - `verbose::Bool`: if `true`, warnings and info logs would be echoed.

### Returns

  - `parameters::NamedTuple`: a named tuple of parameters
"""
function data2para(
    data::AbstractArray,
    seglen::Integer,
    segshift::Integer,
    eplen::Integer,
    freqlist::AbstractArray{Int},
    method::String,
    subave::Bool,
    verbose::Bool
    )
    # data dimension
    if ndims(data) != 2
        data = squeeze(data)
        ndims(data) != 2 && throw(DimensionMismatch("Data must be a 2D-array!"))
        verbose && @info "Data is squeezed to a 2D-array)"
    end
    if size(data, 1) < size(data, 2)
        verbose && @info "Data is transposed to (#samples, #channels)"
        data = reshape(data, size(data, 2), size(data, 1))
    end
    if size(data, 1) < seglen
        throw(DimensionMismatch("seglen must be smaller than number of samples!"))
    end

    # number of samples per channel and number of channels
    nsamples, nchannels = size(data)

    # method shall always be lowercase
    method = lowercase(method)
    if eplen == 0
        verbose && @warn "Epoch length = 0 => No estimation of standard deviation."
        method = "none"
    end

    # if eplen = nsamples: continuous recording
    eplen == 0 && (eplen = nsamples)
    nep = int(nsamples / eplen)  # number of epochs

    if nep == 1 && subave == true
        subave = false
        verbose && @warn "subave is set to false for continuous data (nep = 1)"
    end

    segshift == 0 && (segshift = int(seglen / 2))
    nseg = int((eplen - seglen) / segshift) + 1

    # size(freqlist) = (freqs, nfbands)
    if length(freqlist) == 0
        freqlist = reshape(1:((int(seglen / 2)) - 1), :, 1)
    elseif ndims(freqlist) == 1
        freqlist = reshape(freqlist, :, 1)
    elseif ndims(freqlist) > 2
        throw("freqlist must be a UnitRange or a 2D-Array!")
    elseif size(freqlist, 1) < size(freqlist, 2)
        freqlist = freqlist'
        verbose && @info "freqlist is transposed to (#freq, #nfbands)"
    end
    if maximum(freqlist) >= int(seglen / 2)
        throw("Maximum frequency for freqlist is larger than the Nyquist frequency!")
    elseif minimum(freqlist) < 1
        throw("Minimum frequency for freqlist is 0 or negative!")
    end
    maxfreq = maximum(freqlist)  # max frequency of all frequency bands
    nfbands = size(freqlist, 2)  # number of frequency bands

    parameters = (
        data=data,
        nsamples=nsamples,
        nchan=nchannels,
        eplen=eplen,
        nep=nep,
        method=method,
        subave=subave,
        segshift=segshift,
        nseg=nseg,
        freqlist=freqlist,
        maxfreq=maxfreq,
        nfbands=nfbands,
    )

    return parameters
end

"""
    make_eposeg(data, seglen, nep, nseg, nchan, segshift)

Partitioning data into epochs and segments

### Arguments

  - `data::AbstractArray`: NxM array for N data points in M channels.
  - `seglen::Integer`: segment length
  - `nep::Integer`: number of epochs
  - `nseg::Integer`: number of segments per epoch
  - `nchan::Integer`: number of channels
  - `segshift::Integer`: number of bins by which neighboring segments are shifted.

### Returns

  - `epseg::AbstractArray`: partitioned data into shape `(seglen, nep, nseg, nchan)`

**Note**: returned Array may have more data entries than input data.
"""
function make_eposeg(
    data::AbstractArray,
    seglen::Integer,
    eplen::Integer,
    nep::Integer,
    nseg::Integer,
    nchan::Integer,
    segshift::Integer,
    )::AbstractArray

    # preallocation
    epseg = Array{Float64}(undef, seglen, nep, nseg, nchan)

    for (i, e) in zip(1:nep, 1:eplen:(nep * eplen))
        for (j, s) in zip(1:nseg, 1:segshift:(nseg * seglen))
            @views epseg[:, i, j, :] = data[e:(e + eplen - 1), :][s:(s + seglen - 1), :]
        end
    end

    return epseg
end

"""
    cs2ps(cs)

Cross Spectra to Phase Slope as defined in Eq. 3 [Nolte et al. 2008]:

```math
\\tilde\\Psi_{ij}=\\mathfrak{I}\\left(\\sum_{f \\in F} C_{ij}^*(f)~C_{ij}(f + \\delta f)\\right)
```

where:

```math
C_{ij}(f) = \\frac{S_{ij}(f)}{\\sqrt{S_{ii}(f)~S_{jj}(f)}}
```

is the complext coherency, and \$S\$ is the cross spectral matrix
(returnd by `data2cs` function), \$\\delta f\$ is the frequency resolution,
and \$\\mathfrak{I}\$ denotes taking the imaginary part.

### Arguments

  - `cs::AbstractArray`: Cross Spectral array with shape `(seglen, nchan, nchan)`

### Returns

  - phase slope index (::AbstractArray) with shape `(nchan, nchan)`

**Note**: frequency resolution is assumed to be the resolution of freq. band!
"""
function cs2ps(cs::AbstractArray)

    # size(cs) = (seglen, nchan, nchan)
    # we don't use df but we slice the fband

    # complex coherency
    @einsum coh[f, i, j] := cs[f, i, j] / sqrt(cs[f, i, i] * cs[f, j, j])

    # phase slope (Eq. 3)
    @views imag.(sum(conj(coh[1:(end - 1), :, :]) .* coh[2:end, :, :]; dims=1))
end

"""
    data2ps(data)

Segmented data to Cross Spectra based on Eq. 2 [Nolte et al. 2008]:

```math
S_{ij}(f) = \\langle \\hat{y}_i(f) \\hat{y}_i^*(f)\\rangle .
```

### Arguments

  - `data::AbstractArray`: Segmented data of shape `(maxfreq, nep, nseg, nchan)`

### Return

  - `cs::AbstractArray{Complex}`: Cross Spectral of shape `(maxfreq, nep, nseg, nchan, nchan)`
"""
function data2cs(data::AbstractArray)
    # cs: cross-spectral matrix

    # Eq. 2 size(csepseg) = (maxfreq, nep, nseg, nchan, nchan)
    @einsum cs[f, e, s, i, j] := data[f, e, s, i] * conj(data[f, e, s, j])
end

"""
    _cs_ = cs2cs_(data, cs, fband, nep, segave, subave, method)

preparing Cross Spectra for Phase Slope by segment averaging and subtraction

### Arguments

  - `data::AbstractArray`: Fourier-transformed detrended epoched segmented data.
  - `cs::AbstractArray`: Cross Spectra of data
  - `fband::AbstractArray`: 1D-Array of frequency range for PSI calculation
  - `nep::Integer`: number of epochs
  - `segave::Bool`: if `true`, averages across CS segments
  - `subave::Bool`: if `true`, subtract average across CS segments
  - `method::String`: standard deviation estimation method

### Returns

  - `_cs_::AbstractArray`: segment averaged and subtracted Cross Spectra
"""
function cs2cs_(
    data::AbstractArray,
    cs::AbstractArray,
    fband::AbstractArray,
    nep::Integer,
    segave::Bool,
    subave::Bool,
    method::String,
    )
    if segave
        if method == "bootstrap"
            randboot = rand(1:nep, nep)
            cs_ = dropmean(view(cs, :, randboot, :, :, :), (2, 3))
            av_ = dropmean(view(data, fband, randboot, :, :), (2, 3))
        elseif method == "psi"
            cs_ = dropmean(cs, (2, 3))
            av_ = dropmean(view(data, fband, :, :, :), (2, 3))
        elseif method == "jackknife"
            cs_ = dropmean(cs, 3)
            av_ = dropmean(view(data, fband, :, :, :), 3)
        end
    else
        if method == "bootstrap"
            randboot = rand(1:nep, nep)
            cs_ = dropmean(view(cs, :, randboot, 1, :, :), 2)
            av_ = dropmean(view(data, fband, randboot, 1, :), 2)
        elseif method == "psi"
            cs_ = dropmean(view(cs, :, :, 1, :, :), 2)
            av_ = dropmean(view(data, fband, :, 1, :), 2)
        elseif method == "jackknife"
            cs_ = view(cs, :, :, 1, :, :)
            av_ = view(data, fband, :, 1, :)
        end
    end

    if subave
        if method == "bootstrap" || method == "psi"
            @einsum av_ij[f, i, j] := av_[f, i] * conj(av_[f, j])
        elseif method == "jackknife"
            @einsum av_ij[f, e, i, j] := av_[f, e, i] * conj(av_[f, e, j])
        end
        _cs_ = squeeze(cs_ - av_ij)
    else
        if method == "bootstrap" || method == "psi" || method == "jackknife"
            _cs_ = squeeze(cs_)
        end
    end

    return _cs_
end

"""
    data2psi(data, seglen ; segshift, eplen, freqlist, method,
                            nboot, segave, subave, detrend)

calculates phase slope index (PSI)

### Arguments

  - `data::AbstractArray`: NxM array for N data points in M channels
  - `seglen::Integer`: segment length (determines the frequency resolution).
    If defining `freqlist`, `seglen` must be the same as sampling frequency.

*optional arguments*

  - `segshift::Integer`: number of bins by which neighboring segments are shifted
    (default is `seglen / 2`)
  - `eplen::Integer`: length of epochs (default is number of samples)
  - `freqlist::AbstractArray`: a UnitRange or 2D-Array where each column is a frequency band
    (default is full range). Note that the DC component (0th frequency of FFT)
    is discarded, and `freqlist` shall only include integer values.
  - `method::String`: standard deviation estimation method (default is `"jackknife"`)
  - `subave::Bool`: if `true`, subtract average across Cross Spectra segments (default is `false`)
  - `segave::Bool`: if `true`, average across Cross Spectra segments (default is `true`)
  - `nboot::Integer`: number of bootstrap resampling iterations (default is `100`)
  - `detrend::Bool`: if `true`, performs a 0th-order detrend across raw segments (default is `false`)
  - `window::Function`: window function with interval length as sole necessary argument (default is `Hanning`)
  - `verbose::Bool`: if `true`, warnings and info logs would be echoed. (default is `false`)

### Returns

  - `psi::AbstractArray`: Phase Slope Index with shape `(channel, channel, frequency bands)`
  - `psi_std::AbstractArray`: PSI estimated standard deviation with shape `(channel, channel, frequency bands)`
"""
function data2psi(
    data::AbstractArray,
    seglen::Integer;
    segshift::Integer=0,
    eplen::Integer=0,
    freqlist::AbstractArray{Int}=Int[],
    method::String="jackknife",
    subave::Bool=false,
    segave::Bool=true,
    nboot::Integer=100,
    detrend::Bool=false,
    window::Function=hanning_fun,
    verbose::Bool=false
    )

    (data, nsamples, nchan, eplen, nep, method, subave, segshift, nseg, freqlist, maxfreq, nfbands) = data2para(
        data, seglen, segshift, eplen, freqlist, method, subave, verbose
    )

    eposeg = make_eposeg(data, seglen, eplen, nep, nseg, nchan, segshift)

    if detrend
        detrend!(eposeg, 0)
    end

    eposeg .*= window(seglen)

    eposeg = view(fft(eposeg, 1), 2:(maxfreq + 1), :, :, :)

    # preallocation
    psi = Array{Float64}(undef, nchan, nchan, nfbands)
    if method == "jackknife"
        psi_est = Array{Float64}(undef, nchan, nchan, nfbands, nep)
    elseif method == "bootstrap"
        psi_est = Array{Float64}(undef, nchan, nchan, nfbands, nboot)
    end
    for (f, fband) in enumerate(eachrow(freqlist'))
        cs_full = data2cs(view(eposeg, fband, :, :, :))

        cs_psi = cs2cs_(eposeg, cs_full, fband, nep, segave, subave, "psi")
        psi[:, :, f] = cs2ps(cs_psi)

        if method == "jackknife"
            cs_jack = cs2cs_(eposeg, cs_full, fband, nep, segave, subave, "jackknife")
            for e in 1:nep
                cs_jack_se = (nep * cs_psi - view(cs_jack, :, e, :, :)) / (nep + 1)
                psi_est[:, :, f, e] = cs2ps(cs_jack_se)
            end
        elseif method == "bootstrap"
            for e in 1:nboot
                cs_boot = cs2cs_(eposeg, cs_full, fband, nep, segave, subave, "bootstrap")
                psi_est[:, :, f, e] = cs2ps(cs_boot)
            end
        end
    end

    if method == "jackknife"
        psi_std = sqrt(nep) * squeeze(std(psi_est; corrected=true, dims=4))
    elseif method == "bootstrap"
        psi_std = squeeze(std(psi_est; corrected=true, dims=4))
    else
        psi_std = fill(NaN, (nchan, nchan, nfbands))
    end

    psi = squeeze(psi)
    psi_std = squeeze(psi_std)

    return psi, psi_std
end

end  # module
