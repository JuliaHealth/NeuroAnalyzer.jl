export ftransform
export htransform
export transform
export hanalytic

"""
    ftransform(s; <keyword arguments>)

Calculate Fourier transform (FFT / rFFT).

# Arguments

- `s::AbstractVector`: signal vector
- `pad::Int64=0`: number of zeros to append
- `db::Bool=false`: normalize powers to dB
- `nf::Bool=false`: : if true, return coefficients for negative **and** positive frequencies; otherwise return coefficients for positive frequencies only

# Returns

Named tuple:

- `c::Vector{ComplexF64}`: Fourier coefficients
- `a::Vector{Float64}`: amplitudes
- `p::Vector{Float64}`: powers
- `ph::Vector{Float64}`: phases in radians

# Notes

To obtain the matching frequency vector use `f = freqs(s, fs).f`.
"""
function ftransform(
    s::AbstractVector;
    pad::Int64 = 0,
    db::Bool = false,
    nf::Bool = false,
)::@NamedTuple{
    c::Vector{ComplexF64},
    a::Vector{Float64},
    p::Vector{Float64},
    ph::Vector{Float64}
}

    # number of samples
    n = length(s)

    # compute FFT: full spectrum (positive + negative) or one-sided (positive only)
    ft = nf ? fft0(s, pad) : rfft0(s, pad)

    # amplitudes: normalize by N; compensate for the removed negative-frequency
    # mirror by doubling all bins except DC (index 1) and Nyquist (index end)
    a = abs.(ft) ./ n
    !nf && (a[2:(end - 1)] .*= 2)

    # powers: same normalization and symmetry compensation.
    p = abs2.(ft) ./ n
    !nf && (p[2:(end - 1)] .*= 2)

    # convert powers to dB if requested
    db && (p = pow2db.(p))

    # zero out near-zero coefficients before computing phase to avoid numerical noise polluting the phase angles
    ft[abs.(ft) .< eps()] .= 0
    ph = DSP.angle.(ft)

    return (c = ft, a = a, p = p, ph = ph)

end


"""
    ftransform(s; <keyword arguments>)

Calculate Fourier transform (FFT / rFFT)

# Arguments

- `s::AbstractArray`: signal array (channels, samples, epochs)
- `pad::Int64`: number of zeros to append
- `db::Bool=false`: normalize powers to dB
- `nf::Bool=false`: if true, return Fourier coefficients for negative and positive frequencies, otherwise return Fourier coefficients for positive frequencies only

# Returns

Named tuple:

- `c::Array{ComplexF64, 3}`: Fourier coefficients, shape `(channels, samples, epochs)`
- `a::Array{Float64, 3}`: amplitudes, shape `(channels, samples, epochs)`
- `p::Array{Float64, 3}`: powers, shape `(channels, samples, epochs)`
- `ph::Array{Float64, 3}`: phases in radians, shape `(channels, samples, epochs)`
"""
function ftransform(
    s::AbstractArray;
    pad::Int64 = 0,
    db::Bool = false,
    nf::Bool = false,
)::@NamedTuple{
    c::Array{ComplexF64, 3},
    a::Array{Float64, 3},
    p::Array{Float64, 3},
    ph::Array{Float64, 3}
}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pilot call to determine output length (depends on pad and nf)
    fft_tmp = ftransform(
        @view(s[1, :, 1]),
        pad = pad,
        db = db,
        nf = nf,
    )
    c = zeros(ComplexF64, ch_n, length(fft_tmp.c), ep_n)
    a = zeros(Float64, ch_n, length(fft_tmp.a), ep_n)
    p = zeros(Float64, ch_n, length(fft_tmp.p), ep_n)
    ph = zeros(Float64, ch_n, length(fft_tmp.ph), ep_n)

    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        ftransform_data = ftransform(
            @view(s[ch_idx, :, ep_idx]),
            pad = pad,
            db = db,
            nf = nf
        )
        c[ch_idx, :, ep_idx]  = ftransform_data.c
        a[ch_idx, :, ep_idx]  = ftransform_data.a
        p[ch_idx, :, ep_idx]  = ftransform_data.p
        ph[ch_idx, :, ep_idx] = ftransform_data.ph
    end

    return (c = c, a = a, p = p, ph = ph)

end

"""
    htransform(s; <keyword arguments>)

Calculate Hilbert transform (analytic signal).

# Arguments

- `s::AbstractVector`: signal vector
- `db::Bool=false`: normalize powers to dB

# Returns

Named tuple:

- `c::Vector{ComplexF64}`: analytic (Hilbert) coefficients
- `a::Vector{Float64}`: instantaneous amplitudes
- `p::Vector{Float64}`: instantaneous powers
- `ph::Vector{Float64}`: instantaneous phases in radians
"""
function htransform(
    s::AbstractVector;
    db::Bool = false,
)::@NamedTuple{
    c::Vector{ComplexF64},
    a::Vector{Float64},
    p::Vector{Float64},
    ph::Vector{Float64}
}

    # compute the analytic signal via the Hilbert transform
    ht = DSP.hilbert(s)

    # instantaneous amplitude (envelope)
    a = abs.(ht)

    # instantaneous phase
    ph = DSP.angle.(ht)

    # instantaneous power
    p = abs2.(ht)
    db && (p = pow2db.(p))

    return (c = ht, a = a, p = p, ph = ph)

end

"""
    htransform(s; <keyword arguments>)

Calculate Hilbert transform (analytic signal).

# Arguments

- `s::AbstractArray`: signal array (channels, samples, epochs)
- `db::Bool=false`: normalize powers to dB

# Returns

Named tuple:

- `c::Array{ComplexF64, 3}`: Hilbert coefficients
- `a::Array{Float64, 3}`: amplitudes
- `p::Array{Float64, 3}`: powers
- `ph::Array{Float64, 3}`: phases (in radians)
"""
function htransform(
    s::AbstractArray;
    db::Bool = false,
)::@NamedTuple{
    c::Array{ComplexF64, 3},
    a::Array{Float64, 3},
    p::Array{Float64, 3},
    ph::Array{Float64, 3}
}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # epoch length
    ep_len = size(s, 2)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate outputs
    c = zeros(ComplexF64, ch_n, ep_len, ep_n)
    a = similar(s)
    p = similar(s)
    ph = similar(s)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        htransform_data = htransform(
            @view(s[ch_idx, :, ep_idx]),
            db = db
        )
        c[ch_idx, :, ep_idx] = htransform_data.c
        a[ch_idx, :, ep_idx] = htransform_data.a
        p[ch_idx, :, ep_idx] = htransform_data.p
        ph[ch_idx, :, ep_idx] = htransform_data.ph
    end

    return (c = c, a = a, p = p, ph = ph)

end

"""
    transform(s; <keyword arguments>)

Calculate Fourier or Hilbert transformation.

# Arguments

- `s::AbstractArray`: signal array (channels, samples, epochs)
- `pad::Int64=0`: (FFT only) number of zeros to append
- `h::Bool=false`: if true, use Hilbert transform instead of FFT
- `db::Bool=false`: normalize powers to dB
- `nf::Bool=false`: (FFT only) if true, return negative and positive frequency coefficients; otherwise positive frequencies only

# Returns

Named tuple:

- `c::Array{ComplexF64, 3}`: Fourier or Hilbert coefficients
- `a::Array{Float64, 3}`: amplitudes
- `p::Array{Float64, 3}`: powers
- `ph::Array{Float64, 3}: phases in radians
"""
function transform(
    s::AbstractArray;
    pad::Int64 = 0,
    h::Bool = false,
    db::Bool = false,
    nf::Bool = false,
)::@NamedTuple{
    c::Array{ComplexF64, 3},
    a::Array{Float64, 3},
    p::Array{Float64, 3},
    ph::Array{Float64, 3}
}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    if h
        _warn("htransform() uses Hilbert transform, the signal should be narrowband for best results.")
        return htransform(s, db = db)
    else
        return ftransform(s, pad = pad, db = db, nf = nf)
    end

end

"""
    transform(obj; <keyword arguments>)

Calculate Fourier/Hilbert transform.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pad::Int64=0`: (FFT only) number of zeros to append
- `h::Bool=false`: if true, use Hilbert transform instead of FFT
- `db::Bool=false`: normalize powers to dB
- `nf::Bool=false`: (FFT only) if true, return negative and positive frequency coefficients; otherwise positive frequencies only
# Returns

Named tuple:

- `c::Array{ComplexF64, 3}`: Fourier or Hilbert coefficients
- `a::Array{Float64, 3}`: amplitudes
- `p::Array{Float64, 3}`: powers
- `ph::Array{Float64, 3}: phases in radians
"""
function transform(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    pad::Int64 = 0,
    h::Bool = false,
    db::Bool = false,
    nf::Bool = false,
)::@NamedTuple{
    c::Array{ComplexF64, 3},
    a::Array{Float64, 3},
    p::Array{Float64, 3},
    ph::Array{Float64, 3}
}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    return transform(
        obj.data[ch, :, :],
        pad = pad,
        h = h,
        db = db,
        nf = nf,
    )

end

"""
    hanalytic(s)

Calculate complex analytic signal (`s + i·H(s)`) using Hilbert transformation.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `ha::Vector{ComplexF64}`: complex analytic signal
"""
function hanalytic(s::AbstractVector)::Vector{ComplexF64}

    ha = DSP.hilbert(s)

    return ha

end

"""
    hanalytic(s; <keyword arguments>)

Calculate complex analytic signal (`s + i·H(s)`) using Hilbert transformation.

# Arguments

- `s::AbstractArray`: signal array (channels, samples, epochs)
- `pad::Int64`: number of zeros to append

# Returns

- `ha::Vector{ComplexF64}`: complex analytic signal, shape `(channels, samples, epochs)`
"""
function hanalytic(s::AbstractArray)::Array{ComplexF64, 3}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # epoch length
    ep_len = size(s, 2)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate outputs
    ha = zeros(ComplexF64, ch_n, ep_len, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        ha[ch_idx, :, ep_idx] = DSP.hilbert(@view(s[ch_idx, :, ep_idx]))
    end

    return ha

end

"""
    hanalytic(obj; <keyword arguments>)

Calculate complex analytic signal (`s + i·H(s)`) using Hilbert transformation.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pad::Int64=0`: number of zeros to append

# Returns

- `ha::Vector{ComplexF64}`: complex analytic signal, shape `(channels, samples, epochs)`
"""
function hanalytic(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex}
)::Array{ComplexF64, 3}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    ha = NeuroAnalyzer.hanalytic(@view(obj.data[ch, :, :]))

    return ha

end
