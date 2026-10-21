export filter_g
export filter_g!

"""
    _filter_g_analytic(s; <keyword arguments>)

Band-limited analytic signal: the FFT spectrum of the signal multiplied by a Gaussian kernel centered at frequency `f` (unit gain at its peak, positive frequencies only), scaled so that its absolute value is the amplitude envelope of the band-limited signal and its angle is the instantaneous phase.

# Arguments

- `s::AbstractVector`: signal vector
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `pad::Int64=0`: number of zeros to append before the FFT; must be ≥ 0
- `f::Real`: center frequency of the Gaussian kernel in Hz; must be ≥ 0 and ≤ `fs/2`
- `gw::Real=5`: Gaussian width in Hz (full width parameter); must be > 0

# Returns

- `Vector{ComplexF64}`: analytic band-limited signal, length `length(s)`
"""
function _filter_g_analytic(
    s::AbstractVector;
    fs::Int64,
    pad::Int64 = 0,
    f::Real,
    gw::Real = 5,
)::Vector{ComplexF64}
    # validate
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))
    pad >= 0 || throw(ArgumentError("pad must be ≥ 0."))
    f >= 0 || throw(ArgumentError("f must be ≥ 0."))
    f <= fs / 2 || throw(ArgumentError("f must be ≤ $(fs / 2) Hz (Nyquist)."))
    gw > 0 || throw(ArgumentError("gw must be > 0."))

    # frequency axis matching the FFT output length (including zero-padding)
    n = length(s) + pad
    gf = collect(0:(n - 1)) .* (fs / n)

    # Gaussian standard deviation derived from the width parameter
    gs = (gw * (2π - 1)) / (4π)

    # Build the Gaussian kernel centered at f, then normalize to unit peak gain
    gf .-= f
    g = @. exp(-0.5 * (gf / gs)^2)
    g ./= maximum(abs, g)

    # multiply spectrum by kernel and invert: one-sided spectrum, doubled
    return 2 .* ifft0(fft0(s, pad) .* g, pad)
end

"""
    _filter_g(s; <keyword arguments>)

Filter a signal by multiplying its FFT spectrum with a Gaussian kernel centered at frequency `f`.

The Gaussian is normalized to unit gain at its peak and covers the positive frequencies only, so the inverse FFT is the analytic signal of the band (see `_filter_g_analytic()`); its absolute value is the amplitude envelope of the filtered signal.

# Arguments

- `s::AbstractVector`: signal vector
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `pad::Int64=0`: number of zeros to append before the FFT; must be ≥ 0
- `f::Real`: center frequency of the Gaussian kernel in Hz; must be ≥ 0 and ≤ `fs/2`
- `gw::Real=5`: Gaussian width in Hz (full width parameter); must be > 0

# Returns

- `Vector{Float64}`: amplitude envelope of the filtered signal, length `length(s)`
"""
function _filter_g(
    s::AbstractVector;
    fs::Int64,
    pad::Int64 = 0,
    f::Real,
    gw::Real = 5,
)::Vector{Float64}
    return abs.(_filter_g_analytic(s; fs = fs, pad = pad, f = f, gw = gw))
end

"""
    _filter_g(s; <keyword arguments>)

Filter a 3-dimensional signal array using a Gaussian kernel in the frequency domain for a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `pad::Int64=0`: number of zeros to append; must be ≥ 0
- `f::Real`: center frequency of the Gaussian kernel in Hz; must be ≥ 0 and < `fs/2`
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`) (full width parameter); must be > 0

# Returns

- `Array{Float64, 3}`: filtered array of the same shape as `s`
"""
function _filter_g(
    s::AbstractArray;
    fs::Int64,
    pad::Int64 = 0,
    f::Real,
    gw::Real = 5,
)::Array{Float64, 3}
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    s_new = similar(s, Float64)

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        s_new[ch_idx, :, ep_idx] = _filter_g(
            @views(s[ch_idx, :, ep_idx]),
            fs = fs,
            pad = pad,
            f = f,
            gw = gw,
        )
    end

    return s_new
end

"""
    filter_g(obj; <keyword arguments>)

Replace selected channels of a NEURO object with the amplitude envelope of the band around `f`: the FFT of the signal is multiplied by a Gaussian kernel centered at `f` (positive frequencies only) and transformed back; the absolute value of this band-limited analytic signal is the envelope.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pad::Int64=0`: number of zeros to append; must be ≥ 0
- `f::Real`: center frequency of the Gaussian kernel in Hz; must be ≥ 0 and < `fs/2`
- `gw::Real=5`: Gaussian width in Hz (full width parameter); must be > 0

# Returns

- `NeuroAnalyzer.NEURO`: new object; the selected channels contain the amplitude envelope
"""
function filter_g(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    pad::Int64 = 0,
    f::Real,
    gw::Real = 5,
)::NeuroAnalyzer.NEURO
    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # create new dataset
    obj_tmp = deepcopy(obj)

    obj_tmp.data[ch, :, :] = _filter_g(
        @view(obj.data[ch, :, :]);
        fs = sr(obj),
        pad = pad,
        f = f,
        gw = gw,
    )
    push!(obj_tmp.history, "filter_g(obj; ch=$ch, pad=$pad, f=$f)")

    return obj_tmp
end

"""
    filter_g!(obj; <keyword arguments>)

Replace selected channels of a NEURO object in-place with the amplitude envelope of the band around `f` (see `filter_g()`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pad::Int64=0`: number of zeros to append; must be ≥ 0
- `f::Real`: center frequency of the Gaussian kernel in Hz; must be ≥ 0 and < `fs/2`
- `gw::Real=5`: Gaussian width in Hz (used by `:gh`) (full width parameter); must be > 0

# Returns

- `Nothing`
"""
function filter_g!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    pad::Int64 = 0,
    f::Real,
    gw::Real = 5,
)::Nothing
    obj_tmp = filter_g(obj; ch = ch, pad = pad, f = f, gw = gw)
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj_tmp = nothing

    return nothing
end
