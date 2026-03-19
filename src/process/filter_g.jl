export filter_g

"""
    filter_g(s; <keyword arguments>)

Filter a signal by multiplying its FFT spectrum with a Gaussian kernel centered at frequency `f`.

The Gaussian is normalized to unit gain at its peak. Taking the absolute value of the inverse FFT after multiplication yields the filtered signal envelope.

# Arguments

- `s::AbstractVector`: signal vector
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `pad::Int64=0`: number of zeros to append before the FFT; must be ≥ 0
- `f::Real`: center frequency of the Gaussian kernel in Hz; must be ≥ 0 and < `fs/2`
- `gw::Real=5`: Gaussian width in Hz (full width parameter); must be > 0

# Returns

- `Vector{Float64}`: filtered signal of length `length(s)`

# Throws

- `ArgumentError`: if `fs < 1`, `pad < 0`, `f < 0`, `f ≥ fs/2`, or `gw ≤ 0`

# See also

[`filter_g(::AbstractArray)`](@ref), [`filter_g(::NeuroAnalyzer.NEURO)`](@ref)
"""
function filter_g(
    s::AbstractVector;
    fs::Int64,
    pad::Int64 = 0,
    f::Real,
    gw::Real = 5
)::Vector{Float64}

    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))
    pad >= 0 || throw(ArgumentError("pad must be ≥ 0."))
    f >= 0 || throw(ArgumentError("f must be ≥ 0."))
    f < fs / 2 || throw(ArgumentError("f must be < $(fs/2) Hz (Nyquist)."))
    gw > 0 || throw(ArgumentError("gw must be > 0."))

    # frequency axis matching the FFT output length (including zero-padding)
    n   = length(s) + pad
    gf  = collect(range(0, fs; length=n))

    # Gaussian standard deviation derived from the width parameter
    gs  = (gw * (2π - 1)) / (4π)

    # Build the Gaussian kernel centered at f, then normalize to unit peak gain
    gf .-= f
    g = @. exp(-0.5 * (gf / gs)^2)
    g ./= maximum(abs, g)

    # multiply spectrum by kernel and invert; abs gives the signal envelope
    return abs.(ifft0(fft0(s, pad) .* g, pad))

end

"""
    filter_g(s; <keyword arguments>)

Filter a 3-dimensional signal array using a Gaussian kernel in the frequency domain. Applies [`filter_g(::AbstractVector)`](@ref) to each channel × epoch slice in parallel.

# Arguments

- `s::AbstractArray`: 3-D signal array (channels, samples, epochs)
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `pad::Int64=0`: number of zeros to append; must be ≥ 0
- `f::Real`: center frequency of the Gaussian kernel in Hz; must be ≥ 0 and < `fs/2`
- `gw::Real=5`: Gaussian width in Hz (full width parameter); must be > 0

# Returns

- `Array{Float64, 3}`: filtered array of the same shape as `s`

# Throws

- `ArgumentError`: propagated from [`filter_g(::AbstractVector)`](@ref)

# See also

[`filter_g(::AbstractVector)`](@ref), [`filter_g(::NeuroAnalyzer.NEURO)`](@ref)
"""
function filter_g(
    s::AbstractArray;
    fs::Int64,
    pad::Int64 = 0,
    f::Real,
    gw::Real = 5
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
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        s_new[ch_idx, :, ep_idx] = filter_g(@views(s[ch_idx, :, ep_idx]),
                                            fs = fs,
                                            pad = pad,
                                            f = f,
                                            gw = gw)
    end

    return s_new

end

"""
    filter_g(obj; <keyword arguments>)

Filter selected channels of a NEURO object using a Gaussian kernel in the frequency domain.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pad::Int64=0`: number of zeros to append; must be ≥ 0
- `f::Real`: center frequency of the Gaussian kernel in Hz; must be ≥ 0 and < `fs/2`
- `gw::Real=5`: Gaussian width in Hz (full width parameter); must be > 0

# Returns

- `NeuroAnalyzer.NEURO`: new object with filtered channels

# See also

[`filter_g!`](@ref), [`filter_g(::AbstractArray)`](@ref)
"""
function filter_g(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    pad::Int64 = 0,
    f::Real,
    gw::Real = 5
)::NeuroAnalyzer.NEURO

    ch = get_channel(obj, ch = ch)
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = @views filter_g(obj.data[ch, :, :], fs = sr(obj), pad = pad, f = f, gw = gw)
    push!(obj_new.history, "filter_g(OBJ, ch=$ch, pad=$pad, f=$f)")

    return obj_new

end

"""
    filter_g!(obj; <keyword arguments>)

Filter selected channels of a NEURO object in-place using a Gaussian kernel in the frequency domain.

Mutates `obj.data` and `obj.history`. Delegates to [`filter_g`](@ref) and copies the result back.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pad::Int64=0`: number of zeros to append; must be ≥ 0
- `f::Real`: center frequency of the Gaussian kernel in Hz; must be ≥ 0 and < `fs/2`
- `gw::Real=5`: Gaussian width in Hz (full width parameter); must be > 0

# Returns

- `Nothing`

# See also

[`filter_g`](@ref)
"""
function filter_g!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    pad::Int64 = 0,
    f::Real,
    gw::Real = 5
)::Nothing

    obj_new = filter_g(obj, ch = ch, pad = pad, f = f, gw = gw)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing

end
