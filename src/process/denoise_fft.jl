export denoise_fft
export denoise_fft!

"""
    denoise_fft(s; <keyword arguments>)

Perform FFT-based denoising on a signal vector.

Computes the FFT, zeros all frequency components whose power exceeds the threshold `t`, then reconstructs the signal via the inverse FFT. This keeps low-power (signal) components and suppresses high-power artifact components.

# Arguments

- `s::AbstractVector`: signal vector
- `pad::Int64=0`: number of zeros to append before FFT; must be ≥ 0
- `t::Real=0`: power spectral density threshold; components with power > `t` are zeroed; if `t = 0`, the mean power across all frequency bins is used

# Returns

Named tuple:

- `s::Vector{Float64}`: denoised signal of the same length as the input
- `f_idx::BitVector`: boolean mask; `true` at each frequency index that was zeroed

# Throws

- `ArgumentError`: if `pad < 0`

# See also

[`denoise_fft(::AbstractArray)`](@ref), [`denoise_fft(::NeuroAnalyzer.NEURO)`](@ref)
"""
function denoise_fft(
    s::AbstractVector;
    pad::Int64 = 0,
    t::Real = 0
)::@NamedTuple{
    s::Vector{Float64},
    f_idx::BitVector
}

    pad >= 0 || throw(ArgumentError("pad must be ≥ 0."))

    # compute FFT and power spectrum
    s_fft = fft0(s, pad)
    s_pow = @. abs2(s_fft) / length(s)

    # set threshold to mean power if t=0
    t == 0 && (t = mean(s_pow))

    # zero coefficients above threshold
    f_idx = s_pow .> t
    s_fft[f_idx] .= 0

    return (s = abs.(ifft0(s_fft, pad)), f_idx = f_idx)

end

"""
    denoise_fft(s; <keyword arguments>)

Perform FFT-based denoising on every channel × epoch slice of a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape `(channels, samples, epochs)`
- `pad::Int64=0`: number of zeros to append before FFT; must be ≥ 0
- `t::Real=0`: power spectral density threshold; components with power > `t` are zeroed; if `t = 0`, the mean power across all frequency bins is used

# Returns

- `Array{Float64, 3}`: denoised 3D signal array

# Throws

- `ArgumentError`: if `s` is not a 3D array
"""
function denoise_fft(
    s::AbstractArray;
    pad::Int64 = 0,
    t::Real = 0
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
        dfft_data = denoise_fft(@view(s[ch_idx, :, ep_idx]), pad = pad, t = t)
        s_new[ch_idx, :, ep_idx] = dfft_data.s
    end

    return s_new

end

"""
    denoise_fft(obj; <keyword arguments>)

Perform FFT-based denoising on selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pad::Int64=0`: number of zeros to append before FFT; must be ≥ 0
- `t::Real=0`: power spectral density threshold; components with power > `t` are zeroed; if `t = 0`, the mean power across all frequency bins is used

# Returns

- `NeuroAnalyzer.NEURO`: new object with denoised channels

# See also

[`denoise_fft!`](@ref), [`denoise_fft(::AbstractArray)`](@ref)
"""
function denoise_fft(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    pad::Int64 = 0,
    t::Int64 = 0
)::NeuroAnalyzer.NEURO

    # resolve channel names to integer indices
    ch = get_channel(obj, ch = ch)

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = denoise_fft(@view(obj.data[ch, :, :]), pad = pad, t = t)
    push!(obj_new.history, "denoise_fft(OBJ, ch=$ch, pad=$pad, t=$t)")

    return obj_new

end

"""
    denoise_fft!(obj; <keyword arguments>)

Perform FFT-based denoising in-place on selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pad::Int64=0`: number of zeros to append before FFT; must be ≥ 0
- `t::Real=0`: power spectral density threshold; components with power > `t` are zeroed; if `t = 0`, the mean power across all frequency bins is used

# Returns

- `Nothing`

# See also

[`denoise_fft`](@ref)
"""
function denoise_fft!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    pad::Int64 = 0,
    t::Int64 = 0
)::Nothing

    obj_new = denoise_fft(obj, ch = ch, pad = pad, t = t)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing

end
