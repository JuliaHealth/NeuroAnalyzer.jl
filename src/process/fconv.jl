export fconv

"""
    fconv(s; <keyword arguments>)

Perform convolution in the frequency domain.

Both `s` and `kernel` are zero-padded to length `length(s) + length(kernel) - 1` before the FFT so the result is equivalent to linear (non-circular) convolution.

# Arguments

- `s::AbstractVector`: signal vector
- `kernel::AbstractVector`: convolution kernel (must be non-empty)
- `norm::Bool=true`: normalize the kernel FFT by its maximum magnitude to keep post-convolution amplitudes on the same scale as the input

# Returns

- `Vector{ComplexF64}`: convolved signal (same length as `s`)

# Throws

- `ArgumentError` if `kernel` is empty
"""
function fconv(
    s::AbstractVector;
    kernel::AbstractVector,
    norm::Bool = true
)::Vector{ComplexF64}

    isempty(kernel) &&
        throw(ArgumentError("kernel must be non-empty."))

    s_fft = fft0(s, length(kernel) - 1)
    kernel_fft = fft0(kernel, length(s) - 1)

    if norm
        km = cmax(kernel_fft)
        # guard against a zero kernel (all elements zero → cmax = 0)
        iszero(km) && throw(ArgumentError(
            "kernel is all-zero; convolution would produce NaN output."))
        kernel_fft ./= km
    end

    s_conv = ifft0(s_fft .* kernel_fft)
    s_new = _remove_kernel(s_conv, kernel)

    return s_new

end

"""
    fconv(s; <keyword arguments>)

Perform convolution in the frequency domain.

Both `s` and `kernel` are zero-padded to length `length(s) + length(kernel) - 1` before the FFT so the result is equivalent to linear (non-circular) convolution.

# Arguments
- `s::AbstractArray`: signal array, shape `(channels, samples, epochs)`
- `kernel::AbstractVector`: convolution kernel (must be non-empty)
- `norm::Bool=true`: normalize the kernel FFT by its maximum magnitude to keep post-convolution amplitudes on the same scale as the input

# Returns
- `Array{ComplexF64, 3}`: convolved signal, same shape as `s`

# Throws

- `ArgumentError` if `s` is not 3-dimensional or `kernel` is empty
"""
function fconv(
    s::AbstractArray;
    kernel::AbstractVector,
    norm::Bool = true
)::Array{ComplexF64, 3}

    isempty(kernel) &&
        throw(ArgumentError("kernel must be non-empty."))

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    s_new = zeros(ComplexF64, size(s))

    # initialize progress bar
    progbar = Progress(ep_n * ch_n, dt = 1, barlen = 20, color = :white, enabled = progress_bar)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        s_new[ch_idx, :, ep_idx] = fconv(@view(s[ch_idx, :, ep_idx]), kernel = kernel, norm = norm)

        # update progress bar
        progress_bar && next!(progbar)
    end

    return s_new

end

"""
    fconv(obj; <keyword arguments>)

Perform convolution in the frequency domain on selected channels of a `NeuroAnalyzer.NEURO` object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `kernel::AbstractVector`: convolution kernel (must be non-empty)
- `norm::Bool=true`: normalize the kernel FFT by its maximum magnitude to keep post-convolution amplitudes on the same scale as the input

# Returns

- `Array{ComplexF64, 3}`: convolved signal for the selected channels
"""
function fconv(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    kernel::AbstractVector,
    norm::Bool = true
)::Array{ComplexF64, 3}

    # resolve channel names to integer indices
    ch = get_channel(obj, ch = ch)
    _info("Group delay: $(_group_delay(kernel)) samples")
    
    return fconv(@view(obj.data[ch, :, :]), kernel = kernel, norm = norm)

end
