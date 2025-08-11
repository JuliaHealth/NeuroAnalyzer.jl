export fconv

"""
    fconv(s; <keyword arguments>)

Perform convolution in the frequency domain.

# Arguments

- `s::AbstractVector`
- `kernel::AbstractVector`
- `norm::Bool=true`: normalize kernel to keep the post-convolution results in the same scale as the original data

# Returns

- `s_new::Vector{ComplexF64}`: convoluted signal
"""
function fconv(s::AbstractVector; kernel::AbstractVector, norm::Bool=true)::Vector{ComplexF64}

    s_fft = fft0(s, length(kernel) - 1)
    kernel_fft = fft0(kernel, length(s) - 1)
    norm && (kernel_fft ./= cmax(kernel_fft))
    s_conv = ifft0(s_fft .* kernel_fft)
    s_new = _remove_kernel(s_conv, kernel)

    return s_new

end

"""
    fconv(s; <keyword arguments>)

Perform convolution in the frequency domain.

# Arguments

- `s::AbstractArray`
- `kernel::AbstractVector`: convolution kernel
- `norm::Bool=true`: normalize kernel to keep the post-convolution results in the same scale as the original data

# Returns

- `s_new::Array{ComplexF64, 3}`: convoluted signal
"""
function fconv(s::AbstractArray; kernel::AbstractVector, norm::Bool=true)::Array{ComplexF64, 3}

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = zeros(ComplexF64, size(s))

    # initialize progress bar
    progress_bar && (progbar = Progress(ep_n * ch_n, dt=1, barlen=20, color=:white))

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views fconv(s[ch_idx, :, ep_idx], kernel=kernel, norm=norm)

            # update progress bar
            progress_bar && next!(progbar)
        end
    end

    return s_new

end

"""
    fconv(obj; <keyword arguments>)

Perform convolution in the frequency domain.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `kernel::AbstractVector`: convolution kernel
- `norm::Bool=true`: normalize kernel to keep the post-convolution results in the same scale as the original data

# Returns

- `s_new::Array{ComplexF64, 3}`: convoluted signal
"""
function fconv(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, kernel::AbstractVector, norm::Bool=true)::Array{ComplexF64, 3}

    ch = get_channel(obj, ch=ch)
    s_new = @views fconv(obj.data[ch, :, :], kernel=kernel, norm=norm)
    _info("Group delay: $(_group_delay(kernel)) samples")

    return s_new

end