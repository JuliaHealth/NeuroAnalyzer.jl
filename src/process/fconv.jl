export fconv

"""
    fconv(s; kernel, norm, pad)

Perform convolution in the frequency domain.

# Arguments

- `s::AbstractVector`
- `kernel::AbstractVector`
- `norm::Bool=true`: normalize kernel to keep the post-convolution results in the same scale as the original data
- `pad::Int64=0`: number of zeros to add to the signal for FFT

# Returns

- `s_new::Vector{ComplexF64}`: convoluted signal
"""
function fconv(s::AbstractVector; kernel::AbstractVector, norm::Bool=true, pad::Int64=0)

    n_s = length(s)
    n_k = length(kernel)
    half_k = floor(Int64, n_k / 2)
    s_fft = fft0(s, pad + n_k - 1)
    kernel_fft = fft0(kernel, pad + n_s - 1)
    norm && (kernel_fft ./= cmax(kernel_fft))
    s_conv = ifft0(s_fft .* kernel_fft, pad)

    # remove in- and out- edges
    if mod(n_k, 2) == 0
        s_new = s_conv[half_k:(end - half_k)]
    else
        s_new = s_conv[half_k:(end - half_k - 1)]
    end

    return s_new

end

"""
    fconv(s; kernel, norm, pad)

Perform convolution in the frequency domain.

# Arguments

- `s::AbstractArray`
- `kernel::AbstractVector`: convolution kernel
- `norm::Bool=true`: normalize kernel to keep the post-convolution results in the same scale as the original data
- `pad::Int64=0`: number of zeros to add to the signal for FFT

# Returns

- `s_new::Array{ComplexF64, 3}`: convoluted signal
"""
function fconv(s::AbstractArray; kernel::AbstractVector, norm::Bool=true, pad::Int64=0)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = zeros(ComplexF64, size(s))

    # initialize progress bar
    progress_bar && (progbar = Progress(ep_n * ch_n, dt=1, barlen=20, color=:white))

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views fconv(s[ch_idx, :, ep_idx], kernel=kernel, norm=norm, pad=pad)

            # update progress bar
            progress_bar && next!(progbar)
        end
    end

    return s_new

end

"""
    fconv(obj; ch, kernel, norm, pad)

Perform convolution in the frequency domain.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `kernel::AbstractVector`: convolution kernel
- `norm::Bool=true`: normalize kernel to keep the post-convolution results in the same scale as the original data
- `pad::Int64=0`: number of zeros to add to the signal for FFT

# Returns

- `s_new::Array{ComplexF64, 3}`: convoluted signal
"""
function fconv(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), kernel::AbstractVector, norm::Bool=true, pad::Int64=0)

    _check_channels(obj, ch)

    s_new = @views fconv(obj.data[ch, :, :], kernel=kernel, norm=norm, pad=pad)
    
    return s_new

end