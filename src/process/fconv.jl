export fconv

"""
    fconv(signal; kernel, norm)

Perform convolution in the frequency domain.

# Arguments

- `signal::AbstractArray`
- `kernel::AbstractVector`
- `pad::Int64=0`: number of zeros to add
- `norm::Bool=true`: normalize kernel

# Returns

- `s_conv::Vector{Float64}`
"""
function fconv(signal::AbstractArray; pad::Int64=0, kernel::AbstractVector, norm::Bool=true)

    n_signal = length(signal)
    n_kernel = length(kernel)
    half_kernel = floor(Int64, n_kernel / 2)
    s_fft = fft0(signal, pad + n_kernel - 1)
    kernel_fft = fft0(kernel, pad + n_signal - 1)
    norm == true && (kernel_fft ./= cmax(kernel_fft))
    s_conv = @views ifft0(s_fft .* kernel_fft, pad)

    # remove in- and out- edges
    if mod(n_kernel, 2) == 0 
        return s_conv[half_kernel:(end - half_kernel)]
    else
        return s_conv[(half_kernel + 1):(end - half_kernel)]
    end
end

"""
    fconv(obj; channel, kernel, norm)

Perform convolution in the frequency domain.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all channels
- `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`: kernel for convolution
- `norm::Bool=true`: normalize kernel to keep the post-convolution results in the same scale as the original data
- `pad::Int64=0`: number of zeros to add signal for FFT

# Returns

- `s_convoluted::Array{Float64, 3}`: convoluted signal
"""
function fconv(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), kernel::Union{Vector{<:Real}, Vector{ComplexF64}}, norm::Bool=true, pad::Int64=0)

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    s_convoluted = zeros(ch_n, epoch_len(obj), ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ch_n, 1))

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_convoluted[ch_idx, :, ep_idx] = @views s_fconv(obj.data[channel[ch_idx], :, ep_idx], kernel=kernel, norm=norm, pad=pad)
            
            # update progress bar
            progress_bar == true && next!(p)
        end
    end

    return s_convoluted
end
