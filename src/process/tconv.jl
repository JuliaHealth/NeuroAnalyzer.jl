export tconv

"""
    tconv(signal; kernel)

Performs convolution in the time domain.

# Arguments

- `signal::AbstractVector`
- `kernel::AbstractVector`

# Returns

- `s_conv::Vector{Float64}`
"""
function tconv(signal::AbstractVector; kernel::AbstractVector)

    s_conv = DSP.conv(signal, kernel)

    half_kernel = floor(Int, length(kernel) / 2)

    # remove in- and out- edges
    if mod(length(kernel), 2) == 0 
        return real.(s_conv)[half_kernel:(end - half_kernel)]
    else
        return real.(s_conv)[(half_kernel + 1):(end - half_kernel)]
    end
end

"""
    tconv(obj; channel, kernel)

Perform convolution in the time domain.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all channels
- `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`: kernel used for convolution

# Returns

- `s_convoluted::Array{Float64, 3}`: convoluted signal
"""
function tconv(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), kernel::Union{Vector{<:Real}, Vector{ComplexF64}})

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    s_convoluted = zeros(ch_n, epoch_len(obj), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_convoluted[ch_idx, :, ep_idx] = tconv(obj.data[channel[ch_idx], :, ep_idx], kernel=kernel)
        end
    end

    return s_convoluted
end
