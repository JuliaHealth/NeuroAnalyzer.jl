export denoise_fft
export denoise_fft!

"""
   denoise_fft(signal; pad, threshold) 

Perform FFT denoising.

# Arguments

- `signal::AbstractVector`
- `pad::Int64=0`: number of zeros to add
- `threshold::Real=0`: PSD threshold for keeping frequency components; if 0, use mean signal power value

# Returns

Named tuple containing:
- `s_denoised::Vector{Float64}`
- `frq_idx::BitVector`: index of components zeroed
"""
function denoise_fft(signal::AbstractVector; pad::Int64=0, threshold::Real=0)

    s_fft = fft0(signal, pad)
    s_pow = (real.(s_fft .* conj.(s_fft))) ./ length(signal)
    
    threshold == 0 && (threshold = mean(s_pow))
    
    # zero frequencies with power above threshold
    frq_idx = s_pow .> threshold
    s_fft[frq_idx] .= Complex(0, 0)

    return (s_denoised=ifft0(s_fft, pad), frq_idx=frq_idx)
end

"""
    denoise_fft(obj; channel, pad, threshold)

Perform FFT denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `pad::Int64=0`: number of zeros to add signal for FFT
- `threshold::Int64=100`: PSD threshold for keeping frequency components

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function denoise_fft(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), pad::Int64=0, threshold::Int64=100)

    _check_channels(obj, channel)

    ep_n = epoch_n(obj)

    obj_new = deepcopy(obj)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in eachindex(channel)
                obj_new.data[channel[ch_idx], :, ep_idx], _ = @views denoise_fft(obj_new.data[channel[ch_idx], :, ep_idx], pad=pad, threshold=threshold)
        end
    end

    reset_components!(obj_new)
    push!(obj_new.header.history, "denoise_fft(OBJ, channel=$channel, pad=$pad, threshold=$threshold)")

    return obj_new
end

"""
    denoise_fft!(obj; channel, pad, threshold)

Perform FFT denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `pad::Int64=0`: number of zeros to add signal for FFT
- `threshold::Int64=100`: PSD threshold for keeping frequency components
"""
function denoise_fft!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), pad::Int64=0, threshold::Int64=100)

    obj_tmp = denoise_fft(obj, channel=channel, pad=pad, threshold=threshold)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    reset_components!(obj)

    return nothing
end

