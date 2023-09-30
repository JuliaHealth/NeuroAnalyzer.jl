export denoise_fft
export denoise_fft!

"""
   denoise_fft(s; pad, t)

Perform FFT denoising.

# Arguments

- `s::AbstractVector`
- `pad::Int64=0`: number of zeros to add
- `t::Real=0`: PSD threshold for keeping frequency components; if 0, use mean signal power value

# Returns

Named tuple containing:
- `s::Vector{Float64}`
- `f_idx::BitVector`: index of components zeroed
"""
function denoise_fft(s::AbstractVector; pad::Int64=0, t::Real=0)

    s_fft = fft0(s, pad)
    s_pow = (real.(s_fft .* conj.(s_fft))) ./ length(s)
    
    t == 0 && (t = mean(s_pow))

    # zero frequencies with power above threshold
    f_idx = s_pow .> t
    s_fft[f_idx] .= Complex(0, 0)

    return (s=abs.(ifft0(s_fft, pad)), f_idx=f_idx)

end

"""
    denoise_fft(s; pad, t)

Perform FFT denoising.

# Arguments

- `s::AbstractArray`
- `pad::Int64=0`: number of zeros to add
- `t::Real=0`: PSD threshold for keeping frequency components; if 0, use mean signal power value

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function denoise_fft(s::AbstractArray; pad::Int64=0, t::Real=0)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx], _ = @views denoise_fft(s[ch_idx, :, ep_idx], pad=pad, t=t)
        end
    end

    return s_new
end

"""
    denoise_fft(obj; ch, pad, t)

Perform FFT denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
- `pad::Int64=0`: number of zeros to add signal for FFT
- `t::Int64=100`: PSD threshold for keeping frequency components

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function denoise_fft(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)), pad::Int64=0, t::Int64=100)

    _check_channels(obj, ch)

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = @views denoise_fft(obj.data[ch, :, :], pad=pad, t=t)
    reset_components!(obj_new)
    push!(obj_new.history, "denoise_fft(OBJ, ch=$ch, pad=$pad, t=$t)")

    return obj_new

end

"""
    denoise_fft!(obj; ch, pad, t)

Perform FFT denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
- `pad::Int64=0`: number of zeros to add signal for FFT
- `t::Int64=100`: PSD threshold for keeping frequency components
"""
function denoise_fft!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)), pad::Int64=0, t::Int64=100)

    obj_new = denoise_fft(obj, ch=ch, pad=pad, t=t)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end

