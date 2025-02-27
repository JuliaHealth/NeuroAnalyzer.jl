export denoise_fft
export denoise_fft!

"""
    denoise_fft(s; <keyword arguments>)

Perform FFT denoising.

# Arguments

- `s::AbstractVector`
- `pad::Int64=0`: number of zeros to add
- `t::Real=0`: PSD threshold for keeping frequency components; if 0, use mean signal power value

# Returns

Named tuple containing:
- `s::Vector{Float64}`: denoised signal
- `f_idx::BitVector`: index of components zeroed
"""
function denoise_fft(s::AbstractVector; pad::Int64=0, t::Real=0)::@NamedTuple{s::Vector{Float64}, f_idx::BitVector}

    s_fft = fft0(s, pad)
    s_pow = @. (abs(s_fft * conj(s_fft))) / length(s)

    t == 0 && (t = mean(s_pow))

    # zero frequencies with power above threshold
    f_idx = s_pow .> t
    s_fft[f_idx] .= Complex(0, 0)

    return (s=abs.(ifft0(s_fft, pad)), f_idx=f_idx)

end

"""
    denoise_fft(s; <keyword arguments>)

Perform FFT denoising.

# Arguments

- `s::AbstractArray`
- `pad::Int64=0`: number of zeros to add
- `t::Real=0`: PSD threshold for keeping frequency components; if 0, use mean signal power value

# Returns

- `s_new::Array{Float64, 3}`
"""
function denoise_fft(s::AbstractArray; pad::Int64=0, t::Real=0)::Array{Float64, 3}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)

    @inbounds for ep_idx in 1:ep_n
        if use_cuda
            CUDA.synchronize()
            for ch_idx in 1:ch_n
                s_new[ch_idx, :, ep_idx], _ = @views denoise_fft(s[ch_idx, :, ep_idx], pad=pad, t=t)
            end
            CUDA.synchronize()
        else
            Threads.@threads :static for ch_idx in 1:ch_n
                s_new[ch_idx, :, ep_idx], _ = @views denoise_fft(s[ch_idx, :, ep_idx], pad=pad, t=t)
            end
        end
    end

    return s_new
end

"""
    denoise_fft(obj; <keyword arguments>)

Perform FFT denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `pad::Int64=0`: number of zeros to add signal for FFT
- `t::Int64=100`: PSD threshold for keeping frequency components

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function denoise_fft(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, pad::Int64=0, t::Int64=100)::NeuroAnalyzer.NEURO

    ch = get_channel(obj, ch=ch)
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = @views denoise_fft(obj.data[ch, :, :], pad=pad, t=t)
    reset_components!(obj_new)
    push!(obj_new.history, "denoise_fft(OBJ, ch=$ch, pad=$pad, t=$t)")

    return obj_new

end

"""
    denoise_fft!(obj; <keyword arguments>)

Perform FFT denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `pad::Int64=0`: number of zeros to add signal for FFT
- `t::Int64=100`: PSD threshold for keeping frequency components

# Returns

Nothing
"""
function denoise_fft!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, pad::Int64=0, t::Int64=100)::Nothing

    obj_new = denoise_fft(obj, ch=ch, pad=pad, t=t)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end

