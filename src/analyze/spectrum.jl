export spectrum
export hspectrum

"""
    spectrum(signal; pad)

Calculate FFT, amplitudes, powers and phases.

# Arguments

- `signal::AbstractArray`
- `pad::Int64=0`: number of zeros to add
- `norm::Bool=false`: normalize do dB

# Returns

Named tuple containing:
- `s_fft::Vector{ComplexF64}`
- `s_amp::Vector{Float64}`
- `s_pow::Vector{Float64}`
- `s_pha::Vector{Float64}`
"""
function spectrum(signal::AbstractArray; pad::Int64=0, norm::Bool=false)

    s_fft = fft0(signal, pad)

    # amplitudes
    s_amp = abs.(s_fft) ./ length(signal)       # normalize
    s_amp = s_amp[1:(length(s_amp) ÷ 2)]        # remove negative frequencies
    s_amp[2:end] .*= 2                          # double positive frequencies
    # power
    s_pow = s_amp.^2
    norm == true && (s_pow = pow2db.(s_pow))
    # phases
    s_pha = angle.(s_fft)

    return (s_fft=s_fft, s_amp=s_amp, s_pow=s_pow, s_pha=s_pha)
end

"""
    hspectrum(signal; pad=0)

Calculate amplitudes, powers and phases using Hilbert transform.

# Arguments

- `signal::AbstractArray`
- `pad::Int64`: pad the `signal` with `pad` zeros
- `norm::Bool=true`: normalize do dB

# Returns

Named tuple containing:
- `h::Vector(ComplexF64}`: Hilbert components
- `h_amp::Vector{Float64}`
- `h_pow::Vector{Float64}`
- `h_pha::Vector{Float64}`
"""
function hspectrum(signal::AbstractArray; pad::Int64=0, norm::Bool=true)

    h = hilbert(pad0(signal, pad))

    # amplitudes
    h_amp = @. abs(h)
    # power
    h_pow = h_amp.^2
    norm == true && (h_pow = pow2db.(h_pow))
    # phases
    h_pha = angle.(h)

    return (h=h, h_amp=h_amp, h_pow=h_pow, h_pha=h_pha)
end

"""
    spectrum(obj; channel, pad, h)

Calculate FFT/Hilbert transformation components, amplitudes, powers and phases.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all channels
- `pad::Int64=0`: number of zeros to add signal for FFT
- `h::Bool=false`: use Hilbert transform for calculations instead of FFT
- `norm::Bool=false`: normalize do dB

# Returns

Named tuple containing:
- `c::Array{ComplexF64, 3}`: Fourier or Hilbert components
- `amp::Array{Float64, 3}`: amplitudes
- `pow::Array{Float64, 3}`: powers
- `pha::Array{Float64, 3}: phase angles
"""
function spectrum(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), pad::Int64=0, h::Bool=false, norm::Bool=false)

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    fft_size = epoch_len(obj) + pad

    s_c = zeros(ComplexF64, ch_n, fft_size, ep_n)
    s_pha = zeros(ch_n, fft_size, ep_n)
    if h == true
        s_amp = zeros(ch_n, fft_size, ep_n)
        s_pow = zeros(ch_n, fft_size, ep_n)
    else
        s_amp = zeros(ch_n, fft_size ÷ 2, ep_n)
        s_pow = zeros(ch_n, fft_size ÷ 2, ep_n)
    end        

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            if h == true
                s_c[ch_idx, :, ep_idx], s_amp[ch_idx, :, ep_idx], s_pow[ch_idx, :, ep_idx], s_pha[ch_idx, :, ep_idx] = @views hspectrum(obj.data[channel[ch_idx], :, ep_idx], pad=pad, norm=norm)
            else
                s_c[ch_idx, :, ep_idx], s_amp[ch_idx, :, ep_idx], s_pow[ch_idx, :, ep_idx], s_pha[ch_idx, :, ep_idx] = @views spectrum(obj.data[channel[ch_idx], :, ep_idx], pad=pad, norm=norm)
            end
        end
    end

    return (c=s_c, amp=s_amp, pow=s_pow, pha=s_pha)
end
