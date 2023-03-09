export fft0
export fft2
export ifft0
export nextpow2
export dft

"""
    fft0(x, n)

Zeros-padded FFT.

# Arguments

- `x::AbstractVector`
- `n::Int64`: number of zeros to add

# Returns

- `fft0::Vector{ComplexF64}`
"""
function fft0(x::AbstractVector, n::Int64=0)

    n < 0 && throw(ArgumentError("n must be ≥ 0."))

    if CUDA.functional() && use_cuda
        # _free_gpumem()
        CUDA.memory_status()
        if n == 0
            cx = CuArray(x)
        else
            cx = CuArray(pad0(x, n))
        end
        return Vector(fft(cx))
    else
        if n == 0
            return fft(x)
        else
            return fft(pad0(x, n))
        end
    end
end

"""
    ifft0(x, n)

IFFT of zero-padded vector.

# Arguments

- `x::AbstractVector`
- `n::Int64`: number of zeros added to `x`

# Returns

- `signal::Vector{Float64}`: real part of the signal trimmed to original length
"""
function ifft0(x::AbstractVector, n::Int64=0)

    n < 0 && throw(ArgumentError("n must be ≥ 0."))

    if CUDA.functional() && use_cuda
        # _free_gpumem()
        x = Vector(ifft(CuArray(x)))
    else
        x = ifft(x)
    end
    return round.(real.(x[1:(length(x) - n)]))
end

"""
    fft2(x)

Zeros-padded FFT, so the length of padded vector is a power of 2.

# Arguments

- `x::AbstractVector`

# Returns

- `fft2::Vector{ComplexF64}`
"""
function fft2(x::AbstractVector)
    n = nextpow2(length(x)) - length(x)
    fft0(x, n)
end

"""
    nextpow2(x)

Return the next power of 2 for given number `x`.

# Argument

- `x::Int64`

# Returns

- `nextpow::Int64`
"""
function nextpow2(x::Int64)
    # return x == 0 ? 1 : (2 ^ ndigits(x - 1, base=2))
    return nextpow(2, x)
end

"""
    dft(signal; fs, pad)

Return FFT and DFT sample frequencies for a DFT.

# Arguments

- `signal::AbstractVector`
- `fs::Int64`: sampling rate
- `pad::Int64=0`: number of zeros to add

# Returns

Named tuple containing:
- `s_fft::Vector{ComplexF64}`: FFT
- `s_sf::Vector{Float64}`: sample frequencies
"""
function dft(signal::AbstractVector; fs::Int64, pad::Int64=0)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))

    s_fft = fft0(signal, pad)
    # number of samples
    n = length(signal)
    # time between samples
    d = 1 / fs
    s_sf = Vector(fftfreq(n, d))

    return (s_fft=s_fft, s_sf=s_sf)
end

"""
    dft(obj; channel)

Return FFT and DFT sample frequencies for a DFT.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all channels
- `pad::Int64=0`: number of zeros to add signal for FFT

# Returns

Named tuple containing:
- `s_fft::Array{ComplexF64, 3}`: FFT
- `s_sf::Vector{Float64}`: sample frequencies
"""
function dft(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), pad::Int64=0)

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    fs = sr(obj)
    s_fft = zeros(ComplexF64, ch_n, epoch_len(obj), ep_n)
    s_sf = nothing

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_fft[ch_idx, :, ep_idx], s_sf = @views s_dft(obj.data[channel[ch_idx], :, ep_idx], fs=fs, pad=pad)
        end
    end

    return (s_fft=s_fft, s_sf=s_sf)
end

