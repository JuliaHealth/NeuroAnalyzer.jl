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

    @assert n >=0 "n must be ≥ 0."

    if CUDA.functional() && use_cuda
        # CUDA.memory_status()
        # _free_gpumem()
        CUDA.synchronize()
        if n == 0
            return Vector(rfft(CuVector(x)))
        else
            return Vector(rfft(CuVector(pad0(x, n))))
        end
        CUDA.synchronize()
    else
        if n == 0
            return rfft(x)
        else
            return rfft(pad0(x, n))
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

- `ifft0::Vector{ComplexF64}`: reconstructed signal trimmed to original length
"""
function ifft0(x::AbstractVector, n::Int64=0)

    @assert n >= 0 "n must be ≥ 0."

    if CUDA.functional() && use_cuda
        CUDA.synchronize()
        x = Vector(ifft(CuVector(x)))
        CUDA.synchronize()
    else
        x = ifft(x)
    end

    return x[1:(length(x) - n)]

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
    
    return fft0(x, n)

end

"""
    nextpow2(x)

Return the next power of 2 for given number `x`.

# Argument

- `x::Int64`

# Returns

- `nextpow2::Int64`
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
- `ft::Vector{ComplexF64}`: FFT
- `f::Vector{Float64}`: sample frequencies
"""
function dft(signal::AbstractVector; fs::Int64, pad::Int64=0)

    @assert fs >= 1 "fs must be ≥ 1."

    ft = fft0(signal, pad)

    # number of samples
    n = length(signal)
    
    # time between samples
    d = 1 / fs
    f = Vector(fftfreq(n, d))

    return (ft=ft, f=f)

end

"""
    dft(obj; channel)

Return FFT and DFT sample frequencies for a DFT.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `pad::Int64=0`: number of zeros to add signal for FFT

# Returns

Named tuple containing:
- `ft::Array{ComplexF64, 3}`: FFT
- `f::Vector{Float64}`: sample frequencies
"""
function dft(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), pad::Int64=0)

    _check_channels(obj, ch)
    ch_n = length(ch)
    ep_n = nepochs(obj)

    fs = sr(obj)
    ft = zeros(ComplexF64, ch_n, epoch_len(obj), ep_n)
    f = nothing

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ft[ch_idx, :, ep_idx], f = @views dft(obj.data[ch[ch_idx], :, ep_idx], fs=fs, pad=pad)
        end
    end

    return (ft=ft, f=f)
    
end
