export fft0
export fft2
export ifft0
export nextpow2
export rfft0
export rfft2
export fft_transform

"""
    fft0(x, n)

Perform zeros-padded FFT.

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
            return Vector(fft(CuVector(x)))
        else
            return Vector(fft(CuVector(pad0(x, n))))
        end
        CUDA.synchronize()
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

Perform IFFT of zero-padded vector.

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

Perform zeros-padded FFT, so the length of padded vector is a power of 2.

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

Return the next power of 2 for a given number.

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
    rfft0(x, n)

Perform zeros-padded single-sided FFT.

# Arguments

- `x::AbstractVector`
- `n::Int64`: number of zeros to add

# Returns

- `rfft0::Vector{ComplexF64}`
"""
function rfft0(x::AbstractVector, n::Int64=0)

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
    rfft2(x)

Perform zeros-padded single-sided FFT, so the length of padded vector is a power of 2.

# Arguments

- `x::AbstractVector`

# Returns

- `rfft2::Vector{ComplexF64}`
"""
function rfft2(x::AbstractVector)

    n = nextpow2(length(x)) - length(x)

    return rfft0(x, n)

end

"""
    fft_transform(x; fs, wlen, woverlap, w, demean, pad, mode)

Perform FFT transformation.

# Arguments

- `x::AbstractVector`
- `fs::Int64`: sampling rate
- `wlen::Int64=fs`: window length
- `woverlap::Int64=round(Int64, wlen * 0.97)`:
- `w::Bool=false`: if true, apply Hanning window per segment
- `demean::Bool=false`: if true, demean each segment
- `nfft::Int64=0`: length of input vector to the FFT; if nfft > n_samples, then the input signal will be zero-padded until it is of length nfft
- `mode::Symbol=:r`:
    - `:r`: use one-sided FFT (rfft)
    - `:f`: use two-sided FFT (fft)

# Returns

- `mf::Vector{ComplexF64}`: Fourier coefficients
- `f::Vector{Float64}`: frequencies
"""
function fft_transform(x::AbstractVector; fs::Int64, wlen::Int64=fs, woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=false, nfft::Int64=nextpow(2, wlen), demean::Bool=false, pad::Int64=0, mode::Symbol=:r)

    _check_var(mode, [:r, :f], "mode")
    @assert fs > 0 "fs must be > 0."
    @assert pad >= 0 "pad must be ≥ 0."

    # split into segments
    m = vec2mat(x, wlen=wlen, woverlap=woverlap)

    # remove mean by segments
    demean && (m = delmean(m, dims=2))

    # pad each segment with zeros
    n = size(m, 2)
    nfft > n && (m = pad0(m, nfft - n))

    # apply window per segment
    w = w ? hanning(size(m, 2)) : ones(size(m, 2))
    for idx in 1:size(m, 1)
        m[idx, :] = @views m[idx, :] .* w
    end

    # calculate FFT
    mf = nothing
    if mode === :r
        mf = vec(reshape(rfft(m, 2), 1, :))
        f = round.(rfftfreq(((length(mf) * 2) - 1), fs), digits=3)
        mf .*= 2
    elseif mode === :f
        mf = vec(reshape(fft(m, 2), 1, :))
        f = round.(fftfreq(length(m), fs), digits=3)
    end

    # scale
    # s = 1.0 / length(x)
    # mf .*= s

    return mf, f

end
