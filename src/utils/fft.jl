export fft0
export fft2
export ifft0
export nextpow2
export rfft0
export rfft2

"""
    fft0(x, n)

Perform zeros-padded FFT.

# Arguments

- `x::AbstractVector`
- `n::Int64`: number of zeros to add

# Returns

- `fft0::Vector{ComplexF64}`
"""
function fft0(x::AbstractVector, n::Int64=0)::Vector{ComplexF64}

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
function ifft0(x::AbstractVector, n::Int64=0, cnorm::Bool=false)::Vector{ComplexF64}

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
function fft2(x::AbstractVector)::Vector{ComplexF64}

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
function nextpow2(x::Int64)::Int64

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
function rfft0(x::AbstractVector, n::Int64=0)::Vector{ComplexF64}

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
function rfft2(x::AbstractVector)::Vector{ComplexF64}

    n = nextpow2(length(x)) - length(x)

    return rfft0(x, n)

end
