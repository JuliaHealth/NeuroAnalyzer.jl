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

- `x::AbstractVector`: signal vector
- `n::Int64`: number of zeros to append before transforming

# Returns

- `fft0::Vector{ComplexF64}`: full (two-sided) Fourier coefficients
"""
function fft0(x::AbstractVector, n::Int64 = 0)::Vector{ComplexF64}

    @assert n >= 0 "n must be ≥ 0."

    # when n == 0, skip pad0() to avoid an unnecessary copy.
    return n == 0 ? fft(x) : fft(pad0(x, n))

end

"""
    ifft0(x, n)

Perform IFFT of a zero-padded spectrum and trim the result to the original length.

If a signal of length `L` was zero-padded by `n` samples before the forward FFT, the padded spectrum has length `L + n`. This function recovers the original `L`-sample signal by discarding the trailing `n` samples after the IFFT.

# Arguments

- `x::AbstractVector`: zero-padded spectrum (length `L + n`)
- `n::Int64`: number of zeros that were appended to the original signal

# Returns

- `ifft0::Vector{ComplexF64}`: reconstructed signal of length `length(x) - n`
"""
function ifft0(x::AbstractVector, n::Int64 = 0)::Vector{ComplexF64}

    @assert n >= 0 "n must be ≥ 0."

    # when n == 0 no trimming is required; return the full IFFT
    n == 0 && return ifft(x)

    # trim the IFFT output back to the original signal length
    return ifft(x)[1:(length(x) - n)]

end

"""
    fft2(x)

Perform zero-padded FFT, padding the input to the next power of 2 in length.

# Arguments

- `x::AbstractVector`: signal vector

# Returns

- `fft2::Vector{ComplexF64}`: full (two-sided) Fourier coefficients
"""
function fft2(x::AbstractVector)::Vector{ComplexF64}

    # compute the number of zeros needed to reach the next power of 2
    n = nextpow2(length(x)) - length(x)

    return fft0(x, n)

end

"""
    nextpow2(x)

Return the smallest power of 2 that is ≥ `x`.

# Arguments

- `x::Int64`

# Returns

- `nextpow2::Int64`
"""
function nextpow2(x::Int64)::Int64

    return nextpow(2, x)

end

"""
    rfft0(x, n)

Perform zero-padded one-sided FFT (positive frequencies only).

# Arguments

- `x::AbstractVector`: signal vector
- `n::Int64`: number of zeros to append before transforming

# Returns

- `rfft0::Vector{ComplexF64}`: one-sided Fourier coefficients
"""
function rfft0(x::AbstractVector, n::Int64 = 0)::Vector{ComplexF64}

    @assert n >= 0 "n must be ≥ 0."

    # when n == 0, skip pad0() to avoid an unnecessary copy
    return n == 0 ? rfft(x) : rfft(pad0(x, n))

end

"""
    rfft2(x)

Perform zero-padded one-sided FFT, padding the input to the next power of 2 in length.

# Arguments

- `x::AbstractVector`: signal vector

# Returns

- `rfft2::Vector{ComplexF64}`: one-sided Fourier coefficients
"""
function rfft2(x::AbstractVector)::Vector{ComplexF64}

    n = nextpow2(length(x)) - length(x)

    return rfft0(x, n)

end
