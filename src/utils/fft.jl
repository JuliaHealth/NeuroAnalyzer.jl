export fft0
export fft2
export ifft0
export nextpow2
export rfft0
export rfft2

"""
    fft0(x, n)

Perform a zero-padded full (two-sided) FFT.

Appends `n` zeros to `x` before computing the FFT. When `n = 0` the input is transformed directly, avoiding an unnecessary copy.

# Arguments

- `x::AbstractVector`: signal vector of length `L`
- `n::Int64`: number of zeros to append; must be ≥ 0.

# Returns

- `fft0::Vector{ComplexF64}`: two-sided Fourier coefficients of length `L + n`.

# Throws
- `ArgumentError`: If `n < 0`.

# See also
[`fft2`](@ref), [`rfft0`](@ref)
"""
function fft0(x::AbstractVector, n::Int64 = 0)::Vector{ComplexF64}

    @assert n >= 0 "n must be ≥ 0."

    # when n == 0, skip pad0() to avoid an unnecessary copy.
    return n == 0 ? fft(x) : fft(pad0(x, n))

end

"""
    ifft0(x, n)

Perform an IFFT on a zero-padded spectrum and trim the result to the original signal length.

If a signal of length `L` was zero-padded by `n` samples before the forward FFT, the padded spectrum has length `L + n`. This function recovers the original `L`-sample signal by discarding the trailing `n` samples after the IFFT.

# Arguments

- `x::AbstractVector`: zero-padded spectrum (length `L + n`)
- `n::Int64`: number of zeros that were appended during the forward FFT; must be ≥ 0

# Returns

- `ifft0::Vector{ComplexF64}`: reconstructed signal of length `length(x) - n`

# Throws

- `ArgumentError`: If `n < 0` or `n ≥ length(x)`.

# See also

[`fft0`](@ref)
"""
function ifft0(x::AbstractVector, n::Int64 = 0)::Vector{ComplexF64}

    @assert n >= 0 "n must be ≥ 0."
    @assert n < length(x)  "n must be < length(x); got n=$n, length(x)=$(length(x))."

    # when n == 0 no trimming is needed; return the full IFFT directly
    n == 0 && return ifft(x)

    # trim trailing n samples to recover the original signal length L = length(x) - n
    return ifft(x)[1:(length(x) - n)]

end

"""
    fft2(x)

Perform a full (two-sided) FFT after zero-padding the input to the next power of 2.

Zero-padding to a power-of-2 length maximizes FFT efficiency (radix-2 algorithm). If `length(x)` is already a power of 2, no padding is added.

# Arguments

- `x::AbstractVector`: signal vector of length `L`

# Returns

- `fft2::Vector{ComplexF64}`: two-sided Fourier coefficients of length `nextpow2(L)`

# See also

[`fft0`](@ref), [`rfft2`](@ref), [`nextpow2`](@ref)
"""
function fft2(x::AbstractVector)::Vector{ComplexF64}

    # compute the number of zeros needed to reach the next power of 2
    n = nextpow2(length(x)) - length(x)

    return fft0(x, n)

end

"""
    nextpow2(x)

Return the smallest power of 2 that is ≥ `x`.

Thin wrapper around `Base.nextpow(2, x)` with an explicit positivity guard.

# Arguments

- `x::Int64`: input value; must be > 0

# Returns

- `nextpow2::Int64`: smallest integer `p` such that `p = 2^k ≥ x` for some `k ≥ 0`

# Throws

- `ArgumentError`: if `x ≤ 0`

# Examples

```julia
nextpow2(5)   # → 8
nextpow2(8)   # → 8
nextpow2(9)   # → 16
```
"""
function nextpow2(x::Int64)::Int64

    return nextpow(2, x)

end

"""
    rfft0(x, n)

Perform a zero-padded one-sided (real) FFT.

Appends `n` zeros to `x` before computing `rfft`, returning only the positive-frequency coefficients. When `n = 0` the input is transformed directly, avoiding an unnecessary copy.

# Arguments

- `x::AbstractVector`: signal vector of length `L`
- `n::Int64`: number of zeros to append; must be ≥ 0.

# Returns

- `rfft0::Vector{ComplexF64}`: one-sided Fourier coefficients of length `(L + n) ÷ 2 + 1`

# Throws

- `ArgumentError`: if `n < 0`

# See also

[`rfft2`](@ref), [`fft0`](@ref)
"""
function rfft0(x::AbstractVector, n::Int64 = 0)::Vector{ComplexF64}

    @assert n >= 0 "n must be ≥ 0."

    # when n == 0, skip pad0() to avoid an unnecessary copy
    return n == 0 ? rfft(x) : rfft(pad0(x, n))

end

"""
    rfft2(x)

Perform a one-sided (real) FFT after zero-padding the input to the next power of 2.

Zero-padding to a power-of-2 length maximizes FFT efficiency (radix-2 algorithm). If `length(x)` is already a power of 2, no padding is added.

# Arguments

- `x::AbstractVector`: signal vector of length `L`

# Returns

- `rfft2::Vector{ComplexF64}`: one-sided Fourier coefficients of length `nextpow2(L) ÷ 2 + 1`

# See also

[`rfft0`](@ref), [`fft2`](@ref), [`nextpow2`](@ref)
"""
function rfft2(x::AbstractVector)::Vector{ComplexF64}

    # number of zeros needed to reach the next power-of-2 length
    n = nextpow2(length(x)) - length(x)

    return rfft0(x, n)

end
