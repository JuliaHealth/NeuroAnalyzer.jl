"""
    _fir_response(f, w)

Compute the complex frequency response of a FIR filter.

Evaluates `H(ω) = Σₖ f[k] × exp(−im × ω × (k−1))` for each frequency in `w`, where `k` is 1-based and the first coefficient `f[1]` corresponds to delay 0.

# Arguments

- `f::Vector{<:Real}`: FIR filter coefficients
- `w=range(0, π; length=1024)`: frequency grid in radians; default is 1024 equally spaced points in `[0, π]`

# Returns

- `Vector{ComplexF64}`: complex frequency response, one value per frequency in `w`

# References

Based on Matti Pastell, "FIR filter design with Julia".
"""
function _fir_response(
    f::Vector{<:Real},
    w = range(0, π; length=1024)
)::Vector{ComplexF64}

    n = length(w)
    h = Vector{ComplexF64}(undef, n)

    for i in 1:n
        h[i] = sum(f[j] * exp(-im * w[i] * (j - 1)) for j in eachindex(f))
    end

    return h
end