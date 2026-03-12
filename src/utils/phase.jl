export phases

"""
    phases(s)

Compute the instantaneous phase of a signal via the Hilbert transform. The analytic signal is obtained as `z = DSP.hilbert(s)`, then the wrapped instantaneous phase is extracted as `φ(t) = angle(z(t)) = atan(imag(z), real(z))`.

For the unwrapped phase (suitable for differentiating to obtain instantaneous frequency) apply `DSP.unwrap` to the output.

# Arguments

- `s::AbstractVector`: signal vector; must contain at least 1 element

# Returns

- `Vector{Float64}`: instantaneous phase in radians ∈ (−π, π].

# Throws

- `ArgumentError`: if `s` is empty

# See also

[`DSP.hilbert`](https://docs.juliadsp.org), [`DSP.unwrap`](https://docs.juliadsp.org)
"""
function phases(s::AbstractVector)::Vector{Float64}

    @assert length(s) > 0 "s must not be empty."
    # DSP.hilbert() returns the analytic signal z = s + i·H(s)
    # Base.angle(z) = atan(imag(z), real(z)) gives the wrapped instantaneous phase
    return Base.angle.(DSP.hilbert(s))

end
