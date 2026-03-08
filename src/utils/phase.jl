export phases

"""
    phases(s)

Computes the instantaneous phase of a signal via the Hilbert transform: φ(t) = angle( H(s)(t) )  ∈ (−π, π]

The result is the wrapped instantaneous phase in radians. For the unwrapped phase (suitable for differentiation to obtain instantaneous frequency) use DSP.unwrap() on the output.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `phases::Vector{Float64}`: instantaneous phase in radians ∈ (−π, π]
"""
function phases(s::AbstractVector)::Vector{Float64}

    # DSP.hilbert() returns the analytic signal z = s + i·H(s)
    # angle(z) = atan(imag(z), real(z)) is the instantaneous phase
    return DSP.angle.(DSP.hilbert(s))

end
