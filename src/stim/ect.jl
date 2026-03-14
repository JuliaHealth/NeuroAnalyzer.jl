export ect_charge

"""
    ect_charge(; <keyword arguments>)

Calculate the total charge delivered during an ECT stimulation session.

Each pulse contributes `pw × pint` of charge (pulse width × current = charge per pulse). Multiplying by the number of pulses (`pf × duration`) gives the total charge. The milli-prefixes cancel as follows:

    ms × mA × Hz × s = 10⁻³ s × 10⁻³ A × s⁻¹ × s = 10⁻⁶ C = 10⁻³ mC

so a factor of `10⁻³` is applied to convert the raw product to mC.

# Arguments

- `pw::Real`: pulse width in ms; must be > 0
- `pint::Real`: pulse intensity (current) in mA; must be > 0
- `pf::Real`: pulse frequency in Hz; must be > 0
- `duration::Real`: stimulation duration in seconds; must be > 0

# Returns

- `Float64`: total charge in mC

# Throws
- `ArgumentError`: if any argument is ≤ 0
"""
function ect_charge(; pw::Real, pint::Real, pf::Real, duration::Real)::Float64

    !(pw > 0) && throw(ArgumentError("pw must be > 0."))
    !(pint > 0) && throw(ArgumentError("pint must be > 0."))
    !(pf > 0) && throw(ArgumentError("pf must be > 0."))
    !(duration > 0) && throw(ArgumentError("duration must be > 0."))

    # unit conversion: ms × mA × Hz × s → mC requires a factor of 10⁻³
    return pw * pint * pf * duration * 1e-3

end
