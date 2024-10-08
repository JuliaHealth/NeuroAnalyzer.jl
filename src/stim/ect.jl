export ect_charge

"""
    ect_charge(; <keyword arguments>)

Calculate charge administered during ECT.

# Arguments

- `pw::Real`: pulse width [ms]
- `pint::Real`: pulse intensity [mA]
- `pf::Real`: pulse frequency [Hz]
- `duration::Real`: stimulation duration [s]

# Returns

- `charge::Float64`: charge [mC]
"""
function ect_charge(; pw::Real, pint::Real, pf::Real, duration::Real)::Float64

    return pw * pint * pf * duration

end

