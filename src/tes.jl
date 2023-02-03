"""
    tes_dose(current, pad_area, duration)

Convert `current`, `pad_area` and stimulation `duration` into `charge`, `current_density` and `charge_ density`.

# Arguments

- `current::Real`: stimulation current [mA]
- `pad_area::Real`: electrode pad area [cm²]
- `duration::Int64`: stimulation duration [s]

# Returns

Named tuple containing:
- `charge::Float64`: charge [C]
- `current_density::Float64`: current density [A/m²]
- `charge_density::Float64`: delibvered charge density [kC/m²]

# Source

Chhatbar PY, George MS, Kautz SA, Feng W. Quantitative reassessment of safety limits of tDCS for two animal studies. Brain Stimulation. 2017;10(5):1011–2.
"""
function tes_dose(;current::Real, pad_area::Real, duration::Int64)
    
    charge = (current / 1_000) * duration
    current_density = (current / 1_000) / (pad_area / 1_000)
    charge_density = (charge / 1_000) / (pad_area / 1_000)
    
    return (charge=charge, current_density=current_density, charge_density=charge_density)
end

"""
    ect_charge(; pw, pint, pf, duration)

Calculate charge administered during ECT.

# Arguments

- `pw::Real`: pulse width [ms]
- `pint::Real`: pulse intensity [mA]
- `pf::Real`: pulse frequency [Hz]
- `duration::Real`: stimulation duration [s]

# Returns

- `charge::Float64`: charge [mC]
"""
function ect_charge(; pw::Real, pint::Real, pf::Real, duration::Real)
    return pw * pint * pf * duration
end