################################
#                              #
# Low-level internal functions #
#                              #
################################

################################

"""
    tes_dose(current, pad_area, duration)

Converts `current`, `pad_area` and stimulation `duration` into `charge`, `current_density` and `charge_ density`.

# Arguments

- `current::Real`: stimulation current [mA]
- `pad_area::Real`: electrode pad area [cm^2]
- `duration::Int64`: stimulation duration [s]

# Returns

- `charge::Float64`: charge [C]
- `current_density::Float64`: current density [A/m^2]
- `charge_density::Float64`: delibvered charge density [kC/m^2]

# Source

Chhatbar PY, George MS, Kautz SA, Feng W. Quantitative reassessment of safety limits of tDCS for two animal studies. Brain Stimulation. 2017;10(5):1011–2.
"""
function tes_dose(current::Real, pad_area::Real, duration::Int64)
    
    charge = (current / 1_000) * duration
    current_density = (current / 1_000) / (pad_area / 1_000)
    charge_density = (charge / 1_000) / (pad_area / 1_000)
    
    return charge, current_density, charge_density
end