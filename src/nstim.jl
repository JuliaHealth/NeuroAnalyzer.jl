"""
    tes_dose(current, pad_area, duration)

Converts `current`, `pad_area` and stimulation `duration` into:
- charge (C)
- current density (A/m^2)
- charge density (kC/m^2)

# Source

Chhatbar PY, George MS, Kautz SA, Feng W. Quantitative reassessment of safety limits of tDCS for two animal studies. Brain Stimulation. 2017;10(5):1011–2.
"""

function tes_dose(current, pad_area, duration)
    charge = (current / 1_000) * duration
    current_density = (current / 1_000) / (pad_area / 1_000)
    charge_density = (current / 1_000) / (pad_area / 10_000)
    return charge, current_density, charge_density
end