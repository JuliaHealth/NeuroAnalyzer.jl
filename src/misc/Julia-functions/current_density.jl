"""
    current_density(current, area)

Calculates TES current density for `current` [mA] over `area` [cm^2].
"""

function current_density(current, area)
    # convert mA to A
    current = current/1000
    # convert cm^2 to m^2
    area = area/10000
    current_density = current / area
    return(current_density)
end
