"""
    rms(vector)

Calculates root mean square of the vector `x`.
"""

function rms(v)
    rms = sqrt(mean(x.^2))
    return rms
end