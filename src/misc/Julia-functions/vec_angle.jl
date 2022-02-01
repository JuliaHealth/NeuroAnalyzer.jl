using LinearAlgebra

"""
    vec_angle(x, y)

Calculates the angle between two vectors `x` and `y`.
"""

function vec_angle(x, y)
    dotxy = dot(x, y)
    normx = norm(x)
    normy = norm(y)
    angle_r = dotxy / (normx * normy)
    angle_d = (angle_r * 180) / pi
    return angle_r, angle_d
end
