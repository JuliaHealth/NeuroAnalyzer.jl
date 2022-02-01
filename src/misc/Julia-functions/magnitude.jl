"""
    magnitude(vector)

Calculates `vector` magnitude

||v|| = sqrt(v_1^2 + v_2^2 + ... + v_n^2)
"""

function magnitude(v)
    m = sqrt(sum(v.^2))
    return m
end