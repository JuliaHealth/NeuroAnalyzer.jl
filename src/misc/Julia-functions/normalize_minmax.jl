"""
    normalize_minmax(x)

Normalize the vector `x` to 0..1
"""

function normalize_minmax()
    normalize_minmax = (x - min(x)) / (max(x) - min(x));
    return normalize_minmax
end