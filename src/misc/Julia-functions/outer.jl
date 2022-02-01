"""
    outer(x, y)

Calculates outer product of vectors `x` and `y`.
"""

function outer(x, y)
    op = x * y';
    return op
end