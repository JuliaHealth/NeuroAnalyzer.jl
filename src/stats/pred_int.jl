export pred_int

"""
    pred_int(n)

Return the 95 % prediction interval multiplier adjusted for sample size `n`.

Values for `n ≤ 19` are taken from a lookup table. Values for `n ≥ 20` are approximated from published stepped ranges; a warning is issued because accuracy decreases for larger samples. For very large samples (`n > 200`) the value converges to the normal-distribution limit of 1.96.

Returns `NaN` for `n = 1` (prediction interval is undefined for a single observation).

# Arguments

- `n::Int64`: sample size; must be ≥ 1

# Returns

- `Float64`: prediction interval multiplier

# Throws

- `ArgumentError`: if `n < 1`

# Notes

For `n > 20` the result is approximate and a diagnostic warning is issued.
"""
function pred_int(n::Int64)::Float64

    !(n >= 1) && throw(ArgumentError("n must be ≥ 1."))

    # exact lookup table for n = 1 … 19 (index = n)
    # n = 1: undefined → NaN; n = 2: 15.56; … ; n = 19: 2.10
    table = [
        NaN, 15.56, 4.97, 3.56, 3.04, 2.78, 2.62, 2.51,
        2.43, 2.37, 2.33, 2.29, 2.26, 2.24, 2.22, 2.18,
        2.17, 2.16, 2.10
    ]
    n <= 19 && return table[n]

    _warn("For n > 19 the prediction interval multiplier is approximate.")

    # stepped approximation for n ≥ 20
    n <= 24  && return 2.10
    n <= 29  && return 2.08
    n <= 34  && return 2.06
    n <= 39  && return 2.05
    n <= 49  && return 2.03
    n <= 59  && return 2.02
    n <= 69  && return 2.01
    n <= 89  && return 2.00
    n <= 199 && return 1.99
    n <= 200 && return 1.98
    return 1.96

end
