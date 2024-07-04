export pred_int

"""
    pred_int(n)

Calculate the prediction interval (95% CI adjusted for sample size)

# Arguments

- `n::Int64`: sample size

# Returns

- `pred_int::Float64`
"""
function pred_int(n::Int64)

    @assert n >= 1 "n must be ≥ 1."

    n in 1:19 && return [NaN, 15.56, 4.97, 3.56, 3.04, 2.78, 2.62, 2.51, 2.43, 2.37, 2.33, 2.29, 2.26, 2.24, 2.22, 2.18, 2.17, 2.16, 2.10][n]

    @warn "For n > 20 result may not be accurate."
    n in 20:25 && return 2.10
    n in 25:30 && return 2.08
    n in 31:35 && return 2.06
    n in 35:40 && return 2.05
    n in 41:50 && return 2.03
    n in 51:60 && return 2.02
    n in 61:70 && return 2.01
    n in 71:80 && return 2.00
    n in 81:90 && return 2.00
    n in 91:199 && return 1.99
    n in 200:200 && return 1.98
    n > 200 && return 1.96

end
