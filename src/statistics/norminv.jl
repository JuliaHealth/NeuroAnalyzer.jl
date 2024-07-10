export norminv

"""
    norminv(x::Real)

Convert probability to a normal distribution with a peak at 0.5.

# Arguments

- `x::Real`

# Returns

- `n::Float64`
"""
function norminv(x::Real)

    n = quantile(Distributions.Normal(), x)

    return n

end
