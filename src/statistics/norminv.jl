export norminv

"""
    norminv(x::Real)

Convert probability to a normal distribution with a peak at 0.5.

# Arguments

- `x::Real`

# Returns

- `norminv::Float64`
"""
function norminv(x::Real)

    return quantile(Distributions.Normal(), x)

end
