export norminv

"""
    norminv(x::Real)

Convert probability to a normal distribution with a peak at 0.5.

# Arguments

- `x::Real`

# Returns

- `z::Float64`
"""
function norminv(x::Real)
    return quantile(Normal(), x)
end
