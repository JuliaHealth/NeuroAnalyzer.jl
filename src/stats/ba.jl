export ba

"""
    ba(x, y; <keyword arguments>)

Calculate Bland-Altman comparison between two clinical measurements.

# Arguments

- `x::AbstractVector`
- `y::AbstractVector`
- `la::FLoat64=0.95`: 95% limits of agreement for each comparison (average difference ± 1.96 standard deviation of the difference)

# Returns

Named tuple containing:
- `m::Float64`: mean of the difference
- `s_u::Float64`: +1.96 (for default la) × standard deviation of the difference
- `s_d::Float64`: -1.96 (for default la) × standard deviation of the difference

# Notes

To plot the Bland-Altman plot: `p = hline([0, m]); hline!([0, s_u]); hline!([0, s-d])`
"""
function ba(x::AbstractVector, y::AbstractVector; la::Float64=0.95)::@NamedTuple{m::Float64, s_u::Float64, s_d::Float64}

    @assert la > 0 "la must be > 0.0."
    @assert la < 1 "la must be < 1.0."

    z = p2z(1 - la, twotailed=true)
    d = x .- y

    m = mean(d)
    s_u = z * std(d)
    s_d = -z * std(d)

    return (m=m, s_u=s_u, s_d = s_d)

end