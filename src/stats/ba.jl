export ba

"""
    ba(x, y; <keyword arguments>)

Calculate a Bland-Altman comparison between two sets of clinical measurements.

Computes the mean difference and the upper and lower limits of agreement (LoA), defined as `mean(d) ± z × std(d)`, where `d = x − y` and `z` is the two-tailed Z-score corresponding to `la`.

# Arguments

- `x::AbstractVector`: first measurement vector
- `y::AbstractVector`: second measurement vector; must be the same length as `x`
- `la::Float64=0.95`: confidence level for the limits of agreement; must be in `(0, 1)`; default gives the standard 95 % LoA (mean ± 1.96 SD)

# Returns

Named tuple:

- `m::Float64`: mean of the differences `x − y`
- `s_u::Float64`: upper limit of agreement (`m + z × SD`)
- `s_d::Float64`: lower limit of agreement (`m − z × SD`)

# Throws

- `ArgumentError`: if `la ∉ (0, 1)`, vectors are empty, or lengths differ

# Notes

To produce a Bland-Altman plot:

```julia
means = (x .+ y) ./ 2
diffs = x .- y
scatter(means, diffs)
hline!([m, s_u, s_d])
```

# See also

[`p2z`](@ref)
"""
function ba(
    x::AbstractVector,
    y::AbstractVector;
    la::Float64 = 0.95
)::@NamedTuple{
    m::Float64,
    s_u::Float64,
    s_d::Float64
}

    !(length(x) > 0) && throw(ArgumentError("x must not be empty."))
    !(length(x) == length(y)) && throw(ArgumentError("x and y must have the same length."))
    _in(la, (0, 1), "la")

    # two-tailed Z-score for the requested confidence level
    z = p2z(1 - la; twotailed=true)

    d = x .- y
    m = mean(d)
    sd = std(d)

    s_u = z * sd
    s_d = -z * sd

    return (m=m, s_u=s_u, s_d=s_d)

end
