export res_norm

"""
    res_norm(x, g)

Test whether residuals follow a normal distribution, per group and overall.

For each group (and for the whole sample), residuals are computed as `x .- mean(x_group)`. Normality is then assessed with:

- Anderson–Darling one-sample test against `Normal(0, 1)`.
- Kolmogorov–Smirnov one-sample exact test against `Normal(0, 1)`.

# Arguments
- `x::AbstractVector`: data values; must not be empty
- `g::Vector{Int64}=repeat([1], length(x))`: group membership for each element of `x`; must have the same length as `x`

# Returns

Named tuple:
- `adt_p::Vector{Float64}`: Anderson–Darling p-values; one per group (in sorted group order) plus one for the whole sample at the last index
- `ks_p::Vector{Float64}`: Kolmogorov–Smirnov p-values; same layout as `adt_p`

If there is only one group, both vectors have length 1 (whole-sample result only).

# Throws
- `ArgumentError`: if `x` is empty, `length(x) ≠ length(g)`, or any group has fewer than 3 observations

# Notes

- Results are non-random: both tests compare residuals against the theoretical `Normal(0, 1)` distribution rather than a random reference sample.
- For large groups `ExactOneSampleKSTest` may be slow; consider wrapping in `ApproximateOneSampleKSTest` for `n > 1000`.
"""
function res_norm(
    x::AbstractVector,
    g::Vector{Int64} = repeat([1], length(x))
)::@NamedTuple{
    adt_p::Vector{Float64},
    ks_p::Vector{Float64}
}

    !(length(x) > 0) && throw(ArgumentError("x must not be empty."))
    !(length(x) == length(g)) && throw(ArgumentError("x and g must have the same length."))

    groups = sort(unique(g))
    n_out = length(groups) > 1 ? length(groups) + 1 : 1
    adt_p = zeros(n_out)
    ks_p = zeros(n_out)

    ref = Distributions.Normal(0, 1)

    if length(groups) > 1
        # per-group residual normality tests
        for (i, grp) in enumerate(groups)
            x_grp = x[g .== grp]
            !(length(x_grp) >= 3) && throw(ArgumentError("Group $grp must have at least 3 observations."))
            res = x_grp .- mean(x_grp)
            # use OneSampleADTest against the theoretical distribution
            adt_p[i] = pvalue(OneSampleADTest(res, ref))
            ks_p[i]  = pvalue(ExactOneSampleKSTest(res, ref))
        end
    end

    # whole-sample residual normality test
    res_all = x .- mean(x)
    adt_p[end] = pvalue(OneSampleADTest(res_all, ref))
    ks_p[end] = pvalue(ExactOneSampleKSTest(res_all, ref))

    return (adt_p=adt_p, ks_p=ks_p)

end
