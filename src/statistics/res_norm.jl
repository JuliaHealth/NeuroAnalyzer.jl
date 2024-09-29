export res_norm

"""
    res_norm(x, g)

Test normal distribution of residuals.

# Arguments

- `x::AbstractVector`: data values
- `g::Vector{Int64}`: group(s) to which each data value belongs

# Returns

Named tuple containing:
- `adt_p::Vector{Float64}`: p-values for k-sample Anderson–Darling test vs normal distribution
- `ks_p::Vector{Float64}`: p-values for one-sample exact Kolmogorov–Smirnov test vs normal distribution

# Notes

p-values are reported for each group and for the whole sample. If there is only one group, p-values are returned only for the whole sample p-values are reported.
"""
function res_norm(x::AbstractVector, g::Vector{Int64}=repeat([1], length(x)))::NamedTuple{adt_p::Vector{Float64}, ks_p::Vector{Float64}}

    groups = sort(unique(g))

    if length(groups) > 1
        adt_p = zeros(length(groups) + 1)
        ks_p = zeros(length(groups) + 1)
    else
        adt_p = zeros(1)
        ks_p = zeros(1)
    end

    if length(groups) > 1
        # check residuals normality per groups
        for group_idx in eachindex(groups)
            m = mean(x[g .== groups[group_idx]])
            res = x[g .== groups[group_idx]] .- m
            adt_p[group_idx] = pvalue(KSampleADTest(res, rand(Distributions.Normal(0, 1), length(res))))
            ks_p[group_idx] = pvalue(ExactOneSampleKSTest(res, Distributions.Normal(0, 1)))
        end
    end

    # check residuals normality for the whole sample
    m = mean(x)
    res = x .- m
    adt_p[end] = pvalue(KSampleADTest(res, rand(Distributions.Normal(0, 1), length(res))))
    ks_p[end] = pvalue(ExactOneSampleKSTest(res, Distributions.Normal(0, 1)))

    return (adt_p=adt_p, ks_p=ks_p)

end
