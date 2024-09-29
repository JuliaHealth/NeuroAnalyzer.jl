export cor_test

"""
    cor_test(seg1, seg2)

Calculate correlation between two vectors.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`

# Returns

Named tuple containing:
- `t::CorrelationTest{Float64}`
- `r::Float64`: correlation coefficient
- `rc::Tuple{Float64, Float64}`: correlation coefficient confidence interval
- `ts::Tuple{Float64, String}`: t-statistics
- `df::Int64`: degrees of freedom
- `p::Float64`: p-value
"""
function cor_test(s1::AbstractVector, s2::AbstractVector)::NamedTuple{t::CorrelationTest{Float64}, r::Float64, rc::Tuple{Float64, Float64}, ts::Tuple{Float64, String}, df::Int64, p::Float64}

    @assert length(s1) == length(s2) "Both vectors must have the same length."
    t = CorrelationTest(s1, s2)
    p = pvalue(t)
    p < eps() && (p = eps())
    df = length(s1) + length(s2) - 2

    return (t=t, r=t.r, rc=confint(t), ts=(t.t, "t"), df=df, p=p)

end
