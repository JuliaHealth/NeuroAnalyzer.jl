export rfz
export r1r2_zscore
export cor_test

"""
    rfz(r)

Calculate Fischer's z transformations of correlation coefficient.

# Arguments

- `r::Float64`: correlation coefficient

# Returns

- `z::Float64`: z score

# Notes

Formula: `z = atanh(r)`, gives ∞ for r = 1.0 and -∞ for r = -1.0. This functions returns 18.36840028483855 (atanh(1 - eps())) for r = 1.0 and -18.36840028483855 (atanh(1 - eps())) for r = -1.0.
"""
function rfz(r::Float64)::Float64

    _in(r, (-1.0, 1.0), "r")

    if r == 1.0
        z = atanh(r - eps())
    elseif r == -1.0
        z = atanh(r + eps())
    else
        z = atanh(r)
    end

    return z

end

"""
    r1r2_zscore(; <keyword arguments>)

Calculate z score for the difference of two correlation coefficients.

# Arguments

- `r1::Float64`: correlation coefficient, group 1
- `r2::Float64`: correlation coefficient, group 2
- `n1::Int64`: number of observations, group 1
- `n2::Int64`: number of observations, group 2

# Returns

- `z::Float64`: z score
"""
function r1r2_zscore(; r1::Float64, r2::Float64, n1::Int64, n2::Int64)::Float64

    _in(r1, (-1.0, 1.0), "r1")
    _in(r2, (-1.0, 1.0), "r2")
    @assert n1 > 0 "n1 must be > 0."
    @assert n2 > 0 "n2 must be > 0."

    z1 = rfz(r1)
    z2 = rfz(r2)

    z = (z1 - z2) / sqrt((1 / (n1 - 3)) + (1 / (n2 - 3)))

    return z

end

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
- `p::Float64`: p value
"""
function cor_test(s1::AbstractVector, s2::AbstractVector)::@NamedTuple{t::CorrelationTest{Float64}, r::Float64, rc::Tuple{Float64, Float64}, ts::Tuple{Float64, String}, df::Int64, p::Float64}

    @assert length(s1) == length(s2) "Lengths of s1 and s2 must be equal."

    t = CorrelationTest(s1, s2)
    p = pvalue(t)
    p < eps() && (p = eps())
    df = length(s1) + length(s2) - 2

    return (t=t, r=t.r, rc=confint(t), ts=(t.t, "t"), df=df, p=p)

end
