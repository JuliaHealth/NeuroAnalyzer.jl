export rfz
export r1r2_zscore
export cor_test

"""
    rfz(r)

Compute Fisher's Z transformation of a Pearson correlation coefficient.

Applies `z = atanh(r)`. Because `atanh(±1)` is infinite, the boundary values `r = ±1` are clamped to `±atanh(1 − ε)` where `ε = eps()`, giving `z ≈ ±18.368` in practice.

# Arguments

- `r::Float64`: Pearson correlation coefficient; must be in `[−1, 1]`

# Returns

- `z::Float64`: Fisher z-transformed value

# Throws

- `ArgumentError`: if `r ∉ [−1, 1]`

# Notes

- Boundary clamping: `rfz(1.0) ≈ 18.368`, `rfz(-1.0) ≈ -18.368`.
- For values strictly inside `(−1, 1)` the result equals `atanh(r)` exactly.

# See also

[`r1r2_zscore`](@ref), [`cor_test`](@ref)
"""
function rfz(r::Float64)::Float64

    _in(r, (-1.0, 1.0), "r")
    # clamp boundary values to avoid ±Inf from atanh(±1)
    r == 1.0  && return atanh(1.0 - eps())
    r == -1.0 && return atanh(-1.0 + eps())

    return atanh(r)

end

"""
    r1r2_zscore(; <keyword arguments>)

Calculate the Z-score for the difference between two independent Pearson correlation coefficients using Fisher's Z transformation.

The test statistic is: `z = (atanh(r1) − atanh(r2)) / √(1/(n1 − 3) + 1/(n2 − 3))`

# Arguments

- `r1::Float64`: correlation coefficient for group 1; must be in `(−1, 1)`
- `r2::Float64`: correlation coefficient for group 2; must be in `(−1, 1)`
- `n1::Int64`: number of observations in group 1; must be > 3
- `n2::Int64`: number of observations in group 2; must be > 3

# Returns

- `Float64`: z-score for the difference `r1 − r2`

# Throws

- `ArgumentError`: if `r1` or `r2` are outside `(−1, 1)`, or `n1`/`n2` ≤ 3

# Notes

Both samples must have `n > 3` for the Fisher Z standard error `1/√(n − 3)` to be defined. The original guards (`n > 0`) were insufficient.

# See also

[`rfz`](@ref), [`cor_test`](@ref)
"""
function r1r2_zscore(; r1::Float64, r2::Float64, n1::Int64, n2::Int64)::Float64

    _in(r1, (-1.0, 1.0), "r1")
    _in(r2, (-1.0, 1.0), "r2")
    @assert n1 > 3 "n1 must be > 3 (required for Fisher's Z standard error)."  # was: > 0
    @assert n2 > 3 "n2 must be > 3 (required for Fisher's Z standard error)."  # was: > 0

    # Fisher z-transform both coefficients, then standardize the difference
    z = (rfz(r1) - rfz(r2)) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))

    return z

end

"""
    cor_test(seg1, seg2)

Calculate the Pearson correlation between two vectors and return associated test statistics.

# Arguments

- `s1::AbstractVector`: first signal vector; must have the same length as `s2` and `length > 3` (required for the Fisher Z confidence interval)
- `s2::AbstractVector`: second signal vector

# Returns

Named tuple:

- `t::CorrelationTest{Float64}`: full HypothesisTests.jl result object
- `r::Float64`: Pearson correlation coefficient
- `rc::Tuple{Float64, Float64}`: 95 % confidence interval for `r`
- `ts::Tuple{Float64, String}`: t-statistic and label `"t"`
- `df::Int64`: degrees of freedom (`n1 + n2 − 2`)
- `p::Float64`: two-tailed p-value (clamped to `eps()` if below machine epsilon)

# Throws

- `ArgumentError`: if `s1` and `s2` have different lengths or `length < 4`

# See also

[`rfz`](@ref), [`r1r2_zscore`](@ref)
"""
function cor_test(
    s1::AbstractVector,
    s2::AbstractVector
)::@NamedTuple{
    t::CorrelationTest{Float64},
    r::Float64,
    rc::Tuple{Float64, Float64},
    ts::Tuple{Float64, String},
    df::Int64,
    p::Float64,
}

    @assert length(s1) == length(s2) "s1 and s2 must have the same length."
    @assert length(s1) > 3 "length(s1) must be > 3 (required for confidence interval)."

    t  = CorrelationTest(s1, s2)
    p  = pvalue(t)
    p < eps() && (p = eps())
    df = length(s1) + length(s2) - 2

    return (t=t, r=t.r, rc=confint(t), ts=(t.t, "t"), df=df, p=p)

end
