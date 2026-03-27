export binom_prob
export binom_test

"""
    binom_prob(p, r, n)

Calculate the probability of exactly `r` successes in `n` independent trials.

Uses the binomial probability mass function: `P(X = r) = C(n, r) × pʳ × (1 − p)^(n − r)`

# Arguments

- `p::Float64`: success probability per trial; must be in `[0, 1]`
- `r::Int64`: number of successes; must satisfy `0 ≤ r ≤ n`
- `n::Int64`: number of successes; must satisfy `0 ≤ r ≤ n`

# Returns

- `Float64`: probability of exactly `r` successes
"""
function binom_prob(p::Float64, r::Int64, n::Int64)::Float64

    # validate
    _in(p, (0.0, 1.0), "p")
    n >= 1 || throw(ArgumentError("n must be ≥ 1."))
    r >= 0 || throw(ArgumentError("r must be ≥ 0."))
    r <= n || throw(ArgumentError("r must be ≤ n."))

    return binomial(n, r) * (p^r) * (1 - p)^(n - r)
end

"""
    binom_test(p, n; verbose)

Test whether a proportion differs significantly from 0.5 (two-sided binomial test).

# Arguments

- `prop::Float64`: proportion of level-1 observations; must be in `[0, 1]`
- `n::Int64`: total number of observations; must be ≥ 1
- `verbose::Bool=true`: if `true`, print counts, proportions, 95 % CIs, and p-value

# Returns

Named tuple:

- `x0::Int64`: level-0 count
- `x1::Int64`: level-1 count
- `p0::Float64`: level-0 proportion
- `p1::Float64`: level-1 proportion
- `ci0::Tuple{Float64, Float64}`: level-0 proportion 95 % CI
- `ci1::Tuple{Float64, Float64}`: level-1 proportion 95 % CI
- `p::Float64`: Two-sided binomial test p-value (clamped to `eps()` if below machine epsilon)
"""
function binom_test(
    prop::Float64,
    n::Int64;
    verbose::Bool = true,
)::@NamedTuple{
    x0::Int64,
    x1::Int64,
    p0::Float64,
    p1::Float64,
    ci0::Tuple{Float64, Float64},
    ci1::Tuple{Float64, Float64},
    p::Float64,
}

    # validate
    _in(prop, (0.0, 1.0), "prop")
    n >= 1 || throw(ArgumentError("n must be ≥ 1."))

    x1 = round(Int64, prop * n)
    x0 = n - x1

    ci0 = cip(1 - prop, n)
    ci1 = cip(prop, n)

    if verbose
        println(
            "Level 0: counts: $x0\t proportion: $(round(1 - prop, digits = 3))\t" *
            " 95%CI: $(round(ci0[1], digits = 3)), $(round(ci0[2], digits = 3))",
        )
        println(
            "Level 1: counts: $x1\t proportion: $(round(prop, digits = 3))\t" *
            " 95%CI: $(round(ci1[1], digits = 3)), $(round(ci1[2], digits = 3))",
        )
    end

    p = pvalue(BinomialTest(x1, n, 0.5))

    if verbose
        pv_str = round(p; digits = 3) == 0.0 ? "<0.001" : string(round(p; digits = 3))
        println("Binomial test p value: $pv_str")
        println("Note: level-1 proportion is tested against 0.5")
    end

    # clamp p-value to machine epsilon to avoid exact zero in downstream log transforms
    p < eps() && (p = eps())

    p0 = 1 - prop
    p1 = prop

    return (; x0, x1, p0, p1, ci0, ci1, p)
end

"""
    binom_test(x; verbose)

Test whether the proportion of `true` values in a boolean vector differs from 0.5.

Delegates to [`binom_test(::Float64, ::Int64)`](@ref) after computing the proportion of `true` observations.

# Arguments

- `x::Vector{Bool}`: binary observation vector; must not be empty
- `verbose::Bool=true`: if `true`, print counts, proportions, 95 % CIs, and p-value

# Returns

Named tuple:

- `x0::Int64`: level-0 count
- `x1::Int64`: level-1 count
- `p0::Float64`: level-0 proportion
- `p1::Float64`: level-1 proportion
- `ci0::Tuple{Float64, Float64}`: level-0 proportion 95 % CI
- `ci1::Tuple{Float64, Float64}`: level-1 proportion 95 % CI
- `p::Float64`: Two-sided binomial test p-value (clamped to `eps()` if below machine epsilon)
"""
function binom_test(
    x::Vector{Bool};
    verbose::Bool = true,
)::@NamedTuple{
    x0::Int64,
    x1::Int64,
    p0::Float64,
    p1::Float64,
    ci0::Tuple{Float64, Float64},
    ci1::Tuple{Float64, Float64},
    p::Float64,
}

    # validate
    n = length(x)
    n > 0 || throw(ArgumentError("x must not be empty."))

    n1 = count(identity, x)

    return binom_test(
        n1 / n,
        n;
        verbose = verbose,
    )
end

"""
    binom_test(x, n; verbose)

Test whether the proportion `x / n` differs from 0.5 (two-sided binomial test).

Delegates to [`binom_test(::Float64, ::Int64)`](@ref).

# Arguments

- `x::Int64`: number of level 1 observations; must satisfy `0 ≤ x ≤ n`
- `n::Int64`: number of observations; must be ≥ 1
- `verbose::Bool=true`: if `true`, print counts, proportions, 95 % CIs, and p-value

# Returns

Named tuple:

- `x0::Int64`: level-0 count
- `x1::Int64`: level-1 count
- `p0::Float64`: level-0 proportion
- `p1::Float64`: level-1 proportion
- `ci0::Tuple{Float64, Float64}`: level-0 proportion 95 % CI
- `ci1::Tuple{Float64, Float64}`: level-1 proportion 95 % CI
- `p::Float64`: Two-sided binomial test p-value (clamped to `eps()` if below machine epsilon)
"""
function binom_test(
    x::Int64,
    n::Int64;
    verbose::Bool = true,
)::@NamedTuple{
    x0::Int64,
    x1::Int64,
    p0::Float64,
    p1::Float64,
    ci0::Tuple{Float64, Float64},
    ci1::Tuple{Float64, Float64},
    p::Float64,
}

    # validate
    n >= 1 || throw(ArgumentError("n must be ≥ 1."))
    _in(x, (0, n), 'x')

    return binom_test(
        x / n,
        n;
        verbose = verbose,
    )
end
