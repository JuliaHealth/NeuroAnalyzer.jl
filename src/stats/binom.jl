export binom_prob
export binom_test

"""
    binom_prob(p, r, n)

Calculate probability of exactly `r` successes in `n` trials.

# Arguments

- `p::Float64`: proportion of successes
- `r::Int64`: number of successes
- `n::Int64`: number of trials

# Returns

- `bp::Float64`: probability
"""
function binom_prob(p::Float64, r::Int64, n::Int64)::Float64

    _in(p, (0.0, 1.0), "p")
    @assert r > 0 "r must be > 0.0."
    @assert n > 0 "n must be > 0."

    bp = binomial(n, r) * (p^r) * (1 - p)^(n - r)

    return bp

end

"""
    binom_test(p, n; verbose)

Test proportion against 0.5.

# Arguments

- `p::Float64`: proportion of level 1 observations
- `n::Int64`: number of observations
- `verbose::Bool=true`: print detailed output

# Returns

Named tuple containing:
- `x0::Int64`: level 0 counts
- `x1::Int64`: level 1 counts
- `p0::Float64`: level 0 proportion
- `p1::Float64`: level 1 proportion
- `ci0::Tuple{Float64, Float64}`: level 1 proportion 95%CI
- `ci1::Tuple{Float64, Float64}`: level 1 proportion 95%CI
- `p::Float64`: Binomial test p value
"""
function binom_test(p::Float64, n::Int64; verbose::Bool=true)::@NamedTuple{x0::Int64, x1::Int64, p0::Float64, p1::Float64, ci0::Tuple{Float64, Float64}, ci1::Tuple{Float64, Float64}, p::Float64}

    _in(p, (0.0, 1.0), "p")
    @assert n > 0 "n must be > 0."
    x = round(Int64, p * n)
    ci0 = cip(1 - p, n)
    verbose && println("Level 0: counts: $(n - x)\t proportion: $(round(1 - p, digits=3))\t 95%CI: $(round.(ci0[1], digits=3)), $(round.(ci0[2], digits=3))")
    ci1 = cip(p, n)
    verbose && println("Level 1: counts: $(x)\t proportion: $(round(p, digits=3))\t 95%CI: $(round.(ci1[1], digits=3)), $(round.(ci1[2], digits=3))")

    pv = pvalue(BinomialTest(x, n, 0.5))
    if verbose
        if round(pv, digits=3) == 0.0
            println("Binomial test p value: <0.001")
        else
            println("Binomial test p value: $(round(pv, digits=3))")
        end
        println("Note: level 1 proportion is tested against 0.5")
    end

    pv < eps() && (pv = eps())

    return (x0=(n - x), x1=x, p0=(1 - p), p1=p, ci0=ci0, ci1=ci1, p=pv)

end

"""
    binom_test(x; verbose)

Test proportion against 0.5.

# Arguments

- `x::Vector{Bool}`: vector of level 0 and 1 categories
- `verbose::Bool=true`: print detailed output

# Returns

Named tuple containing:
- `x0::Int64`: level 0 counts
- `x1::Int64`: level 1 counts
- `p0::Float64`: level 0 proportion
- `p1::Float64`: level 1 proportion
- `ci0::Tuple{Float64, Float64}`: level 1 proportion 95%CI
- `ci1::Tuple{Float64, Float64}`: level 1 proportion 95%CI
- `p::Float64`: Binomial test p value
"""
function binom_test(x::Vector{Bool}; verbose::Bool=true)::@NamedTuple{x0::Int64, x1::Int64, p0::Float64, p1::Float64, ci0::Tuple{Float64, Float64}, ci1::Tuple{Float64, Float64}, p::Float64}

    n = length(x)
    @assert n > 0 "Number of observations must be > 0."
    x = length(findall(isequal(true), x))
    p = x / n

    return binom_test(p, n, verbose=verbose)

end

"""
    binom_test(x, n; verbose)

Test proportion against 0.5.

# Arguments

- `x::Int64`: number of level 1 observations
- `n::Int64`: number of observations
- `verbose::Bool=true`: print detailed output

# Returns

Named tuple containing:
- `x0::Int64`: level 0 counts
- `x1::Int64`: level 1 counts
- `p0::Float64`: level 0 proportion
- `p1::Float64`: level 1 proportion
- `ci0::Tuple{Float64, Float64}`: level 1 proportion 95%CI
- `ci1::Tuple{Float64, Float64}`: level 1 proportion 95%CI
- `p::Float64`: Binomial test p value
"""
function binom_test(x::Int64, n::Int64; verbose::Bool=true)::@NamedTuple{x0::Int64, x1::Int64, p0::Float64, p1::Float64, ci0::Tuple{Float64, Float64}, ci1::Tuple{Float64, Float64}, p::Float64}

    _in(x, (0, n), "x")

    return binom_test(x/n, n, verbose=verbose)

end
