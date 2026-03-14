export size_c2g
export size_c1g
export size_p2g
export size_p1g
export power_c2g
export power_c1g
export power_p2g
export power_p1g
export size_c1diff
export size_p1diff
export mde
export size_p
export size_m

"""
    size_c2g(; <keyword arguments>)

Calculate the required sample size for a two-group continuous outcome study.

Assumes equal group standard deviations (both equal to `s1`).

# Arguments
- `m1::Real`: group 1 mean
- `s1::Real`: group 1 standard deviation (assumed equal for both groups); must be > 0
- `m2::Real`: group 2 mean (expected)
- `r::Int64=1`: enrollment ratio (group 2 / group 1); must be ≥ 1
- `alpha::Float64=0.05`: type I error probability; must be in `(0, 1)`
- `power::Float64=0.8`: desired power; must be in `(0, 1)`

# Returns

Named tuple:

- `n1::Int64`: group 1 sample size
- `n2::Int64`: group 2 sample size

# Throws

- `ArgumentError`: if `s1 ≤ 0`, `r < 1`, `m1 == m2`, or `alpha`/`power` out of range

# See also
[`power_c2g`](@ref), [`size_c1g`](@ref)
"""
function size_c2g(;
    m1::Real,
    s1::Real,
    m2::Real,
    r::Int64 = 1,
    alpha::Float64 = 0.05,
    power::Float64 = 0.8
)::@NamedTuple{n1::Int64, n2::Int64}

    _in(alpha, (0, 1.0), "alpha")
    _in(power, (0, 1.0), "power")
    !(s1 > 0) && throw(ArgumentError("s1 must be > 0."))
    !(r  >= 1) && throw(ArgumentError("r must be ≥ 1."))
    !(m1 != m2) && throw(ArgumentError("m1 and m2 must differ (zero effect size)."))

    beta  = 1 - power
    delta = abs(m2 - m1)
    # equal-variance formula: n1 = s1²(1 + 1/r) × (z_α + z_β)² / δ²
    n1 = ceil(Int64, (s1^2 * (1 + 1/r) * (cl2z(1 - alpha) + cl2z(1 - beta))^2) / delta^2)
    n2 = n1 * r

    return (n1=n1, n2=n2)

end

"""
    size_c1g(; <keyword arguments>)

Calculate the required sample size for a one-group continuous outcome study (group vs population).

# Arguments

- `m::Real`: population mean
- `s::Real`: population standard deviation; must be > 0
- `xbar::Real`: expected study group mean; must differ from `m`
- `alpha::Float64=0.05`: type I error probability; must be in `(0, 1)`
- `power::Float64=0.8`: desired power; must be in `(0, 1)`
- `iter::Bool=false`: if `true`, use an iterative power search (n from 2 to 10,000)

# Returns

- `Int64`: required sample size

# Throws

- `ArgumentError`: if `s ≤ 0`, `m == xbar`, or `alpha`/`power` out of range

# See also

[`power_c1g`](@ref), [`size_c2g`](@ref)
"""
function size_c1g(;
    m::Real,
    s::Real,
    xbar::Real,
    alpha::Float64 = 0.05,
    power::Float64 = 0.8,
    iter::Bool = false
)::Int64

    _in(alpha, (0, 1.0), "alpha")
    _in(power, (0, 1.0), "power")
    !(s > 0) && throw(ArgumentError("s must be > 0."))
    !(m != xbar) && throw(ArgumentError("m and xbar must differ."))

    if iter
        powers = [power_c1g(m=m, s=s, xbar=xbar, n=n, alpha=alpha) for n in 2:10_000]
        # +1 because search starts at n=2
        return vsearch(power, powers) + 1
    else
        beta = 1 - power
        return ceil(Int64, (s^2 * (cl2z(1 - beta) + cl2z(1 - alpha))^2) / (m - xbar)^2)
    end

end

"""
    size_p2g(; <keyword arguments>)

Calculate the required sample size for a two-group proportion study.

# Arguments

- `p1::Float64`: group 1 proportion; must be in `(0, 1)`
- `p2::Float64`: group 2 proportion; must be in `(0, 1)` and ≠ `p1`
- `r::Int64=1`: enrollment ratio (group 2 / group 1); must be ≥ 1
- `alpha::Float64=0.05`: type I error probability; must be in `(0, 1)`
- `power::Float64=0.8`: desired power; must be in `(0, 1)`

# Returns

Named tuple:

- `n1::Int64`: group 1 sample size
- `n2::Int64`: group 2 sample size

# Throws

- `ArgumentError`: if proportions are out of range or equal, or `r < 1`

# See also

[`power_p2g`](@ref), [`size_p1g`](@ref)
"""
function size_p2g(;
    p1::Float64,
    p2::Float64,
    r::Int64 = 1,
    alpha::Float64 = 0.05,
    power::Float64 = 0.8
)::@NamedTuple{n1::Int64, n2::Int64}

    _in(alpha, (0, 1.0), "alpha")
    _in(power, (0, 1.0), "power")
    _in(p1, (0.0, 1.0), "p1")
    _in(p2, (0.0, 1.0), "p2")
    !(p1 != p2) && throw(ArgumentError("p1 and p2 must differ."))
    !(r >= 1) && throw(ArgumentError("r must be ≥ 1."))

    beta = 1 - power
    delta = abs(p2 - p1)
    q1 = 1 - p1
    q2 = 1 - p2
    p_dash = (p1 + r * p2) / (1 + r)
    q_dash = 1 - p_dash

    n1 = ceil(
        Int64,
        (cl2z(1 - alpha) * sqrt(p_dash * q_dash * (1 + 1/r)) +
         cl2z(1 - beta)  * sqrt(p1 * q1 + p2 * q2 / r))^2 / delta^2,
    )
    n2 = n1 * r

    return (n1=n1, n2=n2)

end

"""
    size_p1g(; <keyword arguments>)

Calculate the required sample size for a one-group proportion study (group vs population).

# Arguments

- `p1::Float64`: population proportion; must be in `(0, 1)`
- `p2::Float64`: study group proportion; must be in `(0, 1)` and ≠ `p1`
- `alpha::Float64=0.05`: type I error probability; must be in `(0, 1)`
- `power::Float64=0.8`: desired power; must be in `(0, 1)`

# Returns

- `Int64`: required sample size

# Throws

- `ArgumentError`: if proportions are out of range or equal

# See also

[`power_p1g`](@ref), [`size_p2g`](@ref)
"""
function size_p1g(;
    p1::Float64,
    p2::Float64,
    alpha::Float64 = 0.05,
    power::Float64 = 0.8
)::Int64

    _in(alpha, (0, 1.0), "alpha")
    _in(power, (0, 1.0), "power")
    _in(p1, (0.0, 1.0), "p1")
    _in(p2, (0.0, 1.0), "p2")
    !(p1 != p2) && throw(ArgumentError("p1 and p2 must differ."))

    beta = 1 - power
    q0 = 1 - p1
    q1 = 1 - p2

    return ceil(Int64,
        (p1 * q0 * (cl2z(1 - alpha) + cl2z(1 - beta) * sqrt((p2 * q1) / (p1 * q0)))^2) /
        (p2 - p1)^2
    )

end

"""
    power_c2g(; <keyword arguments>)

Calculate study power for a two-group continuous outcome comparison.

# Arguments

- `m1::Real`: group 1 mean
- `s1::Real`: group 1 standard deviation; must be > 0
- `n1::Int64`: group 1 sample size; must be ≥ 1
- `m2::Real`: group 2 mean
- `s2::Real`: group 2 standard deviation; must be > 0
- `n2::Int64`: group 2 sample size; must be ≥ 1
- `alpha::Float64=0.05`: type I error probability; must be in `(0, 1)`

# Returns

- `Float64`: estimated study power

# Throws

- `ArgumentError`: if SDs ≤ 0, sample sizes < 1, or `alpha` out of range

# See also

[`size_c2g`](@ref), [`power_c1g`](@ref)
"""
function power_c2g(;
    m1::Real,
    s1::Real,
    n1::Int64,
    m2::Real,
    s2::Real,
    n2::Int64,
    alpha::Float64 = 0.05
)::Float64

    _in(alpha, (0, 1.0), "alpha")
    !(s1 > 0) && throw(ArgumentError("s1 must be > 0."))
    !(s2 > 0) && throw(ArgumentError("s2 must be > 0."))
    !(n1 >= 1) && throw(ArgumentError("n1 must be ≥ 1."))
    !(n2 >= 1) && throw(ArgumentError("n2 must be ≥ 1."))

    delta = abs(m2 - m1)
    z = -cl2z(1 - alpha) + delta / sqrt(s1^2/n1 + s2^2/n2)

    return z2p(abs(z))

end

"""
    power_c1g(; <keyword arguments>)

Calculate study power for a one-group continuous outcome comparison (group vs population), using the t-distribution.

# Arguments

- `m::Real`: population mean
- `s::Real`: population standard deviation; must be > 0
- `xbar::Real`: study group mean
- `n::Int64`: group sample size; must be ≥ 2
- `alpha::Float64=0.05`: type I error probability; must be in `(0, 1)`

# Returns

- `Float64`: estimated study power

# Throws

- `ArgumentError`: if `s ≤ 0`, `n < 2`, or `alpha` out of range

# See also

[`size_c1g`](@ref), [`power_c2g`](@ref)
"""
function power_c1g(;
    m::Real,
    s::Real,
    xbar::Real,
    n::Int64,
    alpha::Float64=0.05
)::Float64

    _in(alpha, (0, 1.0), "alpha")
    !(s > 0 ) && throw(ArgumentError("s must be > 0."))
    !(n >= 2) && throw(ArgumentError("n must be ≥ 2."))

    t_crit  = crit_t(n - 1, alpha)
    t_stat  = (xbar - m) / (s / sqrt(n))
    # two-tailed power: sum of left and right tail probabilities
    power_l = cdf(TDist(n - 1), -t_crit + t_stat)
    power_r = 1 - cdf(TDist(n - 1), t_crit + t_stat)

    return power_l + power_r

end

"""
    power_p2g(; <keyword arguments>)

Calculate study power for a two-proportion comparison.

# Arguments

- `p1::Float64`: group 1 proportion; must be in `(0, 1)`
- `p2::Float64`: group 2 proportion; must be in `(0, 1)`
- `n1::Int64`: group 1 sample size; must be ≥ 1
- `n2::Int64`: group 2 sample size; must be ≥ 1
- `alpha::Float64=0.05`: type I error probability; must be in `(0, 1)`

# Returns

- `Float64`: estimated study power

# Throws

- `ArgumentError`: if proportions out of range or sample sizes < 1

# See also

[`size_p2g`](@ref), [`power_p1g`](@ref)
"""
function power_p2g(;
    p1::Float64,
    p2::Float64,
    n1::Int64,
    n2::Int64,
    alpha::Float64 = 0.05
)::Float64

    _in(alpha, (0, 1.0), "alpha")
    _in(p1, (0.0, 1.0), "p1")
    _in(p2, (0.0, 1.0), "p2")
    !(n1 >= 1) && throw(ArgumentError("n1 must be ≥ 1."))
    !(n2 >= 1) && throw(ArgumentError("n2 must be ≥ 1."))

    delta = abs(p2 - p1)
    q1 = 1 - p1
    q2 = 1 - p2
    r = n2 / n1
    p_dash = (p1 + r * p2) / (1 + r)
    q_dash = 1 - p_dash

    se_alt = sqrt(p1*q1/n1 + p2*q2/n2)
    se_null = sqrt(p_dash * q_dash * (1/n1 + 1/n2))

    z = delta / se_alt - cl2z(1 - alpha) * (se_null / se_alt)

    return z2p(abs(z))

end

"""
    power_p1g(; <keyword arguments>)

Calculate study power for a one-proportion comparison (group vs population).

# Arguments

- `p1::Float64`: study group proportion; must be in `(0, 1)`
- `p2::Float64`: population proportion; must be in `(0, 1)`
- `n1::Int64`: study group sample size; must be ≥ 1
- `alpha::Float64=0.05`: type I error probability; must be in `(0, 1)`

# Returns

- `Float64`: estimated study power

# Throws

- `ArgumentError`: if proportions out of range or `n1 < 1`

# See also

[`size_p1g`](@ref), [`power_p2g`](@ref)
"""
function power_p1g(;
    p1::Float64,
    p2::Float64,
    n1::Int64,
    alpha::Float64 = 0.05
)::Float64

    _in(alpha, (0, 1.0), "alpha")
    _in(p1, (0.0, 1.0), "p1")
    _in(p2, (0.0, 1.0), "p2")
    !(n1 >= 1) && throw(ArgumentError("n1 must be ≥ 1."))

    q0 = 1 - p2
    q1 = 1 - p1
    z  = (sqrt(n1 * (p1 - p2)^2 / (p2 * q0)) - cl2z(1 - alpha)) /
         sqrt(p1 * q1 / (p2 * q0))

    return z2p(abs(z))

end

"""
    size_c1diff(; <keyword arguments>)

Calculate required sample size for detecting a difference in variance (study SD vs population SD) using a lookup table.

# Arguments

- `s1::Real`: study SD to detect; ratio `s1/s2` must be in `[0.1, 1.5]`
- `s2::Real`: population SD; must be > 0
- `twotailed::Bool=true`: if `true`, return twice the one-sided table value
- `power::Float64=0.8`: desired power; must be in `[0.8, 0.99]`

# Returns

- `Int64`: required sample size (doubled if `twotailed=true`)

# Throws

- `ArgumentError`: if `s2 == 0` or `power ∉ (0, 1)`

# Notes

Values outside the table range are clamped to the nearest boundary and a warning is issued.

# See also

[`size_p1diff`](@ref), [`mde`](@ref)
"""
function size_c1diff(;
    s1::Real,
    s2::Real,
    twotailed::Bool = true,
    power::Float64 = 0.8
)::Int64

    _in(power, (0, 1.0), "power")
    !(s2 != 0) && throw(ArgumentError("s2 must not be zero."))

    sdiff_values = [
        0.1, 0.2, 0.3, 0.4, 0.5,
        0.6, 0.7, 0.8, 0.9, 1.0,
        1.1, 1.2, 1.3, 1.4, 1.5
    ]
    power_values = [0.99, 0.95, 0.9, 0.8]

    sdiff = s1 / s2
    !(sdiff in sdiff_values) && _warn("sdiff=$sdiff not in table; result will be estimated.")
    !(power in power_values) && _warn("power=$power not in table; result will be estimated.")

    sdiff = clamp(sdiff, 0.1, 1.5)
    power = clamp(power, 0.8, 0.99)

    sdiff_idx = vsearch(sdiff, sdiff_values)
    power_idx = vsearch(power, power_values)

    table = [
        3676 2600 2103 1571;
         920  651  527  394;
         410  290  235  176;
         231  164  133  100;
         148  105   86   64;
         104   74   60   45;
          76   54   44   33;
          59   42   34   26;
          47   34   27   21;
          38   27   22   17;
          32   23   19   14;
          27   20   16   12;
          23   17   14   11;
          20   15   12    9;
          18   13   11    8
    ]

    n = table[sdiff_idx, power_idx]
    return twotailed ? 2 * n : n

end

"""
    size_p1diff(; <keyword arguments>)

Calculate required sample size for detecting a difference in proportions (study vs population) using a lookup table.

# Arguments

- `p1::Float64`: study group proportion; must be in `(0, 1)`
- `p2::Float64`: population proportion; must be in `(0, 1)`
- `power::Float64=0.8`: desired power; must be in `[0.8, 0.99]`

# Returns

- `Int64`: total required sample size (both groups combined)

# Throws

- `ArgumentError`: if proportions out of range or `power ∉ (0, 1)`

# Notes

Values outside the table range are clamped and a warning is issued.

# See also

[`size_c1diff`](@ref), [`size_p1g`](@ref)
"""
function size_p1diff(; p1::Float64, p2::Float64, power::Float64 = 0.8)::Int64

    _in(power, (0, 1.0), "power")
    _in(p1, (0.0, 1.0), "p1")
    _in(p2, (0.0, 1.0), "p2")

    p_avg = (p1 + p2) / 2
    sdiff = round((p2 - p1) / sqrt(p_avg * (1 - p_avg)), digits=1)

    sdiff_values = [
        0.1, 0.2, 0.3, 0.4, 0.5,
        0.6, 0.7, 0.8, 0.9, 1.0,
        1.1, 1.2, 1.3, 1.4, 1.5
    ]
    power_values = [0.99, 0.95, 0.9, 0.8]

    !(sdiff in sdiff_values) && _warn("sdiff=$sdiff not in table; result will be estimated.")
    !(power in power_values) && _warn("power=$power not in table; result will be estimated.")

    sdiff = clamp(sdiff, 0.1, 1.5)
    power = clamp(power, 0.8, 0.99)

    sdiff_idx = vsearch(sdiff, sdiff_values)
    power_idx = vsearch(power, power_values)

    table = [
        3676 2600 2103 1571;
         920  651  527  394;
         410  290  235  176;
         231  164  133  100;
         148  105   86   64;
         104   74   60   45;
          76   54   44   33;
          59   42   34   26;
          47   34   27   21;
          38   27   22   17;
          32   23   19   14;
          27   20   16   12;
          23   17   14   11;
          20   15   12    9;
          18   13   11    8
    ]

    return 2 * table[sdiff_idx, power_idx]
end

"""
    mde(; <keyword arguments>)

Calculate the minimum detectable effect (MDE) for a given sample size.

Computed as `MDE = (z_α + z_β)² × s² / n`.

# Arguments

- `n::Int64`: number of subjects per group; must be ≥ 1
- `s::Real`: standard deviation; must be > 0
- `alpha::Float64=0.05`: type I error probability; must be in `(0, 1)`
- `beta::Float64=0.20`: type II error probability; must be in `(0, 1)`
- `verbose::Bool=true`: if `true`, print the critical Z-scores

# Returns

- `Float64`: minimum detectable effect size

# Throws

- `ArgumentError`: if `n < 1`, `s ≤ 0`, or `alpha`/`beta` out of range

# See also

[`size_c2g`](@ref), [`size_c1g`](@ref)
"""
function mde(;
    n::Int64,
    s::Real,
    alpha::Float64=0.05,
    beta::Float64=0.2,
    verbose::Bool=true
)::Float64

    _in(alpha, (0, 1.0), "alpha")
    _in(beta, (0, 1.0), "beta")
    !(n >= 1) && throw(ArgumentError("n must be ≥ 1."))
    !(s >  0) && throw(ArgumentError("s must be > 0."))

    z_alpha = crit_z(alpha)
    z_beta = cl2z(1 - beta)

    verbose && println("z_α = $z_alpha")
    verbose && println("z_β = $z_beta")

    return (z_alpha + z_beta)^2 * s^2 / n

end

"""
    size_p(; <keyword arguments>)

Calculate the required sample size for estimating a proportion within a margin of error `E`.

# Arguments

- `p::Union{Float64, Nothing}=nothing`: expected proportion; if `nothing`, uses the conservative estimate `p = 0.5`
- `alpha::Float64=0.05`: type I error probability; must be in `(0, 1)`
- `E::Float64`: desired margin of error; must be in `(0, 1)`

# Returns

- `Int64`: required sample size

# Throws

- `ArgumentError`: if `E` or `alpha` are out of range

# See also

[`size_m`](@ref), [`size_p1g`](@ref)
"""
function size_p(;
    p::Union{Float64, Nothing} = nothing,
    alpha::Float64 = 0.05,
    E::Float64
)::Int64

    _in(E, (0, 1.0), "E")
    _in(alpha, (0, 1.0), "alpha")

    z = crit_z(alpha / 2; twotailed=false)
    if isnothing(p)
        # conservative: p=0.5 maximizes p(1-p)
        return ceil(Int64, z^2 * 0.25 / E^2)
    else
        _in(p, (0.0, 1.0), "p")
        return ceil(Int64, z^2 * p * (1 - p) / E^2)
    end

end

"""
    size_m(; <keyword arguments>)

Calculate the required sample size for estimating a population mean within a margin of error `E`.

# Arguments

- `sigma::Real`: population standard deviation; must be > 0
- `alpha::Float64=0.05`: type I error probability; must be in `(0, 1)`
- `E::Real`: desired margin of error; must be > 0

# Returns

- `Int64`: required sample size

# Throws

- `ArgumentError`: if `sigma ≤ 0`, `E ≤ 0`, or `alpha` out of range

# See also

[`size_p`](@ref), [`size_c1g`](@ref)
"""
function size_m(;
    sigma::Real,
    alpha::Float64 = 0.05,
    E::Real
)::Int64

    _in(alpha, (0, 1.0), "alpha")
    !(sigma > 0) && throw(ArgumentError("sigma must be > 0."))
    !(E > 0) && throw(ArgumentError("E must be > 0."))

    return ceil(Int64, (crit_z(alpha / 2; twotailed=false) * sigma / E)^2)

end
