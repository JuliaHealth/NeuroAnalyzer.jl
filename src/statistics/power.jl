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

"""
    size_c2g(; m1, s1, m2, r, alpha, power)

Calculate required sample size for a continuous variable (group 1 vs group 2).

# Arguments

- `m1::Real`: group 1 mean
- `s1::Real`: group 1 standard deviation
- `m2::Real`: group 2 mean (expected)
- `r::Int64=1`: enrollment ratio -- the ratio of group 2 to group 1 enrollment
- `alpha::Float64=0.05`: the probability of type I error
- `power::Float64=0.8`: the ability to detect a difference between groups (power = 1 - beta, where beta is the probability of type II error)

# Returns

Named tuple containing:
- `n1::Int64`: group 1 sample size
- `n2::Int64`: group 2 sample size
"""
function size_c2g(; m1::Real, s1::Real, m2::Real, r::Int64=1, alpha::Float64=0.05, power::Float64=0.8)

    beta = 1 - power
    delta = abs(m2 - m1)

    n1 = round(Int64, ((s1^2 + s1^2 / r) * (ci2z(1 - alpha / 2) + ci2z(1 - beta))^2) / delta^2)
    n2 = n1 * r

    return (n1=n1, n2=n2)

end

"""
    size_c1g(; m, s, xbar, alpha, power, iter)

Calculate required sample size for a continuous variable (group 1 vs population).

# Arguments

- `m::Real`: population mean
- `s::Real`: population standard deviation
- `xbar::Real`: expected group mean
- `alpha::Float64=0.05`: the probability of type I error
- `power::Float64=0.8`: the ability to detect a difference between groups (power = 1 - beta, where beta is the probability of type II error)
- `iter::Bool=false`: use iterative method

# Returns

- `n::Int64`: group sample size
"""
function size_c1g(; m::Real, s::Real, xbar::Real, alpha::Float64=0.05, power::Float64=0.8, iter::Bool=false)

    if iter
        n = zeros(length(2:10_000))
        n = [power_c1g(m=m, s=s, xbar=xbar, n=idx, alpha=alpha) for idx in 2:10_000]
        n = vsearch(power, n)
    else
        beta = 1 - power
        n = round(Int64, (s^2 * (ci2z(1 - beta) + ci2z(1 - alpha / 2))^2) / ((m - xbar)^2))
    end

    return n

end

"""
    size_p2g(; p1, p2, r, alpha, power)

Calculate required sample size for a proportion (group 1 vs group 2).

# Arguments

- `p1::Float64`: group 1 anticipated incidence
- `p2::Float64`: group 2 anticipated incidence
- `r::Int64=1`: enrollment ratio -- the ratio of group 2 to group 1 enrollment
- `alpha::Float64=0.05`: the probability of type I error
- `power::Float64=0.8`: the ability to detect a difference between groups (power = 1 - beta, where beta is the probability of type II error)

# Returns

Named tuple containing:
- `n1::Int64`: group 1 sample size
- `n2::Int64`: group 2 sample size
"""
function size_p2g(; p1::Float64, p2::Float64, r::Int64=1, alpha::Float64=0.05, power::Float64=0.8)

    beta = 1 - power
    delta = abs(p2 - p1)
    q1 = 1 - p1
    q2 = 1 - p2
    p_dash = (p1 + r * p2) / (1 + r)
    q_dash = 1 - p_dash

    n1 = round(Int64, ((ci2z(1 - alpha / 2) * sqrt(p_dash * q_dash * (1 + 1 / r)) + ci2z(1 - beta) * sqrt(p1 * q1 + ((p2 * q2) / r))))^2 / delta^2)
    n2 = n1 * r

    return (n1=n1, n2=n2)

end

"""
    size_p1g(; p1, p2, alpha, power)

Calculate required sample size for a proportion (group 1 vs population).

# Arguments

- `p1::Float64`: population incidence
- `p2::Float64`: group anticipated incidence
- `alpha::Float64=0.05`: the probability of type I error
- `power::Float64=0.8`: the ability to detect a difference between groups (power = 1 - beta, where beta is the probability of type II error)

# Returns

- `n::Int64`: group 1 sample size
"""
function size_p1g(; p1::Float64, p2::Float64, alpha::Float64=0.05, power::Float64=0.8)

    beta = 1 - power
    q0 = 1 - p1
    q1 = 1 - p2

    n = round(Int64, (p1 * q0 * (ci2z(1 - alpha / 2) + ci2z(1 -beta) * sqrt((p2 * q1) / (p1 * q0)))^2) / (p2 - p1)^2)

    return n

end

"""
    power_c2g(; m1, s1, n1, m2, s2, n2, alpha)

Calculate study power for a continuous variable (group 1 vs group 2).

# Arguments

- `m1::Real`: group 1 mean
- `s1::Real`: group 1 standard deviation
- `n1::Int64`: group 1 sample size
- `m2::Real`: group 2 mean
- `s2::Real`: group 2 standard deviation
- `n2::Int64`: group 2 sample size
- `alpha::Float64=0.05`: the probability of type I error

# Returns

- `p::Float64`: study power
"""
function power_c2g(; m1::Real, s1::Real, n1::Int64, m2::Real, s2::Real, n2::Int64, alpha::Float64=0.05)

    delta = abs(m2 - m1)

    z = -ci2z(1 - alpha / 2) + (delta / sqrt((s1^2 / n1) + (s2^2 / n2)))
    p = z2p(z)

    return p

end


"""
    power_c1g(; m, s, xbar, n, alpha)

Calculate study power for a continuous variable (group 1 vs population).

# Arguments

- `m::Real`: population mean
- `s::Real`: population standard deviation
- `xbar::Real`: group mean
- `n::Int64`: group sample size
- `alpha::Float64=0.05`: the probability of type I error

# Returns

- `p::Float64`: study power
"""
function power_c1g(; m::Real, s::Real, xbar::Real, n::Int64, alpha::Float64=0.05)

    # delta = abs(m - xbar)
    # z = -ci2z(1 - alpha / 2) + (delta * sqrt(n1) / s)
    # return z2p(z)

    t_r = crit_t(n - 1, alpha, twosided=true)
    t_l = -t_r
    t = (xbar - m) / (s / sqrt(n))
    power_l = cdf(TDist(n - 1), t_l + t)
    power_r = 1 - cdf(TDist(n - 1), t_r + t)
    p = power_l + power_r

    return p

end

"""
    power_p2g(; p1, p2, n1, n2, alpha)

Calculate required sample size for a proportion (group 1 vs group 2).

# Arguments

- `p1::Float64`: group 1 incidence
- `p2::Float64`: group 2 incidence
- `n1::Int64`: group 1 sample size
- `n2::Int64`: group 2 sample size
- `alpha::Float64=0.05`: the probability of type I error

# Returns

- `p::Float64`: study power
"""
function power_p2g(; p1::Float64, p2::Float64, n1::Int64, n2::Int64, alpha::Float64=0.05)

    delta = abs(p2 - p1)
    k = n2 / n1
    q1 = 1 - p1
    q2 = 1 - p2
    r = n2 / n1
    p_dash = (p1 + r * p2) / (1 + r)
    q_dash = 1 - p_dash
    z = (delta / (sqrt(((p1 * q1) / n1) + (p2 * q2) / n2))) - ci2z(1 - alpha / 2) * ((sqrt(p_dash * q_dash * ((1 / n1) + (1 / n2)))) / (sqrt(((p1 * q1)  / n1) + ((p2 * q2) / n2))))
    p = z2p(z)

    return p

end

"""
    power_p1g(; p1, p2, n1, alpha)

Calculate required sample size for a proportion (group 1 vs population).

# Arguments

- `p1::Float64`: population incidence
- `p2::Float64`: group 1 anticipated incidence
- `n1::Int64`: group 1 sample size
- `alpha::Float64=0.05`: the probability of type I error

# Returns

- `p::Float64`: study power
"""
function power_p1g(; p1::Float64, p2::Float64, n1::Int64, alpha::Float64=0.05)

    q0 = 1 - p1
    q1 = 1 - p2
    z = (sqrt(n1 * ((p2 - p1)^2 / (p1 * q0))) - ci2z(1 - alpha / 2)) / (sqrt((p2 * q1)/(p1 * q0)))
    p = z2p(z)

    return p

end

"""
    size_c1diff(; s1, s2, twosided, power)

Calculate required sample size for detecting a difference in a continuous variable (group 1 vs population).

# Arguments

- `s1::Real`: population standard deviation
- `s2::Real`: study standard deviation that we want to detect
- `twosided::Bool=true`: if true, the estimation is for two-sided difference
- `power::Float64=0.8`: the ability to detect a difference between groups (power = 1 - beta, where beta is the probability of type II error)

# Returns

- `n::Int64`: study sample size
"""
function size_c1diff(; s1::Real, s2::Real, twosided::Bool=true, power::Float64=0.8)

    sdiff = s2 / s1

    sdiff_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5]
    power_values = [0.99, 0.95, 0.9, 0.8]

    !(sdiff in sdiff_values) && _warn("Sample size will be estimated.")
    !(power in power_values) && _warn("Sample size will be estimated.")

    if sdiff < 0.1
        sdiff = 0.1
        _warn("Sample size will be estimated.")
    elseif sdiff > 1.5
        sdiff = 1.5
        _warn("Sample size will be estimated.")
    end

    if power < 0.8
        power = 0.8
        _warn("Sample size will be estimated.")
    elseif power > 0.99
        power = 0.99
        _warn("Sample size will be estimated.")
    end

    sdiff_idx = vsearch(sdiff, sdiff_values)
    power_idx = vsearch(power, power_values)

    table = [3676 2600 2103 1571;
             920 651 527 394;
             410 290 235 176;
             231 164 133 100;
             148 105 86 64;
             104 74 60 45;
             76 54 44 33;
             59 42 34 26;
             47 34 27 21;
             38 27 22 17;
             32 23 19 14;
             27 20 16 12;
             23 17 14 11;
             20 15 12 9;
             18 13 11 8]

    n = table[sdiff_idx, power_idx]

    return twosided ? 2 * n : n

end

"""
    size_p1diff(; p1, p2, power)

Calculate required sample size for detecting a difference in a proportion (group 1 vs population).

# Arguments

- `p1::Real`: population proportion
- `p2::Real`: study proportion that we want to detect
- `power::Float64=0.8`: the ability to detect a difference between groups (power = 1 - beta, where beta is the probability of type II error)

# Returns

- `n::Int64`: study sample size (for both study groups)
"""
function size_p1diff(; p1::Real, p2::Real, power::Float64=0.8)

    p = (p1 + p2) / 2
    sdiff = round((p1 - p2) / sqrt(p * (1 - p)), digits=1)

    sdiff_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5]
    power_values = [0.99, 0.95, 0.9, 0.8]

    !(sdiff in sdiff_values) && _warn("Sample size will be estimated.")
    !(power in power_values) && _warn("Sample size will be estimated.")

    if sdiff < 0.1
        sdiff = 0.1
        _warn("Sample size will be estimated.")
    elseif sdiff > 1.5
        sdiff = 1.5
        _warn("Sample size will be estimated.")
    end

    if power < 0.8
        power = 0.8
        _warn("Sample size will be estimated.")
    elseif power > 0.99
        power = 0.99
        _warn("Sample size will be estimated.")
    end

    sdiff_idx = vsearch(sdiff, sdiff_values)
    power_idx = vsearch(sdiff, sdiff_values)

    table = [3676 2600 2103 1571;
             920 651 527 394;
             410 290 235 176;
             231 164 133 100;
             148 105 86 64;
             104 74 60 45;
             76 54 44 33;
             59 42 34 26;
             47 34 27 21;
             38 27 22 17;
             32 23 19 14;
             27 20 16 12;
             23 17 14 11;
             20 15 12 9;
             18 13 11 8]

    n = 2 * table[sdiff_idx, power_idx]

    return n

end