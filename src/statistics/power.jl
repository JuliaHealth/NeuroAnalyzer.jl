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
    size_c2g(m1, s1, m2, r, alpha, power)

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

    n1 = round(Int64, ((s1^2 + s1^2 / r) * (crit_z(1 - alpha / 2) + crit_z(1 - beta))^2) / delta^2)
    n2 = n1 * r

    return (n1=n1, n2=n2)

end

"""
    size_c1g(m0, s0, m1, alpha, power)

Calculate required sample size for a continuous variable (group 1 vs population).

# Arguments

- `m0::Real`: population mean
- `s0::Real`: population standard deviation
- `m1::Real`: group 1 mean (expected)
- `alpha::Float64=0.05`: the probability of type I error
- `power::Float64=0.8`: the ability to detect a difference between groups (power = 1 - beta, where beta is the probability of type II error)

# Returns

Named tuple containing:
- `n::Int64`: group 1 sample size
"""
function size_c1g(; m0::Real, s0::Real, m1::Real, alpha::Float64=0.05, power::Float64=0.8)

    beta = 1 - power

    n = round(Int64, (s0^2 * (crit_z(1 - beta) + crit_z(1 - alpha / 2))^2) / ((m0 - m1)^2))
    
    return n

end

"""
    size_p2g(x)

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

    n1 = round(Int64, ((crit_z(1 - alpha / 2) * sqrt(p_dash * q_dash * (1 + 1 / r)) + crit_z(1 - beta) * sqrt(p1 * q1 + ((p2 * q2) / r))))^2 / delta^2)
    n2 = n1 * r

    return (n1=n1, n2=n2)

end

"""
    size_p1g(x)

Calculate required sample size for a proportion (group 1 vs population).

# Arguments

- `p0::Float64`: population incidence
- `p1::Float64`: group anticipated incidence
- `alpha::Float64=0.05`: the probability of type I error
- `power::Float64=0.8`: the ability to detect a difference between groups (power = 1 - beta, where beta is the probability of type II error)

# Returns

- `n::Int64`: group 1 sample size
"""
function size_p1g(; p0::Float64, p1::Float64, r::Int64=1, alpha::Float64=0.05, power::Float64=0.8)

    beta = 1 - power
    q0 = 1 - p0
    q1 = 1 - p1

    n = round(Int64, (p0 * q0 * (crit_z(1 - alpha / 2) + crit_z(1 -beta) * sqrt((p1 * q1) / (p0 * q0)))^2) / (p1 - p0)^2)

    return n

end

"""
    power_c2g(x)

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
function power_c2g(; m1::Real, s1::Real, m2::Real, s2::Real, n1::Int64, n2::Int64, alpha::Float64=0.05)

    delta = abs(m2 - m1)

    z = -crit_z(1 - alpha / 2) + (delta / sqrt((s1^2 / n1) + (s2^2 / n2)))

    return z2pow(z)

end


"""
    power_c1g(x)

Calculate study power for a continuous variable (group 1 vs population).

# Arguments

- `m0::Real`: population mean
- `s0::Real`: population standard deviation
- `m1::Real`: group 1 mean
- `n1::Int64`: group 1 sample size
- `alpha::Float64=0.05`: the probability of type I error

# Returns

- `p::Float64`: study power
"""
function power_c1g(; m0::Real, s0::Real, m1::Real, n1::Int64, alpha::Float64=0.05)

    delta = abs(m0 - m1)

    z = -crit_z(1 - alpha / 2) + (delta * sqrt(n1) / s0)

    return z2pow(z)

end

"""
    power_p2g(x)

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

    z = (delta / (sqrt(((p1 * q1) / n1) + (p2 * q2) / n2))) - crit_z(1 - alpha / 2) * ((sqrt(p_dash * q_dash * ((1 / n1) + (1 / n2)))) / (sqrt(((p1 * q1)  / n1) + ((p2 * q2) / n2))))

    return z2pow(z)

end

"""
    power_p1g(x)

Calculate required sample size for a proportion (group 1 vs population).

# Arguments

- `p0::Float64`: population incidence
- `p1::Float64`: group 1 anticipated incidence
- `n1::Int64`: group 1 sample size
- `alpha::Float64=0.05`: the probability of type I error

# Returns

- `p::Float64`: study power
"""
function power_p1g(; p0::Float64, p1::Float64, n1::Int64, alpha::Float64=0.05)

    q0 = 1 - p0
    q1 = 1 - p1

    z = (sqrt(n1 * ((p1 - p0)^2 / (p0 * q0))) - crit_z(1 - alpha / 2)) / (sqrt((p1 * q1)/(p0 * q0)))

    return z2pow(z)

end

"""
    size_c1diff(s1, s2, alpha, power)

Calculate required sample size for detecting a difference in a continuous variable (group 1 vs population).

# Arguments

- `s0::Real`: population standard deviation
- `s1::Real`: study standard deviation that we want to detect
- `two_sided::Bool=true`: if true, the estimation is for two-sided difference
- `power::Float64=0.8`: the ability to detect a difference between groups (power = 1 - beta, where beta is the probability of type II error)

# Returns

- `n::Int64`: study sample size
"""
function size_c1diff(; s0::Real, s1::Real, two_sided::Bool=true, power::Float64=0.8)

    sdiff = s1 / s0

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

    return two_sided ? 2 * n : n

end

"""
    size_p1diff(p0, p1, power)

Calculate required sample size for detecting a difference in a proportion (group 1 vs population).

# Arguments

- `p0::Real`: population proportion
- `p1::Real`: study proportion that we want to detect
- `power::Float64=0.8`: the ability to detect a difference between groups (power = 1 - beta, where beta is the probability of type II error)

# Returns

- `n::Int64`: study sample size (for both study groups)
"""
function size_p1diff(; p0::Real, p1::Real, power::Float64=0.8)

    p = (p0 + p1) / 2
    sdiff = round((p0 - p1) / sqrt(p * (1 - p)), digits=1)

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

    n = table[sdiff_idx, power_idx]

    return 2 * n

end