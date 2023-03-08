export s_normalize
export s_normalize_zscore
export s_normalize_minmax
export s_normalize_n
export s_normalize_log
export s_normalize_gauss
export s_normalize_log10
export s_normalize_neglog
export s_normalize_neglog10
export s_normalize_neg
export s_normalize_pos
export s_normalize_perc

"""
    s_normalize(signal, n; method)

Normalize.

# Arguments

- `signal::AbstractArray`
- `n::Real`
- `method::Symbol`: :zscore, :minmax, :log, :log10, :neglog, :neglog10, :neg, :pos, :perc, :gauss, :invroot, :n, :none

# Returns

- `s_normalized::Vector{Float64}`
"""
function s_normalize(signal::AbstractArray, n::Real=1.0; method::Symbol)

    _check_var(method, [:zscore, :minmax, :log, :log10, :neglog, :neglog10, :neg, :pos, :perc, :gauss, :invroot, :n, :none], "method")

    if method === :zscore
        return s_normalize_zscore(signal)
    elseif method === :minmax
        return s_normalize_minmax(signal)
    elseif method === :log
        return s_normalize_log(signal)
    elseif method === :log10
        return s_normalize_log10(signal)
    elseif method === :neglog
        return s_normalize_neglog(signal)
    elseif method === :neglog10
        return s_normalize_neglog10(signal)
    elseif method === :neg
        return s_normalize_neg(signal)
    elseif method === :pos
        return s_normalize_pos(signal)
    elseif method === :perc
        return s_normalize_perc(signal)
    elseif method === :gauss
        return s_normalize_gauss(signal)
    elseif method === :invroot
        return s_normalize_invroot(signal)
    elseif method === :n
        return s_normalize_n(signal, n)
    elseif method === :none
        return signal
    end
end

"""
    s_normalize_zscore(signal)

Normalize (by z-score).

# Arguments

- `signal::AbstractArray`

# Returns

- `s_normalized::Vector{Float64}`
"""
function s_normalize_zscore(signal::AbstractArray)
    m = mean(signal)
    s = std(signal)
    return @. (signal - m) / s
end

"""
    s_normalize_minmax(signal)

Normalize in [-1, +1].

# Arguments

- `signal::AbstractArray`

# Returns

- `s_normalized::AbstractArray`
"""
function s_normalize_minmax(signal::AbstractArray)
    mi = minimum(signal)
    mx = maximum(signal)
    mxi = mx - mi
    return @. (2 * (signal - mi) / mxi) - 1
end

"""
    s_normalize_n(signal, n)

Normalize in [0, n], default is [0, +1].

# Arguments

- `signal::AbstractArray`
- `n::Real=1.0`

# Returns

- `s_normalized::AbstractArray`
"""
function s_normalize_n(signal::AbstractArray, n::Real=1.0)
    max_x = maximum(signal)
    min_x = minimum(signal)
    return n .* (signal .- min_x) ./ (max_x - min_x)
end

"""
    s_normalize_log(signal)

Normalize using log-transformation.

# Arguments

- `signal::AbstractArray`

# Returns

- `s_normalized::AbstractArray`
"""
function s_normalize_log(signal::AbstractArray)
    m = abs(minimum(signal))
    return @. log(1 + signal + m)
end

"""
    s_normalize_gauss(signal, dims)

Normalize to Gaussian.

# Arguments

- `signal::AbstractArray`
- `dims::Int64=1`: dimension for cumsum()

# Returns

- `s_normalized::Vector{Float64}`
"""
function s_normalize_gauss(signal::AbstractArray, dims::Int64=1)
    dims in 1:ndims(signal) || throw(ArgumentError("dims must be in: 1:$(ndims(signal)).")) 
    l = length(signal) + 1
    return atanh.((tiedrank(cumsum(signal, dims=dims)) ./ l .- 0.5) .* 2)
end

"""
    s_normalize_log10(signal)

Normalize using log10-transformation.

# Arguments

- `signal::AbstractArray`

# Returns

- `s_normalized::Vector{Float64}`
"""
function s_normalize_log10(signal::AbstractArray)
    m = 1 + abs(minimum(signal))
    return @. log10(signal + m)
end

"""
    s_normalize_neglog(signal)

Normalize to using -log-transformation.

# Arguments

- `signal::AbstractArray`

# Returns

- `s_normalized::Vector{Float64}`
"""
function s_normalize_neglog(signal::AbstractArray)
    return @. -log(signal)
end

"""
    s_normalize_neglog10(signal)

Normalize using -log10-transformation.

# Arguments

- `signal::AbstractArray`

# Returns

- `s_normalized::Vector{Float64}`
"""
function s_normalize_neglog10(signal::AbstractArray)
    return @. -log10(signal)
end

"""
    s_normalize_neg(signal)

Normalize in [0, -∞].

# Arguments

- `signal::AbstractArray`

# Returns

- `s_normalized::Vector{Float64}`
"""
function s_normalize_neg(signal::AbstractArray)
    m = maximum(signal)
    return @. signal - m
end

"""
    s_normalize_pos(signal)

Normalize in [0, +∞].

# Arguments

- `signal::AbstractArray`

# Returns

- `s_normalized::Vector{Float64}`
"""
function s_normalize_pos(signal::AbstractArray)
    m = abs(minimum(signal))
    return @. signal + m
end

"""
    s_normalize_perc(signal)

Normalize in percentages.

# Arguments

- `signal::AbstractArray`

# Returns

- `s_normalized::Vector{Float64}`
"""
function s_normalize_perc(signal::AbstractArray)
    m1 = minimum(signal)
    m2 = maximum(signal)
    m = m2 - m1
    return (signal .- m1) ./ m
end
