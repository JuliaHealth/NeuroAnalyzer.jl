"""
    s_normalize_neglog(signal)

Normalize `signal` to -log.

# Arguments

- `signal::AbstractArray`

# Returns

- `s_normalized::Vector{Float64}`
"""
function s_normalize_neglog(signal::AbstractArray)

    s_normalized = @. -log(signal)

    return s_normalized
end

"""
    s_normalize_neglog10(signal)

Normalize `signal` to -log10.

# Arguments

- `signal::AbstractArray`

# Returns

- `s_normalized::Vector{Float64}`
"""
function s_normalize_neglog10(signal::AbstractArray)

    s_normalized = @. -log10(signal)

    return s_normalized
end

"""
    s_normalize_neg(signal)

Normalize `signal` in [0, -∞].

# Arguments

- `signal::AbstractArray`

# Returns

- `s_normalized::Vector{Float64}`
"""
function s_normalize_neg(signal::AbstractArray)

    m = maximum(signal)
    s_normalized = @. signal - m

    return s_normalized
end

"""
    s_normalize_pos(signal)

Normalize `signal` in [0, +∞].

# Arguments

- `signal::AbstractArray`

# Returns

- `s_normalized::Vector{Float64}`
"""
function s_normalize_neg(signal::AbstractArray)

    m = minimum(signal)
    s_normalized = @. signal + abs(m)

    return s_normalized
end

"""
    s_normalize_perc(signal)

Normalize `signal` in percentages.

# Arguments

- `signal::AbstractArray`

# Returns

- `s_normalized::Vector{Float64}`
"""
function s_normalize_neg(signal::AbstractArray)

    s_normalized = (signal .- minimum(signal)) ./ (maximum(signal) .- minimum(signal))
    
    return s_normalized
end

"""
    s_normalize_absmin(signal)

Normalize `signal` in [-abs(min), abs(min)].

# Arguments

- `signal::AbstractArray`

# Returns

- `s_normalized::Vector{Float64}`
"""
function s_normalize_absmin(signal::AbstractArray)

    m = minimum(signal)
    s_normalized = @. signal + abs(m)

    return s_normalized
end