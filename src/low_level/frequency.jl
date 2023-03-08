export hz2rads
export rads2hz
export s2t
export t2s

"""
    hz2rads(f)

Convert frequency `f` in Hz to rad/s.

# Arguments

- `f::Real`

# Returns

- `f_rads::Float64`
"""
function hz2rads(f::Real)
    return 2 * pi * f
end

"""
    rads2hz(f)

Convert frequency `f` in rad/s to Hz.

# Arguments

- `f::Real`

# Returns

- `f_rads::Float64`
"""
function rads2hz(f::Real)
    return f / 2 * pi
end

"""
    t2f(t)

Convert cycle length in ms `t` to frequency.

# Arguments

- `t::Real`: cycle length in ms

# Returns

- `f::Float64`: frequency in Hz
"""
function t2f(t::Real)
    t <= 0 && throw(ArgumentError("t must be > 0."))
    return round(1000 / t, digits=2)
end

"""
    f2t(f)

Convert frequency `f` to cycle length in ms.

# Arguments

- `f::Real`: frequency in Hz

# Returns

- `f::Float64`: cycle length in ms
"""
function f2t(f::Real)
    f <= 0 && throw(ArgumentError("f must be > 0."))
    return round(1000 / f, digits=2)
end
