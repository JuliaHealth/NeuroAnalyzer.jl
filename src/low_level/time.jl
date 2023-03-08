export t2f
export f2t

"""
    t2s(t, fs)

Convert time to sample number.

# Arguments

- `t::Real`: time in s
- `fs::Int64`: sampling rate

# Returns

- `s::Int64`: sample number
"""
function t2s(t::Real, fs::Int64)
    t < 0 && throw(ArgumentError("t must be ≥ 0."))
    if t == 0
        return 1
    else
        return round(Int64, t * fs)
    end
end

"""
    s2t(s, fs)

Convert sample number to time.

# Arguments

- `t::Int64`: sample number
- `fs::Int64`: sampling rate

# Returns

- `t::Float64`: time in s
"""
function s2t(s::Int64, fs::Int64)
    s < 0 && throw(ArgumentError("s must be > 0."))
    return s / fs
end

