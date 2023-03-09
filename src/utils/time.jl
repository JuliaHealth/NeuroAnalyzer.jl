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

"""
    t2s(obj; t)

Convert time in seconds to samples.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `t::Real`: time in seconds

# Returns

- `t_s::Int64`: time in samples
"""
function t2s(obj::NeuroAnalyzer.NEURO; t::Real)
    return floor(Int64, t * sr(obj)) + 1
end

"""
    s2t(obj; t)

Convert time in samples to seconds.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `t::Int64`: time in samples

# Returns

- `t_s::Float64`: time in seconds
"""
function s2t(obj::NeuroAnalyzer.NEURO; t::Int64)
    return round(t / sr(obj), digits=2)
end
