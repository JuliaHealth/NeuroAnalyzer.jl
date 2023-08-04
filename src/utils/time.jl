export t2s
export s2t

"""
    t2s(t, fs)

Convert time to sample number.

# Arguments

- `t::T`: time in s
- `fs::Int64`: sampling rate

# Returns

- `t2s::Int64`: sample number
"""
function t2s(t::T, fs::Int64) where {T<:Real}
    
    @assert t >= 0 "t must be ≥ 0."
    return t == 0 ? 1 : round(Int64, t * fs)

end

"""
    s2t(s, fs)

Convert sample number to time.

# Arguments

- `t::Int64`: sample number
- `fs::Int64`: sampling rate

# Returns

- `s2t::Float64`: time in s
"""
function s2t(s::Int64, fs::Int64)

    @assert s > 0 "s must be > 0."
    
    return s / fs

end

"""
    t2s(obj; t)

Convert time in seconds to samples.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `t::T`: time in seconds

# Returns

- `t2s::Int64`: time in samples
"""
function t2s(obj::NeuroAnalyzer.NEURO; t::T) where {T<:Real}

    return floor(Int64, t * sr(obj)) + 1

end

"""
    s2t(obj; s)

Convert time in samples to seconds.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `s::Int64`: time in samples

# Returns

- `s2t::Float64`: time in seconds
"""
function s2t(obj::NeuroAnalyzer.NEURO; s::Int64)

    return round(s / sr(obj), digits=2)

end
