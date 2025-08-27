export t2s
export s2t
export markers_s2t
export markers_s2t!

"""
    t2s(t, fs)

Convert time in seconds to sample number.

# Arguments

- `t::Real`: time in s
- `fs::Int64`: sampling rate

# Returns

- `t2s::Int64`: sample number
"""
function t2s(t::Real, fs::Int64)::Int64

    @assert t >= 0 "t must be ≥ 0."

    return t == 0 ? 1 : ceil(Int64, t * fs)

end

"""
    s2t(s, fs)

Convert sample number to time in seconds.

# Arguments

- `s::Real`: sample number
- `fs::Int64`: sampling rate

# Returns

- `s2t::Float64`: time in s
"""
function s2t(s::Real, fs::Int64)::Float64

    if s == 0
        _info("Sample number 0 replaced with 1")
        s = 1
    end

    return round(s / fs - (1/fs), digits=3)

end

"""
    t2s(obj; <keyword arguments>)

Convert time in seconds to sample number.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `t::Real`: time in seconds

# Returns

- `t2s::Int64`: sample number
"""
function t2s(obj::NeuroAnalyzer.NEURO; t::Real)::Int64

    return t2s(t, sr(obj))

end

"""
    s2t(obj; <keyword arguments>)

Convert time in samples to seconds.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `s::Int64`: time in samples

# Returns

- `s2t::Float64`: time in seconds
"""
function s2t(obj::NeuroAnalyzer.NEURO; s::Int64)::Float64

    return s2t(s, sr(obj))

end

"""
    markers_s2t(m; <keyword arguments>)

Convert markers start and length from samples to seconds.

# Arguments

- `m::DataFrame`: markers
- `fs::Int64`: sampling rate

# Returns

- `m_new::DataFrame`
"""
function markers_s2t(m::DataFrame; fs::Int64)::DataFrame

    m_new = deepcopy(m)
    m_new[:, :start] = s2t.(m[:, :start], fs)
    m_new[:, :length] = s2t.(m[:, :length], fs)

    return m_new

end

"""
    markers_s2t!(m; <keyword arguments>)

Convert markers start and length from samples to seconds.

# Arguments

- `m::DataFrame`: markers
- `fs::Int64`: sampling rate

# Returns

Nothing
"""
function markers_s2t!(m::DataFrame; fs::Int64)::Nothing

    m[:, :start] = s2t.(m[:, :start], fs)
    m[:, :length] = s2t.(m[:, :length], fs)

    return nothing

end

"""
    markers_s2t(obj; <keyword arguments>)

Convert markers start and length from samples to seconds.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `m_new::DataFrame`
"""
function markers_s2t(obj::NeuroAnalyzer.NEURO)::DataFrame

    fs = sr(obj)
    m_new = deepcopy(obj.markers)
    m_new[:, :start] = s2t.(m_new[:, :start], fs)
    m_new[:, :length] = s2t.(m_new[:, :length], fs)

    return m_new

end

"""
    markers_s2t!(obj; <keyword arguments>)

Convert markers start and length from samples to seconds.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

Nothing
"""
function markers_s2t!(obj::NeuroAnalyzer.NEURO)::Nothing

    fs = sr(obj)
    obj.markers[:, :start] = s2t.(obj.markers[:, :start], fs)
    obj.markers[:, :length] = s2t.(obj.markers[:, :length], fs)

    return nothing

end
