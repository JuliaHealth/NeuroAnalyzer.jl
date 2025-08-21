export t2s
export s2t
export markers_s2t
export markers_s2t!

"""
    t2s(t, fs)

Convert time in seconds to sample number.

# Arguments

- `t::T`: time in s
- `fs::Int64`: sampling rate

# Returns

- `t2s::Int64`: sample number
"""
function t2s(t::T, fs::Int64)::Int64 where {T<:Real}

    @assert t >= 0 "t must be ≥ 0."

    return t == 0 ? 1 : round(Int64, t * fs)

end

"""
    s2t(s, fs)

Convert sample number to time in seconds.

# Arguments

- `t::Int64`: sample number
- `fs::Int64`: sampling rate

# Returns

- `s2t::Float64`: time in s
"""
function s2t(s::Int64, fs::Int64)::Float64

    @assert s > 0 "s must be > 0."

    return s / fs

end

"""
    t2s(obj; <keyword arguments>)

Convert time in seconds to sample number.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `t::T`: time in seconds

# Returns

- `t2s::Int64`: time in samples
"""
function t2s(obj::NeuroAnalyzer.NEURO; t::T)::Int64 where {T<:Real}

    return floor(Int64, t * sr(obj)) + 1

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

    return round(s / sr(obj), digits=3)

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
    m_new[:, :start] = round.(m[:, :start] ./ fs, digits=3)
    m_new[:, :length] = round.(m[:, :length] ./ fs, digits=3)

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

    m[:, :start] = round.(m[:, :start] ./ fs, digits=3)
    m[:, :length] = round.(m[:, :length] ./ fs, digits=3)

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

    m_new = deepcopy(obj.markers)
    m_new[:, :start] = round.(m_new[:, :start] ./ sr(obj), digits=2)
    m_new[:, :length] = round.(m_new[:, :length] ./ sr(obj), digits=2)

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

    obj.markers[:, :start] = round.(obj.markers[:, :start] ./ sr(obj), digits=3)
    obj.markers[:, :length] = round.(obj.markers[:, :length] ./ sr(obj), digits=3)

    return nothing

end
