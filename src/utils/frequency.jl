export hz2rads
export rads2hz
export t2f
export f2t
export freqs

"""
    hz2rads(f)

Convert frequency in Hz to rad/s.

# Arguments

- `f::Real`

# Returns

- `f_rads::Float64`
"""
function hz2rads(f::Real)

    return 2pi * f

end

"""
    rads2hz(f)

Convert frequency in rad/s to Hz.

# Arguments

- `f::Real`

# Returns

- `f_rads::Float64`
"""
function rads2hz(f::Real)

    return f / 2pi

end

"""
    t2f(t)

Convert cycle length in ms to frequency.

# Arguments

- `t::Real`: cycle length in ms

# Returns

- `f::Float64`: frequency in Hz
"""
function t2f(t::Real)

    @assert t > 0 "t must be > 0."

    return round(1000 / t, digits=2)

end

"""
    f2t(f)

Convert frequency in Hz to cycle length in ms.

# Arguments

- `f::Real`: frequency in Hz

# Returns

- `f::Float64`: cycle length in ms
"""
function f2t(f::Real)

    @assert f > 0 "f must be > 0."

    return round(1000 / f, digits=2)

end

"""
    freqs(t)

Return vector of frequencies and Nyquist frequency for time vector.

# Arguments

- `t::AbstractVector, AbstractRange}`: time vector

# Returns

- `hz::Vector{Float64}`
- `nf::Float64`
"""
function freqs(t::Union{AbstractVector, AbstractRange})

    typeof(t) <: AbstractRange && (t = collect(t))

    # sampling interval
    dt = t[2] - t[1]

    # sampling rate
    fs = round(Int64, 1 / dt)

    # Nyquist frequency
    nf = fs / 2

    # frequency array
    # hz = linspace(0, nf, floor(Int64, length(t) / 2))
    hz = Vector(rfftfreq(length(t), fs))

    return hz, nf

end

"""
    freqs(s, fs)

Return vector of frequencies and Nyquist frequency for signal.

# Arguments

- `s::AbstractVector`
- `fs::Int64`

# Returns

- `hz::Vector{Float64`: signal vector
- `nf::Float64`
"""
function freqs(s::AbstractVector, fs::Int64)

    @assert fs >= 1 "fs must be ≥ 1."

    # Nyquist frequency
    nf = fs / 2
    # frequency array
    # hz = linspace(0, nf, floor(Int64, length(s) / 2))
    hz = round.(Vector(rfftfreq(length(s), fs)), digits=3)

    return hz, nf

end

"""
    freqs(obj)

Return vector of frequencies and Nyquist frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

Named tuple containing:
- `hz::Vector{Float64}`
- `nf::Float64`
"""
function freqs(obj::NeuroAnalyzer.NEURO)

    hz, nf = freqs(obj.data[1, :, 1], sr(obj))

    return (hz=hz, nf=nf)

end
