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
function hz2rads(f::Real)::Float64

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
function rads2hz(f::Real)::Float64

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
function t2f(t::Real)::Float64

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
function f2t(f::Real)::Float64

    @assert f > 0 "f must be > 0."

    return round(1000 / f, digits=2)

end

"""
    freqs(t; nf)

Return vector of frequencies and Nyquist frequency for time vector.

# Arguments

- `t::AbstractVector, AbstractRange}`: time vector
- `nf::Bool=false`: if true, return negative and positive frequencies, otherwise return positive frequencies only

# Returns

- `hz::Vector{Float64}`
- `nqf::Float64`
"""
function freqs(t::Union{AbstractVector, AbstractRange}; nf::Bool=false)::Tuple{Vector{Float64}, Float64}

    typeof(t) <: AbstractRange && (t = collect(t))

    # sampling interval
    dt = t[2] - t[1]

    # sampling rate
    fs = round(Int64, 1 / dt)

    # Nyquist frequency
    nqf = fs / 2

    # frequency array
    # hz = linspace(0, nf, floor(Int64, length(t) / 2))
    if nf
        hz = fftshift(round.(Vector(fftfreq(length(t), fs)), digits=3))
    else
        hz = round.(Vector(rfftfreq(length(t), fs)), digits=3)
    end

    return hz, nqf

end

"""
    freqs(s, fs; nf)

Return vector of frequencies and Nyquist frequency for signal.

# Arguments

- `s::AbstractVector`
- `fs::Int64`
- `nf::Bool=false`: if true, return negative and positive frequencies, otherwise return positive frequencies only

# Returns

- `hz::Vector{Float64`: signal vector
- `nqf::Float64`
"""
function freqs(s::AbstractVector, fs::Int64; nf::Bool=false)::Tuple{Vector{Float64}, Float64}

    @assert fs >= 1 "fs must be ≥ 1."

    # Nyquist frequency
    nqf = fs / 2
    # frequency array
    # hz = linspace(0, nf, floor(Int64, length(s) / 2) + 1)
    if nf
        hz = fftshift(round.(Vector(fftfreq(length(s), fs)), digits=3))
    else
        hz = round.(Vector(rfftfreq(length(s), fs)), digits=3)
    end

    return hz, nqf

end

"""
    freqs(obj; nf)

Return vector of frequencies and Nyquist frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `nf::Bool=false`: if true, return negative and positive frequencies, otherwise return positive frequencies only

# Returns

Named tuple containing:
- `hz::Vector{Float64}`
- `nqf::Float64`
"""
function freqs(obj::NeuroAnalyzer.NEURO; nf::Bool=false)::@NamedTuple{hz::Vector{Float64}, nqf::Float64}

    hz, nqf = freqs(obj.data[1, :, 1], sr(obj), nf=nf)

    return (hz=hz, nqf=nqf)

end
