export hz2rads
export rads2hz
export t2f
export f2t
export freqs

"""
    hz2rads(f)

Convert a frequency from Hz to rad/s.

Computes `2π × f`.

# Arguments

- `f::Real`: frequency in Hz

# Returns

- `Float64`: frequency in rad/s

# See also

[`rads2hz`](@ref)
"""
function hz2rads(f::Real)::Float64

    return 2pi * f

end

"""
    rads2hz(f)

Convert a frequency from rad/s to Hz.

Computes `f / 2π`.

# Arguments

- `f::Real`: frequency in rad/s

# Returns

- `Float64`: frequency in Hz

# See also

[`hz2rads`](@ref)
"""
function rads2hz(f::Real)::Float64

    return f / 2pi

end

"""
    t2f(t)

Convert a cycle length in milliseconds to a frequency in Hz.

Computes `1000 / t`, rounded to 2 decimal places.

# Arguments

- `t::Real`: cycle length in ms; must be > 0

# Returns

- `f::Float64`: frequency in Hz

# Throws

- `ArgumentError`: if `t ≤ 0`.

# See also

[`f2t`](@ref)
"""
function t2f(t::Real)::Float64

    @assert t > 0 "t must be > 0."

    return round(1000 / t, digits = 2)

end

"""
    f2t(f)

Convert a frequency in Hz to a cycle length in milliseconds.

Computes `1000 / f`, rounded to 2 decimal places.

# Arguments

- `f::Real`: frequency in Hz; must be > 0

# Returns

- `Float64`: cycle length in ms

# Throws

- `ArgumentError`: if `f ≤ 0`

# See also

[`t2f`](@ref)
"""
function f2t(f::Real)::Float64

    @assert f > 0 "f must be > 0."

    return round(1000 / f, digits = 2)

end

"""
    freqs(t; nf)

Return the frequency vector and Nyquist frequency for a time vector.

The sampling rate is inferred as `1 / (t[2] - t[1])`.

# Arguments

- `t::AbstractVector, AbstractRange}`: time vector; must contain at least 2 elements
- `nf::Bool=false`: if `true`, return both negative and positive frequencies (via `fftfreq` + `fftshift`); if `false`, return positive frequencies only (via `rfftfreq`)

# Returns

- `hz::Vector{Float64}`: frequency vector in Hz, rounded to 3 decimal places
- `nqf::Float64`: Nyquist frequency in Hz

# Throws
- `ArgumentError`: if `length(t) < 2`

# See also

[`freqs(::AbstractVector, ::Int64)`](@ref),
[`freqs(::Int64, ::Int64)`](@ref)
"""
function freqs(
        t::Union{AbstractVector, AbstractRange}; nf::Bool = false
    )::Tuple{Vector{Float64}, Float64}

    @assert length(t) >= 2 "t must contain at least 2 elements."
    
    # materialize ranges so indexing is always valid
    t = collect(t)

    # infer sampling rate from the interval between the first two samples
    dt = t[2] - t[1]
    fs = round(Int64, 1 / dt)
    # Nyquist frequency
    nqf = fs / 2
    # frequency vector
    hz = nf ? fftshift(round.(Vector(fftfreq(length(t), fs)), digits=3)) :
              round.(Vector(rfftfreq(length(t), fs)), digits=3)

    return hz, nqf

end

"""
    freqs(s, fs; <keyword arguments>)

Return the frequency vector and Nyquist frequency for a signal vector.

# Arguments

- `s::AbstractVector`: signal vector
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `nf::Bool=false`: if `true`, return both negative and positive frequencies (via `fftfreq` + `fftshift`); if `false`, return positive frequencies only (via `rfftfreq`)

# Returns

- `hz::Vector{Float64}`: frequency vector in Hz, rounded to 3 decimal places
- `nqf::Float64`: Nyquist frequency in Hz

# Throws
- `ArgumentError`: if `fs < 1`

# See also
[`freqs(::Union{AbstractVector, AbstractRange})`](@ref),
[`freqs(::Int64, ::Int64)`](@ref)
"""
function freqs(
    s::AbstractVector,
    fs::Int64;
    nf::Bool = false
)::Tuple{Vector{Float64}, Float64}

    @assert fs >= 1 "fs must be ≥ 1."
    # Nyquist frequency
    nqf = fs / 2
    # frequency vector
    hz = nf ? fftshift(round.(Vector(fftfreq(length(s), fs)), digits=3)) :
              round.(Vector(rfftfreq(length(s), fs)), digits=3)

    return hz, nqf

end

"""
    freqs(n, fs; <keyword arguments>)

Return the frequency vector and Nyquist frequency for a signal of `n` samples.

# Arguments

- `n::Int64`: number of samples; must be ≥ 1
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `nf::Bool=false`: if `true`, return both negative and positive frequencies (via `fftfreq` + `fftshift`); if `false`, return positive frequencies only (via `rfftfreq`)

# Returns

- `hz::Vector{Float64}`: frequency vector in Hz, rounded to 3 decimal places
- `nqf::Float64`: Nyquist frequency in Hz

# Throws

- `ArgumentError`: if `n < 1` or `fs < 1`.

# See also

[`freqs(::Union{AbstractVector, AbstractRange})`](@ref),
[`freqs(::AbstractVector, ::Int64)`](@ref)
"""
function freqs(
    n::Int64,
    fs::Int64;
    nf::Bool = false
)::Tuple{Vector{Float64}, Float64}

    @assert fs >= 1 "fs must be ≥ 1."

    # Nyquist frequency
    nqf = fs / 2
    # frequency vector
    hz = nf ? fftshift(round.(Vector(fftfreq(n, fs)), digits=3)) :
              round.(Vector(rfftfreq(n, fs)), digits=3)

    return hz, nqf

end

"""
    freqs(obj; <keyword arguments>)

Return the frequency vector and Nyquist frequency for a NEURO object.

Uses the first channel and first epoch of `obj` to infer signal length, and reads the sampling rate via `sr(obj)`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `nf::Bool=false`: if `true`, return both negative and positive frequencies; if `false`, return positive frequencies only

# Returns

Named tuple:

- `hz::Vector{Float64}`: frequency vector in Hz, rounded to 3 decimal places
- `nqf::Float64`: Nyquist frequency in Hz

# See also

[`freqs(::AbstractVector, ::Int64)`](@ref)
"""
function freqs(
    obj::NeuroAnalyzer.NEURO; nf::Bool = false
)::@NamedTuple{
    hz::Vector{Float64},
    nqf::Float64
}

    hz, nqf = freqs(obj.data[1, :, 1], sr(obj); nf = nf)

    return (hz = hz, nqf = nqf)

end
