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

- `Float64`: frequency in Hz
"""
function t2f(t::Real)::Float64
    # validate
    t > 0 || throw(ArgumentError("t must be > 0."))

    return round(1000 / t; digits = 2)
end

"""
    f2t(f)

Convert a frequency in Hz to a cycle length in milliseconds.

Computes `1000 / f`, rounded to 2 decimal places.

# Arguments

- `f::Real`: frequency in Hz; must be > 0

# Returns

- `Float64`: cycle length in ms
"""
function f2t(f::Real)::Float64
    # validate
    f > 0 || throw(ArgumentError("f must be > 0."))

    return round(1000 / f; digits = 2)
end

"""
    freqs(t; <keyword arguments>)

Return the frequency vector and Nyquist frequency for a time vector.

The sampling rate is inferred from the mean sampling interval, `(t[end] - t[1]) / (length(t) - 1)`, and rounded to an integer.

# Arguments

- `t::AbstractVector`: uniformly sampled, increasing time vector in seconds; must contain at least 2 elements
- `nf::Bool=false`: if `true`, return both negative and positive frequencies (via `fftfreq` + `fftshift`); if `false`, return positive frequencies only (via `rfftfreq`)

# Returns

- `Vector{Float64}`: frequency vector in Hz, rounded to 3 decimal places
- `Float64`: Nyquist frequency in Hz
"""
function freqs(t::AbstractVector; nf::Bool = false)::Tuple{Vector{Float64}, Float64}
    # validate
    n = length(t)
    n >= 2 || throw(ArgumentError("t must contain at least 2 elements."))

    # mean sampling interval (more robust to floating-point steps than t[2] - t[1])
    dt = (t[end] - t[1]) / (n - 1)
    dt > 0 || throw(ArgumentError("t must be increasing."))
    fs = round(Int64, 1 / dt)
    fs >= 1 || throw(ArgumentError("Sampling rate inferred from t is < 1 Hz."))

    return freqs(n, fs; nf = nf)
end

"""
    freqs(s, fs; <keyword arguments>)

Return the frequency vector and Nyquist frequency for a signal vector.

# Arguments

- `s::AbstractVector`: signal vector
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `nf::Bool=false`: if `true`, return both negative and positive frequencies (via `fftfreq` + `fftshift`); if `false`, return positive frequencies only (via `rfftfreq`)

# Returns

- `Vector{Float64}`: frequency vector in Hz, rounded to 3 decimal places
- `Float64`: Nyquist frequency in Hz
"""
function freqs(
    s::AbstractVector,
    fs::Int64;
    nf::Bool = false,
)::Tuple{Vector{Float64}, Float64}
    return freqs(length(s), fs; nf = nf)
end

"""
    freqs(n, fs; <keyword arguments>)

Return the frequency vector and Nyquist frequency for a signal of `n` samples.

# Arguments

- `n::Int64`: number of samples; must be ≥ 1
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `nf::Bool=false`: if `true`, return both negative and positive frequencies (via `fftfreq` + `fftshift`); if `false`, return positive frequencies only (via `rfftfreq`)

# Returns

- `Vector{Float64}`: frequency vector in Hz, rounded to 3 decimal places
- `Float64`: Nyquist frequency in Hz

# Notes

Bins are at `k × fs / n`. For odd `n` there is no bin at the Nyquist frequency; the highest positive bin is `fs/2 − fs/(2n)`.
"""
function freqs(
    n::Int64,
    fs::Int64;
    nf::Bool = false,
)::Tuple{Vector{Float64}, Float64}
    # validate
    n >= 1 || throw(ArgumentError("n must be ≥ 1."))
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))

    # Nyquist frequency
    nqf = fs / 2
    # frequency vector
    hz =
        nf ? fftshift(round.(Vector(fftfreq(n, fs)); digits = 3)) :
        round.(Vector(rfftfreq(n, fs)); digits = 3)

    return hz, nqf
end

"""
    freqs(obj; <keyword arguments>)

Return the frequency vector and Nyquist frequency for a NEURO object.

Uses the epoch length and sampling rate of `obj`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `nf::Bool=false`: if `true`, return both negative and positive frequencies; if `false`, return positive frequencies only

# Returns

- `Vector{Float64}`: frequency vector in Hz, rounded to 3 decimal places
- `Float64`: Nyquist frequency in Hz
"""
function freqs(obj::NeuroAnalyzer.NEURO; nf::Bool = false)::Tuple{Vector{Float64}, Float64}
    return freqs(epoch_len(obj), sr(obj); nf = nf)
end