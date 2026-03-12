export cwtfrq

"""
    cwtfrq(s; <keyword arguments>)

Return the mean frequencies of a collection of analytic or real wavelets for a signal of a given length.

The first frequency bin is set to `0.0` Hz because `getMeanFreq` returns a non-physical low-frequency artifact for the lowest scale.

# Arguments

- `s::AbstractVector`: signal vector; used only for its length
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`: wavelet to use; see the [ContinuousWavelets.jl documentation](https://github.com/dsweber2/ContinuousWavelets.jl) for the full list of available wavelets

# Returns

- `Vector{Float64}`: center frequencies in Hz (rounded to 2 decimal places), length determined by the number of wavelet scales

# Throws
- `ArgumentError`: if `fs < 1`

# See also

[`cwtfrq(::AbstractArray)`](@ref), [`cwtfrq(::NeuroAnalyzer.NEURO)`](@ref)
"""
function cwtfrq(
    s::AbstractVector;
    fs::Int64,
    wt::T = wavelet(Morlet(2π), β = 2)
) where {T <: CWT}

    @assert fs >= 1 "fs must be ≥ 1."

    f = round.(ContinuousWavelets.getMeanFreq(length(s), wt, fs), digits=2)
    # lowest scale returns a non-physical frequency; replace with DC (0 Hz)
    f[1] = 0.0

    return f

end

"""
    cwtfrq(s; <keyword arguments>)

Return the mean frequencies of a collection of analytic or real wavelets for a 3-dimensional signal array.

Delegates to the vector method using the first channel and first epoch `s[1, :, 1]` to determine the wavelet frequency grid (all channels and epochs share the same grid for a fixed signal length).

# Arguments

- `s::AbstractArray`: 3-dimensional signal array `(channels, samples, epochs)`
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`: wavelet to use; see the [ContinuousWavelets.jl documentation](https://github.com/dsweber2/ContinuousWavelets.jl) for the full list of available wavelets

# Returns

- `Vector{Float64}`: center frequencies in Hz (rounded to 2 decimal places)

# Throws

- `ArgumentError`: if `s` is not 3-dimensional or `fs < 1`

# See also

[`cwtfrq(::AbstractVector)`](@ref), [`cwtfrq(::NeuroAnalyzer.NEURO)`](@ref)
"""
function cwtfrq(
    s::AbstractArray;
    fs::Int64,
    wt::T = wavelet(Morlet(2π), β = 2)
) where {T <: CWT}

    _chk3d(s)

    # all epochs/channels share the same length → any single slice is representative
    return cwtfrq(@view(s[1, :, 1]); fs=fs, wt=wt)

end

"""
    cwtfrq(obj; <keyword arguments>)

Return the mean frequencies of a collection of analytic or real wavelets for a NEURO object.

Uses the first channel and first epoch to determine the wavelet frequency grid.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`: wavelet to use; see the [ContinuousWavelets.jl documentation](https://github.com/dsweber2/ContinuousWavelets.jl) for the full list of available wavelets

# Returns

- `Vector{Float64}`: center frequencies in Hz (rounded to 2 decimal places)

# See also
[`cwtfrq(::AbstractVector)`](@ref), [`cwtfrq(::AbstractArray)`](@ref)
"""
function cwtfrq(
    obj::NeuroAnalyzer.NEURO;
    wt::T=wavelet(Morlet(2π), β=2)
) where {T <: CWT}

    _log_off()
    f = @views cwtfrq(obj.data[1, :, 1], fs = sr(obj), wt = wt)
    _log_on()

    return f

end
