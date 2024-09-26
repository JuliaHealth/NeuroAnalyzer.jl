export cwtfrq

"""
    cwtfrq(s; <keyword arguments>)

Return mean frequencies of a collection of analytic or real wavelets for a given signal.

# Arguments

- `s::AbstractVector`
- `fs::Int64`: sampling rate
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

- `f::Vector{Float64}`: frequencies
"""
function cwtfrq(s::AbstractVector; fs::Int64, wt::T=wavelet(Morlet(2π), β=32, Q=128)) where {T <: CWT}

    @assert fs >= 1 "fs must be ≥ 1."
    f = round.(ContinuousWavelets.getMeanFreq(length(s), wt, fs), digits=2)
    f[1] = 0

    return f

end

"""
    cwtfrq(s; <keyword arguments>)

Return mean frequencies of a collection of analytic or real wavelets for a given signal.

# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

- `f::Vector{Float64}`: frequencies
"""
function cwtfrq(s::AbstractArray; fs::Int64, wt::T=wavelet(Morlet(2π), β=32, Q=128)) where {T <: CWT}

    _chk3d(s)
    @assert fs >= 1 "fs must be ≥ 1."
    f = round.(ContinuousWavelets.getMeanFreq(size(s, 2), wt, fs), digits=2)
    f[1] = 0

    return f

end

"""
    cwtfrq(obj; <keyword arguments>)

Return mean frequencies of a collection of analytic or real wavelets for a given signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

- `f::Vector{Float64}`: frequencies
"""
function cwtfrq(obj::NeuroAnalyzer.NEURO; wt::T=wavelet(Morlet())) where {T <: CWT}

    f = @views cwtfrq(obj.data[1, :, 1], fs=sr(obj), wt=wt)

    return f

end