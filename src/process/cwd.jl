export cwd
export icwd

"""
    cwd(s; <keyword arguments>)

Perform continuous wavelet decomposition (CWD) on a signal vector.

Returns the real part of the CWT coefficient matrix. Each row corresponds to one wavelet scale (frequency band); each column corresponds to one time point.

# Arguments

- `s::AbstractVector`: signal vector
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

- `Matrix{Float64}`: CWT coefficient matrix of shape `(n_scales,  length(s))`

# See also

[`icwd`](@ref), [`cwd(::AbstractArray)`](@ref)
"""
function cwd(s::AbstractVector; wt::T = wavelet(Morlet(2π), β = 2))::Matrix{Float64} where {T <: CWT}

    # ContinuousWavelets.cwt returns (samples × scales); transpose to (scales × samples)
    return Matrix(real.(ContinuousWavelets.cwt(s, wt))')

    return ct

end

"""
    cwd(s; <keyword arguments>)

Perform continuous wavelet decomposition on a 3-D signal array.

Applies [`cwd(::AbstractVector)`](@ref) to every channel × epoch slice in parallel and returns a 4-D coefficient array.

# Arguments

- `s::AbstractArray`: 3-D signal array `(channels, samples, epochs)`
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

- `Array{Float64, 4}`: CWT coefficients of shape `(channels, n_scales, samples, epochs)`

# See also

[`cwd(::AbstractVector)`](@ref), [`cwd(::NeuroAnalyzer.NEURO)`](@ref)
"""
function cwd(s::AbstractArray; wt::T = wavelet(Morlet(2π), β = 2))::Array{Float64, 4} where {T <: CWT}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels, epoch length and number of epochs
    ch_n, ep_len, ep_n = size(s)

    # dry run on the first slice to determine the number of wavelet scales
    _log_off()
    n_scales = size(ContinuousWavelets.cwt(s[1, :, 1], wt), 2)

    # pre-allocate output
    ct = zeros(ch_n, n_scales, ep_len, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        ct[ch_idx, :, :, ep_idx] = cwd(@view(s[ch_idx, :, ep_idx]), wt=wt)
    end

    _log_on()

    return ct

end

"""
    cwd(obj; <keyword arguments>)

Perform continuous wavelet decomposition on selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

- `Array{Float64, 4}`: CWT coefficients of shape `(channels, n_scales, samples, epochs)`

# See also

[`icwd`](@ref), [`cwd(::AbstractArray)`](@ref)
"""
function cwd(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    wt::T=wavelet(Morlet(2π), β=2)
)::Array{Float64, 4} where {T <: CWT}

    # resolve channel names to integer indices
    ch = get_channel(obj, ch=ch)

    return cwd(@view(obj.data[ch, :, :]); wt=wt)

end

"""
    icwd(ct; <keyword arguments>)

Perform inverse continuous wavelet transformation (iCWT).

Reconstructs the original signal from a CWT coefficient matrix produced by [`cwd`](@ref).

# Arguments

- `ct::Matrix{Float64}`: CWT coefficient matrix of shape `(n_scales, samples)` (as returned by `cwd`)
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`: continuous wavelet; must match the one used in the forward transform
- `type::Symbol=:pd`: reconstruction method:
    - `:pd`: PenroseDelta (default; generally most accurate)
    - `:nd`: NaiveDelta
    - `:df`: DualFrames

# Returns

- `Vector{Float64}`: reconstructed signal

# Throws

- `ArgumentError`: if `type` is not `:pd`, `:nd`, or `:df`

# See also

[`cwd`](@ref)
"""
function icwd(
    ct::Matrix{Float64};
    wt::T = wavelet(Morlet(2π), β = 2),
    type::Symbol = :pd
)::Vector{Float64} where {T <: CWT}

    _check_var(type, [:nd, :pd, :df], "type")

    # transpose back from (scales × samples) to the (samples × scales) layout expected by ContinuousWavelets.icwt
    ct_t = Matrix(ct')

    s = if type === :pd
        ContinuousWavelets.icwt(ct_t, wt, PenroseDelta())
    elseif type === :nd
        ContinuousWavelets.icwt(ct_t, wt, NaiveDelta())
    elseif type === :df
        real.(ContinuousWavelets.icwt(ct_t, wt, DualFrames()))
    end

    return vec(s)

end
