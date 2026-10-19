export cwd
export icwd

"""
    _cwd(s; <keyword arguments>)

Perform continuous wavelet decomposition (CWD) on a signal vector.

Returns the CWT coefficient matrix (real part, or complex coefficients). Each row corresponds to one wavelet scale (frequency band); each column corresponds to one time point.

# Arguments

- `s::AbstractVector`: signal vector
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `complex::Bool=false`: return complex coefficients (required for an accurate inverse transform with complex wavelets, e.g. Morlet)

# Returns

- `Matrix{Float64}` or `Matrix{ComplexF64}`: CWT coefficient matrix of shape `(n_scales, length(s))`
"""
function _cwd(
    s::AbstractVector;
    wt::T = wavelet(Morlet(2π); β = 2),
    complex::Bool = false,
)::Union{Matrix{Float64}, Matrix{ComplexF64}} where {T <: CWT}
    # ContinuousWavelets.cwt returns (samples × scales); transpose to (scales × samples)
    c = ContinuousWavelets.cwt(s, wt)
    return complex ? permutedims(ComplexF64.(c)) : permutedims(real.(c))
end

"""
    _cwd(s; <keyword arguments>)

Perform continuous wavelet decomposition on a 3-D signal array.

Applies `_cwd(::AbstractVector)` to every channel × epoch slice in parallel and returns a 4-D coefficient array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `complex::Bool=false`: return complex coefficients

# Returns

- `Array{Float64, 4}` or `Array{ComplexF64, 4}`: CWT coefficients of shape `(channels, n_scales, samples, epochs)`
"""
function _cwd(
    s::AbstractArray;
    wt::T = wavelet(Morlet(2π); β = 2),
    complex::Bool = false,
)::Union{Array{Float64, 4}, Array{ComplexF64, 4}} where {T <: CWT}
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels, epoch length and number of epochs
    ch_n, ep_len, ep_n = size(s)

    # ContinuousWavelets.jl logs warnings for every slice; logging is restored on errors
    _log_off()
    try
        # dry run on the first slice to determine the number of wavelet scales
        n_scales = size(ContinuousWavelets.cwt(s[1, :, 1], wt), 2)

        # pre-allocate output
        ct = zeros(complex ? ComplexF64 : Float64, ch_n, n_scales, ep_len, ep_n)

        # calculate over channel and epochs
        @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
            ch_idx, ep_idx = idx[1], idx[2]
            ct[ch_idx, :, :, ep_idx] =
                NeuroAnalyzer._cwd(@view(s[ch_idx, :, ep_idx]); wt = wt, complex = complex)
        end

        return ct
    finally
        _log_on()
    end
end

"""
    cwd(obj; <keyword arguments>)

Perform continuous wavelet decomposition on selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `wt<:CWT=wavelet(Morlet(2π), β=2)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets
- `complex::Bool=false`: return complex coefficients (required for an accurate inverse transform `icwd()` with complex wavelets, e.g. Morlet); `false`: real part

# Returns

- `Array{Float64, 4}` or `Array{ComplexF64, 4}`: CWT coefficients of shape `(channels, n_scales, samples, epochs)`
"""
function cwd(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    wt::T = wavelet(Morlet(2π); β = 2),
    complex::Bool = false,
)::Union{Array{Float64, 4}, Array{ComplexF64, 4}} where {T <: CWT}
    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    return _cwd(@view(obj.data[ch, :, :]); wt = wt, complex = complex)
end

"""
    icwd(ct; <keyword arguments>)

Perform inverse continuous wavelet transformation (iCWT).

Reconstructs the original signal from a CWT coefficient matrix produced by [`cwd`](@ref).

# Arguments

- `ct::AbstractMatrix{<:Number}`: CWT coefficient matrix of shape `(n_scales, samples)` (as returned by `cwd`); with complex wavelets (e.g. Morlet) the complex coefficients (`cwd(...; complex = true)`) are required for an accurate reconstruction
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`: continuous wavelet; must match the one used in the forward transform
- `type::Symbol=:pd`: reconstruction method:
    - `:pd`: PenroseDelta (default; generally most accurate)
    - `:nd`: NaiveDelta
    - `:df`: DualFrames

# Returns

- `Vector{Float64}`: reconstructed signal
"""
function icwd(
    ct::AbstractMatrix{<:Number};
    wt::T = wavelet(Morlet(2π); β = 2),
    type::Symbol = :pd,
)::Vector{Float64} where {T <: CWT}
    # validate
    _check_var(type, [:nd, :pd, :df], "type")

    # transpose back from (scales × samples) to the (samples × scales) layout expected by ContinuousWavelets.icwt
    ct_t = permutedims(ct)

    s = if type === :pd
        ContinuousWavelets.icwt(ct_t, wt, PenroseDelta())
    elseif type === :nd
        ContinuousWavelets.icwt(ct_t, wt, NaiveDelta())
    elseif type === :df
        ContinuousWavelets.icwt(ct_t, wt, DualFrames())
    end

    return real.(vec(s))
end
