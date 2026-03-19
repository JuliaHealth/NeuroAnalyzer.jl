export dwd
export idwd

"""
    dwd(s; <keyword arguments>)

Perform discrete wavelet decomposition (DWD).

Returns the decomposition coefficient matrix. Each row corresponds to one subspace node; each column corresponds to one time sample.

# Arguments

- `s::AbstractVector`: signal vector
- `wt<:DiscreteWavelet=wavelet(WT.haar)`: discrete wavelet; see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: decomposition type:
    - `:sdwt`: stationary discrete wavelet transform
    - `:acdwt`: discrete autocorrelation wavelet transform
- `l::Int64=maxtransformlevels(s)`: number of decomposition levels; must be ≥ 1 and ≤ `maxtransformlevels(s)`

# Returns

- `Matrix{Float64}`: DWD coefficient matrix of shape `((1 + Σ 2^k for k=1..l), length(s))`

# Throws

- `ArgumentError`: if `type` is invalid or `l > maxtransformlevels(s)`

# See also

[`idwd`](@ref), [`dwd(::AbstractArray)`](@ref)
"""
function dwd(
    s::AbstractVector;
    wt::T=wavelet(WT.haar),
    type::Symbol,
    l::Int64=maxtransformlevels(s),
)::Matrix{Float64} where {T <: DiscreteWavelet}

    _check_var(type, [:sdwt, :acdwt], "type")

    l <= maxtransformlevels(s) || throw(ArgumentError("l must be ≤ $(maxtransformlevels(s))."))

    if type === :sdwt
        dc = swpd(s, wt, l)
    elseif type === :acdwt
        dc = acwpd(s, wt, l)
    end

    return Matrix(dc')

end

"""
    dwd(s; <keyword arguments>)

Perform discrete wavelet decomposition (DWD).

Returns the decomposition coefficient matrix. Each row corresponds to one subspace node; each column corresponds to one time sample.

# Arguments

- `s::AbstractArray`: signal array, shape `(channels, samples, epochs)`
- `wt<:DiscreteWavelet=wavelet(WT.haar)`: discrete wavelet; see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: transformation type:
    - `:sdwt`: stationary discrete wavelet transform
    - `:acdwt`: discrete autocorrelation wavelet transform
- `l::Int64=maxtransformlevels(s)`: number of decomposition levels; must be ≥ 1 and ≤ `maxtransformlevels(s)`

# Returns

- `Array{Float64, 4}`: DWD coefficients of shape `(channels, n_nodes, samples, epochs)`

# Throws

- `ArgumentError`: if `s` is not 3-dimensional or memory would be exceeded

# See also

[`dwd(::AbstractVector)`](@ref), [`dwd(::NeuroAnalyzer.NEURO)`](@ref)
"""
function dwd(
    s::AbstractArray;
    wt::T=wavelet(WT.haar),
    type::Symbol,
    l::Int64=maxtransformlevels(s[1, :, 1])
)::Array{Float64, 4} where {T <: DiscreteWavelet}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels, epoch length and number of epochs
    ch_n, ep_len, ep_n = size(s)

    n_nodes = 1 + sum(2 .^ (1:l))
    # 8 bytes per Float64
    mem_mb  = (ch_n * n_nodes * ep_len * ep_n * 8) / 1_048_576
    mem_mb < _fmem() || throw(ArgumentError("Insufficient memory: need ≈ $(round(mem_mb, digits=1)) MB."))

    dc = zeros(ch_n, n_nodes, ep_len, ep_n)

    _log_off()

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        dc[ch_idx, :, :, ep_idx] = dwd(
            @view(s[ch_idx, :, ep_idx]),
            wt=wt,
            type=type,
            l=l
        )
    end

    _log_on()

    return dc

end

"""
    dwd(obj; <keyword arguments>)

Perform discrete wavelet decomposition on selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `wt<:DiscreteWavelet=wavelet(WT.haar)`: discrete wavelet; see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: transformation type:
    - `:sdwt`: stationary discrete wavelet transform
    - `:acdwt`: discrete autocorrelation wavelet transform
- `l::Int64=maxtransformlevels(s)`: number of decomposition levels; must be ≥ 1 and ≤ `maxtransformlevels(s)`

# Returns

- `Array{Float64, 4}`: DWD coefficients of shape `(channels, n_nodes, samples, epochs)`

# See also

[`idwd`](@ref), [`dwd(::AbstractArray)`](@ref)
"""
function dwd(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    wt::T = wavelet(WT.haar),
    type::Symbol,
    l::Int64 = 0,
)::Array{Float64, 4} where {T <: DiscreteWavelet}

    if l == 0
        l = maxtransformlevels(obj.data[1, :, 1])
        _info("Calculating DWD using maximum level: $l")
    end

    # resolve channel names to integer indices
    ch = get_channel(obj, ch = ch)

    return dwd(@view(obj.data[ch, :, :]); wt=wt, type=type, l=l)

end

"""
    idwd(dc; <keyword arguments>)

Perform inverse discrete wavelet decomposition (iDWD).

Reconstructs a signal from a subset (or all) of the DWD coefficient rows.

# Arguments

- `dc::Matrix{Float64}`: DWD coefficient matrix of shape `(n_nodes, samples)` as returned by [`dwd`](@ref)
- `wt<:DiscreteWavelet=wavelet(WT.haar)`: discrete wavelet; see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: transformation type:
    - `:sdwt`: average-based stationary discrete wavelet transform
    - `:acdwt`: discrete autocorrelation wavelet transform
- `c::Union{Int64, Vector{Int64}, AbstractRange}=axes(dc, 1)`: row indices of coefficients to use for reconstruction; default uses all rows; indices must be in `[1, size(dc, 1)]`

# Returns

- `Vector{Float64}`: reconstructed signal

# Throws

- `ArgumentError`: if `type` is invalid or any index in `c` is out of range

# See also

[`dwd`](@ref)
"""
function idwd(
    dc::Matrix{Float64};
    wt::T = wavelet(WT.haar),
    type::Symbol,
    c::Union{Int64, Vector{Int64}, AbstractRange} = axes(dc, 1),
)::Vector{Float64} where {T <: DiscreteWavelet}

    _check_var(type, [:sdwt, :acdwt], "type")

    # validate and normalize coefficient index selection
    c_vec = c isa Int64 ? [c] : sort(collect(c))
    n = size(dc, 1)
    for idx in c_vec
        (idx >= 1 && idx <= n) || throw(ArgumentError("c index $idx must be in [1, $n]."))
    end

    s = if type === :sdwt
        isdwt(dc[c_vec, :]', wt)
    elseif type === :acdwt
        iacdwt(dc[c_vec, :]', wt)
    end

    return vec(s)
end
