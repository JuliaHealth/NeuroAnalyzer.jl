export dw_trans
export idw_trans

"""
    dw_trans(s; <keyword arguments>)

Perform discrete wavelet transformation (DWT).

# Arguments

- `s::AbstractVector`
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: transformation type:
    - `:sdwt`: Stationary Wavelet Transforms
    - `:acdwt`: Autocorrelation Wavelet Transforms
- `l::Int64=0`: number of levels, default maximum number of levels available or total transformation

# Returns

- `dt::Array{Float64, 2}`: DWT coefficients cAl, cD1, ..., cDl (by rows)
"""
function dw_trans(s::AbstractVector; wt::T, type::Symbol, l::Int64=0) where {T <: DiscreteWavelet}
    _check_var(type, [:sdwt, :acdwt], "type")

    @assert l <= maxtransformlevels(s) "l must be ≤ $(maxtransformlevels(s))."

    if l == 0
        l = maxtransformlevels(s)
        _info("Calculating DWT using maximum level: $l")
    end

    if type === :sdwt
        dwt_coefs = sdwt(s, wt, l)
    elseif type === :acdwt
        dwt_coefs = acdwt(s, wt, l)
    end

    dwt_c = zeros(size(dwt_coefs, 2), size(dwt_coefs, 1))
    dwt_c[1, :] = @view dwt_coefs[:, 1]
    @inbounds for idx in 2:(l + 1)
        dwt_c[idx, :] = @views dwt_coefs[:, (end - idx + 2)]
    end

    return dwt_c

end

"""
    dw_trans(s; <keyword arguments>)

Perform discrete wavelet transformation (DWT).

# Arguments

- `s::AbstractArray`
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: transformation type:
    - `:sdwt`: Stationary Wavelet Transforms
    - `:acdwt`: Autocorrelation Wavelet Transforms
- `l::Int64=0`: number of levels, default is maximum number of levels available or total transformation

# Returns

- `dt::Array{Float64, 4}`: DWT coefficients cAl, cD1, ..., cDl (by rows)
"""
function dw_trans(s::AbstractArray; wt::T, type::Symbol, l::Int64=0) where {T <: DiscreteWavelet}

    if l == 0
        l = maxtransformlevels(s[1, :, 1])
        _info("Calculating DWT using maximum level: $l")
    end

    ch_n, ep_len, ep_n = size(s)

    dt = zeros(ch_n, (l + 1), ep_len, ep_n)
    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            dt[ch_idx, :, :, ep_idx] = @views dw_trans(s[ch_idx, :, ep_idx], wt=wt, type=type, l=l)
        end
    end

    return dt

end

"""
    dw_trans(obj; <keyword arguments>)

Perform discrete wavelet transformation (DWT).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: list of channels
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: transformation type:
    - `:sdwt`: Stationary Wavelet Transforms
    - `:acdwt`: Autocorrelation Wavelet Transforms
- `l::Int64=0`: number of levels, default is maximum number of levels available or total transformation

# Returns

- `dt::Array{Float64, 4}`: DWT coefficients cAl, cD1, ..., cDl (by rows)
"""
function dw_trans(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, wt::T, type::Symbol, l::Int64=0) where {T <: DiscreteWavelet}

    ch = _ch_idx(obj, ch)
    dt = @views dw_trans(obj.data[ch, :, :], wt=wt, type=type, l=l)

    return dt

end

"""
    idw_trans(dwt_coefs; <keyword arguments>)

Perform inverse discrete wavelet transformation (iDWT) of the `dwt_coefs`.

# Arguments

- `dwt_coefs::AbstractArray`: DWT coefficients cAl, cD1, ..., cDl (by rows)
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`, see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: transformation type:
    - `:sdwt`: Stationary Wavelet Transforms
    - `:acdwt`: Autocorrelation Wavelet Transforms

# Returns

- `s_new::Vector{Float64}`: reconstructed signal
"""
function idw_trans(dwt_coefs::AbstractArray; wt::T, type::Symbol) where {T <: DiscreteWavelet}

    _check_var(type, [:sdwt, :acdwt], "type")

    # reconstruct array of DWT coefficients as returned by Wavelets.jl functions
    dwt_c = zeros(size(dwt_coefs, 2), size(dwt_coefs, 1))
    dwt_c[:, 1] = @view dwt_coefs[1, :]
    @inbounds for idx in 2:size(dwt_coefs, 1)
        dwt_c[:, idx] = @views dwt_coefs[(end - idx + 2), :]
    end

    if type === :sdwt
        return isdwt(dwt_c, wt)
    elseif type === :acdwt
        return iacdwt(dwt_c, wt)
    end

end
