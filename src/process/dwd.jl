export dwd
export idwd

"""
    dwd(s; <keyword arguments>)

Perform discrete wavelet decomposition (DWD).

# Arguments

- `s::AbstractVector`
- `wt<:DiscreteWavelet=wavelet(WT.haar)`: discrete wavelet, see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: decomposition type:
    - `:sdwt`: stationary discrete wavelet transform
    - `:acdwt`: discrete autocorrelation wavelet transform
- `l::Int64=maxtransformlevels(s)`: number of levels, default maximum number of levels available or total transformation

# Returns

- `dc::Matrix{Float64}`: DWD coefficients (by rows)
"""
function dwd(s::AbstractVector; wt::T=wavelet(WT.haar), type::Symbol, l::Int64=maxtransformlevels(s))::Matrix{Float64} where {T <: DiscreteWavelet}

    _check_var(type, [:sdwt, :acdwt], "type")

    @assert l <= maxtransformlevels(s) "l must be ≤ $(maxtransformlevels(s))."

    if type === :sdwt
        dc = swpd(s, wt, l)
    elseif type === :acdwt
        dc = acwpd(s, wt, l)
    end

    return Matrix(dc')

end

"""
    dwd(s; <keyword arguments>)

Perform discrete wavelet transformation (DWT).

# Arguments

- `s::AbstractArray`
- `wt<:DiscreteWavelet=wavelet(WT.haar)`: discrete wavelet, see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: transformation type:
    - `:sdwt`: stationary discrete wavelet transform
    - `:acdwt`: discrete autocorrelation wavelet transform
- `l::Int64=maxtransformlevels(s[1, :, 1])`: number of levels, default is the maximum number of levels available or total transformation

# Returns

- `dc::Array{Float64, 4}`: DWD coefficients
"""
function dwd(s::AbstractArray; wt::T=wavelet(WT.haar), type::Symbol, l::Int64=maxtransformlevels(s[1, :, 1]))::Array{Float64, 4} where {T <: DiscreteWavelet}

    _chk3d(s)

    ch_n, ep_len, ep_n = size(s)

    @assert (ch_n * (1 + sum(2 .^ (1:l))) * ep_len * ep_n) / 1048576 < _fmem() "Not enough memory."
    dc = zeros(ch_n, (1 + sum(2 .^ (1:l))), ep_len, ep_n)

    _log_off()
    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            dc[ch_idx, :, :, ep_idx] = @views dwd(s[ch_idx, :, ep_idx], wt=wt, type=type, l=l)
        end
    end
    _log_on()

    return dc

end

"""
    dwd(obj; <keyword arguments>)

Perform discrete wavelet transformation (DWT).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `wt<:DiscreteWavelet=wavelet(WT.haar)`: discrete wavelet, see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: transformation type:
    - `:sdwt`: stationary discrete wavelet transform
    - `:acdwt`: discrete autocorrelation wavelet transform
- `l::Int64=0`: number of levels, default is the maximum number of levels available or total transformation

# Returns

- `dc::Array{Float64, 4}`: DWD coefficients
"""
function dwd(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, wt::T=wavelet(WT.haar), type::Symbol, l::Int64=0)::Array{Float64, 4} where {T <: DiscreteWavelet}

    if l == 0
        l = maxtransformlevels(obj.data[1, :, 1])
        _info("Calculating DWD using maximum level: $l")
    end

    ch = get_channel(obj, ch=ch)
    dc = @views dwd(obj.data[ch, :, :], wt=wt, type=type, l=l)

    return dc

end

"""
    idwd(dc; <keyword arguments>)

Perform inverse discrete wavelet transformation (iDWT).

# Arguments

- `dc::Matrix{Float64}`: DWD coefficients (by rows)
- `wt<:DiscreteWavelet=wavelet(WT.haar)`: discrete wavelet, see Wavelets.jl documentation for the list of available wavelets
- `type::Symbol`: transformation type:
    - `:sdwt`: average-based stationary discrete wavelet transform
    - `:acdwt`: discrete autocorrelation wavelet transform
- `c::Union{Int64, Vector{Int64}, AbstractRange}=axes(dc, 1)`: which coefficients are used for reconstruction, default is all coefficients

# Returns

- `s::AbstractArray`: reconstructed signal
"""
function idwd(dc::Matrix{Float64}; wt::T=wavelet(WT.haar), type::Symbol, c::Union{Int64, Vector{Int64}, AbstractRange}=axes(dc, 1))::AbstractArray where {T <: DiscreteWavelet}

    _check_var(type, [:sdwt, :acdwt], "type")

    if length(c) > 1
        c = sort(c)
        for idx in c
            _in(idx, (1, size(dc, 1)))
        end
    else
        _in(c[1], (1, size(dc, 1)), "c")
    end

    if type === :sdwt
        s = @views isdwt(dc[c, :]', wt)
    elseif type === :acdwt
        s = @views iacdwt(dc[c, :]', wt)
    end

    return s

end
