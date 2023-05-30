export corm

"""
   corm(s; norm=true)

Calculate correlation matrix of `s * s'`.

# Arguments

- `s::AbstractVector`
- `norm::Bool`: normalize correlation matrix

# Returns

- `cm::Matrix{Float64}`
"""
function corm(s::AbstractVector; norm::Bool=false)
    # channels-vs-channels
    if CUDA.functional() && use_cuda
        cm = Matrix(cor(CuVector(s) * CuVector(s)'))
    else
        cm = cor(s * s')
    end

    # normalize
    norm == true && (cm = m_norm(cm))

    return cm
    
end

"""
   corm(s1, s2; norm=true)

Calculate correlation matrix of `s1 * s2'`.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`
- `norm::Bool`: normalize correlation matrix

# Returns

- `cm::Matrix{Float64}`
"""
function corm(s1::AbstractVector, s2::AbstractVector; norm::Bool=false)

    length(s1) == length(s2) || throw(ArgumentError("s1 and s2 must have the same length."))

    # channels-vs-channels
    if CUDA.functional() && use_cuda
        cm = Matrix(cor(CuVector(s1) * CuVector(s2)'))
    else
        cm = cor(s1 * s2')
    end

    # normalize
    norm == true && (cm = m_norm(cm))

    return cm

end

"""
   corm(s; norm=true)

Calculate correlation matrix.

# Arguments

- `s::AbstractArray`
- `norm::Bool=false`: normalize covariance

# Returns

- `cm::Array{Float64, 4}`
"""
function corm(s::AbstractArray; norm::Bool=false)

    ch_n = size(s, 1)
    ep_len = size(s, 2)
    ep_n = size(s, 3)

    # initialize progress bar
    progress_bar == true && (pb = Progress(ep_len * ep_n, 1))

    cm = zeros(ch_n, ch_n, ep_len, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for s_idx in 1:ep_len
            @views @inbounds cm[:, :, s_idx, ep_idx] = corm(s[:, s_idx, ep_idx], norm=norm)

            # update progress bar
            progress_bar == true && next!(pb)
        end
    end
    
    return cm
    
end

"""
    corm(obj; ch, norm)

Calculate correlation matrix.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `norm::Bool=true`: normalize matrix

# Returns

- `cm::Array{Float64, 3}`: correlation matrix for each epoch
"""
function corm(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), norm::Bool=false)

    _check_channels(obj, ch)

    cm = corm(obj.data[ch, :, :], norm=norm)

    return cm
    
end