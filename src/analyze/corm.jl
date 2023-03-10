export corm

"""
   corm(signal; norm=true)

Calculate correlation matrix of `signal * signal'`.

# Arguments

- `signal::AbstractVector`
- `norm::Bool`: normalize correlation matrix

# Returns

- `cor_mat::Matrix{Float64}`
"""
function corm(signal::AbstractVector; norm::Bool=false)
    # channels-vs-channels
    if CUDA.functional() && use_cuda
        cor_mat = Matrix(cor(CuVector(signal) * CuVector(signal)'))
    else
        cor_mat = cor(signal * signal')
    end

    # normalize
    norm == true && (cor_mat = m_norm(cor_mat))

    return cor_mat
end

"""
   corm(signal1, signal2; norm=true)

Calculate correlation matrix between `signal1` and `signal2`.

# Arguments

- `signal1::AbstractVector`
- `signal2::AbstractVector`
- `norm::Bool`: normalize correlation matrix

# Returns

- `cor_mat::Matrix{Float64}`
"""
function corm(signal1::AbstractVector, signal2::AbstractVector; norm::Bool=false)

    length(signal1) == length(signal2) || throw(ArgumentError("Both signals must be of the same length."))

    # channels-vs-channels
    if CUDA.functional() && use_cuda
        cor_mat = Matrix(cor(CuVector(signal1) * CuVector(signal2)'))
    else
        cor_mat = cor(signal1 * signal2')
    end

    # normalize
    norm == true && (cor_mat = m_norm(cor_mat))

    return cor_mat
end

"""
   corm(signal; norm=true)

Calculate correlation matrix.

# Arguments

- `signal::AbstractArray`
- `norm::Bool=false`: normalize covariance

# Returns

- `cor_mat::Array{Float64, 4}`
"""
function corm(signal::AbstractArray; norm::Bool=false)

    ch_n = size(signal, 1)
    ep_len = size(signal, 2)
    ep_n = size(signal, 3)

    # initialize progress bar
    progress_bar == true && (pb = Progress(ep_len * ep_n, 1))

    cor_mat = zeros(ch_n, ch_n, ep_len, ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for signal_idx in 1:ep_len
            @views @inbounds cor_mat[:, :, signal_idx, ep_idx] = corm(signal[:, signal_idx, ep_idx], norm=norm)

            # update progress bar
            progress_bar == true && next!(pb)
        end
    end
    
    return cor_mat
end

"""
    corm(obj; channel, norm)

Calculate correlation matrix.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `norm::Bool=true`: normalize matrix

# Returns

- `cor_mat::Array{Float64, 3}`: correlation matrix for each epoch
"""
function corm(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), norm::Bool=false)

    _check_channels(obj, channel)

    return corm(obj.data[channel, :, :])

end