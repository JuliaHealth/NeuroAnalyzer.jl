export covm

"""
   covm(s; norm=true)

Calculate covariance matrix of `s * s'`.

# Arguments

- `s::AbstractVector`
- `norm::Bool=false`: normalize covariance

# Returns

- `cm::Matrix{Float64}`
"""
function covm(s::AbstractVector; norm::Bool=false)

    # channels-vs-channels
    if CUDA.functional() && use_cuda
        cm = Matrix(cov(CuVector(s) * CuVector(s)'))
    else
        cm = cov(s * s')
    end

    # normalize
    norm == true && (cm = m_norm(cm))

    return cm

end

"""
   covm(s1, s2; norm=true)

Calculate covariance matrix of `s1 * s2'`.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`
- `norm::Bool=false`: normalize covariance

# Returns

- `cm::Matrix{Float64}`
"""
function covm(s1::AbstractVector, s2::AbstractVector; norm::Bool=false)

    @assert length(s1) == length(s2) "s1 and s2 must have the same length."

    # channels-vs-channels
    if CUDA.functional() && use_cuda
        cm = Matrix(cov(CuVector(s1) * CuVector(s2)'))
    else
        cm = cov(s1 * s2')
    end

    # normalize
    norm == true && (cm = m_norm(cm))

    return cm

end

"""
   covm(s; norm=true)

Calculate covariance matrix.

# Arguments

- `s::AbstractArray`
- `norm::Bool=false`: normalize covariance

# Returns

- `cm::Matrix{Float64}`
"""
function covm(s::AbstractArray; norm::Bool=false)

    ch_n = size(s, 1)
    ep_len = size(s, 2)
    ep_n = size(s, 3)

    # initialize progress bar
    progress_bar == true && (progbar = Progress(ep_len * ep_n, dt=1, barlen=20, color=:white))

    cm = zeros(ch_n, ch_n, ep_len, ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for s_idx in 1:ep_len
            @views @inbounds cm[:, :, s_idx, ep_idx] = covm(s[:, s_idx, ep_idx], norm=norm)

            # update progress bar
            progress_bar == true && next!(progbar)
        end
    end

    return cm

end

"""
    covm(obj; ch, norm)

Calculate covariance matrix of `signal * signal'`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `norm::Bool=false`: normalize matrix

# Returns

- `cm::Array{Float64, 3}`: covariance matrix for each epoch
"""
function covm(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), norm::Bool=false)

    _check_channels(obj, ch)

    cm = covm(obj.data[ch, :, :], norm=norm)

    return cm

end