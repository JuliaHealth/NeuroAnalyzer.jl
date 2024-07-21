export covm

"""
    covm(s; <keyword arguments>)

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
    norm && (cm = m_norm(cm))

    return cm

end

"""
    covm(s1, s2; <keyword arguments>)

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
    norm && (cm = m_norm(cm))

    return cm

end

"""
    covm(s; <keyword arguments>)

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
    progress_bar && (progbar = Progress(ep_len * ep_n, dt=1, barlen=20, color=:white))

    cm = zeros(ch_n, ch_n, ep_len, ep_n)
    @inbounds for ep_idx in 1:ep_n
        if use_cuda
            CUDA.synchronize()
            for s_idx in 1:ep_len
                @views @inbounds cm[:, :, s_idx, ep_idx] = covm(s[:, s_idx, ep_idx], norm=norm)

                # update progress bar
                progress_bar && next!(progbar)
            end
            CUDA.synchronize()
        else
            Threads.@threads for s_idx in 1:ep_len
                @views @inbounds cm[:, :, s_idx, ep_idx] = covm(s[:, s_idx, ep_idx], norm=norm)

                # update progress bar
                progress_bar && next!(progbar)
            end
        end
    end

    return cm

end

"""
    covm(obj; <keyword arguments>)

Calculate covariance matrix of `signal * signal'`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}: list of channels
- `norm::Bool=false`: normalize matrix

# Returns

- `cm::Array{Float64, 3}`: covariance matrix for each epoch
"""
function covm(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, norm::Bool=false)

    ch = _ch_idx(obj, ch)

    cm = covm(obj.data[ch, :, :], norm=norm)

    return cm

end