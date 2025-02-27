export corm

"""
    corm(s; <keyword arguments>)

Calculate correlation matrix of `s * s'`.

# Arguments

- `s::AbstractVector`
- `norm::Bool`: normalize correlation matrix

# Returns

- `cm::Matrix{Float64}`
"""
function corm(s::AbstractVector; norm::Bool=false)::Matrix{Float64}

    # channels-vs-channels
    if CUDA.functional() && use_cuda
        cm = Matrix(cor(CuVector(s) * CuVector(s)'))
    else
        cm = cor(s * s')
    end

    # normalize
    norm && (cm = m_norm(cm))

    return cm

end

"""
    corm(s1, s2; <keyword arguments>)

Calculate correlation matrix of `s1 * s2'`.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`
- `norm::Bool`: normalize correlation matrix

# Returns

- `cm::Matrix{Float64}`
"""
function corm(s1::AbstractVector, s2::AbstractVector; norm::Bool=false)::Matrix{Float64}

    @assert length(s1) == length(s2) "s1 and s2 must have the same length."

    # channels-vs-channels
    if CUDA.functional() && use_cuda
        cm = Matrix(cor(CuVector(s1) * CuVector(s2)'))
    else
        cm = cor(s1 * s2')
    end

    # normalize
    norm && (cm = m_norm(cm))

    return cm

end

"""
    corm(s; <keyword arguments>)

Calculate correlation matrix.

# Arguments

- `s::AbstractArray`
- `norm::Bool=false`: normalize covariance

# Returns

- `cm::Array{Float64, 4}`
"""
function corm(s::AbstractArray; norm::Bool=false)::Array{Float64, 4}

    _chk3d(s)
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
                @views @inbounds cm[:, :, s_idx, ep_idx] = corm(s[:, s_idx, ep_idx], norm=norm)

                # update progress bar
                progress_bar && next!(progbar)
            end
            CUDA.synchronize()
        else
            Threads.@threads :static for s_idx in 1:ep_len
                @views @inbounds cm[:, :, s_idx, ep_idx] = corm(s[:, s_idx, ep_idx], norm=norm)

                # update progress bar
                progress_bar && next!(progbar)
            end
        end
    end

    return cm

end

"""
     corm(obj; <keyword arguments>)

Calculate correlation matrix.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}: list of channels
- `norm::Bool=true`: normalize matrix

# Returns

- `cm::Array{Float64, 4}`: correlation matrix for each epoch
"""
function corm(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, norm::Bool=false)::Array{Float64, 4}

    ch = get_channel(obj, ch=ch)

    cm = corm(obj.data[ch, :, :], norm=norm)

    return cm

end