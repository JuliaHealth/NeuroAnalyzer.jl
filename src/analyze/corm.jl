export corm

"""
    corm(s; <keyword arguments>)

Calculate correlation matrix of `s * s'`.

# Arguments

- `s::AbstractVector`
- `norm::Bool`: normalize correlation matrix

# Returns

- `cm::Matrix{Float64}`: correlation matrix
"""
function corm(s::AbstractVector; norm::Bool=false)::Matrix{Float64}

    # channels-vs-channels
    if na_gpu === :cuda
        cm = cor(Matrix(CuVector(s) * CuVector(s)'))
    elseif na_gpu === :amdgpu
        cm = cor(Matrix(ROCVector(s) * ROCVector(s)'))
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

- `cm::Matrix{Float64}`: correlation matrix
"""
function corm(s1::AbstractVector, s2::AbstractVector; norm::Bool=false)::Matrix{Float64}

    @assert length(s1) == length(s2) "s1 and s2 must have the same length."

    # channels-vs-channels
    if na_gpu === :cuda
        cm = Matrix(cor(CuVector(s1) * CuVector(s2)'))
    elseif na_gpu === :amdgpu
        cm = cor(Matrix(ROCVector(s1) * ROCVector(s2)'))
    else
        cm = cor(s1 * s2')
    end

    # normalize
    norm && (cm = m_norm(cm))

    return cm

end


"""
    corm(s; <keyword arguments>)

Calculate corelation matrix of channels × time points matrix.

# Arguments

- `s::AbstractMatrix`
- `norm::Bool=false`: normalize correlation

# Returns

- `cm::Matrix{Float64}`: corelation matrix
"""
function corm(s::AbstractMatrix; norm::Bool=false)::Matrix{Float64}

    # channels-vs-channels
    cm = Statistics.cor(s')

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

- `cm::Array{Float64, 3}`: correlation matrix for each epoch
"""
function corm(s::AbstractArray; norm::Bool=false)::Array{Float64, 3}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    cm = zeros(ch_n, ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        @views @inbounds cm[:, :, ep_idx] = corm(s[:, :, ep_idx], norm=norm)
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

- `cm::Array{Float64, 3}`: correlation matrix for each epoch
"""
function corm(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, norm::Bool=false)::Array{Float64, 3}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")

    cm = corm(obj.data[ch, :, :], norm=norm)

    return cm

end