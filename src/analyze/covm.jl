export covm

"""
    covm(s; <keyword arguments>)

A single channel has no cross-channel covariance; the 1×1 result is simply the signal's own variance.

# Arguments

- `s::AbstractVector`: signal vector
- `norm::Bool=false`: normalize covariance matrix matrix

# Returns

- `cm::Matrix{Float64}`: 1×1 covariance matrix (`[var(s)]`)
"""
function covm(s::AbstractVector; norm::Bool = false)::Matrix{Float64}

    # 1×1 variance matrix
    cm = fill(Statistics.var(s), 1, 1)

    # normalize if requested
    norm && (cm = m_norm(cm))

    return cm

end

"""
    covm(s1, s2; <keyword arguments>)

Calculate covariance matrix of two signals.

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector
- `norm::Bool=false`: normalize covariance matrix

# Returns

- `cm::Matrix{Float64}`: 2×2 covariance matrix
"""
function covm(s1::AbstractVector, s2::AbstractVector; norm::Bool = false)::Matrix{Float64}

    @assert length(s1) == length(s2) "s1 and s2 must have the same length."

    # compute the 2×2 channels-vs-channels covariance matrix
    # hcat → n×2; cor → 2×2
    cm = Statistics.cov(hcat(s1, s2))

    # normalize if requested
    norm && (cm = m_norm(cm))

    return cm

end

"""
    covm(s; <keyword arguments>)

Calculate covariance matrix of a matrix.

# Arguments

- `s::AbstractMatrix`: signal matrix (channels × samples)
- `norm::Bool=false`: normalize covariance matrix

# Returns

- `cm::Matrix{Float64}`: covariance matrix of shape (channels × channels)
"""
function covm(s::AbstractMatrix; norm::Bool = false)::Matrix{Float64}

    # transpose so columns are channels
    # Statistics.cov then returns the channels × channels covariance matrix
    cm = Statistics.cov(s')

    # normalize if requested
    norm && (cm = m_norm(cm))

    return cm

end

"""
    covm(s; <keyword arguments>)

Calculate covariance matrix of an array.

# Arguments

- `s::AbstractArray`: signal array (channels × samples × epochs)
- `norm::Bool=false`: normalize covariance matrix

# Returns

- `cm::Array{Float64, 3}`: covariance matrix of shape (channels × channels × epochs)
"""
function covm(s::AbstractArray; norm::Bool = false)::Array{Float64, 3}

    # validate that the input is a proper 3-D array (channels × samples × epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of vectors
    ep_n = size(s, 3)

    # pre-allocate output
    cm = zeros(ch_n, ch_n, ep_n)

    # calculate over epochs
    @inbounds Threads.@threads :dynamic for ep_idx in 1:ep_n
        cm[:, :, ep_idx] = covm(@view(s[:, :, ep_idx]), norm = norm)
    end

    return cm

end

"""
    covm(obj; <keyword arguments>)

Calculate covariance matrix between all channel pairs within a single object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}: channel name(s)
- `norm::Bool=false`: normalize matrix

# Returns

- `cm::Array{Float64, 3}`: covariance matrix of shape (channels × channels × epochs)
"""
function covm(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, norm::Bool = false)::Array{Float64, 3}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    cm = covm(@view(obj.data[ch, :, :]), norm = norm)

    return cm

end
