export corm

"""
    corm(s1, s2; <keyword arguments>)

Computes the channel × channel Pearson correlation matrix.

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector
- `norm::Bool = false`: normalize correlation matrix

# Returns

- `cm::Matrix{Float64}`: 2×2 correlation matrix
"""
function corm(s1::AbstractVector, s2::AbstractVector; norm::Bool = false)::Matrix{Float64}

    @assert length(s1) == length(s2) "s1 and s2 must have the same length."

    # compute the 2×2 channels-vs-channels correlation matrix
    # hcat → n×2; cor → 2×2
    cm = Statistics.cor(hcat(s1, s2))

    # normalize if requested
    norm && (cm = m_norm(cm))

    return cm

end


"""
    corm(s; <keyword arguments>)

Calculate corelation matrix of channels × time points matrix.

# Arguments

- `s::AbstractMatrix`: signal matrix (channels × samples)
- `norm::Bool=false`: normalize correlation matrix

# Returns

- `cm::Matrix{Float64}`: correlation matrix (channels × channels)
"""
function corm(s::AbstractMatrix; norm::Bool = false)::Matrix{Float64}

    # transpose so columns are channels
    # Statistics.cor then returns the channels × channels Pearson correlation matrix
    cm = Statistics.cor(s')

    # normalize if requested
    norm && (cm = m_norm(cm))

    return cm

end

"""
    corm(s; <keyword arguments>)

Calculate correlation matrix for each epoch of a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array (channels × samples × epochs)
- `norm::Bool=false`: normalize correlation matrix

# Returns

- `cm::Array{Float64, 3}`: correlation matrix for each epoch, shape `(channels, channels, epochs)`
"""
function corm(s::AbstractArray; norm::Bool = false)::Array{Float64, 3}

    # validate that the input is a proper 3-D array (channels × samples × epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    cm = zeros(ch_n, ch_n, ep_n)

    # calculate over epochs
    @inbounds Threads.@threads for ep_idx in 1:ep_n
        cm[:, :, ep_idx] = corm(@view(s[:, :, ep_idx]), norm = norm)
    end

    return cm

end

"""
     corm(obj; <keyword arguments>)

Calculate correlation matrix.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}: channel name(s)
- `norm::Bool=true`: normalize matrix

# Returns

- `cm::Array{Float64, 3}`: correlation matrix for each epoch, shape `(channels, channels, epochs)`
"""
function corm(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, norm::Bool = false)::Array{Float64, 3}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    cm = corm(@view(obj.data[ch, :, :]), norm = norm)

    return cm

end
