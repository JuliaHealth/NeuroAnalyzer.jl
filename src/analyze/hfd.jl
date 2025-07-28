export hfd

"""
    hfd(s)

Calculate the Higuchi fractal dimension (Higuchi, 1988).

# Arguments

- `s::AbstractVector`

# Returns

- `hd::Float64`
"""
function hfd(s::AbstractVector)::Float64

    hd = higuchi_dim(s)

    return hd

end

"""
    hfd(s)

Calculate the Higuchi fractal dimension (Higuchi, 1988).

# Arguments

- `s::AbstractArray`

# Returns

- `hd::Matrix{Float64}`
"""
function hfd(s::AbstractArray)::Matrix{Float64}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    hd = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            hd[ch_idx, ep_idx] = @views hfd(s[ch_idx, :, ep_idx])
        end
    end

    return hd

end

"""
    hfd(obj; <keyword arguments>)

Calculate the Higuchi fractal dimension (Higuchi, 1988).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns

- `hd::Matrix{Float64}`
"""
function hfd(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Matrix{Float64}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")

    hd = @views hfd(obj.data[ch, :, :])

    return hd

end
