export dirinrg

"""
    dirinrg(s)

Calculate Dirichlet energy.

# Arguments

- `s::AbstractVector`

# Returns

- `dn::Float64`
"""
function dirinrg(s::AbstractVector)::Float64

    dn = norm(diff(s), 2)^2

    return dn

end

"""
    dirinrg(s)

Calculate Dirichlet energy.

# Arguments

- `s::AbstractArray`

# Returns

- `dn::Matrix{Float64}`
"""
function dirinrg(s::AbstractArray)::Matrix{Float64}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    hd = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            hd[ch_idx, ep_idx] = @views dirinrg(s[ch_idx, :, ep_idx])
        end
    end

    return hd

end

"""
    dirinrg(obj; <keyword arguments>)

Calculate Dirichlet energy.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns

- `dn::Matrix{Float64}`
"""
function dirinrg(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Matrix{Float64}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")

    dn = @views dirinrg(obj.data[ch, :, :])

    return dn

end
