export dirinrg

"""
    dirinrg(s)

Calculate Dirichlet energy. It measures the "roughness" of a discrete signal by summing the squared differences between consecutive samples: E(s) = Σ (s[i+1] - s[i])²  =  ‖diff(s)‖²

A smooth, slowly-varying signal has low Dirichlet energy; a noisy or rapidly oscillating signal has high energy.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `dn::Float64`: Dirichlet energy
"""
function dirinrg(s::AbstractVector)::Float64

    # sum(abs2, diff(s)) = Σ(Δs)²
    # equivalent to norm(diff(s), 2)^2 but avoids the sqrt in norm() immediately followed by squaring it back
    dn = sum(abs2, diff(s))

    return dn

end

"""
    dirinrg(s)

Calculate Dirichlet energy. It measures the "roughness" of a discrete signal by summing the squared differences between consecutive samples: E(s) = Σ (s[i+1] - s[i])²  =  ‖diff(s)‖²

A smooth, slowly-varying signal has low Dirichlet energy; a noisy or rapidly oscillating signal has high energy.

# Arguments

- `s::AbstractArray`: signal array (channels × samples × epochs)

# Returns

- `dn::Matrix{Float64}`: Dirichlet energy of shape `(channels, epochs)`
"""
function dirinrg(s::AbstractArray)::Matrix{Float64}

    # validate that the input is a proper 3-D array (channels × samples × epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    dn = zeros(ch_n, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        dn[ch_idx, ep_idx] = dirinrg(@view(s[ch_idx, :, ep_idx]))
    end

    return dn

end

"""
    dirinrg(obj; <keyword arguments>)

Calculate Dirichlet energy. It measures the "roughness" of a discrete signal by summing the squared differences between consecutive samples: E(s) = Σ (s[i+1] - s[i])²  =  ‖diff(s)‖²

A smooth, slowly-varying signal has low Dirichlet energy; a noisy or rapidly oscillating signal has high energy.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `dn::Matrix{Float64}`: Dirichlet energy of shape `(channels, epochs)`
"""
function dirinrg(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Matrix{Float64}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    dn = dirinrg(@view(obj.data[ch, :, :]))

    return dn

end
