export tkeo

"""
    tkeo(s)

Calculate Teager-Kaiser energy-tracking operator: y(t) = x(t)^2 - x(t-1) × x(t+1)

# Arguments

- `s::AbstractVector`

# Returns

- `t::Vector{Float64}`
"""
function tkeo(s::AbstractVector)

    t = zeros(length(s))
    t[1] = s[1]
    t[end] = s[end]

    @inbounds @simd for idx in 2:(length(s) - 1)
        t[idx] = s[idx]^2 - (s[idx - 1] * s[idx + 1])
    end

    return t

end

"""
    tkeo(s; channel)

Calculate Teager-Kaiser energy-tracking operator: y(t) = x(t)^2 - x(t-1) × x(t+1)

# Arguments

- `s::AbstractArray`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `t::Array{Float64, 3}`
"""
function tkeo(s::AbstractArray)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    t = similar(s)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            t[ch_idx, :, ep_idx] = @views tkeo(s[ch_idx, :, ep_idx])
        end
    end

    return t

end


"""
    tkeo(obj; channel)

Calculate Teager-Kaiser energy-tracking operator: y(t) = x(t)^2 - x(t-1) × x(t+1)

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `t::Array{Float64, 3}`
"""
function tkeo(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, ch)

    t = @views tkeo(obj.data[ch, :, :])

    return t

end
