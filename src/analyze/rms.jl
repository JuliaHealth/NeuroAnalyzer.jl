export rms

"""
    rms(s)

Calculate Root Mean Square.

# Arguments

- `s::AbstractVector`

# Returns

- `rms::Float64`
"""
function rms(s::AbstractVector)

    # rms = sqrt(mean(s.^2))    
    return norm(s) / sqrt(length(s))

end

"""
    rms(s)

Calculate Root Mean Square.

# Arguments

- `s::AbstractArray`

# Returns

- `r::Array{Float64, 2}`
"""
function rms(s::AbstractArray)

    ch_n = size(s, 1)
    ep_n = size(s, 3)
    
    r = zeros(ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            r[ch_idx, ep_idx] = @views rms(s[ch_idx, :, ep_idx])
        end
    end

    return r

end

"""
    rms(obj; channel)

Calculate Root Mean Square.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `r::Array{Float64, 2}`
"""
function rms(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj))

    _check_channels(obj, ch)

    r = @views rms(obj.data[ch, :, :])

    return r

end
