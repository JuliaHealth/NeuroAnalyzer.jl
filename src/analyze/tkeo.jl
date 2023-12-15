export tkeo

"""
    tkeo(s; method)

Calculate Teager-Kaiser energy-tracking operator.

# Arguments

- `s::AbstractVector`
- `method::Symbol=:pow`:
    - `:pow`: y(t) = x(t)^2 - x(t-1) × x(t+1)
    - `:der`: y(t) = f'(t) - f(t) × f''(t)

# Returns

- `t::Vector{Float64}`
"""
function tkeo(s::AbstractVector; method::Symbol=:pow)

    _check_var(method, [:pow, :der], "method")

    if method === :pow
        t = zeros(length(s))
        t[1] = s[1]
        t[end] = s[end]

        @inbounds @simd for idx in 2:(length(s) - 1)
            t[idx] = s[idx]^2 - (s[idx - 1] * s[idx + 1])
        end
        
        return t
    
    else

        d1 = derivative(s)
        d2 = derivative(d1)
        t = @. d1 - s * d2

        return t

    end

end

"""
    tkeo(s; method)

Calculate Teager-Kaiser energy-tracking operator

# Arguments

- `s::AbstractArray`
- `method::Symbol=:pow`:
    - `:pow`: y(t) = x(t)^2 - x(t-1) × x(t+1)
    - `:der`: y(t) = f'(t) - f(t) × f''(t)

# Returns

- `t::Array{Float64, 3}`
"""
function tkeo(s::AbstractArray; method::Symbol=:pow)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    t = similar(s)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            t[ch_idx, :, ep_idx] = @views tkeo(s[ch_idx, :, ep_idx], method=method)
        end
    end

    return t

end


"""
    tkeo(obj; channel, method)

Calculate Teager-Kaiser energy-tracking operator.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `method::Symbol=:pow`:
    - `:pow`: y(t) = x(t)^2 - x(t-1) × x(t+1)
    - `:der`: y(t) = f'(t) - f(t) × f''(t)

# Returns

- `t::Array{Float64, 3}`
"""
function tkeo(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), method::Symbol=:pow)

    _check_channels(obj, ch)

    t = @views tkeo(obj.data[ch, :, :], method=method)

    return t

end
