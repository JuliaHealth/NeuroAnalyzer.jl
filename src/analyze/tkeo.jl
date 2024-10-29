export tkeo

"""
    tkeo(s, t; <keyword arguments>)

Calculate Teager-Kaiser energy-tracking operator.

# Arguments

- `s::AbstractVector`: signal
- `t::AbstractVector=collect(1:length(s))`: time points
- `method::Symbol=:pow`:
    - `:pow`: TKEO = x(t)^2 - x(t-1) × x(t+1)
    - `:der`: TKEO = f'(t) - f(t) × f''(t)
    - `:amp`: TKEO = envelope(amplitude)^2

# Returns

- `tk::Vector{Float64}`
"""
function tkeo(s::AbstractVector, t::AbstractVector=collect(1:length(s)); method::Symbol=:pow)::Vector{Float64}

    _check_var(method, [:pow, :der, :amp], "method")

    tk = nothing

    if method === :pow
        tk = zeros(length(s))
        tk[1] = s[1]
        tk[end] = s[end]
        @inbounds for idx in 2:(length(s) - 1)
            tk[idx] = s[idx]^2 - (s[idx - 1] * s[idx + 1])
        end
    elseif method === :der
        d1 = derivative(s)
        d2 = derivative(d1)
        tk = @. d1 - s * d2
    else
        tk = env_up(s, t, d=8).^2
    end

    return tk

end

"""
    tkeo(s, t; <keyword arguments>)

Calculate Teager-Kaiser energy-tracking operator

# Arguments

- `s::AbstractArray`: signal
- `t::AbstractArray=collect(1:length(s))`: time points
- `method::Symbol=:pow`:
    - `:pow`: TKEO = x(t)^2 - x(t-1) × x(t+1)
    - `:der`: TKEO = f'(t) - f(t) × f''(t)
    - `:amp`: TKEO = envelope(amplitude)^2

# Returns

- `tk::Array{Float64, 3}`
"""
function tkeo(s::AbstractArray, t::AbstractVector=collect(1:length(s)); method::Symbol=:pow)::Array{Float64, 3}

    _chk3d(s)

    tk = similar(s)
    @inbounds for ep_idx in axes(s, 3)
        Threads.@threads for ch_idx in axes(s, 1)
            tk[ch_idx, :, ep_idx] = @views tkeo(s[ch_idx, :, ep_idx], t, method=method)
        end
    end

    return tk

end

"""
    tkeo(obj; <keyword arguments>)

Calculate Teager-Kaiser energy-tracking operator.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}: list of channels
- `method::Symbol=:pow`:
    - `:pow`: TKEO = x(t)^2 - x(t-1) × x(t+1)
    - `:der`: TKEO = f'(t) - f(t) × f''(t)
    - `:amp`: TKEO = envelope(amplitude)^2

# Returns

- `tk::Array{Float64, 3}`
"""
function tkeo(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, method::Symbol=:pow)::Array{Float64, 3}

    ch = get_channel(obj, ch=ch)
    tk = @views tkeo(obj.data[ch, :, :], obj.epoch_time, method=method)

    return tk

end
