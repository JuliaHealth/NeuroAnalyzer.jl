export tkeo

"""
    tkeo(s, t; <keyword arguments>)

Calculate Teager-Kaiser energy-tracking operator for a 1-D signal vector.

# Arguments

- `s::AbstractVector`: signal vector
- `t::AbstractVector=collect(1:length(s))`: time points
- `method::Symbol=:pow`:
    - `:pow`: TKEO = x(t)^2 - x(t-1) × x(t+1)
    - `:der`: TKEO = f'(t) - f(t) × f''(t)
    - `:amp`: TKEO = envelope(amplitude)^2

# Returns

- `Vector{Float64}`
"""
function tkeo(
        s::AbstractVector,
        t::AbstractVector = collect(1:length(s));
        method::Symbol = :pow,
    )::Vector{Float64}
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
        tk = env_up(s, t; d = 8) .^ 2
    end

    return tk
end

"""
    tkeo(s, t; <keyword arguments>)

Calculate Teager-Kaiser energy-tracking operator for a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `t::AbstractArray=collect(1:length(s))`: time points
- `method::Symbol=:pow`:
    - `:pow`: TKEO = x(t)^2 - x(t-1) × x(t+1)
    - `:der`: TKEO = f'(t) - f(t) × f''(t)
    - `:amp`: TKEO = envelope(amplitude)^2

# Returns

- `Array{Float64, 3}`
"""
function tkeo(
        s::AbstractArray,
        t::AbstractVector = collect(1:length(s));
        method::Symbol = :pow,
    )::Array{Float64, 3}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    # pre-allocate output
    tk = similar(s, Float64)

    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        tk[ch_idx, :, ep_idx] = tkeo(@view(s[ch_idx, :, ep_idx]), t, method = method)
    end

    return tk
end

"""
    tkeo(obj; <keyword arguments>)

Calculate Teager-Kaiser energy-tracking operator for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}: list of channels
- `method::Symbol=:pow`:
    - `:pow`: TKEO = x(t)^2 - x(t-1) × x(t+1)
    - `:der`: TKEO = f'(t) - f(t) × f''(t)
    - `:amp`: TKEO = envelope(amplitude)^2

# Returns

- `Array{Float64, 3}`
"""
function tkeo(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex},
        method::Symbol = :pow,
    )::Array{Float64, 3}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ? get_channel(obj; ch = ch, exclude = "bad") :
                       get_channel(obj; ch = ch, exclude = "")

    return tkeo(@view(obj.data[ch, :, :]), obj.epoch_time; method = method)
end
