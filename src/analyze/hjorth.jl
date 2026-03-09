export hjorth

# internal helper: Hjorth mobility = sqrt( var(s') / var(s) )
_h_mob(s::AbstractVector)::Float64 = sqrt(var(derivative(s)) / var(s))

"""
    hjorth(s)

Calculate Hjorths parameters: three features derived from a signal and its successive derivatives:

- Activity - total power of the signal (`var(s)`)
- Mobility - estimate of the mean frequency (`sqrt(var(s') / var(s))`)
- Complexity - similarity of the signal shape to a pure sine wave (`Mobility(s') / Mobility(s)`)
where s' = derivative(s).

# Arguments

- `s::AbstractVector`: signal vector

# Returns

Named tuple containing:

- `h_act::Float64`: activity
- `h_mob::Float64`: mobility
- `h_comp::Float64`: complexity
"""
function hjorth(s::AbstractVector)::@NamedTuple{h_act::Float64, h_mob::Float64, h_comp::Float64}

    h_act  = var(s)
    h_mob  = _h_mob(s)
    h_comp = _h_mob(derivative(s)) / h_mob

    return (h_act = h_act, h_mob = h_mob, h_comp = h_comp)

end

"""
    hjorth(s)

Calculate Hjorths parameters: three features derived from a signal and its successive derivatives:

- Activity - total power of the signal (`var(s)`)
- Mobility - estimate of the mean frequency (`sqrt(var(s') / var(s))`)
- Complexity - similarity of the signal shape to a pure sine wave (`Mobility(s') / Mobility(s)`)
where s' = derivative(s).

# Arguments

- `s::AbstractArray`: signal array (channels × samples × epochs)

# Returns

Named tuple containing:

- `h_act::Matrix{Float64}`: activity, shape `(channels, epochs)`
- `h_mob::Matrix{Float64}`: mobility, shape `(channels, epochs)`
- `h_comp::Matrix{Float64}`: complexity, shape `(channels, epochs)`
"""
function hjorth(s::AbstractArray)::@NamedTuple{h_act::Matrix{Float64}, h_mob::Matrix{Float64}, h_comp::Matrix{Float64}}

    # validate that the input is a proper 3-D array (channels × samples × epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    h_act = zeros(ch_n, ep_n)
    h_mob = zeros(ch_n, ep_n)
    h_comp = zeros(ch_n, ep_n)

    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        result = hjorth(@view(s[ch_idx, :, ep_idx]))
        h_act[ch_idx, ep_idx]  = result.h_act
        h_mob[ch_idx, ep_idx]  = result.h_mob
        h_comp[ch_idx, ep_idx] = result.h_comp
    end

    return (h_act = h_act, h_mob = h_mob, h_comp = h_comp)

end

"""
    hjorth(obj; <keyword arguments>)

Calculate Hjorths parameters: three features derived from a signal and its successive derivatives:

- Activity - total power of the signal (`var(s)`)
- Mobility - estimate of the mean frequency (`sqrt(var(s') / var(s))`)
- Complexity - similarity of the signal shape to a pure sine wave (`Mobility(s') / Mobility(s)`)
where s' = derivative(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

Named tuple containing:

- `h_act::Matrix{Float64}`: activity, shape `(channels, epochs)`
- `h_mob::Matrix{Float64}`: mobility, shape `(channels, epochs)`
- `h_comp::Matrix{Float64}`: complexity, shape `(channels, epochs)`
"""
function hjorth(
        obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}
    )::@NamedTuple{h_act::Matrix{Float64}, h_mob::Matrix{Float64}, h_comp::Matrix{Float64}}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    h = hjorth(@view(obj.data[ch, :, :]))

    return h

end
