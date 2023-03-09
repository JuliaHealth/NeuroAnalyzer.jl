export normalize
export normalize_zscore
export normalize_minmax
export normalize_n
export normalize_log
export normalize_gauss
export normalize_log10
export normalize_neglog
export normalize_neglog10
export normalize_neg
export normalize_pos
export normalize_perc
export normalize_invroot

"""
    normalize(signal, n; method)

Normalize.

# Arguments

- `signal::AbstractArray`
- `n::Real`
- `method::Symbol`: :zscore, :minmax, :log, :log10, :neglog, :neglog10, :neg, :pos, :perc, :gauss, :invroot, :n, :none

# Returns

- `normalized::Vector{Float64}`
"""
function normalize(signal::AbstractArray, n::Real=1.0; method::Symbol)

    _check_var(method, [:zscore, :minmax, :log, :log10, :neglog, :neglog10, :neg, :pos, :perc, :gauss, :invroot, :n, :none], "method")

    if method === :zscore
        return normalize_zscore(signal)
    elseif method === :minmax
        return normalize_minmax(signal)
    elseif method === :log
        return normalize_log(signal)
    elseif method === :log10
        return normalize_log10(signal)
    elseif method === :neglog
        return normalize_neglog(signal)
    elseif method === :neglog10
        return normalize_neglog10(signal)
    elseif method === :neg
        return normalize_neg(signal)
    elseif method === :pos
        return normalize_pos(signal)
    elseif method === :perc
        return normalize_perc(signal)
    elseif method === :gauss
        return normalize_gauss(signal)
    elseif method === :invroot
        return normalize_invroot(signal)
    elseif method === :n
        return normalize_n(signal, n)
    elseif method === :none
        return signal
    end
end

"""
    normalize_zscore(signal)

Normalize (by z-score).

# Arguments

- `signal::AbstractArray`

# Returns

- `normalized::Vector{Float64}`
"""
function normalize_zscore(signal::AbstractArray)
    m = mean(signal)
    s = std(signal)
    return @. (signal - m) / s
end

"""
    normalize_minmax(signal)

Normalize in [-1, +1].

# Arguments

- `signal::AbstractArray`

# Returns

- `normalized::AbstractArray`
"""
function normalize_minmax(signal::AbstractArray)
    mi = minimum(signal)
    mx = maximum(signal)
    mxi = mx - mi
    return @. (2 * (signal - mi) / mxi) - 1
end

"""
    normalize_n(signal, n)

Normalize in [0, n], default is [0, +1].

# Arguments

- `signal::AbstractArray`
- `n::Real=1.0`

# Returns

- `normalized::AbstractArray`
"""
function normalize_n(signal::AbstractArray, n::Real=1.0)
    max_x = maximum(signal)
    min_x = minimum(signal)
    return n .* (signal .- min_x) ./ (max_x - min_x)
end

"""
    normalize_log(signal)

Normalize using log-transformation.

# Arguments

- `signal::AbstractArray`

# Returns

- `normalized::AbstractArray`
"""
function normalize_log(signal::AbstractArray)
    m = abs(minimum(signal))
    return @. log(1 + signal + m)
end

"""
    normalize_gauss(signal, dims)

Normalize to Gaussian.

# Arguments

- `signal::AbstractArray`
- `dims::Int64=1`: dimension for cumsum()

# Returns

- `normalized::Vector{Float64}`
"""
function normalize_gauss(signal::AbstractArray, dims::Int64=1)
    dims in 1:ndims(signal) || throw(ArgumentError("dims must be in: 1:$(ndims(signal)).")) 
    l = length(signal) + 1
    return atanh.((tiedrank(cumsum(signal, dims=dims)) ./ l .- 0.5) .* 2)
end

"""
    normalize_log10(signal)

Normalize using log10-transformation.

# Arguments

- `signal::AbstractArray`

# Returns

- `normalized::Vector{Float64}`
"""
function normalize_log10(signal::AbstractArray)
    m = 1 + abs(minimum(signal))
    return @. log10(signal + m)
end

"""
    normalize_neglog(signal)

Normalize to using -log-transformation.

# Arguments

- `signal::AbstractArray`

# Returns

- `normalized::Vector{Float64}`
"""
function normalize_neglog(signal::AbstractArray)
    return @. -log(signal)
end

"""
    normalize_neglog10(signal)

Normalize using -log10-transformation.

# Arguments

- `signal::AbstractArray`

# Returns

- `normalized::Vector{Float64}`
"""
function normalize_neglog10(signal::AbstractArray)
    return @. -log10(signal)
end

"""
    normalize_neg(signal)

Normalize in [0, -∞].

# Arguments

- `signal::AbstractArray`

# Returns

- `normalized::Vector{Float64}`
"""
function normalize_neg(signal::AbstractArray)
    m = maximum(signal)
    return @. signal - m
end

"""
    normalize_pos(signal)

Normalize in [0, +∞].

# Arguments

- `signal::AbstractArray`

# Returns

- `normalized::Vector{Float64}`
"""
function normalize_pos(signal::AbstractArray)
    m = abs(minimum(signal))
    return @. signal + m
end

"""
    normalize_perc(signal)

Normalize in percentages.

# Arguments

- `signal::AbstractArray`

# Returns

- `normalized::Vector{Float64}`
"""
function normalize_perc(signal::AbstractArray)
    m1 = minimum(signal)
    m2 = maximum(signal)
    m = m2 - m1
    return (signal .- m1) ./ m
end

"""
    normalize_invroot(signal)

Normalize in inverse root (1/sqrt(x)).

# Arguments

- `signal::AbstractArray`

# Returns

- `normalized::Vector{Float64}`
"""
function normalize_invroot(signal::AbstractArray)
    # make signal > 0
    return 1 ./ (sqrt.(signal .+ abs(minimum(signal)) .+ eps()))
end

"""
    normalize(obj; channel, method)

Normalize channel(s)

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `method::Symbol`: method for normalization, see `normalize()` for details

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function normalize(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), method::Symbol)

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    obj_new = deepcopy(obj)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            @views obj_new.data[channel[ch_idx], :, ep_idx] = normalize(obj_new.data[channel[ch_idx], :, ep_idx], method=method)
        end
    end

    reset_components!(obj_new)
    push!(obj_new.header.history, "normalize(OBJ, channel=$channel, method=$method)")

    return obj_new
end

"""
    normalize!(obj; channel, method)

Normalize channel(s)

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `method::Symbol`: method for normalization, see `normalize()` for details
"""
function normalize!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), method::Symbol)

    obj_tmp = normalize(obj, channel=channel, method=method)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    reset_components!(obj)

    return nothing
end