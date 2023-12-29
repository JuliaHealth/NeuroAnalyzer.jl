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
export normalize_softmax

"""
    normalize(s, n; method)

Normalize.

# Arguments

- `s::AbstractArray`
- `m::Real=0.0`
- `n::Real=1.0`
- `method::Symbol`:
    - `:zscore`: by z-score
    - `:minmax`: in [-1, +1]
    - `:log`: using log-transformation
    - `:log10`: using log10-transformation
    - `:neglog`: using -log-transformation
    - `:neglog10`: using -log10-transformation
    - `:neg`: in [0, -∞]
    - `:pos`: in [0, +∞]
    - `:perc`: in percentages
    - `:gauss`: to Gaussian
    - `:invroot`: in inverse root (1/sqrt(x))
    - `:softmax`: exp(x_i) / sum(exp(x))
    - `:n`: in [0, n1], default is [0, 1]; to normalize to [n1, n2], use `normalize_n(s) .* (n2 - n1) .+ n1`
    - `:none`

# Returns

- `normalized::Vector{Float64}`
"""
function normalize(s::AbstractArray, n::Real=1.0; method::Symbol)

    _check_var(method, [:zscore, :minmax, :log, :log10, :neglog, :neglog10, :neg, :pos, :perc, :gauss, :invroot, :n, :mn, :none], "method")

    if method === :zscore
        return normalize_zscore(s)
    elseif method === :minmax
        return normalize_minmax(s)
    elseif method === :log
        return normalize_log(s)
    elseif method === :log10
        return normalize_log10(s)
    elseif method === :neglog
        return normalize_neglog(s)
    elseif method === :neglog10
        return normalize_neglog10(s)
    elseif method === :neg
        return normalize_neg(s)
    elseif method === :pos
        return normalize_pos(s)
    elseif method === :perc
        return normalize_perc(s)
    elseif method === :gauss
        return normalize_gauss(s)
    elseif method === :invroot
        return normalize_invroot(s)
    elseif method === :softmax
        return normalize_softmax(s)
    elseif method === :n
        return normalize_n(s, n)
    elseif method === :none
        return s
    end

end

"""
    normalize_zscore(s)

Normalize (by z-score).

# Arguments

- `s::AbstractArray`

# Returns

- `normalize_zscore::Vector{Float64}`
"""
function normalize_zscore(s::AbstractArray)

    m = mean(s)
    sd = std(s)

    if sd != 0
        return @. (s - m) / sd
    else
        return @. (s - m) / sd
    end

end

"""
    normalize_minmax(s)

Normalize in [-1, +1].

# Arguments

- `s::AbstractArray`

# Returns

- `normalize_minmax::AbstractArray`
"""
function normalize_minmax(s::AbstractArray)

    if length(unique(s)) == 1
        return ones(length(s))
    else
        mi = minimum(s)
        mx = maximum(s)
        mxi = mx - mi
        return @. (2 * (s - mi) / mxi) - 1
    end

end

"""
    normalize_n(s, n)

Normalize in [0, n], default is [0, +1].

# Arguments

- `s::AbstractArray`
- `n::Real=1.0`

# Returns

- `normalize_n::AbstractArray`
"""
function normalize_n(s::AbstractArray, n::Real=1.0)

    if length(unique(s)) == 1
        return ones(length(s)) .* n
    else
        smin, smax = extrema(s)
        return @. n * (s - smin) / (smax - smin)
    end

end

"""
    normalize_log(s)

Normalize using log-transformation.

# Arguments

- `s::AbstractArray`

# Returns

- `normalize_log::AbstractArray`
"""
function normalize_log(s::AbstractArray)

    m = abs(minimum(s))

    return @. log(1 + s + m)

end

"""
    normalize_gauss(s, dims)

Normalize to Gaussian.

# Arguments

- `s::AbstractArray`
- `dims::Int64=1`: dimension for cumsum()

# Returns

- `normalize_gauss::Vector{Float64}`
"""
function normalize_gauss(s::AbstractArray, dims::Int64=1)

    @assert dims in 1:ndims(s) "dims must be in: 1:$(ndims(s))."

    l = length(s) + 1

    return atanh.((tiedrank(cumsum(s, dims=dims)) ./ l .- 0.5) .* 2)

end

"""
    normalize_log10(s)

Normalize using log10-transformation.

# Arguments

- `s::AbstractArray`

# Returns

- `normalize_log10::Vector{Float64}`
"""
function normalize_log10(s::AbstractArray)

    m = 1 + abs(minimum(s))

    return @. log10(s + m)

end

"""
    normalize_neglog(s)

Normalize to using -log-transformation.

# Arguments

- `s::AbstractArray`

# Returns

- `normalize_neglog::Vector{Float64}`
"""
function normalize_neglog(s::AbstractArray)

    return @. -log(s)

end

"""
    normalize_neglog10(s)

Normalize using -log10-transformation.

# Arguments

- `s::AbstractArray`

# Returns

- `normalize_neglog::Vector{Float64}`
"""
function normalize_neglog10(s::AbstractArray)

    return @. -log10(s)

end

"""
    normalize_neg(s)

Normalize in [0, -∞].

# Arguments

- `s::AbstractArray`

# Returns

- `normalize_neg::Vector{Float64}`
"""
function normalize_neg(s::AbstractArray)

    m = maximum(s)

    return @. s - m

end

"""
    normalize_pos(s)

Normalize in [0, +∞].

# Arguments

- `s::AbstractArray`

# Returns

- `normalize_pos::Vector{Float64}`
"""
function normalize_pos(s::AbstractArray)

    m = abs(minimum(s))

    return @. s + m

end

"""
    normalize_perc(s)

Normalize in percentages.

# Arguments

- `s::AbstractArray`

# Returns

- `normalize_perc::Vector{Float64}`
"""
function normalize_perc(s::AbstractArray)

    m1 = minimum(s)
    m2 = maximum(s)
    m = m2 - m1

    if m != 0
        return (s .- m1) ./ m
    else
        return (s .- m1)
    end

end

"""
    normalize_invroot(s)

Normalize in inverse root (1/sqrt(x)).

# Arguments

- `s::AbstractArray`

# Returns

- `normalize_invroot::Vector{Float64}`
"""
function normalize_invroot(s::AbstractArray)

    # make s > 0
    return 1 ./ (sqrt.(s .+ abs(minimum(s)) .+ eps()))

end

"""
    normalize_softmax(s)

Softmax normalize: `exp(x_i) / sum(exp(x))`

# Arguments

- `s::AbstractArray`

# Returns

- `normalize_softmax::Vector{Float64}`
"""
function normalize_softmax(s::AbstractArray)

    return exp.(s) ./ sum(exp.(s))

end

"""
    normalize(obj; ch, method)

Normalize channel(s)

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
- `method::Symbol`: method for normalization, see `normalize()` for details

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function normalize(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)), method::Symbol)

    _check_channels(obj, ch)
    ch_n = length(ch)
    ep_n = nepochs(obj)

    obj_new = deepcopy(obj)
    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            @views obj_new.data[ch[ch_idx], :, ep_idx] = normalize(obj_new.data[ch[ch_idx], :, ep_idx], method=method)
        end
    end

    reset_components!(obj_new)
    push!(obj_new.history, "normalize(OBJ, ch=$ch, method=$method)")

    return obj_new
end

"""
    normalize!(obj; ch, method)

Normalize channel(s)

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
- `method::Symbol`: method for normalization, see `normalize()` for details
"""
function normalize!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)), method::Symbol)

    obj_new = normalize(obj, ch=ch, method=method)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end