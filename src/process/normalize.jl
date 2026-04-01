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
export normalize_sigmoid
export normalize_mad
export normalize_rank
export normalize_fisher

# ------------------------------------------------------------------ #
# normalize() — meta-function dispatcher                             #
# ------------------------------------------------------------------ #

"""
    normalize(s, n; <keyword arguments>)

Normalize.

# Arguments

- `s::AbstractVector`: signal vector
- `n::Real=1`: scaling parameter used by `:minmax`, `:n`
- `method::Symbol`: normalization method:
    - `:zscore`: by z-score
    - `:minmax`: in [-n, +n]
    - `:log`: using log-transformation
    - `:log10`: using log10-transformation
    - `:neglog`: using -log-transformation
    - `:neglog10`: using -log10-transformation
    - `:neg`: in [-∞, 0]
    - `:pos`: in [0, +∞]
    - `:perc`: in percentages
    - `:gauss`: to Gaussian
    - `:invroot`: to inverse root: 1/sqrt(x)
    - `:n`: in [0, n], default is [0, 1]
    - `:softmax`: using softmax function: exp(x_i) / sum(exp(x))
    - `:sigmoid`: using sigmoid function: 1 /  1 + exp(-x_i)
    - `:mad`: by MAD
    - `:rank`: using tiedranks
    - `:none`

# Returns

- `normalized::AbstractVector`
"""
function normalize(s::AbstractVector, n::Real = 1; method::Symbol)::AbstractVector
    _check_var(
        method,
        [
            :zscore,
            :minmax,
            :log,
            :log10,
            :neglog,
            :neglog10,
            :neg,
            :pos,
            :perc,
            :gauss,
            :invroot,
            :n,
            :softmax,
            :sigmoid,
            :mad,
            :rank,
            :none,
        ],
        "method",
    )

    method === :zscore && return normalize_zscore(s)
    method === :minmax && return normalize_minmax(s, n)
    method === :log && return normalize_log(s)
    method === :log10 && return normalize_log10(s)
    method === :neglog && return normalize_neglog(s)
    method === :neglog10 && return normalize_neglog10(s)
    method === :neg && return normalize_neg(s)
    method === :pos && return normalize_pos(s)
    method === :perc && return normalize_perc(s)
    method === :gauss && return normalize_gauss(s)
    method === :invroot && return normalize_invroot(s)
    method === :n && return normalize_n(s, n)
    method === :softmax && return normalize_softmax(s)
    method === :sigmoid && return normalize_sigmoid(s)
    method === :mad && return normalize_mad(s)
    method === :rank && return normalize_rank(s)
    method === :fisher && return normalize_fisher(s)
    return method === :none && return s
end

"""
    normalize(s, n; <keyword arguments>)

Normalize a signal array using the specified method.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `n::Real=1`: scaling parameter used by `:minmax`, `:n`
- `bych::Bool=false`: if `true`, normalize each channel separately
- `method::Symbol`: normalization method:
    - `:zscore`: by z-score
    - `:minmax`: in [-n, +n]
    - `:log`: using log-transformation
    - `:log10`: using log10-transformation
    - `:neglog`: using -log-transformation
    - `:neglog10`: using -log10-transformation
    - `:neg`: in [-∞, 0]
    - `:pos`: in [0, +∞]
    - `:perc`: in percentages
    - `:gauss`: to Gaussian
    - `:invroot`: to inverse root: 1/sqrt(x)
    - `:n`: in [0, n], default is [0, 1]; <keyword arguments>) .+ n1`
    - `:softmax`: using softmax function: exp(x_i) / sum(exp(x))
    - `:sigmoid`: using sigmoid function: 1 /  1 + exp(-x_i)
    - `:none`

# Returns

- `AbstractArray`: normalized signal, same shape as `s`
"""
function normalize(
    s::AbstractArray,
    n::Real = 1;
    bych::Bool = false,
    method::Symbol,
)::AbstractArray
    _check_var(
        method,
        [
            :zscore,
            :minmax,
            :log,
            :log10,
            :neglog,
            :neglog10,
            :neg,
            :pos,
            :perc,
            :gauss,
            :invroot,
            :n,
            :softmax,
            :sigmoid,
            :none,
        ],
        "method",
    )

    method === :zscore && return normalize_zscore(s; bych = bych)
    method === :minmax && return normalize_minmax(s, n; bych = bych)
    method === :log && return normalize_log(s; bych = bych)
    method === :log10 && return normalize_log10(s; bych = bych)
    method === :neglog && return normalize_neglog(s; bych = bych)
    method === :neglog10 && return normalize_neglog10(s; bych = bych)
    method === :neg && return normalize_neg(s; bych = bych)
    method === :pos && return normalize_pos(s; bych = bych)
    method === :perc && return normalize_perc(s; bych = bych)
    method === :gauss && return normalize_gauss(s; bych = bych)
    method === :invroot && return normalize_invroot(s; bych = bych)
    method === :n && return normalize_n(s, n; bych = bych)
    method === :softmax && return normalize_softmax(s; bych = bych)
    method === :sigmoid && return normalize_sigmoid(s; bych = bych)
    return method === :none && return s
end

"""
    normalize(obj; <keyword arguments>)

Normalize selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `method::Symbol`: normalization method:
    - `:zscore`: by z-score
    - `:minmax`: in [-n, +n]
    - `:log`: using log-transformation
    - `:log10`: using log10-transformation
    - `:neglog`: using -log-transformation
    - `:neglog10`: using -log10-transformation
    - `:neg`: in [-∞, 0]
    - `:pos`: in [0, +∞]
    - `:perc`: in percentages
    - `:gauss`: to Gaussian
    - `:invroot`: to inverse root: 1/sqrt(x)
    - `:n`: in [0, n], default is [0, 1]
    - `:softmax`: using softmax function: exp(x_i) / sum(exp(x))
    - `:sigmoid`: using sigmoid function: 1 /  1 + exp(-x_i)
    - `:mad`: by MAD
    - `:rank`: using tiedranks
    - `:none`
- `bych::Bool=false`: if `true`, normalize each channel separately
- `n::Real=1`: scaling parameter used by `:minmax`, `:n`

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function normalize(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    method::Symbol,
    bych::Bool = false,
    n::Real = 1,
)::NeuroAnalyzer.NEURO

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    ch_n = length(ch)
    ep_n = nepochs(obj)

    # create new dataset
    obj_new = deepcopy(obj)

    if bych
        # normalize each (channel, epoch) slice independently
        @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
            ch_idx, ep_idx = idx[1], idx[2]
            obj_new.data[ch[ch_idx], :, ep_idx] = NeuroAnalyzer.normalize(
                @view(obj_new.data[ch[ch_idx], :, ep_idx]),
                n,
                method = method,
            )
        end
    else
        # normalize the entire selected channel block at once
        obj_new.data[ch, :, :] = NeuroAnalyzer.normalize(
            obj_new.data[ch, :, :],
            n;
            method = method,
            bych = false,
        )
    end

    push!(obj_new.history, "normalize(obj; ch=$ch, method=$method, n=$n)")

    return obj_new
end

"""
    normalize!(obj; <keyword arguments>)

Normalize selected channels in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `method::Symbol`: normalization method:
    - `:zscore`: by z-score
    - `:minmax`: in [-n, +n]
    - `:log`: using log-transformation
    - `:log10`: using log10-transformation
    - `:neglog`: using -log-transformation
    - `:neglog10`: using -log10-transformation
    - `:neg`: in [-∞, 0]
    - `:pos`: in [0, +∞]
    - `:perc`: in percentages
    - `:gauss`: to Gaussian
    - `:invroot`: to inverse root: 1/sqrt(x)
    - `:n`: in [0, n], default is [0, 1]
    - `:softmax`: using softmax function: exp(x_i) / sum(exp(x))
    - `:sigmoid`: using sigmoid function: 1 /  1 + exp(-x_i)
    - `:mad`: by MAD
    - `:rank`: using tiedranks
    - `:none`
- `bych::Bool=false`: if `true`, normalize each channel separately
- `n::Real=1`: scaling parameter used by `:minmax`, `:n`

# Returns

- `Nothing`
"""
function normalize!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    method::Symbol,
    bych::Bool = false,
    n::Real = 1,
)::Nothing
    obj_new = NeuroAnalyzer.normalize(obj; ch = ch, method = method, bych = bych, n = n)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing
end

# ------------------------------------------------------------------ #
# individual normalization methods                                   #
# ------------------------------------------------------------------ #

"""
    normalize_zscore(s)

Normalize by z-score: `(x − x̄) / σ`. If `σ = 0`, subtracts the mean only.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `AbstractVector`: normalized signal, same shape as `s`
"""
function normalize_zscore(s::AbstractVector)::AbstractVector
    m = mean(s)
    sd = std(s)
    if sd != 0
        return @. (s - m) / sd
    else
        _warn("STD is 0; values normalized to (x − x̄).")
        return @. (s - m)
    end
end

"""
    normalize_zscore(s; <keyword arguments>)

Normalize by z-score across the whole array, or per channel when `bych=true`.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `AbstractArray`: normalized signal, same shape as `s`
"""
function normalize_zscore(s::AbstractArray; bych::Bool = false)::AbstractArray

    # validate
    ndims(s) <= 3 ||
        throw(ArgumentError("normalize_zscore() only works for arrays of ≤ 3 dimensions."))

    if !bych
        m = mean(s)
        sd = std(s)
        if sd != 0
            return @. (s - m) / sd
        else
            _warn("STD is 0; values normalized to (x − x̄).")
            return @. (s - m)
        end
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = normalize_zscore(@view(s[idx, :]))
            end
        else
            for ep in axes(s, 3), ch in axes(s, 1)
                sn[ch, :, ep] = normalize_zscore(@view(s[ch, :, ep]))
            end
        end
        return sn
    end
end

"""
    normalize_minmax(s, n)

Normalize in [-n, +n]. Constant signals are mapped to `+n`.

# Arguments

- `s::AbstractVector`: signal vector
- `n::Real=1`: scaling parameter

# Returns

- `AbstractVector`: normalized signal, same shape as `s`
"""
function normalize_minmax(s::AbstractVector, n::Real = 1)::AbstractVector

    # replace negative zero to avoid unexpected behavior with extrema
    s = replace(s, -0.0 => 0.0)

    if length(unique(s)) == 1
        return ones(length(s)) .* n
    end
    mi, mx = extrema(s)

    return @. ((2 * (s - mi) / (mx - mi)) - 1) * n
end

"""
    normalize_minmax(s, n; <keyword arguments>)

Normalize to `[−n, +n]` across the array, or per channel when `bych=true`. Constant signals are mapped to `+n`.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `n::Real=1`: scaling parameter
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `AbstractArray`: normalized signal, same shape as `s`
"""
function normalize_minmax(s::AbstractArray, n::Real = 1; bych::Bool = false)::AbstractArray

    # replace negative zero to avoid unexpected behavior with extrema
    s = replace(s, -0.0 => 0.0)
    length(unique(s)) == 1 && return ones(size(s)) .* n

    ndims(s) <= 3 ||
        throw(ArgumentError("normalize_minmax() only works for arrays of ≤ 3 dimensions."))

    length(unique(s)) == 1 && return ones(size(s)) .* n
    if !bych
        mi, mx = extrema(s)
        return @. ((2 * (s - mi) / (mx - mi)) - 1) * n
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = normalize_minmax(@view(s[idx, :]), n)
            end
        else
            for ep in axes(s, 3), ch in axes(s, 1)
                sn[ch, :, ep] = normalize_minmax(@view(s[ch, :, ep]), n)
            end
        end
        return sn
    end
end

"""
    normalize_n(s, n)

Normalize to `[0, n]` (default `[0, 1]`). Constant signals map to `n`.

# Arguments

- `s::AbstractVector`: signal vector
- `n::Real=1`: scaling parameter

# Returns

- `AbstractVector`: normalized signal, same shape as `s`
"""
function normalize_n(s::AbstractVector, n::Real = 1)::AbstractVector

    # replace negative zero to avoid unexpected behavior with extrema
    s = replace(s, -0.0 => 0.0)

    if length(unique(s)) == 1
        return ones(length(s)) .* n
    end
    smin, smax = extrema(s)

    return @. n * (s - smin) / (smax - smin)
end

"""
    normalize_n(s, n; <keyword arguments>)

Normalize to `[0, n]`, or per channel when `bych=true`.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `n::Real=1`
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `AbstractArray`: normalized signal, same shape as `s`
"""
function normalize_n(s::AbstractArray, n::Real = 1; bych::Bool = false)::AbstractArray
    ndims(s) <= 3 ||
        throw(ArgumentError("normalize_n() only works for arrays of ≤ 3 dimensions."))

    # replace negative zero to avoid unexpected behavior with extrema
    s = replace(s, -0.0 => 0.0)

    if length(unique(s)) == 1
        return ones(size(s)) .* n
    end

    if !bych
        smin, smax = extrema(s)
        return @. n * (s - smin) / (smax - smin)
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = normalize_n(@view(s[idx, :]), n)
            end
        else
            for ep in axes(s, 3), ch in axes(s, 1)
                sn[ch, :, ep] = normalize_n(@view(s[ch, :, ep]), n)
            end
        end
        return sn
    end
end

"""
    normalize_log(s)

Log-normalize: `log(1 + x + |min(x)|)`, ensuring the argument is always ≥ 1.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `AbstractVector`: normalized signal, same shape as `s`
"""
function normalize_log(s::AbstractVector)::AbstractVector
    m = abs(minimum(s))
    sn = @. log(1 + s + m)

    return sn
end

"""
    normalize_log(s; <keyword arguments>)

Normalize using log-transformation.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `AbstractArray`: normalized signal, same shape as `s`
"""
function normalize_log(s::AbstractArray; bych::Bool = false)::AbstractArray

    # validate
    ndims(s) <= 3 ||
        throw(ArgumentError("normalize_log() only works for arrays of ≤ 3 dimensions."))

    if !bych
        m = abs(minimum(s))
        return @. log(1 + s + m)
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = normalize_log(@view(s[idx, :]))
            end
        else
            for ep in axes(s, 3), ch in axes(s, 1)
                sn[ch, :, ep] = normalize_log(@view(s[ch, :, ep]))
            end
        end
        return sn
    end
end

"""
    normalize_gauss(s)

Normalize to Gaussian via rank-based inverse normal transform (Fisher–Yates): `atanh(2 * tiedrank(x) / (n+1) − 1)`.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `AbstractVector`: normalized signal, same shape as `s`
"""
function normalize_gauss(s::AbstractVector)::AbstractVector
    l = length(s) + 1
    sn = (tiedrank(s) ./ l .- 0.5) .* 2
    return atanh.(sn)
end

"""
    normalize_gauss(s; <keyword arguments>)

Normalize to Gaussian via rank-based inverse normal transform (Fisher–Yates): `atanh(2 * tiedrank(x) / (n+1) − 1)`.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `AbstractArray`: normalized signal, same shape as `s`
"""
function normalize_gauss(s::AbstractArray; bych::Bool = false)::AbstractArray

    # validate
    ndims(s) <= 3 ||
        throw(ArgumentError("normalize_gauss() only works for arrays of ≤ 3 dimensions."))

    if !bych
        l = length(s) + 1
        sn = (tiedrank(s) ./ l .- 0.5) .* 2
        return atanh.(sn)
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = normalize_gauss(@view(s[idx, :]))
            end
        else
            for ep in axes(s, 3), ch in axes(s, 1)
                sn[ch, :, ep] = normalize_gauss(@view(s[ch, :, ep]))
            end
        end
        return sn
    end
end

"""
    normalize_log10(s)

Log₁₀-normalize: `log10(x + 1 + |min(x)|)`.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `AbstractVector`: normalized signal, same shape as `s`
"""
function normalize_log10(s::AbstractVector)::AbstractVector
    m = 1 + abs(minimum(s))
    return @. log10(s + m)
end

"""
    normalize_log10(s; <keyword arguments>)

Log₁₀-normalize: `log10(x + 1 + |min(x)|)`.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `AbstractArray`: normalized signal, same shape as `s`
"""
function normalize_log10(s::AbstractArray; bych::Bool = false)::AbstractArray

    # validate
    ndims(s) <= 3 ||
        throw(ArgumentError("normalize_log10() only works for arrays of ≤ 3 dimensions."))

    if !bych
        m = 1 + abs(minimum(s))
        return @. log10(s + m)
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = normalize_log10(@view(s[idx, :]))
            end
        else
            for ep in axes(s, 3), ch in axes(s, 1)
                sn[ch, :, ep] = normalize_log10(@view(s[ch, :, ep]))
            end
        end
        return sn
    end
end

"""
    normalize_neglog(s; <keyword arguments>)

Negative log-normalize: `−log(x)`. `bych` is ignored (element-wise operation).

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `bych::Bool=false`: ignored

# Returns

- `Vector{Float64}`
"""
function normalize_neglog(s::AbstractArray; bych::Bool = false)::AbstractArray
    return @. -log(s)
end

"""
    normalize_neglog10(s; <keyword arguments>)

Negative log-normalize: `−log(x)`. `bych` is ignored (element-wise operation).

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `bych::Bool=false`: ignored

# Returns

- `AbstractArray`: normalized signal, same shape as `s`
"""
function normalize_neglog10(s::AbstractArray; bych::Bool = false)::AbstractArray
    return @. -log10(s)
end

"""
    normalize_neg(s)

Shift signal to `(−∞, 0]` by subtracting the maximum.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `AbstractVector`: normalized signal, same shape as `s`
"""
function normalize_neg(s::AbstractVector)::AbstractVector
    return s .- maximum(s)
end

"""
    normalize_neg(s)

Shift signal to `(−∞, 0]` by subtracting the maximum.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `AbstractArray`: normalized signal, same shape as `s`
"""
function normalize_neg(s::AbstractArray; bych::Bool = false)::AbstractArray

    # validate
    ndims(s) <= 3 ||
        throw(ArgumentError("normalize_neg() only works for arrays of ≤ 3 dimensions."))

    if !bych
        return s .- maximum(s)
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = normalize_neg(@view(s[idx, :]))
            end
        else
            for ep in axes(s, 3), ch in axes(s, 1)
                sn[ch, :, ep] = normalize_neg(@view(s[ch, :, ep]))
            end
        end
        return sn
    end
end

"""
    normalize_pos(s)

Shift signal to `[0, +∞)` by adding `|min(x)|`.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `AbstractVector`: normalized signal, same shape as `s`
"""
function normalize_pos(s::AbstractVector)::AbstractVector
    return s .+ abs(minimum(s))
end

"""
    normalize_pos(s; <keyword arguments>)

Shift signal to `[0, +∞)` by adding `|min(x)|`.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `AbstractArray`: normalized signal, same shape as `s`
"""
function normalize_pos(s::AbstractArray; bych::Bool = false)::AbstractArray

    # validate
    ndims(s) <= 3 ||
        throw(ArgumentError("normalize_pos() only works for arrays of ≤ 3 dimensions."))

    if !bych
        return s .+ abs(minimum(s))
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = normalize_pos(@view(s[idx, :]))
            end
        else
            for ep in axes(s, 3), ch in axes(s, 1)
                sn[ch, :, ep] = normalize_pos(@view(s[ch, :, ep]))
            end
        end
        return sn
    end
end

"""
    normalize_perc(s)

Normalize to percentages: `(x − min) / (max − min)`. Constant signals produce all-zeros.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `AbstractVector`: normalized signal, same shape as `s`
"""
function normalize_perc(s::AbstractVector)::AbstractVector
    m1 = minimum(s)
    m = maximum(s) - m1
    return m != 0 ? (s .- m1) ./ m : (s .- m1)
end

"""
    normalize_perc(s; <keyword arguments>)

Normalize to percentages: `(x − min) / (max − min)`. Constant signals produce all-zeros.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `AbstractArray`: normalized signal, same shape as `s`
"""
function normalize_perc(s::AbstractArray; bych::Bool = false)::AbstractArray

    # validate
    ndims(s) <= 3 ||
        throw(ArgumentError("normalize_perc() only works for arrays of ≤ 3 dimensions."))

    if !bych
        m1 = minimum(s)
        m = maximum(s) - m1
        return m != 0 ? (s .- m1) ./ m : (s .- m1)
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = normalize_perc(@view(s[idx, :]))
            end
        else
            for ep in axes(s, 3), ch in axes(s, 1)
                sn[ch, :, ep] = normalize_perc(@view(s[ch, :, ep]))
            end
        end
        return sn
    end
end

"""
    normalize_invroot(s)

Inverse-root normalize: `1 / √x`. Exact zeros are replaced with `eps()` to avoid division by zero. A copy of `s` is made to avoid mutating the input.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `AbstractVector`: normalized signal, same shape as `s`
"""
function normalize_invroot(s::AbstractVector)::AbstractVector
    sc = copy(s)
    idx = findall(iszero, sc)
    isempty(idx) || (sc[idx] .= eps())

    return 1 ./ sqrt.(sc)
end

"""
    normalize_invroot(s; <keyword arguments>)

Inverse-root normalize: `1 / √x`. Exact zeros are replaced with `eps()` to avoid division by zero. A copy of `s` is made to avoid mutating the input.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `AbstractArray`: normalized signal, same shape as `s`
"""
function normalize_invroot(s::AbstractArray; bych::Bool = false)::AbstractArray

    # validate
    ndims(s) <= 3 ||
        throw(ArgumentError("normalize_invroot() only works for arrays of ≤ 3 dimensions."))

    if !bych
        sc = copy(s)
        idx = findall(iszero, sc)
        isempty(idx) || (sc[idx] .= eps())
        return 1 ./ sqrt.(sc)
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = normalize_invroot(@view(s[idx, :]))
            end
        else
            for ep in axes(s, 3), ch in axes(s, 1)
                sn[ch, :, ep] = normalize_invroot(@view(s[ch, :, ep]))
            end
        end
        return sn
    end
end

"""
    normalize_softmax(s; <keyword arguments>)

Softmax normalization: `exp(xᵢ) / Σexp(x)`. `bych` is ignored.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `bych::Bool=false`: ignored

# Returns

- `AbstractArray`: normalized signal, same shape as `s`
"""
function normalize_softmax(s::AbstractArray; bych::Bool = false)::AbstractArray
    ex = exp.(s)
    return ex ./ sum(ex)
end

"""
    normalize_sigmoid(s; <keyword arguments>)

Sigmoid normalization: `1 / (1 + e^{−xᵢ})`. `bych` is ignored.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `bych::Bool=false`: ignored

# Returns

- `AbstractArray`: normalized signal, same shape as `s`
"""
function normalize_sigmoid(s::AbstractArray; bych::Bool = false)::AbstractArray
    return @. 1 / (1 + exp(-s))
end

"""
    normalize_mad(s)

Normalize by MAD: `(x − x̃) / (1.4826 × MAD)`. If MAD = 0, subtracts the median only.

The `1.4826` factor is the standard consistency constant for a normal distribution (scales MAD to estimate σ).

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `AbstractVector`: normalized signal, same shape as `s`
"""
function normalize_mad(s::AbstractVector)::AbstractVector
    m = median(s)
    md = 1.4826 * mad(s)
    if md != 0
        return @. (s - m) / md
    else
        _warn("MAD is 0; values normalized to (x − x̃).")
        return @. (s - m)
    end
end

"""
    normalize_mad(s; <keyword arguments>)

Normalize by MAD: `(x − x̃) / (1.4826 × MAD)`. If MAD = 0, subtracts the median only.

The `1.4826` factor is the standard consistency constant for a normal distribution (scales MAD to estimate σ).

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `AbstractArray`: normalized signal, same shape as `s`
"""
function normalize_mad(s::AbstractArray; bych::Bool = false)::AbstractArray

    # validate
    ndims(s) <= 3 ||
        throw(ArgumentError("normalize_mad() only works for arrays of ≤ 3 dimensions."))

    if !bych
        m = median(s)
        # FIX: original used `mad(s)` without the 1.4826 consistency factor,
        # making the array method inconsistent with the vector method.
        md = 1.4826 * mad(s)
        if md != 0
            return @. (s - m) / md
        else
            _warn("MAD is 0; values normalized to (x − x̃).")
            return @. (s - m)
        end
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = normalize_mad(@view(s[idx, :]))
            end
        else
            for ep in axes(s, 3), ch in axes(s, 1)
                sn[ch, :, ep] = normalize_mad(@view(s[ch, :, ep]))
            end
        end
        return sn
    end
end

"""
    normalize_rank(s)

Normalize using tied ranks (result is in `[1, n]`).

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `AbstractVector`: normalized signal, same shape as `s`
"""
function normalize_rank(s::AbstractVector)::AbstractVector
    sn = tiedrank(s)

    return sn
end

"""
    normalize_rank(s; <keyword arguments>)

Normalize using tied ranks (result is in `[1, n]`).

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `AbstractArray`: normalized signal, same shape as `s`
"""
function normalize_rank(s::AbstractArray; bych::Bool = false)::AbstractArray

    # validate
    length(unique(s)) == 1 && return ones(length(s))

    ndims(s) <= 3 ||
        throw(ArgumentError("normalize_rank() only works for arrays of ≤ 3 dimensions."))

    if !bych
        return tiedrank(s)
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = normalize_rank(@view(s[idx, :]))
            end
        else
            for ep in axes(s, 3), ch in axes(s, 1)
                sn[ch, :, ep] = normalize_rank(@view(s[ch, :, ep]))
            end
        end
        return sn
    end
end

"""
    normalize_fisher(s)

Fisher z-transform normalization. First normalizes to `[−1, +1]` via `normalize_minmax`, then applies `atanh(x) = 0.5 × ln((1+x)/(1−x))`. Endpoint values ±1 are nudged inward by `eps()` to avoid infinite output.

This is mathematically equivalent to `atanh`, and converts a uniform distribution into a normal distribution.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `AbstractVector`: normalized signal, same shape as `s`
"""
function normalize_fisher(s::AbstractVector)::AbstractVector
    sn = normalize_minmax(s)
    # clamp boundary values to avoid atanh(±1) = ±Inf.
    sn[sn .== -1.0] .= -1.0 + eps()
    sn[sn .== 1.0] .= 1.0 - eps()
    # `log(ℯ, x)` = `atanh`.
    return atanh.(sn)
end

"""
    normalize_fisher(s; <keyword arguments>)

Fisher z-transform normalization. First normalizes to `[−1, +1]` via `normalize_minmax`, then applies `atanh(x) = 0.5 × ln((1+x)/(1−x))`. Endpoint values ±1 are nudged inward by `eps()` to avoid infinite output.

This is mathematically equivalent to `atanh`, and converts a uniform distribution into a normal distribution.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `bych::Bool=false`: ignored

# Returns

- `AbstractArray`: normalized signal, same shape as `s`
"""
function normalize_fisher(s::AbstractArray; bych::Bool = false)::AbstractArray
    sn = normalize_minmax(s)
    sn[sn .== -1.0] .= -1.0 + eps()
    sn[sn .== 1.0] .= 1.0 - eps()
    # `log(ℯ, x)` = `atanh`.
    return atanh.(sn)
end
