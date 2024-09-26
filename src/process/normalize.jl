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

"""
    normalize(s, n; <keyword arguments>)

Normalize.

# Arguments

- `s::AbstractVector`
- `n::Real=1.0`
- `method::Symbol`:
    - `:zscore`: by z-score
    - `:minmax`: in [-1, +1]
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
    - `:mad`: by MAD
    - `:rank`: using tiedranks
    - `:none`

# Returns

- `normalized::AbstractVector`
"""
function normalize(s::AbstractVector, n::Float64=1.0; method::Symbol)

    _check_var(method, [:zscore, :minmax, :log, :log10, :neglog, :neglog10, :neg, :pos, :perc, :gauss, :invroot, :n, :softmax, :sigmoid, :mad, :rank, :none], "method")

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
    elseif method === :n
        return normalize_n(s, n)
    elseif method === :softmax
        return normalize_softmax(s)
    elseif method === :sigmoid
        return normalize_sigmoid(s)
    elseif method === :mad
        return normalize_mad(s)
    elseif method === :rank
        return normalize_rank(s)
    elseif method === :none
        return s
    end

end

"""
    normalize(s, n; <keyword arguments>)

Normalize.

# Arguments

- `s::AbstractArray`
- `n::Real=1.0`
- `bych::Bool=false`: if true, normalize each channel separately
- `method::Symbol`:
    - `:zscore`: by z-score
    - `:minmax`: in [-1, +1]
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

- `sn::AbstractArray`
"""
function normalize(s::AbstractArray, n::Float64=1.0; bych::Bool=false, method::Symbol)

    _check_var(method, [:zscore, :minmax, :log, :log10, :neglog, :neglog10, :neg, :pos, :perc, :gauss, :invroot, :n, :softmax, :sigmoid, :none], "method")

    if method === :zscore
        return normalize_zscore(s, bych=bych)
    elseif method === :minmax
        return normalize_minmax(s, bych=bych)
    elseif method === :log
        return normalize_log(s, bych=bych)
    elseif method === :log10
        return normalize_log10(s, bych=bych)
    elseif method === :neglog
        return normalize_neglog(s, bych=bych)
    elseif method === :neglog10
        return normalize_neglog10(s, bych=bych)
    elseif method === :neg
        return normalize_neg(s, bych=bych)
    elseif method === :pos
        return normalize_pos(s, bych=bych)
    elseif method === :perc
        return normalize_perc(s, bych=bych)
    elseif method === :gauss
        return normalize_gauss(s, bych=bych)
    elseif method === :invroot
        return normalize_invroot(s, bych=bych)
    elseif method === :n
        return normalize_n(s, n, bych=bych)
    elseif method === :softmax
        return normalize_softmax(s, bych=bych)
    elseif method === :sigmoid
        return normalize_sigmoid(s, bych=bych)
    elseif method === :none
        return s
    end

end

"""
    normalize_zscore(s)

Normalize by z-score.

# Arguments

- `s::AbstractVector`

# Returns

- `sn::AbstractVector`
"""
function normalize_zscore(s::AbstractVector)

    m = mean(s)
    sd = std(s)
    if sd != 0
        sn = @. (s - m) / sd
    else
        _warn("STD is 0, all values will be normalized to (x - x̄).")
        sn = @. (s - m)
    end

    return sn

end

"""
    normalize_zscore(s; <keyword arguments>)

# Arguments

- `s::AbstractArray`
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `sn::AbstractArray`
"""
function normalize_zscore(s::AbstractArray; bych::Bool=false)

    @assert ndims(s) <= 3 "normalize_zscore() only works for arrays of ≤ 3 dimensions."

    if bych == false
        m = mean(s)
        sd = std(s)
        if sd != 0
            sn = @. (s - m) / sd
        else
            _warn("STD is 0, all values will be normalized to (x - x̄).")
            sn = @. (s - m)
        end
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = @views normalize_zscore(s[idx, :])
            end
        elseif ndims(s) == 3
            for idx1 in axes(s, 3)
                for idx2 in axes(s, 1)
                    sn[idx2, :, idx1] = @views normalize_zscore(s[idx2, :, idx1])
                end
            end
        end
    end

    return sn

end

"""
    normalize_minmax(s)

Normalize in [-1, +1]. If all elements are the same, they are normalized to +1.0.

# Arguments

- `s::AbstractVector`

# Returns

- `sn::AbstractVector`
"""
function normalize_minmax(s::AbstractVector)

    if length(unique(s)) == 1
        sn = ones(length(s))
    else
        mi, mx = extrema(s)
        mxi = mx - mi
        sn = @. (2 * (s - mi) / mxi) - 1
    end

    return sn

end

"""
    normalize_minmax(s; <keyword arguments>)

Normalize in [-1, +1]. If all elements are the same, they are normalized to +1.0.

# Arguments

- `s::AbstractArray`
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `sn::AbstractArray`
"""

function normalize_minmax(s::AbstractArray; bych::Bool=false)

    length(unique(s)) == 1 && return ones(length(s))

    @assert ndims(s) <= 3 "normalize_minmax() only works for arrays of ≤ 3 dimensions."

    if bych == false
        mi, mx = extrema(s)
        mxi = mx - mi
        sn = @. (2 * (s - mi) / mxi) - 1
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = @views normalize_minmax(s[idx, :])
            end
        elseif ndims(s) == 3
            for idx1 in axes(s, 3)
                for idx2 in axes(s, 1)
                    sn[idx2, :, idx1] = @views normalize_minmax(s[idx2, :, idx1])
                end
            end
        end
    end

    return sn

end

"""
    normalize_n(s, n)

Normalize in [0, n], default is [0, +1].

# Arguments

- `s::AbstractVector`
- `n::Real=1.0`

# Returns

- `sn::AbstractVector`
"""
function normalize_n(s::AbstractVector, n::Real=1.0)

    if length(unique(s)) == 1
        sn = ones(length(s)) .* n
    else
        smin, smax = extrema(s)
        sn = @. n * (s - smin) / (smax - smin)
    end

    return sn

end

"""
    normalize_n(s, n; <keyword arguments>)

Normalize in [0, n], default is [0, +1].

# Arguments

- `s::AbstractArray`
- `n::Real=1.0`
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `sn::AbstractArray`
"""
function normalize_n(s::AbstractArray, n::Real=1.0; bych::Bool=false)

    @assert ndims(s) <= 3 "normalize_n() only works for arrays of ≤ 3 dimensions."

    if length(unique(s)) == 1
        sn = zeros(length(s))
    else
        sn = zeros(size(s))
        if bych == false
            smin, smax = extrema(s)
            sn = @. n * (s - smin) / (smax - smin)
        else
            if ndims(s) == 2
                for idx in axes(s, 1)
                    sn[idx, :] = @views normalize_n(s[idx, :], n)
                end
            elseif ndims(s) == 3
                for idx1 in axes(s, 3)
                    for idx2 in axes(s, 1)
                        sn[idx2, :, idx1] = @views normalize_n(s[idx2, :, idx1], n)
                    end
                end
            end
        end
    end

    return sn

end

"""
    normalize_log(s)

Normalize using log-transformation.

# Arguments

- `s::AbstractVector`

# Returns

- `sn::AbstractVector`
"""
function normalize_log(s::AbstractVector)

    m = abs(minimum(s))
    sn = @. log(1 + s + m)

    return sn

end

"""
    normalize_log(s; <keyword arguments>)

Normalize using log-transformation.

# Arguments

- `s::AbstractArray`
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `sn::AbstractArray`
"""
function normalize_log(s::AbstractArray; bych::Bool=false)

    @assert ndims(s) <= 3 "normalize_log() only works for arrays of ≤ 3 dimensions."

    if bych == false
        m = abs(minimum(s))
        sn = @. log(1 + s + m)
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = @views normalize_log(s[idx, :])
            end
        elseif ndims(s) == 3
            for idx1 in axes(s, 3)
                for idx2 in axes(s, 1)
                    sn[idx2, :, idx1] = @views normalize_log(s[idx2, :, idx1])
                end
            end
        end
    end

    return sn

end

"""
    normalize_gauss(s)

Normalize to Gaussian.

# Arguments

- `s::AbstractVector`

# Returns

- `sn::AbstractVector`
"""
function normalize_gauss(s::AbstractVector)

    l = length(s) + 1
    sn = (tiedrank(s) ./ l .- 0.5) .* 2
    sn = atanh.(sn)

    # sn = tiedrank(s)
    # sn = normalize_minmax(sn)
    # sn[sn .== -1] .= -1 + eps()
    # sn[sn .== 1] .= 1 - eps()
    # sn = atanh.(sn)

    return sn

end

"""
    normalize_gauss(s; <keyword arguments>)

Normalize to Gaussian.

# Arguments

- `s::AbstractArray`
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `sn::AbstractArray`
"""
function normalize_gauss(s::AbstractArray; bych::Bool=false)

    @assert ndims(s) <= 3 "normalize_gauss() only works for arrays of ≤ 3 dimensions."

    if bych == false
        l = length(s) + 1
        sn = (tiedrank(s) ./ l .- 0.5) .* 2
        sn = atanh.(sn)
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = @views normalize_gauss(s[idx, :])
            end
        elseif ndims(s) == 3
            for idx1 in axes(s, 3)
                for idx2 in axes(s, 1)
                    sn[idx2, :, idx1] = @views normalize_gauss(s[idx2, :, idx1])
                end
            end
        end
    end

    return sn

end

"""
    normalize_log10(s)

Normalize using log10-transformation.

# Arguments

- `s::AbstractVector`

# Returns

- `sn::AbstractVector`
"""
function normalize_log10(s::AbstractVector)

    m = 1 + abs(minimum(s))
    sn = @. log10(s + m)

    return sn

end

"""
    normalize_log10(s; <keyword arguments>)

Normalize using log10-transformation.

# Arguments

- `s::AbstractArray`
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `sn::AbstractArray`
"""
function normalize_log10(s::AbstractArray; bych::Bool=false)

    @assert ndims(s) <= 3 "normalize_log10() only works for arrays of ≤ 3 dimensions."

    if bych == false
        m = 1 + abs(minimum(s))
        sn = @. log10(s + m)
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = @views normalize_log10(s[idx, :])
            end
        elseif ndims(s) == 3
            for idx1 in axes(s, 3)
                for idx2 in axes(s, 1)
                    sn[idx2, :, idx1] = @views normalize_log10(s[idx2, :, idx1])
                end
            end
        end
    end

    return sn

end

"""
    normalize_neglog(s; <keyword arguments>)

Normalize to using -log-transformation.

# Arguments

- `s::AbstractArray`
- `bych::Bool=false`: ignored

# Returns

- `sn::Vector{Float64}`
"""
function normalize_neglog(s::AbstractArray; bych::Bool=false)

    sn = @. -log(s)

    return sn

end

"""
    normalize_neglog10(s; <keyword arguments>)

Normalize using -log10-transformation.

# Arguments

- `s::AbstractArray`
- `bych::Bool=false`: ignored

# Returns

- `sn::AbstractArray`
"""
function normalize_neglog10(s::AbstractArray; bych::Bool=false)

    sn = @. -log10(s)

    return sn

end

"""
    normalize_neg(s)

Normalize in [-∞, 0].

# Arguments

- `s::AbstractVector`

# Returns

- `sn::AbstractVector`
"""
function normalize_neg(s::AbstractVector)

    m = maximum(s)
    sn = @. s - m

    return sn

end

"""
    normalize_neg(s)

Normalize in [-∞, 0].

# Arguments

- `s::AbstractArray`
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `sn::AbstractArray`
"""
function normalize_neg(s::AbstractArray; bych::Bool=false)

    @assert ndims(s) <= 3 "normalize_neg() only works for arrays of ≤ 3 dimensions."

    if bych == false
        m = maximum(s)
        sn = @. s - m
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = @views normalize_neg(s[idx, :])
            end
        elseif ndims(s) == 3
            for idx1 in axes(s, 3)
                for idx2 in axes(s, 1)
                    sn[idx2, :, idx1] = @views normalize_neg(s[idx2, :, idx1])
                end
            end
        end
    end

    return sn

end

"""
    normalize_pos(s)

Normalize in [0, +∞].

# Arguments

- `s::AbstractVector`

# Returns

- `sn::AbstractVector`
"""
function normalize_pos(s::AbstractVector)

    m = abs(minimum(s))
    sn = @. s + m

    return sn

end

"""
    normalize_pos(s; <keyword arguments>)

Normalize in [0, +∞].

# Arguments

- `s::AbstractArray`
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `sn::AbstractArray`
"""
function normalize_pos(s::AbstractArray; bych::Bool=false)

    @assert ndims(s) <= 3 "normalize_pos() only works for arrays of ≤ 3 dimensions."

    if bych == false
        m = abs(minimum(s))
        sn = @. s + m
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = @views normalize_pos(s[idx, :])
            end
        elseif ndims(s) == 3
            for idx1 in axes(s, 3)
                for idx2 in axes(s, 1)
                    sn[idx2, :, idx1] = @views normalize_pos(s[idx2, :, idx1])
                end
            end
        end
    end

    return sn

end

"""
    normalize_perc(s)

Normalize in percentages.

# Arguments

- `s::AbstractVector`

# Returns

- `sn::AbstractVector`
"""
function normalize_perc(s::AbstractVector)

    m1 = minimum(s)
    m2 = maximum(s)
    m = m2 - m1

    if m != 0
        sn = (s .- m1) ./ m
    else
        sn = (s .- m1)
    end

    return sn

end

"""
    normalize_perc(s; <keyword arguments>)

Normalize in percentages.

# Arguments

- `s::AbstractArray`
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `sn::AbstractArray`
"""
function normalize_perc(s::AbstractArray; bych::Bool=false)

    @assert ndims(s) <= 3 "normalize_perc() only works for arrays of ≤ 3 dimensions."

    if bych == false
        m1 = minimum(s)
        m2 = maximum(s)
        m = m2 - m1
        if m != 0
            sn = (s .- m1) ./ m
        else
            sn = (s .- m1)
        end
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = @views normalize_perc(s[idx, :])
            end
        elseif ndims(s) == 3
            for idx1 in axes(s, 3)
                for idx2 in axes(s, 1)
                    sn[idx2, :, idx1] = @views normalize_perc(s[idx2, :, idx1])
                end
            end
        end
    end

    return sn

end

"""
    normalize_invroot(s)

Normalize in inverse root (1/sqrt(x)).

# Arguments

- `s::AbstractVector`

# Returns

- `sn::AbstractVector`
"""
function normalize_invroot(s::AbstractVector)

    # make s > 0
    idx = findall(x -> isequal(x, 0), s)
    length(idx) > 0 && (s[idx] .= eps())
    sn = 1 ./ (sqrt.(s))

    return sn

end

"""
    normalize_invroot(s; <keyword arguments>)

Normalize in inverse root (1/sqrt(x)).

# Arguments

- `s::AbstractArray`
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `sn::AbstractArray`
"""
function normalize_invroot(s::AbstractArray; bych::Bool=false)

    @assert ndims(s) <= 3 "normalize_invroot() only works for arrays of ≤ 3 dimensions."

    if bych == false
        idx = findall(x -> isequal(x, 0), s)
        length(idx) > 0 && (s[idx] .= eps())
        sn = 1 ./ (sqrt.(s))
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = @views normalize_invroot(s[idx, :])
            end
        elseif ndims(s) == 3
            for idx1 in axes(s, 3)
                for idx2 in axes(s, 1)
                    sn[idx2, :, idx1] = @views normalize_invroot(s[idx2, :, idx1])
                end
            end
        end
    end

    return sn

end

"""
    normalize_softmax(s; <keyword arguments>)

Softmax normalize: `exp(x_i) / sum(exp(x))`

# Arguments

- `s::AbstractArray`
- `bych::Bool=false`: ignored

# Returns

- `sn::AbstractArray`
"""
function normalize_softmax(s::AbstractArray; bych::Bool=false)

    return exp.(s) ./ sum(exp.(s))

end

"""
    normalize_sigmoid(s; <keyword arguments>)

Normalize using sigmoid function: `1 / (1 + e^-x_i)`

# Arguments

- `s::AbstractArray`
- `bych::Bool=false`: ignored

# Returns

- `sn::Vector{Float64}`
"""
function normalize_sigmoid(s::AbstractArray; bych::Bool=false)

    return @. 1 / (1 + exp(-s))

end

"""
    normalize_mad(s)

Normalize by MAD.

# Arguments

- `s::AbstractVector`

# Returns

- `sn::AbstractVector`
"""
function normalize_mad(s::AbstractVector)

    m = median(s)
    md = 1.4826 * mad(s)
    if md != 0
        sn = @. (s - m) / md
    else
        _warn("MAD is 0, all values will be normalized to (x - x̃).")
        sn = @. (s - m)
    end

    return sn

end

"""
    normalize_mad(s; <keyword arguments>)

# Arguments

- `s::AbstractArray`
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `sn::AbstractArray`
"""
function normalize_mad(s::AbstractArray; bych::Bool=false)

    @assert ndims(s) <= 3 "normalize_mad() only works for arrays of ≤ 3 dimensions."

    if bych == false
        m = median(s)
        md = mad(s)
        if md != 0
            sn = @. (s - m) / md
        else
            _warn("MAD is 0, all values will be normalized to (x - x̃).")
            sn = @. (s - m)
        end
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = @views normalize_mad(s[idx, :])
            end
        elseif ndims(s) == 3
            for idx1 in axes(s, 3)
                for idx2 in axes(s, 1)
                    sn[idx2, :, idx1] = @views normalize_mad(s[idx2, :, idx1])
                end
            end
        end
    end

    return sn

end

"""
    normalize_rank(s)

Normalize using tiedranks.

# Arguments

- `s::AbstractVector`

# Returns

- `sn::AbstractVector`
"""
function normalize_rank(s::AbstractVector)

    sn = tiedrank(s)

    return sn

end

"""
    normalize_rank(s; <keyword arguments>)

Normalize using tiedranks.

# Arguments

- `s::AbstractArray`
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `sn::AbstractArray`
"""

function normalize_rank(s::AbstractArray; bych::Bool=false)

    length(unique(s)) == 1 && return ones(length(s))

    @assert ndims(s) <= 3 "normalize_rank() only works for arrays of ≤ 3 dimensions."

    if bych == false
        sn = tiedrank(s)
    else
        sn = zeros(size(s))
        if ndims(s) == 2
            for idx in axes(s, 1)
                sn[idx, :] = @views normalize_rank(s[idx, :])
            end
        elseif ndims(s) == 3
            for idx1 in axes(s, 3)
                for idx2 in axes(s, 1)
                    sn[idx2, :, idx1] = @views normalize_rank(s[idx2, :, idx1])
                end
            end
        end
    end

    return sn

end

"""
    normalize_fisher(s)

Normalize using Fisher z-transform. Converts uniform distribution into normal distribution.

# Arguments

- `s::AbstractVector`

# Returns

- `sn::AbstractVector`
"""
function normalize_fisher(s::AbstractVector)

    sn = normalize_minmax(s)
    sn[sn .== -1] .= -1 + eps()
    sn[sn .== 1] .= 1 - eps()
    sn = @. 0.5 * log(ℯ, (1 + sn) / (1 - sn))

    return sn

end

"""
    normalize_fisher(s; <keyword arguments>)

Normalize using Fisher z-transform. Converts uniform distribution into normal distribution.

# Arguments

- `s::AbstractArray`
- `bych::Bool=false`: ignored

# Returns

- `sn::AbstractArray`
"""

function normalize_fisher(s::AbstractArray; bych::Bool=false)

    sn = normalize_minmax(s)
    sn[sn .== -1] .= -1 + eps()
    sn[sn .== 1] .= 1 - eps()
    sn = @. 0.5 * log(ℯ, (1 + sn) / (1 - sn))

    return sn

end

"""
    normalize(obj; <keyword arguments>)

Normalize channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channel name or list of channel names
- `method::Symbol`: method for normalization, see `normalize()` for details
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function normalize(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, method::Symbol, bych::Bool=false)

    ch = get_channel(obj, ch=ch)
    ch_n = length(ch)
    ep_n = nepochs(obj)

    obj_new = deepcopy(obj)
    if bych
        @inbounds for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
                @views obj_new.data[ch[ch_idx], :, ep_idx] = normalize(obj_new.data[ch[ch_idx], :, ep_idx], method=method)
            end
        end
    else
        obj_new.data[ch, :, :] = normalize(obj_new.data[ch, :, :], method=method, bych=false)
    end

    reset_components!(obj_new)
    push!(obj_new.history, "normalize(OBJ, ch=$ch, method=$method)")

    return obj_new

end

"""
    normalize!(obj; <keyword arguments>)

Normalize channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channel name or list of channel names
- `method::Symbol`: method for normalization, see `normalize()` for details
- `bych::Bool=false`: if true, normalize each channel separately
"""
function normalize!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, method::Symbol, bych::Bool=false)

    obj_new = normalize(obj, ch=ch, method=method, bych=bych)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end