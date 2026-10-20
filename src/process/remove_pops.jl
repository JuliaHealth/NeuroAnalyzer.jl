export remove_pops
export remove_pops!

"""
    _pop_detect(s; <keyword arguments>)

Detect electrode pops: abrupt steps of the signal (outliers of the first difference) after which the signal stays at the new level for a while.

# Arguments

- `s::AbstractVector`: signal vector
- `z::Real`: detection threshold (robust z-score of the first difference)
- `r::Int64`: minimum distance between pops in samples
- `wlen::Int64`: window length in samples for the robust statistics (the remainder is added to the last window)

# Returns

- `Vector{Int64}`: first sample of each pop (after the jump), sorted
"""
function _pop_detect(s::AbstractVector; z::Real, r::Int64, wlen::Int64)::Vector{Int64}
    n = length(s)
    d = diff(s)
    nd = length(d)

    # robust z-score of the first difference, statistics per window (amplitude may change during the recording)
    zs = zeros(nd)
    w_starts = collect(1:wlen:nd)
    length(w_starts) > 1 && nd - w_starts[end] + 1 < wlen && pop!(w_starts)
    for (i, w1) in enumerate(w_starts)
        w2 = i == length(w_starts) ? nd : w1 + wlen - 1
        dw = @view d[w1:w2]
        md = median(dw)
        sd = 1.4826 * mad(dw; center = md, normalize = false)
        sd > 0 && (zs[w1:w2] = abs.(dw .- md) ./ sd)
    end

    # level before / after the jump is compared over k samples
    k = max(3, r ÷ 10)
    # candidates, largest jumps first; d[i] is the jump between samples i and i + 1
    cand = sort(findall(>(z), zs); by = i -> -zs[i])
    pops = Int64[]
    for i in cand
        p = i + 1
        any(q -> abs(q - p) < r, pops) && continue
        # the levels before and after the jump must be measurable (e.g. not a filter transient at the signal start)
        (p - k < 1 || p + k - 1 > n) && continue
        # a pop is a step: the level after the jump differs from the level before by at least half of the jump
        # (a single-sample spike returns to the previous level)
        pre = @view s[(p - k):(p - 1)]
        post = @view s[p:(p + k - 1)]
        step = mean(post) - mean(pre)
        (sign(step) == sign(d[i]) && abs(step) > 0.5 * abs(d[i])) || continue
        push!(pops, p)
    end

    return sort!(pops)
end

"""
    _pop_fit(x)

Fit an exponential decay `a × exp(-t / τ)`, `t = 0, 1, …`, to a pop (least squares: `a` analytically for a given `τ`, `τ` by a grid search refined by golden-section search).

# Arguments

- `x::AbstractVector`: signal after the pop minus the level before the pop

# Returns

Named tuple:

- `a::Float64`: pop amplitude
- `tau::Float64`: time constant in samples
"""
function _pop_fit(x::AbstractVector)::@NamedTuple{a::Float64, tau::Float64}
    t = 0:(length(x) - 1)
    # least-squares amplitude and residual sum of squares for a given τ
    function _sse(τ)
        e = exp.(-t ./ τ)
        a = dot(x, e) / dot(e, e)
        return sum(abs2, x .- a .* e), a
    end

    # grid search over log τ (0.3 sample … 4 × fitted length)
    lτ = range(-0.5, log10(4 * length(x)); length = 80)
    sse = [_sse(exp10(v))[1] for v in lτ]
    i = argmin(sse)

    # golden-section search between the neighbouring grid points
    lo, hi = lτ[max(1, i - 1)], lτ[min(length(lτ), i + 1)]
    g = (sqrt(5) - 1) / 2
    for _ = 1:40
        m1 = hi - g * (hi - lo)
        m2 = lo + g * (hi - lo)
        _sse(exp10(m1))[1] < _sse(exp10(m2))[1] ? (hi = m2) : (lo = m1)
    end
    tau = exp10((lo + hi) / 2)
    a = _sse(tau)[2]

    return (; a, tau)
end

"""
    _remove_pops(s; <keyword arguments>)

Detect and repair electrode pops in a signal vector.

A pop is modeled as a step followed by an exponential return to the previous level: `a × exp(-t / τ)`. The model is fitted to `r` samples after the jump; candidates that are small (`|a|` < `a_min` × robust SD of the signal) or do not follow the model (R² < `r2_min`) are rejected (R² is calculated over the first 3τ samples). Accepted pops are subtracted from the signal (up to the next pop or the end of the signal); the underlying signal is preserved.

# Arguments

- `s::AbstractVector`: signal vector
- `r::Int64=20`: fitting segment length in samples, also the minimum distance between pops
- `z::Real=10`: detection threshold (robust z-score of the first difference)
- `wlen::Int64=length(s)`: window length in samples for the robust statistics
- `a_min::Real=5`: minimum pop amplitude relative to the robust SD of the signal (within `wlen`)
- `r2_min::Real=0.8`: minimum R² of the exponential model
- `repair::Bool=true`: subtract the fitted pops if `true`

# Returns

Named tuple:

- `s::Vector{Float64}`: repaired signal (a copy of the input if `repair=false`)
- `pop_loc::Vector{Int64}`: first sample of each pop
- `amp::Vector{Float64}`: pop amplitudes (signal units)
- `tau::Vector{Float64}`: time constants of the recovery in samples
"""
function _remove_pops(
    s::AbstractVector;
    r::Int64 = 20,
    z::Real = 10,
    wlen::Int64 = length(s),
    a_min::Real = 5,
    r2_min::Real = 0.8,
    repair::Bool = true,
)::@NamedTuple{s::Vector{Float64}, pop_loc::Vector{Int64}, amp::Vector{Float64}, tau::Vector{Float64}}
    # validate
    r >= 5 || throw(ArgumentError("r must be ≥ 5."))
    z > 0 || throw(ArgumentError("z must be > 0."))
    length(s) >= 2 * r + 1 || throw(ArgumentError("s length must be ≥ $(2 * r + 1)."))
    wlen >= 2 * r + 1 || throw(ArgumentError("wlen must be ≥ $(2 * r + 1)."))
    a_min >= 0 || throw(ArgumentError("a_min must be ≥ 0."))
    0 <= r2_min <= 1 || throw(ArgumentError("r2_min must be in [0, 1]."))

    # work on a copy (the input is not modified)
    s = Vector{Float64}(s)
    n = length(s)
    k = max(3, r ÷ 10)

    # candidates: keep those large enough and following the model (evaluated on the original signal)
    cand = _pop_detect(s; z = z, r = r, wlen = wlen)
    pop_loc = Int64[]
    for (i, p) in enumerate(cand)
        p_end = i < length(cand) ? cand[i + 1] - 1 : n
        x = @view(s[p:min(p + r - 1, p_end)]) .- mean(@view s[max(1, p - k):(p - 1)])
        f = _pop_fit(x)
        # goodness of fit where the pop is present (3τ), not over the following signal
        l = clamp(round(Int64, 3 * f.tau), k, length(x))
        xe = @view x[1:l]
        ss = sum(abs2, xe .- mean(xe))
        r2 = ss > 0 ? 1 - sum(abs2, xe .- f.a .* exp.(-(0:(l - 1)) ./ f.tau)) / ss : 0.0
        w1 = clamp(p - wlen ÷ 2, 1, max(1, n - wlen + 1))
        sd = 1.4826 * mad(@view(s[w1:min(n, w1 + wlen - 1)]); normalize = false)
        abs(f.a) >= a_min * sd || continue
        # τ at the upper limit: a step that does not return to its level (e.g. amplifier saturation, lost electrode)
        if f.tau > 0.99 * 4 * length(x)
            _warn("Step at sample $p does not recover within $(length(x)) samples (not a pop, not repaired); check this segment (e.g. amplifier saturation).")
            continue
        end
        r2 >= r2_min && push!(pop_loc, p)
    end
    amp = zeros(length(pop_loc))
    tau = zeros(length(pop_loc))

    for (i, p) in enumerate(pop_loc)
        # level before the pop (already repaired if there was an earlier pop)
        b = mean(@view s[max(1, p - k):(p - 1)])
        # fit within r samples, subtract up to the next pop or the end of the signal
        p_end = i < length(pop_loc) ? pop_loc[i + 1] - 1 : n
        x = @view(s[p:min(p + r - 1, p_end)]) .- b
        f = _pop_fit(x)
        amp[i], tau[i] = f.a, f.tau
        repair && (s[p:p_end] .-= f.a .* exp.(-(0:(p_end - p)) ./ f.tau))
    end

    return (; s, pop_loc, amp, tau)
end

"""
    remove_pops(obj; <keyword arguments>)

Detect and repair electrode pops: abrupt steps of the signal followed by a slow (exponential) return to the previous level, caused by a sudden change of the electrode-skin contact.

Candidates are outliers of the first difference of the signal (robust z-score > `z`, median and MAD calculated within `window`) after which the signal stays at the new level (single-sample spikes are not pops). Each candidate is modeled as `amplitude × exp(-t / tau)`, fitted to `r` samples after the jump; it is accepted as a pop if `|amplitude|` ≥ `a_min` × robust SD of the signal and the model explains the signal (R² ≥ `r2_min`), so that blinks, muscle and other fast activity are not treated as pops. The fitted pops are subtracted, so the underlying signal is preserved.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object (continuous)
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `repair::Bool=true`: repair the pops if `true`, only detect them if `false`
- `window::Real=10.0`: window length (in seconds) for the robust statistics of the detection
- `r::Int64=sr(obj)÷2`: fitting segment length (in samples) after each pop, also the minimum distance between pops
- `z::Real=10`: detection threshold (robust z-score of the first difference)
- `a_min::Real=5`: minimum pop amplitude relative to the robust SD of the signal (within `window`)
- `r2_min::Real=0.8`: minimum R² of the exponential model (lower values accept pops of other shapes)

# Returns

- `NeuroAnalyzer.NEURO`: repaired NEURO object (only if `repair=true`)
- `DataFrame`: detected pops: `:channel` (label), `:sample`, `:time` (s), `:amplitude` (channel units), `:tau` (recovery time constant, s)

# Notes

It is recommended to remove the DC offset (`remove_dc()`) and slow drifts (high-pass filter) before detecting pops.

Large steps that do not return to the previous level within `r` samples (e.g. amplifier saturation, a lost electrode) are not pops: they are not repaired, a warning is shown and these segments should be checked (and rejected if necessary).
"""
function remove_pops(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    repair::Bool = true,
    window::Real = 10.0,
    r::Int64 = sr(obj) ÷ 2,
    z::Real = 10,
    a_min::Real = 5,
    r2_min::Real = 0.8,
)::Union{Tuple{NeuroAnalyzer.NEURO, DataFrame}, DataFrame}
    # validate
    nepochs(obj) == 1 || throw(ArgumentError("remove_pops() must be applied to continuous object."))
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    wlen = round(Int64, window * sr(obj))
    2 * r + 1 <= wlen <= signal_len(obj) ||
        throw(ArgumentError("window must be ≥ $((2 * r + 1) / sr(obj)) and ≤ $(signal_len(obj) / sr(obj)) s."))

    obj_tmp = repair ? deepcopy(obj) : obj
    pops = DataFrame(
        :channel => String[],
        :sample => Int64[],
        :time => Float64[],
        :amplitude => Float64[],
        :tau => Float64[],
    )

    for c in ch
        p = _remove_pops(@view(obj.data[c, :, 1]); r = r, z = z, wlen = wlen, a_min = a_min, r2_min = r2_min, repair = repair)
        repair && (obj_tmp.data[c, :, 1] = p.s)
        for i in eachindex(p.pop_loc)
            push!(pops, (labels(obj)[c], p.pop_loc[i], obj.time_pts[p.pop_loc[i]], p.amp[i], p.tau[i] / sr(obj)))
        end
    end

    repair || return pops
    push!(obj_tmp.history, "remove_pops(obj; ch=$ch, repair=true, window=$window, r=$r, z=$z, a_min=$a_min, r2_min=$r2_min)")

    return obj_tmp, pops
end

"""
    remove_pops!(obj; <keyword arguments>)

Detect and repair electrode pops in-place (see `remove_pops()`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object (continuous)
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `repair::Bool=true`: repair the pops if `true`, only detect them if `false`
- `window::Real=10.0`: window length (in seconds) for the robust statistics of the detection
- `r::Int64=sr(obj)÷2`: fitting segment length (in samples) after each pop, also the minimum distance between pops
- `z::Real=10`: detection threshold (robust z-score of the first difference)
- `a_min::Real=5`: minimum pop amplitude relative to the robust SD of the signal (within `window`)
- `r2_min::Real=0.8`: minimum R² of the exponential model (lower values accept pops of other shapes)

# Returns

- `DataFrame`: detected pops: `:channel` (label), `:sample`, `:time` (s), `:amplitude` (channel units), `:tau` (recovery time constant, s)
"""
function remove_pops!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    repair::Bool = true,
    window::Real = 10.0,
    r::Int64 = sr(obj) ÷ 2,
    z::Real = 10,
    a_min::Real = 5,
    r2_min::Real = 0.8,
)::DataFrame
    repair || return remove_pops(obj; ch = ch, repair = false, window = window, r = r, z = z, a_min = a_min, r2_min = r2_min)
    obj_tmp, pops = remove_pops(obj; ch = ch, repair = true, window = window, r = r, z = z, a_min = a_min, r2_min = r2_min)
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history

    return pops
end
