export outlier_detect
export grubbs

# internal helper: given a local index in the shortened vector, recover the
# original index in x by accounting for previously deleted positions.
function _original_index(local_idx::Int, removed::Vector{Int}, n::Int)::Int
    pos = 0
    count = 0
    for i in 1:n
        i in removed && continue
        count += 1
        count == local_idx && return i
    end
    # fallback (should not be reached)
    return local_idx
end

"""
    outlier_detect(x; <keyword arguments>)

Detect outliers in a vector using a selected method.

# Arguments

- `x::AbstractVector`: input vector; must contain at least 7 elements when using `method = :g`
- `method::Symbol=:iqr`: detection method:
    - `:iqr`: interquartile range (flags values outside `Q1 − 1.5×IQR` and `Q3 + 1.5×IQR`)
    - `:z`: Z-score (flags values with `|z| > 3`)
    - `:g`: iterative Grubbs test applied to both tails

# Returns

- `Vector{Bool}`: boolean mask; `true` at each index identified as an outlier

# Throws

- `ArgumentError`: if `x` is empty, `method` is invalid, or `method = :g` with `length(x) ≤ 6`

# Notes

- The `:g` method iteratively removes the most extreme value while `length(x_tmp) > 6` (Grubbs requires at least 7 observations), checking both the maximum and minimum tails separately.
- Outlier indices are tracked relative to the **original** vector, so deletion from the working copy does not affect index mapping.

# See also

[`grubbs`](@ref)
"""
function outlier_detect(x::AbstractVector; method::Symbol = :iqr)::Vector{Bool}

    @assert length(x) > 0 "x must not be empty."
    _check_var(method, [:iqr, :z, :g], "method")

    o = zeros(Bool, length(x))

    if method === :iqr

        # compute iqr once
        iqr_val = iqr(x)
        lo = quantile(x, 0.25) - 1.5 * iqr_val
        hi = quantile(x, 0.75) + 1.5 * iqr_val
        o .= (x .< lo) .| (x .> hi)

    elseif method === :z

        z  = z_score(x)
        o .= (z .< -3) .| (z .> 3)

    elseif method === :g

        # iterative Grubbs test, both tails
        @assert length(x) > 6 "method = :g requires length(x) > 6."

        # upper-tail pass: repeatedly test whether the current maximum is an outlier
        # working copy; deletions track removals
        x_tmp = collect(Float64, x)
        # original indices already removed
        removed = Int[]
        for _ in (length(x_tmp)):-1:7
            m_idx_local = argmax(x_tmp)
            if grubbs(x_tmp; t=1)
                # map local index back to the original index
                orig_idx = _original_index(m_idx_local, removed, length(x))
                o[orig_idx] = true
                push!(removed, orig_idx)
                deleteat!(x_tmp, m_idx_local)
            else
                # Grubbs is non-significant; no further removals needed
                break
            end
        end

        # lower-tail pass: reset working copy and repeat for the minimum
        x_tmp   = collect(Float64, x)
        removed = Int[]
        for _ in (length(x_tmp)):-1:7
            m_idx_local = argmin(x_tmp)
            if grubbs(x_tmp; t=-1)
                orig_idx = _original_index(m_idx_local, removed, length(x))
                o[orig_idx] = true
                push!(removed, orig_idx)
                deleteat!(x_tmp, m_idx_local)
            else
                break
            end
        end
    end

    return o

end

"""
    grubbs(x; <keyword arguments>)

Perform the Grubbs test for a single outlier.

Tests whether the most extreme value in `x` (maximum, minimum, or both) is a statistically significant outlier at confidence level `alpha`.

The test statistic is:

- Two-tailed: `G = max|xᵢ − x̄| / s`
- Upper one-tailed: `G = (max(x) − x̄) / s`
- Lower one-tailed: `G = (x̄ − min(x)) / s`

The critical value is derived from the t-distribution with `df = n − 2`.

# Arguments
- `x::AbstractVector`: input vector; must contain at least 7 elements
- `alpha::Float64=0.95`: confidence level (not significance level); must be in `(0, 1)`; default 0.95 corresponds to a 5 % significance level
- `t::Int64=0`: test type:
    - `0`: two-tailed (tests both extremes)
    - `1`: upper one-tailed (tests maximum)
    - `-1`: lower one-tailed (tests minimum)

# Returns

- `Bool`: `true` if an outlier is detected; `false` otherwise

# Throws

- `ArgumentError`: if `length(x) < 7`, `alpha ∉ (0, 1)`, or `t ∉ {−1, 0, 1}`

# References
Grubbs FE. Procedures for detecting outlying observations in samples. Technometrics. 1969;11(1):1–21.

# See also

[`outlier_detect`](@ref)
"""
function grubbs(x::AbstractVector; alpha::Float64 = 0.95, t::Int64 = 0)::Bool

    @assert length(x) >= 7 "x must contain at least 7 elements."
    @assert alpha > 0.0 "alpha must be > 0."
    @assert alpha < 1.0 "alpha must be < 1."
    @assert t in (-1, 0, 1) "t must be -1, 0, or 1."

    n  = length(x)
    df = n - 2
    m  = mean(x)
    s  = std(x)

    two_tailed, g = if t == 0
        true, maximum(abs.(x .- m)) / s
    elseif t == 1
        false, (maximum(x) - m) / s
    elseif t == -1
        false, (m - minimum(x)) / s
    end

    # significance level for the t-distribution lookup
    p = two_tailed ? (1 - alpha) / (2 * n) : (1 - alpha) / n
    t_crit = quantile(TDist(df), 1 - p)
    threshold = (n - 1) / sqrt(n) * sqrt(t_crit^2 / (df + t_crit^2))

    return g >= threshold

end
