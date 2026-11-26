export p_adjust
export mass_ttest
export mass_anova
export mass_anova_rm

"""
    p_adjust(p; <keyword arguments>)

Adjust p-values for multiple comparisons. `NaN` values are ignored (and kept).

# Arguments

- `p::AbstractArray{<:Real}`: p-values (any shape, e.g. channels × time points)
- `method::Symbol=:fdr`: method:
    - `:bonferroni`: Bonferroni (family-wise error rate)
    - `:holm`: Holm step-down (family-wise error rate, more powerful than Bonferroni)
    - `:fdr`: Benjamini-Hochberg (false discovery rate, independent or positively dependent tests)
    - `:by`: Benjamini-Yekutieli (false discovery rate, any dependence)

# Returns

- `Array{Float64}`: adjusted p-values (same shape as `p`), at most 1

# Examples

```julia
r = mass_ttest(x, y)
p_fdr = p_adjust(r.p; method = :fdr)
```
"""
function p_adjust(p::AbstractArray{<:Real}; method::Symbol = :fdr)::Array{Float64}
    _check_var(method, [:bonferroni, :holm, :fdr, :by], "method")
    out = fill(NaN, size(p))
    idx = findall(!isnan, p)
    m = length(idx)
    m == 0 && return out
    v = Float64[p[i] for i in idx]
    all(x -> 0 <= x <= 1, v) || throw(ArgumentError("p-values must be in [0, 1]."))
    adj = similar(v)
    if method === :bonferroni
        adj .= min.(1.0, v .* m)
    elseif method === :holm
        o = sortperm(v)
        running = 0.0
        for k = 1:m
            running = max(running, (m - k + 1) * v[o[k]])
            adj[o[k]] = min(1.0, running)
        end
    else
        # Benjamini-Hochberg; Benjamini-Yekutieli: multiplied by the harmonic sum
        c = method === :by ? sum(1 / k for k = 1:m) : 1.0
        o = sortperm(v)
        running = 1.0
        for k = m:-1:1
            running = min(running, v[o[k]] * m * c / k)
            adj[o[k]] = min(1.0, running)
        end
    end
    for (k, i) in enumerate(idx)
        out[i] = adj[k]
    end
    return out
end

"""
    mass_ttest(x, y; <keyword arguments>)

Mass-univariate t-test: a t-test at every point (e.g. channel × time point) of the data. Observations (subjects, epochs) are along the first dimension; the other dimensions are tested point by point. `NaN` values (e.g. bad channels) are left out of the test of their point (for the paired test: the pair).

- `mass_ttest(x)`: one-sample t-test (`mean(x) = mu`)
- `mass_ttest(x, y)`: independent two-sample t-test (Student, `equal_var = true`, or Welch)
- `mass_ttest(x, y; paired = true)`: paired t-test (one-sample test of `x - y`)

The p-values are not corrected for multiple comparisons: see `p_adjust()` and `cluster_perm()`.

# Arguments

- `x::AbstractArray{<:Real}`: observations × points (e.g. subjects × channels × time points)
- `y::Union{Nothing, AbstractArray{<:Real}}=nothing`: second sample (the same shape, except the number of observations of an independent sample)
- `paired::Bool=false`: if `true`, paired t-test
- `mu::Real=0`: mean under the null hypothesis (one-sample and paired tests)
- `equal_var::Bool=true`: if `false`, Welch's t-test (independent samples)
- `tail::Symbol=:both`: alternative hypothesis: `:both` (two-tailed), `:right` (`x > y`), `:left` (`x < y`)

# Returns

Named tuple:

- `t::Array{Float64}`: t statistics, shape: the dimensions of `x` without the first
- `df::Array{Float64}`: degrees of freedom
- `p::Array{Float64}`: p-values (uncorrected)
"""
function mass_ttest(
    x::AbstractArray{<:Real},
    y::Union{Nothing, AbstractArray{<:Real}} = nothing;
    paired::Bool = false,
    mu::Real = 0,
    equal_var::Bool = true,
    tail::Symbol = :both,
)::@NamedTuple{t::Array{Float64}, df::Array{Float64}, p::Array{Float64}}
    _check_var(tail, [:both, :right, :left], "tail")
    sz = size(x)[2:end]
    if isnothing(y)
        t, df = _mt_t1(_mt_mat(x), mu)
    else
        size(y)[2:end] == sz ||
            throw(ArgumentError("x and y must have the same dimensions (except the first)."))
        if paired
            size(y, 1) == size(x, 1) ||
                throw(ArgumentError("Paired samples must have the same number of observations."))
            t, df = _mt_t1(_mt_mat(x) .- _mt_mat(y), mu)
        else
            mu == 0 || throw(ArgumentError("mu is used only by one-sample and paired tests."))
            t, df = _mt_t2(_mt_mat(x), _mt_mat(y); equal_var = equal_var)
        end
    end
    p = _mt_tp(t, df, tail)
    return (t = reshape(t, sz), df = reshape(df, sz), p = reshape(p, sz))
end

"""
    mass_anova(x1, x2, xs...)

Mass-univariate one-way ANOVA (between groups): an F-test at every point (e.g. channel × time point) of the data. Observations (subjects) are along the first dimension; the other dimensions are tested point by point. `NaN` values are left out of the test of their point.

# Arguments

- `x1::AbstractArray{<:Real}, x2::AbstractArray{<:Real}, xs::AbstractArray{<:Real}...`: groups: observations × points (the same dimensions, except the number of observations)

# Returns

Named tuple:

- `f::Array{Float64}`: F statistics, shape: the dimensions of `x1` without the first
- `df1::Array{Float64}`, `df2::Array{Float64}`: degrees of freedom
- `p::Array{Float64}`: p-values (uncorrected)
"""
function mass_anova(
    x1::AbstractArray{<:Real}, x2::AbstractArray{<:Real}, xs::AbstractArray{<:Real}...,
)::@NamedTuple{f::Array{Float64}, df1::Array{Float64}, df2::Array{Float64}, p::Array{Float64}}
    g = (x1, x2, xs...)
    sz = size(x1)[2:end]
    all(size(x)[2:end] == sz for x in g) ||
        throw(ArgumentError("Groups must have the same dimensions (except the first)."))
    f, df1, df2 = _mt_f1([_mt_mat(x) for x in g])
    p = _mt_fp(f, df1, df2)
    return (f = reshape(f, sz), df1 = reshape(df1, sz), df2 = reshape(df2, sz),
            p = reshape(p, sz))
end

"""
    mass_anova_rm(x)

Mass-univariate one-way repeated-measures ANOVA (within subjects): an F-test of the conditions at every point (e.g. channel × time point) of the data. Subjects with `NaN` in any condition are left out of the test of that point.

The F-test assumes sphericity (no Greenhouse-Geisser correction); `cluster_perm(x; test = :anova_rm)` does not.

# Arguments

- `x::AbstractArray{<:Real}`: subjects × conditions × points (e.g. subjects × conditions × channels × time points)

# Returns

Named tuple:

- `f::Array{Float64}`: F statistics, shape: the dimensions of `x` without the first two
- `df1::Array{Float64}`, `df2::Array{Float64}`: degrees of freedom
- `p::Array{Float64}`: p-values (uncorrected)
"""
function mass_anova_rm(
    x::AbstractArray{<:Real},
)::@NamedTuple{f::Array{Float64}, df1::Array{Float64}, df2::Array{Float64}, p::Array{Float64}}
    ndims(x) >= 2 || throw(ArgumentError("x must have at least 2 dimensions (subjects × conditions)."))
    size(x, 2) >= 2 || throw(ArgumentError("x must have at least 2 conditions."))
    sz = size(x)[3:end]
    f, df1, df2 = _mt_frm(reshape(collect(Float64, x), size(x, 1), size(x, 2), :))
    p = _mt_fp(f, df1, df2)
    return (f = reshape(f, sz), df1 = reshape(df1, sz), df2 = reshape(df2, sz),
            p = reshape(p, sz))
end

# observations × points matrix
_mt_mat(x::AbstractArray{<:Real})::Matrix{Float64} = reshape(collect(Float64, x), size(x, 1), :)

# one-sample t statistic and degrees of freedom of every column
function _mt_t1(M::AbstractMatrix{Float64}, mu::Real)
    n, np = size(M)
    t = fill(NaN, np)
    df = fill(NaN, np)
    if !any(isnan, M)
        n < 2 && return t, df
        m = vec(mean(M; dims = 1))
        s = vec(std(M; dims = 1))
        @. t = ifelse(s > 0, (m - mu) / (s / sqrt(n)), NaN)
        df .= n - 1
        return t, df
    end
    @inbounds for j = 1:np
        v = Base.filter(!isnan, @view M[:, j])
        k = length(v)
        k < 2 && continue
        s = std(v)
        s > 0 || continue
        t[j] = (mean(v) - mu) / (s / sqrt(k))
        df[j] = k - 1
    end
    return t, df
end

# two-sample t statistic and degrees of freedom of every column
function _mt_t2(A::AbstractMatrix{Float64}, B::AbstractMatrix{Float64}; equal_var::Bool = true)
    np = size(A, 2)
    t = fill(NaN, np)
    df = fill(NaN, np)
    function one(a, b)
        n1, n2 = length(a), length(b)
        (n1 < 2 || n2 < 2) && return NaN, NaN
        m1, m2 = mean(a), mean(b)
        v1, v2 = var(a), var(b)
        if equal_var
            sp = ((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2)
            se = sqrt(sp * (1 / n1 + 1 / n2))
            d = n1 + n2 - 2.0
        else
            se = sqrt(v1 / n1 + v2 / n2)
            d = (v1 / n1 + v2 / n2)^2 /
                ((v1 / n1)^2 / (n1 - 1) + (v2 / n2)^2 / (n2 - 1))
        end
        se > 0 || return NaN, NaN
        return (m1 - m2) / se, d
    end
    if !any(isnan, A) && !any(isnan, B) && equal_var
        n1, n2 = size(A, 1), size(B, 1)
        (n1 < 2 || n2 < 2) && return t, df
        m1 = vec(mean(A; dims = 1))
        m2 = vec(mean(B; dims = 1))
        v1 = vec(var(A; dims = 1))
        v2 = vec(var(B; dims = 1))
        se = @. sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2) * (1 / n1 + 1 / n2))
        @. t = ifelse(se > 0, (m1 - m2) / se, NaN)
        df .= n1 + n2 - 2
        return t, df
    end
    @inbounds for j = 1:np
        t[j], df[j] = one(Base.filter(!isnan, @view A[:, j]), Base.filter(!isnan, @view B[:, j]))
    end
    return t, df
end

# one-way ANOVA F statistic and degrees of freedom of every column
function _mt_f1(G::Vector{Matrix{Float64}})
    np = size(G[1], 2)
    k = length(G)
    f = fill(NaN, np)
    df1 = fill(NaN, np)
    df2 = fill(NaN, np)
    if !any(M -> any(isnan, M), G)
        ns = [size(M, 1) for M in G]
        N = sum(ns)
        N > k || return f, df1, df2
        ms = [vec(mean(M; dims = 1)) for M in G]
        gm = sum(ns[i] .* ms[i] for i = 1:k) ./ N
        ssb = sum(ns[i] .* (ms[i] .- gm) .^ 2 for i = 1:k)
        ssw = sum(vec(sum((G[i] .- ms[i]') .^ 2; dims = 1)) for i = 1:k)
        @. f = ifelse(ssw > 0, (ssb / (k - 1)) / (ssw / (N - k)), NaN)
        df1 .= k - 1
        df2 .= N - k
        return f, df1, df2
    end
    @inbounds for j = 1:np
        vs = [Base.filter(!isnan, @view M[:, j]) for M in G]
        vs = [v for v in vs if !isempty(v)]
        kk = length(vs)
        N = sum(length, vs; init = 0)
        (kk < 2 || N <= kk) && continue
        gm = sum(sum, vs) / N
        ssb = sum(length(v) * (mean(v) - gm)^2 for v in vs)
        ssw = sum(sum((v .- mean(v)) .^ 2) for v in vs)
        ssw > 0 || continue
        f[j] = (ssb / (kk - 1)) / (ssw / (N - kk))
        df1[j] = kk - 1
        df2[j] = N - kk
    end
    return f, df1, df2
end

# repeated-measures ANOVA F statistic and degrees of freedom of every column of a
# subjects × conditions × points array
function _mt_frm(X::AbstractArray{Float64, 3})
    n, k, np = size(X)
    f = fill(NaN, np)
    df1 = fill(NaN, np)
    df2 = fill(NaN, np)
    function one(Y::AbstractMatrix)
        nn = size(Y, 1)
        nn < 2 && return NaN, NaN, NaN
        gm = mean(Y)
        ssc = nn * sum((vec(mean(Y; dims = 1)) .- gm) .^ 2)
        sss = k * sum((vec(mean(Y; dims = 2)) .- gm) .^ 2)
        sst = sum((Y .- gm) .^ 2)
        sse = sst - ssc - sss
        d1, d2 = k - 1.0, (k - 1.0) * (nn - 1)
        sse > 0 || return NaN, d1, d2
        return (ssc / d1) / (sse / d2), d1, d2
    end
    @inbounds for j = 1:np
        Y = @view X[:, :, j]
        if any(isnan, Y)
            ok = [!any(isnan, @view Y[s, :]) for s = 1:n]
            f[j], df1[j], df2[j] = one(Y[ok, :])
        else
            f[j], df1[j], df2[j] = one(Y)
        end
    end
    return f, df1, df2
end

# p-values of t statistics
function _mt_tp(t::AbstractVector, df::AbstractVector, tail::Symbol)::Vector{Float64}
    return map(t, df) do x, d
        (isnan(x) || isnan(d)) && return NaN
        dist = Distributions.TDist(d)
        tail === :both ? 2 * ccdf(dist, abs(x)) : tail === :right ? ccdf(dist, x) : cdf(dist, x)
    end
end

# p-values of F statistics
function _mt_fp(f::AbstractVector, df1::AbstractVector, df2::AbstractVector)::Vector{Float64}
    return map(f, df1, df2) do x, d1, d2
        (isnan(x) || isnan(d1) || isnan(d2)) && return NaN
        ccdf(Distributions.FDist(d1, d2), x)
    end
end
