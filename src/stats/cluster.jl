export cluster_perm
export channel_adjacency

"""
    cluster_perm(xs...; <keyword arguments>)

Cluster-based permutation test (Maris & Oostenveld 2007) of data with observations (subjects, epochs) along the first dimension, e.g. subjects × channels × time points. It controls the family-wise error rate over all points without assuming independent tests:

1. A test statistic (t or F, see `mass_ttest()`, `mass_anova()`, `mass_anova_rm()`) is calculated at every point.
2. Points whose statistic exceeds the cluster-forming threshold (the critical value at `cluster_alpha`) are grouped into clusters of adjacent points (separately for positive and negative t values); the cluster mass is the sum of the statistic over its points.
3. The data are permuted `nperm` times (sign flips for one-sample and paired tests, shuffled group labels for independent samples, shuffled conditions within subjects for repeated measures), and the largest cluster mass of each permutation forms the null distribution.
4. The p-value of a cluster is the proportion of permutations with a larger maximal cluster mass.

Adjacency: if the points have one dimension (e.g. time points), neighbouring points are adjacent. If they have more dimensions, the first one is channels: channels are adjacent as given by `ch_adj` (see `channel_adjacency()`; default: none, each channel forms its own clusters); in the other dimensions (time points, frequencies) neighbouring points are adjacent.

Tests:

- `:ttest` with one array: one-sample t-test (`mean = mu`)
- `:ttest` with two arrays: independent t-test, or paired t-test (`paired = true`)
- `:anova` with two or more arrays: one-way ANOVA between groups
- `:anova_rm` with one array subjects × conditions × points: one-way repeated-measures ANOVA

# Arguments

- `xs::AbstractArray{<:Real}...`: data: observations × points
- `test::Symbol=:auto`: test (see above); `:auto`: `:ttest` for one or two arrays, `:anova` for more
- `paired::Bool=false`: paired t-test (`:ttest` with two arrays)
- `mu::Real=0`: mean under the null hypothesis (one-sample t-test)
- `tail::Symbol=:both`: alternative hypothesis of t-tests: `:both`, `:right` (positive clusters), `:left` (negative clusters); F-tests are right-tailed
- `nperm::Int64=1000`: number of permutations
- `cluster_alpha::Real=0.05`: cluster-forming threshold (point-wise p-value; two-tailed t-tests: `cluster_alpha / 2` in each tail)
- `alpha::Real=0.05`: significance level of the clusters
- `ch_adj::Union{Nothing, AbstractMatrix{Bool}}=nothing`: channel adjacency (channels × channels) of the first dimension of the points
- `rng::AbstractRNG=Random.default_rng()`: random number generator (e.g. `Random.Xoshiro(1)` for reproducible results)

# Returns

Named tuple:

- `stat::Array{Float64}`: t or F statistic at every point
- `p::Array{Float64}`: point-wise p-values (uncorrected)
- `threshold::Float64`: cluster-forming threshold of the statistic
- `clusters::DataFrame`: one row per cluster (sorted by p-value): `id`, `sign` (+1, -1), `mass`, `size` (number of points), `p`
- `cluster_id::Array{Int64}`: cluster of every point (0: none)
- `p_cluster::Array{Float64}`: p-value of the cluster of every point (`NaN` outside clusters)
- `mask::BitArray`: points of the significant clusters (`p < alpha`)
- `null::Vector{Float64}`: null distribution of the maximal cluster mass

# Notes

A significant cluster indicates that the conditions or groups differ; its extent in channels and time is a description, not a test of where or when the effect starts (Sassenhagen & Draschkow 2019).

# References

- Maris E, Oostenveld R. Nonparametric statistical testing of EEG- and MEG-data. J Neurosci Methods. 2007;164(1):177-190.
- Sassenhagen J, Draschkow D. Cluster-based permutation tests of MEG/EEG data do not establish significance of effect latency or location. Psychophysiology. 2019;56(6):e13335.
"""
function cluster_perm(
    xs::AbstractArray{<:Real}...;
    test::Symbol = :auto,
    paired::Bool = false,
    mu::Real = 0,
    tail::Symbol = :both,
    nperm::Int64 = 1000,
    cluster_alpha::Real = 0.05,
    alpha::Real = 0.05,
    ch_adj::Union{Nothing, AbstractMatrix{Bool}} = nothing,
    rng::AbstractRNG = Random.default_rng(),
)::NamedTuple
    # validate
    isempty(xs) && throw(ArgumentError("No data."))
    _check_var(test, [:auto, :ttest, :anova, :anova_rm], "test")
    _check_var(tail, [:both, :right, :left], "tail")
    nperm >= 1 || throw(ArgumentError("nperm must be ≥ 1."))
    0 < cluster_alpha < 1 || throw(ArgumentError("cluster_alpha must be in (0, 1)."))
    0 < alpha < 1 || throw(ArgumentError("alpha must be in (0, 1)."))
    test === :auto && (test = length(xs) <= 2 ? :ttest : :anova)
    if test === :ttest
        length(xs) <= 2 || throw(ArgumentError("t-test requires one or two arrays."))
    elseif test === :anova
        length(xs) >= 2 || throw(ArgumentError("ANOVA requires at least two arrays."))
    else
        length(xs) == 1 ||
            throw(ArgumentError("Repeated-measures ANOVA requires one array (subjects × conditions × points)."))
        ndims(xs[1]) >= 2 && size(xs[1], 2) >= 2 ||
            throw(ArgumentError("Repeated-measures ANOVA requires at least 2 conditions."))
    end
    ftest = test !== :ttest
    ftest && tail !== :both && throw(ArgumentError("F-tests are right-tailed: tail must be :both."))
    sz = test === :anova_rm ? size(xs[1])[3:end] : size(xs[1])[2:end]
    for x in xs
        (test === :anova_rm ? size(x)[3:end] : size(x)[2:end]) == sz ||
            throw(ArgumentError("Data must have the same dimensions (except observations)."))
    end
    if !isnothing(ch_adj)
        length(sz) >= 2 || throw(ArgumentError("ch_adj requires points with at least 2 dimensions (channels × ...)."))
        size(ch_adj) == (sz[1], sz[1]) ||
            throw(ArgumentError("ch_adj must be a $(sz[1]) × $(sz[1]) matrix (channels × channels)."))
    end

    # data and statistic of a permutation
    if test === :ttest && (length(xs) == 1 || paired)
        D = length(xs) == 1 ? _mt_mat(xs[1]) .- mu : _mt_mat(xs[1]) .- _mt_mat(xs[2])
        paired && size(xs[1], 1) != size(xs[2], 1) &&
            throw(ArgumentError("Paired samples must have the same number of observations."))
        n = size(D, 1)
        statf = perm -> _mt_t1(isnothing(perm) ? D : D .* perm, 0)
        permf = () -> rand(rng, (-1.0, 1.0), n)
    elseif test === :ttest
        A, B = _mt_mat(xs[1]), _mt_mat(xs[2])
        n1 = size(A, 1)
        AB = vcat(A, B)
        statf = perm -> isnothing(perm) ? _mt_t2(A, B) :
                        _mt_t2(AB[perm[1:n1], :], AB[perm[(n1 + 1):end], :])
        permf = () -> randperm(rng, size(AB, 1))
    elseif test === :anova
        G = [_mt_mat(x) for x in xs]
        ns = cumsum([0; [size(M, 1) for M in G]])
        GG = vcat(G...)
        statf = perm -> isnothing(perm) ? _mt_f1(G)[1:2] :
                        _mt_f1([GG[perm[(ns[i] + 1):ns[i + 1]], :] for i in eachindex(G)])[1:2]
        permf = () -> randperm(rng, size(GG, 1))
    else
        X = reshape(collect(Float64, xs[1]), size(xs[1], 1), size(xs[1], 2), :)
        n, k = size(X, 1), size(X, 2)
        statf = perm -> begin
            isnothing(perm) && return _mt_frm(X)[1:2]
            Y = similar(X)
            for s = 1:n
                Y[s, :, :] = X[s, perm[s], :]
            end
            _mt_frm(Y)[1:2]
        end
        permf = () -> [randperm(rng, k) for _ = 1:n]
    end

    # observed statistic and cluster-forming threshold (degrees of freedom of the
    # complete data)
    if ftest
        st = statf(nothing)[1]
        d1 = test === :anova ? length(xs) - 1 : size(xs[1], 2) - 1
        d2 = test === :anova ? sum(size(x, 1) for x in xs) - length(xs) :
             (size(xs[1], 2) - 1) * (size(xs[1], 1) - 1)
        d2 > 0 || throw(ArgumentError("Not enough observations."))
        thr = quantile(Distributions.FDist(d1, d2), 1 - cluster_alpha)
        p = test === :anova ? vec(mass_anova(xs...).p) : vec(mass_anova_rm(xs[1]).p)
    else
        st, df = statf(nothing)
        d = length(xs) == 1 || paired ? size(xs[1], 1) - 1 : size(xs[1], 1) + size(xs[2], 1) - 2
        d > 0 || throw(ArgumentError("Not enough observations."))
        thr = quantile(Distributions.TDist(d), 1 - (tail === :both ? cluster_alpha / 2 : cluster_alpha))
        p = _mt_tp(st, df, tail)
    end

    # observed clusters
    obs = _cp_clusters(st, sz, thr, ftest, tail, ch_adj)

    # null distribution of the maximal cluster mass (permutations drawn in advance, so
    # that the result does not depend on the number of threads)
    perms = [permf() for _ = 1:nperm]
    null = zeros(nperm)
    Threads.@threads :dynamic for i = 1:nperm
        s = statf(perms[i])[1]
        c = _cp_clusters(s, sz, thr, ftest, tail, ch_adj)
        null[i] = isempty(c.mass) ? 0.0 : maximum(abs, c.mass)
    end

    # cluster p-values
    nc = length(obs.mass)
    cp = [(1 + count(>=(abs(m) - 1e-10), null)) / (nperm + 1) for m in obs.mass]
    p_cluster = fill(NaN, length(st))
    for i in eachindex(st)
        obs.id[i] > 0 && (p_cluster[i] = cp[obs.id[i]])
    end
    clusters = DataFrame(
        :id   => collect(1:nc),
        :sign => obs.sign,
        :mass => obs.mass,
        :size => [count(==(c), obs.id) for c = 1:nc],
        :p    => cp,
    )
    sort!(clusters, [:p, DataFrames.order(:mass; by = abs, rev = true)])

    return (
        stat = reshape(st, sz),
        p = reshape(p, sz),
        threshold = thr,
        clusters = clusters,
        cluster_id = reshape(obs.id, sz),
        p_cluster = reshape(p_cluster, sz),
        mask = reshape(.!isnan.(p_cluster) .& (p_cluster .< alpha), sz),
        null = null,
    )
end

# clusters of the points exceeding the threshold: cluster number of every point, sign and
# mass of every cluster
function _cp_clusters(
    st::AbstractVector{Float64},
    sz::Tuple,
    thr::Real,
    ftest::Bool,
    tail::Symbol,
    ch_adj::Union{Nothing, AbstractMatrix{Bool}},
)
    id = zeros(Int64, length(st))
    sgn = Int64[]
    mass = Float64[]
    signs = ftest ? (1,) : tail === :both ? (1, -1) : tail === :right ? (1,) : (-1,)
    for s in signs
        above = [!isnan(v) && s * v > thr for v in st]
        lab, n = _cp_label(above, sz, ch_adj)
        for c = 1:n
            pts = findall(==(c), lab)
            push!(sgn, s)
            push!(mass, sum(st[pts]))
            id[pts] .= length(mass)
        end
    end
    return (id = id, sign = sgn, mass = mass)
end

# connected components of the points of `mask` (linear indices of an array of shape `sz`)
function _cp_label(
    mask::AbstractVector{Bool}, sz::Tuple, ch_adj::Union{Nothing, AbstractMatrix{Bool}},
)
    lab = zeros(Int64, length(mask))
    any(mask) || return lab, 0
    ci = CartesianIndices(sz)
    li = LinearIndices(sz)
    nd = length(sz)
    # the first dimension is channels if the points have ≥ 2 dimensions
    chdim = nd >= 2
    nbr_ch = chdim && !isnothing(ch_adj) ?
             [Base.filter(!=(c), findall(@view ch_adj[c, :])) for c = 1:sz[1]] : nothing
    n = 0
    stack = Int64[]
    for i in eachindex(mask)
        (mask[i] && lab[i] == 0) || continue
        n += 1
        lab[i] = n
        push!(stack, i)
        while !isempty(stack)
            j = pop!(stack)
            c = ci[j]
            # neighbours along time, frequency (and the only dimension)
            for d in (chdim ? (2:nd) : (1:nd)), step in (-1, 1)
                v = c[d] + step
                1 <= v <= sz[d] || continue
                k = li[CartesianIndex(ntuple(q -> q == d ? v : c[q], nd))]
                if mask[k] && lab[k] == 0
                    lab[k] = n
                    push!(stack, k)
                end
            end
            # neighbouring channels at the same time point (frequency)
            if !isnothing(nbr_ch)
                for ch in nbr_ch[c[1]]
                    k = li[CartesianIndex(ntuple(q -> q == 1 ? ch : c[q], nd))]
                    if mask[k] && lab[k] == 0
                        lab[k] = n
                        push!(stack, k)
                    end
                end
            end
        end
    end
    return lab, n
end

"""
    channel_adjacency(obj; <keyword arguments>)

Channel adjacency for `cluster_perm()`: channels are adjacent if the distance of their locations (on the unit sphere) is at most `dist`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object with channel locations
- `ch::Union{String, Vector{String}, Regex}`: channel name(s), in the order of the data
- `dist::Union{Nothing, Real}=nothing`: largest distance of adjacent channels (unit sphere); default: 1.5 × the median distance of the channels to their nearest channel

# Returns

- `BitMatrix`: channels × channels (`true`: adjacent)
"""
function channel_adjacency(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    dist::Union{Nothing, Real} = nothing,
)::BitMatrix
    _has_locs(obj)
    chi = ch isa Vector{String} ? [findfirst(==(c), labels(obj)) for c in ch] :
          _get_ch(obj, ch, false)
    any(isnothing, chi) && throw(ArgumentError("Channels not found: $(join(ch[isnothing.(chi)], ", "))."))
    lbl = labels(obj)[chi]
    locs = _locs_unit(obj.locs)
    li = [findfirst(==(l), locs.label) for l in lbl]
    any(isnothing, li) &&
        throw(ArgumentError("Channels without locations: $(join(lbl[isnothing.(li)], ", "))."))
    xyz = hcat(locs.loc_x[li], locs.loc_y[li], locs.loc_z[li])
    n = length(li)
    d = [sqrt(sum((xyz[i, :] .- xyz[j, :]) .^ 2)) for i = 1:n, j = 1:n]
    if isnothing(dist)
        n > 1 || return trues(1, 1)
        dist = 1.5 * median([minimum(d[i, j] for j = 1:n if j != i) for i = 1:n])
    end
    dist > 0 || throw(ArgumentError("dist must be > 0."))
    return BitMatrix(d .<= dist)
end
