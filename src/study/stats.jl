export study_stats
export cluster_summary

"""
    study_stats(study, kind; <keyword arguments>)

Group statistics of a study measure (see `study_measure!()`): a test at every point (channel × time point, channel × frequency or channel × frequency × time point), corrected for multiple comparisons. The recordings of one subject are averaged first (see `study_average()`), so each subject is one observation; bad channels (`NaN`) are left out of the test of their points.

Factors:

- a design column, e.g. `:group` (between subjects): independent t-test (2 levels) or one-way ANOVA (more levels), for one condition
- `:condition` (within subjects): paired t-test (2 conditions) or repeated-measures ANOVA (more conditions)
- `:session` (within subjects): paired t-test (2 sessions) or repeated-measures ANOVA (more sessions), for one condition; recordings are not averaged over sessions

# Arguments

- `study::NeuroAnalyzer.STUDY`
- `kind::Symbol`: measure (`:erp`, `:psd` or `:tf`)
- `factor::Symbol=:group`: factor (see above)
- `levels::Union{Nothing, Vector{String}}=nothing`: levels to compare (e.g. two of three groups, conditions or sessions), in this order; `nothing`: all levels
- `condition::Union{Nothing, String}=nothing`: condition (between-subject factors and `:session`); may be omitted if the measure has one condition
- `group::Union{Nothing, String}=nothing`: subjects of one group (within-subject factors); `nothing`: all subjects
- `correction::Symbol=:cluster`: correction for multiple comparisons: `:cluster` (cluster-based permutation test, see `cluster_perm()`), `:fdr`, `:by`, `:holm`, `:bonferroni` (see `p_adjust()`) or `:none`
- `alpha::Real=0.05`: significance level
- `db::Bool=true`: for `:psd` and `:tf`, test the power of each subject in dB
- `tail::Symbol=:both`: alternative hypothesis of t-tests: `:both`, `:right` (first level > second level), `:left`
- `nperm::Int64=1000`: number of permutations (`:cluster`)
- `cluster_alpha::Real=0.05`: cluster-forming threshold (`:cluster`)
- `ch_adj::Union{Nothing, AbstractMatrix{Bool}}=nothing`: channel adjacency (`:cluster`, see `channel_adjacency()`); `nothing`: channels are not adjacent
- `rng::AbstractRNG=Random.default_rng()`: random number generator (`:cluster`)

# Returns

Named tuple:

- `kind::Symbol`, `factor::Symbol`, `levels::Vector{String}`: compared levels
- `test::Symbol`: `:ttest`, `:ttest_paired`, `:anova` or `:anova_rm`
- `n::Vector{Int64}`: number of subjects of each level (between-subject factors) or of subjects with data of all levels (within-subject factors; the others are left out)
- `stat::Array{Float64}`: t or F statistic at every point, shape (channel, time or frequency) or (channel, frequency, time)
- `p::Array{Float64}`: point-wise p-values (uncorrected)
- `p_adj::Array{Float64}`: corrected p-values: adjusted p-values, or for `:cluster` the p-value of the cluster of the point (`NaN` outside clusters)
- `mask::BitArray`: significant points (`p_adj < alpha`)
- `clusters::DataFrame`: clusters (`:cluster`; see `cluster_perm()`), otherwise empty
- `cluster_id::Array{Int64}`: cluster of every point (`:cluster`; 0: none)
- `null::Vector{Float64}`: null distribution of the maximal cluster mass (`:cluster`, see `cluster_perm()`), otherwise empty
- `threshold::Float64`: cluster-forming threshold of the statistic (`:cluster`), otherwise `NaN`
- `labels::Vector{String}`, `t::Vector{Float64}`, `f::Vector{Float64}`: channels, time points, frequencies
- `correction::Symbol`, `db::Bool`
"""
function study_stats(
    study::NeuroAnalyzer.STUDY,
    kind::Symbol;
    factor::Symbol = :group,
    levels::Union{Nothing, Vector{String}} = nothing,
    condition::Union{Nothing, String} = nothing,
    group::Union{Nothing, String} = nothing,
    correction::Symbol = :cluster,
    alpha::Real = 0.05,
    db::Bool = true,
    tail::Symbol = :both,
    nperm::Int64 = 1000,
    cluster_alpha::Real = 0.05,
    ch_adj::Union{Nothing, AbstractMatrix{Bool}} = nothing,
    rng::AbstractRNG = Random.default_rng(),
)::NamedTuple
    # validate
    _check_var(kind, [:erp, :psd, :tf], "kind")
    _check_var(correction, [:cluster, :fdr, :by, :holm, :bonferroni, :none], "correction")
    _check_var(tail, [:both, :right, :left], "tail")
    0 < alpha < 1 || throw(ArgumentError("alpha must be in (0, 1)."))
    m = _study_measure(study, kind)
    d, idx = _study_measure_rows(study, m)
    X = collect(selectdim(m.data, 1, idx))
    db = db && kind !== :erp
    db && (X = map(x -> x > 0 ? 10 * log10(x) : NaN, X))
    rest = size(X)[3:end]
    within = factor in (:condition, :session)
    factor === :subject && throw(ArgumentError("factor must not be :subject."))
    within || factor in propertynames(d) ||
        throw(ArgumentError("Design table has no column $factor."))
    if !within && !isnothing(group)
        throw(ArgumentError("group selects subjects of within-subject factors; use levels for $factor."))
    end

    # one condition: between-subject factors and sessions
    function _one_condition()
        if isnothing(condition)
            length(m.conditions) == 1 || throw(
                ArgumentError(
                    "condition must be given (conditions: $(join(repr.(m.conditions), ", "))).",
                ),
            )
            return 1
        end
        return only(_study_conditions(m, condition))
    end

    if factor === :condition
        !isnothing(condition) &&
            throw(ArgumentError("condition is not used with factor = :condition (see levels)."))
        lv = isnothing(levels) ? m.conditions : levels
        ci = _study_conditions(m, lv)
        sel = isnothing(group) ? trues(DataFrames.nrow(d)) : d.group .== group
        any(sel) || throw(ArgumentError("No recordings of group $(repr(group))."))
        # subjects × conditions × points
        Y = cat([reshape(_study_subject_avg(selectdim(X, 2, c), d, sel, false),
                         :, 1, rest...) for c in ci]...; dims = 2)
        n = [size(Y, 1)]
    elseif factor === :session
        c = _one_condition()
        sel = isnothing(group) ? trues(DataFrames.nrow(d)) : d.group .== group
        any(sel) || throw(ArgumentError("No recordings of group $(repr(group))."))
        lv = isnothing(levels) ? unique(d.session[sel]) : levels
        subj = unique(d.subject[sel])
        Y = fill(NaN, length(subj), length(lv), rest...)
        for (k, s) in enumerate(subj), (l, ses) in enumerate(lv)
            r = findfirst(i -> sel[i] && d.subject[i] == s && d.session[i] == ses,
                          1:DataFrames.nrow(d))
            isnothing(r) && continue
            selectdim(selectdim(Y, 1, k), 1, l) .= selectdim(selectdim(X, 1, r), 1, c)
        end
        n = [length(subj)]
    else
        c = _one_condition()
        col = string.(d[!, factor])
        lv = isnothing(levels) ? unique(col) : levels
        for l in lv
            l in col || throw(ArgumentError("Level $(repr(l)) not found in $factor."))
        end
        # a subject must belong to one level
        for s in unique(d.subject)
            length(unique(col[d.subject .== s])) == 1 || throw(
                ArgumentError("Subject $s has several levels of $factor (a within-subject factor?)."),
            )
        end
        G = [_study_subject_avg(selectdim(X, 2, c), d, col .== l, false) for l in lv]
        n = [size(g, 1) for g in G]
    end
    length(lv) >= 2 || throw(ArgumentError("At least 2 levels are needed (levels: $(join(repr.(lv), ", ")))."))
    if within
        # subjects with data of all levels
        ok = [all(l -> !all(isnan, selectdim(selectdim(Y, 1, k), 1, l)), axes(Y, 2))
              for k in axes(Y, 1)]
        nd = count(.!ok)
        nd > 0 && _info("$nd subject(s) without data of all levels left out.")
        Y = Y[ok, ntuple(_ -> Colon(), ndims(Y) - 1)...]
        n = [size(Y, 1)]
        n[1] >= 2 || throw(
            ArgumentError("Fewer than 2 subjects have data of all levels of $factor."),
        )
    else
        all(>=(1), n) && sum(n) > length(n) ||
            throw(ArgumentError("Not enough subjects (subjects per level: $(join(n, ", ")))."))
    end
    test = within ? (length(lv) == 2 ? :ttest_paired : :anova_rm) :
           (length(lv) == 2 ? :ttest : :anova)
    test in (:anova, :anova_rm) && tail !== :both &&
        throw(ArgumentError("F-tests are right-tailed: tail must be :both."))

    if correction === :cluster
        r = if test === :ttest
            cluster_perm(G[1], G[2]; tail = tail, nperm = nperm, cluster_alpha = cluster_alpha,
                         alpha = alpha, ch_adj = ch_adj, rng = rng)
        elseif test === :ttest_paired
            cluster_perm(selectdim(Y, 2, 1), selectdim(Y, 2, 2); paired = true, tail = tail,
                         nperm = nperm, cluster_alpha = cluster_alpha, alpha = alpha,
                         ch_adj = ch_adj, rng = rng)
        elseif test === :anova
            cluster_perm(G...; test = :anova, nperm = nperm, cluster_alpha = cluster_alpha,
                         alpha = alpha, ch_adj = ch_adj, rng = rng)
        else
            cluster_perm(Y; test = :anova_rm, nperm = nperm, cluster_alpha = cluster_alpha,
                         alpha = alpha, ch_adj = ch_adj, rng = rng)
        end
        stat, p, p_adj, mask = r.stat, r.p, r.p_cluster, r.mask
        clusters, cluster_id = r.clusters, r.cluster_id
        null, threshold = r.null, r.threshold
    else
        r = if test === :ttest
            mass_ttest(G[1], G[2]; tail = tail)
        elseif test === :ttest_paired
            mass_ttest(selectdim(Y, 2, 1), selectdim(Y, 2, 2); paired = true, tail = tail)
        elseif test === :anova
            mass_anova(G...)
        else
            mass_anova_rm(Y)
        end
        stat = test in (:ttest, :ttest_paired) ? r.t : r.f
        p = r.p
        p_adj = correction === :none ? copy(p) : p_adjust(p; method = correction)
        mask = .!isnan.(p_adj) .& (p_adj .< alpha)
        clusters = DataFrame(:id => Int64[], :sign => Int64[], :mass => Float64[],
                             :size => Int64[], :p => Float64[])
        cluster_id = zeros(Int64, size(stat))
        null, threshold = Float64[], NaN
    end

    return (
        kind = kind,
        factor = factor,
        levels = string.(lv),
        test = test,
        n = n,
        stat = stat,
        p = p,
        p_adj = p_adj,
        mask = mask,
        clusters = clusters,
        cluster_id = cluster_id,
        null = null,
        threshold = threshold,
        labels = m.labels,
        t = m.t,
        f = m.f,
        correction = correction,
        db = db,
    )
end

"""
    cluster_summary(st)

Table of the clusters of a cluster-based permutation test of a study measure (`study_stats()` with `correction = :cluster`): their mass, size, p-value, extent in time and/or frequency, channels and peak.

# Arguments

- `st::NamedTuple`: result of `study_stats()` with `correction = :cluster`

# Returns

- `DataFrame`: one row per cluster (in the order of `st.clusters`, sorted by p-value):
    - `id`, `sign` (+1, -1), `mass`, `size` (number of points), `p`: as in `st.clusters`
    - `significant`: the cluster is significant (its points are in `st.mask`)
    - `t_start`, `t_end`: extent in seconds (`:erp`, `:tf`)
    - `f_start`, `f_end`: extent in Hz (`:psd`, `:tf`)
    - `n_channels`, `channels`: channels of the cluster (labels joined by ", ")
    - `peak`: statistic of largest magnitude in the cluster, `peak_channel` and `peak_t` (`:erp`, `:tf`) and/or `peak_f` (`:psd`, `:tf`): its point
"""
function cluster_summary(st::NamedTuple)::DataFrame
    all(k -> haskey(st, k), (:kind, :clusters, :cluster_id, :stat, :mask, :labels, :correction)) ||
        throw(ArgumentError("st must be a result of study_stats()."))
    st.correction === :cluster ||
        throw(ArgumentError("st must be a result of study_stats() with correction = :cluster."))
    has_t = st.kind in (:erp, :tf)
    has_f = st.kind in (:psd, :tf)
    df = DataFrame(
        :id => Int64[], :sign => Int64[], :mass => Float64[], :size => Int64[],
        :p => Float64[], :significant => Bool[],
    )
    has_t && (df.t_start = Float64[]; df.t_end = Float64[])
    has_f && (df.f_start = Float64[]; df.f_end = Float64[])
    df.n_channels = Int64[]
    df.channels = String[]
    df.peak = Float64[]
    df.peak_channel = String[]
    has_t && (df.peak_t = Float64[])
    has_f && (df.peak_f = Float64[])
    for c in eachrow(st.clusters)
        pts = findall(==(c.id), st.cluster_id)
        isempty(pts) && continue
        chs = sort(unique(p[1] for p in pts))
        pk = pts[argmax([abs(st.stat[p]) for p in pts])]
        r = Dict{Symbol, Any}(
            :id => c.id, :sign => c.sign, :mass => c.mass, :size => c.size, :p => c.p,
            :significant => any(st.mask[p] for p in pts),
            :n_channels => length(chs), :channels => join(st.labels[chs], ", "),
            :peak => st.stat[pk], :peak_channel => st.labels[pk[1]],
        )
        # time: last dimension; frequency: second dimension
        if has_t
            ti = [p[ndims(st.stat)] for p in pts]
            r[:t_start], r[:t_end] = st.t[minimum(ti)], st.t[maximum(ti)]
            r[:peak_t] = st.t[pk[ndims(st.stat)]]
        end
        if has_f
            fi = [p[2] for p in pts]
            r[:f_start], r[:f_end] = st.f[minimum(fi)], st.f[maximum(fi)]
            r[:peak_f] = st.f[pk[2]]
        end
        push!(df, r)
    end
    return df
end
