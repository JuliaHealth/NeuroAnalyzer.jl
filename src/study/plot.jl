export plot_study_erp
export plot_study_stats

"""
    plot_study_erp(study; <keyword arguments>)

Plot group ERPs (ERFs for MEG) of a study: the grand average of each level of a factor (e.g. each group or condition) with its standard error or 95% confidence interval across subjects. Optionally the ERPs of the subjects are drawn behind the grand averages (butterfly), the difference wave of two levels is added and significant time ranges (from `study_stats()`) are drawn as bars below the waveforms. Requires the `:erp` measure (see `study_measure!()`).

The recordings of one subject are averaged first (unless `factor = :session`), so each subject counts once. With several channels, the ERP of each subject is the mean of the channels (region of interest); bad channels of a subject (`NaN`) are left out of its mean.

# Arguments

- `study::NeuroAnalyzer.STUDY`
- `ch::Union{String, Vector{String}}`: channel, or channels averaged as a region of interest
- `factor::Symbol=:condition`: factor whose levels are plotted: `:condition` or `:session` (within subjects), or a design column, e.g. `:group` (between subjects)
- `levels::Union{Nothing, Vector{String}}=nothing`: levels to plot, in this order; `nothing`: all levels
- `condition::Union{Nothing, String}=nothing`: condition (between-subject factors and `:session`); may be omitted if the measure has one condition
- `group::Union{Nothing, String}=nothing`: subjects of one group (within-subject factors); `nothing`: all subjects
- `band::Symbol=:ci95`: band around the grand averages: `:ci95` (95% confidence interval across subjects, t distribution), `:se` (± standard error) or `:none`
- `subjects::Bool=false`: if `true`, draw the ERPs of the subjects behind the grand averages (thin lines in the colors of their levels)
- `diff::Bool=false`: if `true`, add the difference wave of the two plotted levels (first − second, dashed black line, with its band): for within-subject factors the mean of the differences of the subjects with data of both levels, for between-subject factors the difference of the grand averages (band: Welch standard error)
- `stats::Union{Nothing, NamedTuple}=nothing`: result of `study_stats()` for the `:erp` measure; its significant time points (`mask`) at the plotted channel(s) (with several channels: at any of them) are drawn as bars below the waveforms
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `title::String="default"`: plot title
- `yrev::Bool=false`: if `true`, reverse the y-axis (negative up)
- `zl::Bool=true`: if `true`, draw a vertical line at t = 0
- `leg::Bool=true`: if `true`, draw a legend (outside the plot, on the right)
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_study_erp(
    study::NeuroAnalyzer.STUDY;
    ch::Union{String, Vector{String}},
    factor::Symbol = :condition,
    levels::Union{Nothing, Vector{String}} = nothing,
    condition::Union{Nothing, String} = nothing,
    group::Union{Nothing, String} = nothing,
    band::Symbol = :ci95,
    subjects::Bool = false,
    diff::Bool = false,
    stats::Union{Nothing, NamedTuple} = nothing,
    xlabel::String = "default",
    ylabel::String = "default",
    title::String = "default",
    yrev::Bool = false,
    zl::Bool = true,
    leg::Bool = true,
    mono::Bool = false,
)::GLMakie.Figure
    # validate
    _check_var(band, [:ci95, :se, :none], "band")
    m = _study_measure(study, :erp)
    d, idx = _study_measure_rows(study, m)
    chs = ch isa String ? [ch] : ch
    isempty(chs) && throw(ArgumentError("No channels selected."))
    chi = indexin(chs, m.labels)
    any(isnothing, chi) && throw(
        ArgumentError(
            "Channels not in measure :erp: $(join(chs[isnothing.(chi)], ", ")).",
        ),
    )
    chi = Int64.(chi)
    length(unique(m.unit[chi])) == 1 ||
        throw(ArgumentError("All channels must have the same unit."))

    # levels of the factor: subject ERPs (subject × time) and subject IDs of each level
    lv, S, subj = _study_level_erps(m, d, idx, chi; factor, levels, condition, group)
    length(lv) >= 1 || throw(ArgumentError("No levels to plot."))
    diff && length(lv) != 2 &&
        throw(ArgumentError("diff requires 2 levels (levels: $(join(lv, ", ")))."))

    t = m.t
    unit = m.unit[chi[1]]
    within = factor in (:condition, :session)

    # grand averages and their bands
    ga = [_ga_band(s, band) for s in S]
    if diff
        if within
            # subjects with data of both levels: the mean of their differences
            common = intersect(subj[1], subj[2])
            length(common) >= 2 ||
                throw(ArgumentError("diff requires ≥ 2 subjects with data of both levels."))
            D = S[1][indexin(common, subj[1]), :] .- S[2][indexin(common, subj[2]), :]
            gd = _ga_band(D, band)
        else
            gd = _ga_band_welch(S[1], S[2], band)
        end
    end

    # significant time points at the plotted channel(s)
    sig = nothing
    if !isnothing(stats)
        (haskey(stats, :kind) && stats.kind === :erp && haskey(stats, :mask)) ||
            throw(ArgumentError("stats must be a result of study_stats() for :erp."))
        length(stats.t) == length(t) ||
            throw(ArgumentError("stats and the :erp measure have different time points."))
        si = indexin(chs, stats.labels)
        any(isnothing, si) && throw(
            ArgumentError("Channels not in stats: $(join(chs[isnothing.(si)], ", "))."),
        )
        sig = vec(any(stats.mask[Int64.(si), :]; dims = 1))
    end

    # prepare plot
    GLMakie.activate!(; title = "plot_study_erp()")
    fig = GLMakie.Figure(; size = (900, 450))
    fl = factor === :condition ? "condition" : string(factor)
    ch_str = length(chs) == 1 ? chs[1] : "mean of $(length(chs)) channels"
    xl, yl, tt = _set_defaults(
        xlabel, ylabel, title,
        "Time [s]", "Amplitude [$unit]",
        "Grand average ERP, $ch_str, by $fl" *
        (band === :none ? "" : band === :ci95 ? " (95% CI)" : " (± SE)"),
    )
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel = xl,
        ylabel = yl,
        title = tt,
        _TIME_AXIS_KWARGS...,
        yticks = LinearTicks(10),
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(10),
        yreversed = yrev,
        xautolimitmargin = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )
    _style_axis!(ax)
    zl && GLMakie.vlines!(ax, 0; color = :gray, linestyle = :dash, linewidth = 2)
    GLMakie.hlines!(ax, 0; color = :gray, linewidth = 0.5)

    # categorical colors of the levels
    colors = mono ? _trace_colors(length(lv); mono = true) :
             GLMakie.Makie.wong_colors()[mod1.(1:length(lv), 7)]

    # subjects behind the grand averages
    if subjects
        for (l, s) in enumerate(S), k in axes(s, 1)
            all(isnan, @view(s[k, :])) && continue
            GLMakie.lines!(ax, t, s[k, :]; color = (colors[l], 0.25), linewidth = 0.75)
        end
    end

    # bands, then the grand averages
    for (l, g) in enumerate(ga)
        band === :none ||
            GLMakie.band!(ax, t, g.ll, g.ul; color = (colors[l], 0.2))
    end
    diff && band !== :none && GLMakie.band!(ax, t, gd.ll, gd.ul; color = (:black, 0.12))
    for (l, g) in enumerate(ga)
        GLMakie.lines!(
            ax, t, g.m;
            color = colors[l], linewidth = 2, label = "$(lv[l]) (n = $(g.n))",
        )
    end
    diff && GLMakie.lines!(
        ax, t, gd.m;
        color = :black, linewidth = 2, linestyle = :dash,
        label = "$(lv[1]) − $(lv[2])",
    )

    # y range of everything drawn; room for the significance bars below the waveforms
    ys = Float64[]
    for g in ga
        append!(ys, g.m)
        band === :none || append!(ys, g.ll, g.ul)
    end
    if diff
        append!(ys, gd.m)
        band === :none || append!(ys, gd.ll, gd.ul)
    end
    subjects && foreach(s -> append!(ys, vec(s)), S)
    ys = Base.filter(isfinite, ys)
    y1, y2 = isempty(ys) ? (-1.0, 1.0) : extrema(ys)
    y1 == y2 && (y1 -= 1; y2 += 1)
    pad = 0.05 * (y2 - y1)
    y_lo = y1 - pad
    if !isnothing(sig)
        # bars at the bottom of the plot (the top if the y-axis is reversed: below the
        # waveforms on the screen)
        y_bar = y1 - 0.08 * (y2 - y1)
        y_lo = y1 - 0.13 * (y2 - y1)
        dt = length(t) > 1 ? (t[2] - t[1]) / 2 : 0.0
        runs = _true_runs(sig)
        for (r, (a, b)) in enumerate(runs)
            GLMakie.lines!(
                ax, [t[a] - dt, t[b] + dt], [y_bar, y_bar];
                color = mono ? :black : :red, linewidth = 6, linecap = :butt,
                label = r == 1 ? "significant" : nothing,
            )
        end
        isempty(runs) && _info("No significant time points at the plotted channel(s).")
    end
    GLMakie.ylims!(ax, yrev ? (y2 + pad, y_lo) : (y_lo, y2 + pad))

    leg && _side_legend!(fig, ax)

    return fig
end

# subject ERPs of the levels of `factor`: level names, (subject × time) matrices (mean of the
# channels `chi`, bad channels left out) and the subject IDs of their rows
function _study_level_erps(
    m, d::DataFrame, idx::Vector{Int64}, chi::Vector{Int64};
    factor::Symbol,
    levels::Union{Nothing, Vector{String}},
    condition::Union{Nothing, String},
    group::Union{Nothing, String},
)
    X = selectdim(m.data, 1, idx)
    within = factor in (:condition, :session)
    factor === :subject && throw(ArgumentError("factor must not be :subject."))
    within || factor in propertynames(d) ||
        throw(ArgumentError("Design table has no column $factor."))
    !within && !isnothing(group) && throw(
        ArgumentError("group selects subjects of within-subject factors; use levels for $factor."),
    )
    gsel = isnothing(group) ? trues(DataFrames.nrow(d)) : d.group .== group
    any(gsel) || throw(ArgumentError("No recordings of group $(repr(group))."))

    # condition of between-subject factors and sessions
    if factor !== :condition
        if isnothing(condition)
            length(m.conditions) == 1 || throw(
                ArgumentError(
                    "condition must be given (conditions: $(join(repr.(m.conditions), ", "))).",
                ),
            )
            condition = m.conditions[1]
        end
        c = only(_study_conditions(m, condition))
    end

    all_lv = factor === :condition ? m.conditions : unique(string.(d[gsel, factor]))
    lv = isnothing(levels) ? all_lv : levels
    miss = setdiff(lv, all_lv)
    isempty(miss) || throw(
        ArgumentError(
            "Levels not found: $(join(repr.(miss), ", ")) (levels: $(join(repr.(all_lv), ", "))).",
        ),
    )

    S = Matrix{Float64}[]
    subj = Vector{String}[]
    for l in lv
        if factor === :condition
            sel = gsel
            Xc = selectdim(X, 2, only(_study_conditions(m, l)))
        else
            sel = gsel .& (string.(d[!, factor]) .== l)
            Xc = selectdim(X, 2, c)
        end
        # (subject, channel, time); with factor = :session each recording is one row
        A = _study_subject_avg(Xc, d, sel, factor === :session)
        rows = findall(sel)
        ids = factor === :session ? string.(d.subject[rows]) : unique(string.(d.subject[rows]))
        # mean of the channels (bad channels: NaN, left out)
        R = fill(NaN, size(A, 1), size(A, 3))
        for k in axes(A, 1), j in axes(A, 3)
            v = Base.filter(!isnan, @view A[k, chi, j])
            isempty(v) || (R[k, j] = mean(v))
        end
        # subjects without data of this level (e.g. no epochs of the condition)
        has = [!all(isnan, @view R[k, :]) for k in axes(R, 1)]
        push!(S, R[has, :])
        push!(subj, ids[has])
    end
    return lv, S, subj
end

# grand average over the rows of `s` (subjects) with a band: 95% CI (t distribution) or SE
function _ga_band(s::AbstractMatrix, band::Symbol)
    mu, sd, n = _nanstats(s)
    se = sd ./ sqrt.(max.(n, 1))
    h = band === :ci95 ? se .* [nn > 1 ? quantile(TDist(nn - 1), 0.975) : NaN for nn in n] :
        band === :se ? se : zero(se)
    return (m = mu, ll = mu .- h, ul = mu .+ h, n = isempty(n) ? 0 : maximum(n))
end

# difference of the grand averages of independent groups `s1` − `s2`; band: Welch standard
# error, 95% CI with the Welch-Satterthwaite degrees of freedom
function _ga_band_welch(s1::AbstractMatrix, s2::AbstractMatrix, band::Symbol)
    m1, sd1, n1 = _nanstats(s1)
    m2, sd2, n2 = _nanstats(s2)
    (minimum(n1) >= 2 && minimum(n2) >= 2) ||
        throw(ArgumentError("diff requires ≥ 2 subjects in each level."))
    v1 = sd1 .^ 2 ./ n1
    v2 = sd2 .^ 2 ./ n2
    se = sqrt.(v1 .+ v2)
    mu = m1 .- m2
    h = if band === :ci95
        df = (v1 .+ v2) .^ 2 ./ (v1 .^ 2 ./ (n1 .- 1) .+ v2 .^ 2 ./ (n2 .- 1))
        se .* [isfinite(f) && f > 0 ? quantile(TDist(f), 0.975) : NaN for f in df]
    elseif band === :se
        se
    else
        zero(se)
    end
    return (m = mu, ll = mu .- h, ul = mu .+ h, n = 0)
end

# (first, last) indices of the runs of `true` values
function _true_runs(x::AbstractVector{Bool})::Vector{Tuple{Int64, Int64}}
    runs = Tuple{Int64, Int64}[]
    a = 0
    for i in eachindex(x)
        if x[i] && a == 0
            a = i
        elseif !x[i] && a != 0
            push!(runs, (a, i - 1))
            a = 0
        end
    end
    a != 0 && push!(runs, (a, lastindex(x)))
    return runs
end

"""
    plot_study_stats(st; <keyword arguments>)

Plot the result of `study_stats()` as an image: the statistic (t or F) or the p-values at every channel × time point (`:erp`) or channel × frequency (`:psd`), or for `:tf` at every frequency × time point of one channel. Significant points (`mask`) are outlined and/or the non-significant points are dimmed; channels can be ordered by region.

# Arguments

- `st::NamedTuple`: result of `study_stats()`
- `value::Symbol=:stat`: plotted value:
    - `:stat`: t (diverging color scale, symmetric around 0) or F statistic
    - `:p`: point-wise p-values, as -log₁₀(p)
    - `:p_adj`: corrected p-values, as -log₁₀(p) (for `correction = :cluster` the p-values of the clusters; points outside clusters are empty)
- `sig::Symbol=:outline`: significant points: `:outline` (outline of the significant regions), `:dim` (non-significant points dimmed), `:both` or `:none`
- `order::Symbol=:region`: order of the channels: `:region` (standard 10-20/10-10 labels grouped into frontal, central, temporal, parietal and occipital regions, from left to right within a region; other channels follow in their order) or `:none` (the order of `st.labels`)
- `ch::Union{Nothing, String}=nothing`: channel of a `:tf` result (frequency × time image)
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `title::String="default"`: plot title
- `cb_title::String="default"`: color bar title
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_study_stats(
    st::NamedTuple;
    value::Symbol = :stat,
    sig::Symbol = :outline,
    order::Symbol = :region,
    ch::Union{Nothing, String} = nothing,
    xlabel::String = "default",
    ylabel::String = "default",
    title::String = "default",
    cb_title::String = "default",
    mono::Bool = false,
)::GLMakie.Figure
    # validate
    all(k -> haskey(st, k), (:kind, :test, :stat, :p, :p_adj, :mask, :labels, :levels)) ||
        throw(ArgumentError("st must be a result of study_stats()."))
    _check_var(value, [:stat, :p, :p_adj], "value")
    _check_var(sig, [:outline, :dim, :both, :none], "sig")
    _check_var(order, [:region, :none], "order")
    tf = st.kind === :tf

    V = value === :stat ? st.stat : value === :p ? st.p : st.p_adj
    M = st.mask
    if tf
        isnothing(ch) && throw(ArgumentError("ch must be given for a :tf result."))
        ci = findfirst(==(ch), st.labels)
        isnothing(ci) && throw(ArgumentError("Channel $ch not in st."))
        # (frequency, time) → x: time, y: frequency
        V = V[ci, :, :]
        M = M[ci, :, :]
        x, y = st.t, st.f
        img = permutedims(V)
        msk = permutedims(M)
        ylab = String[]
    else
        isnothing(ch) || throw(ArgumentError("ch is used only for a :tf result."))
        x = st.kind === :erp ? st.t : st.f
        # channel order (bottom to top: the first channel at the top)
        ord = order === :region ? _region_order(st.labels) : collect(eachindex(st.labels))
        ord = reverse(ord)
        y = 1:length(ord)
        img = permutedims(V[ord, :])
        msk = permutedims(M[ord, :])
        ylab = st.labels[ord]
    end
    value === :stat || (img = map(p -> isnan(p) ? NaN : -log10(max(p, eps())), img))

    # color scale
    is_t = st.test in (:ttest, :ttest_paired)
    fin = Base.filter(isfinite, img)
    vmax = isempty(fin) ? 1.0 : maximum(abs, fin)
    vmax > 0 || (vmax = 1.0)
    if value === :stat && is_t
        crange = (-vmax, vmax)
        pal = _palette(:bluesreds; mono = mono)
    else
        crange = (0.0, vmax)
        pal = _palette(:viridis; mono = mono)
    end

    stat_name = is_t ? "t" : "F"
    vname = value === :stat ? stat_name : value === :p ? "-log₁₀(p)" : "-log₁₀(p adj.)"
    lv = length(st.levels) == 2 ? "$(st.levels[1]) vs $(st.levels[2])" : join(st.levels, ", ")
    xl, yl, tt = _set_defaults(
        xlabel, ylabel, title,
        st.kind === :psd ? "Frequency [Hz]" : "Time [s]",
        tf ? "Frequency [Hz]" : "Channel",
        "$lv: $vname" * (tf ? ", $ch" : "") * " ($(st.test), correction: $(st.correction))",
    )
    cb_title == "default" && (cb_title = vname)

    # prepare plot
    GLMakie.activate!(; title = "plot_study_stats()")
    h = tf ? 450 : clamp(120 + 16 * length(y), 400, 1000)
    fig = GLMakie.Figure(; size = (900, h))
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel = xl,
        ylabel = yl,
        title = tt,
        yticks = tf ? GLMakie.Makie.automatic : (collect(y), ylab),
        yticklabelsize = length(y) > 64 ? 8 : 11,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )
    _style_axis!(ax)
    hm = GLMakie.heatmap!(
        ax, x, y, img;
        colormap = pal, colorrange = crange, nan_color = :lightgray,
    )

    # non-significant points dimmed: a white, semi-transparent layer
    if sig in (:dim, :both)
        dim = map(m -> m ? NaN : 1.0, msk)
        GLMakie.heatmap!(
            ax, x, y, dim;
            colormap = [GLMakie.RGBAf(1, 1, 1, 0.65), GLMakie.RGBAf(1, 1, 1, 0.65)],
            colorrange = (0, 1), nan_color = :transparent,
        )
    end
    # outline of the significant regions
    if sig in (:outline, :both) && any(msk)
        GLMakie.contour!(
            ax, x, y, Float64.(msk);
            levels = [0.5], color = :black, linewidth = 1.5,
        )
    end

    # regions: separator lines and names on the right
    if !tf && order === :region
        reg = [_channel_region(l)[1] for l in ylab]
        bounds = [i + 0.5 for i = 1:(length(reg) - 1) if reg[i] != reg[i + 1]]
        isempty(bounds) || GLMakie.hlines!(ax, bounds; color = :black, linewidth = 1)
        if length(unique(reg)) > 1
            names = unique(reg)
            mids = [mean(findall(==(r), reg)) for r in names]
            ax_r = GLMakie.Axis(
                fig[1, 1];
                yaxisposition = :right,
                yticks = (mids, names),
                yticklabelsize = 11,
                ygridvisible = false,
                xgridvisible = false,
                _AXIS_LOCK_KWARGS...,
            )
            GLMakie.hidexdecorations!(ax_r)
            GLMakie.hidespines!(ax_r)
            GLMakie.linkyaxes!(ax, ax_r)
            GLMakie.linkxaxes!(ax, ax_r)
        end
    end

    # tight limits (half a point around the first and last points)
    dx = length(x) > 1 ? (x[2] - x[1]) / 2 : 0.5
    dy = length(y) > 1 ? (y[2] - y[1]) / 2 : 0.5
    GLMakie.limits!(ax, x[1] - dx, x[end] + dx, y[1] - dy, y[end] + dy)

    GLMakie.Colorbar(fig[1, 2], hm; label = cb_title, ticks = _cb_tickfun(; n = 5))

    return fig
end

# region of a channel label (standard 10-20/10-10 names, a suffix such as "-avg" is
# ignored): (region name, rank of the region, rank of the label prefix, left-right position);
# other labels: ("other", 6, 0, 0)
function _channel_region(label::AbstractString)::Tuple{String, Int64, Int64, Int64}
    l = split(label, r"[-_ ]")[1]
    mt = match(r"^([A-Za-z]+?)(z|Z|\d+)$", l)
    isnothing(mt) && return ("other", 6, 0, 0)
    pre = uppercase(mt.captures[1])
    num = mt.captures[2]
    prefixes = ["FP", "AF", "F", "FC", "C", "FT", "T", "TP", "CP", "P", "PO", "O", "I"]
    regions = Dict(
        "FP" => 1, "AF" => 1, "F" => 1, "FC" => 2, "C" => 2, "FT" => 3, "T" => 3,
        "TP" => 3, "CP" => 4, "P" => 4, "PO" => 5, "O" => 5, "I" => 5,
    )
    names = ["frontal", "central", "temporal", "parietal", "occipital"]
    p_idx = findfirst(==(pre), prefixes)
    isnothing(p_idx) && return ("other", 6, 0, 0)
    # left (odd numbers, the most lateral first), midline (z), right (even numbers)
    pos = lowercase(num) == "z" ? 0 : (n = parse(Int64, num); isodd(n) ? -n : n)
    r = regions[pre]
    return (names[r], r, p_idx, pos)
end

# channel order by region (frontal → occipital), prefix and left-right position; other
# channels follow in their order
function _region_order(labels::AbstractVector{<:AbstractString})::Vector{Int64}
    key = [(_channel_region(l)[2:4]..., i) for (i, l) in enumerate(labels)]
    return sortperm(key)
end
