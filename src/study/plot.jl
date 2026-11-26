export plot_study_erp
export plot_study_psd
export plot_study_stats
export plot_study_topo
export plot_study_tf
export plot_study_clusters
export plot_study_values
export plot_study_covariate
export plot_study_overview

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
- `band::Symbol=:ci95`: band around the grand averages: `:ci95` (95% confidence interval across subjects, t distribution), `:se` (± standard error), `:sd` (± standard deviation, the spread of the subjects) or `:none`
- `subjects::Bool=false`: if `true`, draw the ERPs of the subjects behind the grand averages (thin lines in the colors of their levels)
- `diff::Bool=false`: if `true`, add the difference wave of the two plotted levels (first − second, dashed black line, with its band): for within-subject factors the mean of the differences of the subjects with data of both levels, for between-subject factors the difference of the grand averages (band: Welch standard error; no band for `band = :sd`)
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
    _check_var(band, [:ci95, :se, :sd, :none], "band")
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
        (band === :none ? "" : band === :ci95 ? " (95% CI)" : band === :sd ? " (± SD)" : " (± SE)"),
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

"""
    plot_study_psd(study; <keyword arguments>)

Plot group power spectra of a study: the mean spectrum of each level of a factor (e.g. each group or condition) with its standard error or 95% confidence interval across subjects. Optionally the spectra of the subjects are drawn behind the means, the difference of two levels is added and significant frequencies (from `study_stats()`) are shaded. Requires the `:psd` measure (see `study_measure!()`); band-power topographical maps: `plot_study_topo(study, :psd)`.

The recordings of one subject are averaged first (unless `factor = :session`), so each subject counts once. The power of each subject is converted to dB before averaging (`db = true`). With several channels, the spectrum of each subject is the mean of the channels (region of interest); bad channels of a subject (`NaN`) are left out of its mean.

# Arguments

- `study::NeuroAnalyzer.STUDY`
- `ch::Union{String, Vector{String}}`: channel, or channels averaged as a region of interest
- `factor::Symbol=:condition`: factor whose levels are plotted: `:condition` or `:session` (within subjects), or a design column, e.g. `:group` (between subjects)
- `levels::Union{Nothing, Vector{String}}=nothing`: levels to plot, in this order; `nothing`: all levels
- `condition::Union{Nothing, String}=nothing`: condition (between-subject factors and `:session`); may be omitted if the measure has one condition
- `group::Union{Nothing, String}=nothing`: subjects of one group (within-subject factors); `nothing`: all subjects
- `db::Bool=true`: average the power of the subjects in dB
- `band::Symbol=:se`: band around the means: `:se` (± standard error), `:ci95` (95% confidence interval across subjects, t distribution), `:sd` (± standard deviation, the spread of the subjects) or `:none`
- `subjects::Bool=false`: if `true`, draw the spectra of the subjects behind the means (thin lines in the colors of their levels)
- `diff::Bool=false`: if `true`, add the difference of the two plotted levels (first − second, dashed black line, with its band): for within-subject factors the mean of the differences of the subjects with data of both levels, for between-subject factors the difference of the means (band: Welch standard error; no band for `band = :sd`)
- `stats::Union{Nothing, NamedTuple}=nothing`: result of `study_stats()` for the `:psd` measure; its significant frequencies (`mask`) at the plotted channel(s) (with several channels: at any of them) are shaded
- `flim::Union{Nothing, Tuple{Real, Real}}=nothing`: frequency limits; `nothing`: all frequencies
- `frq::Symbol=:lin`: frequency scaling: `:lin` (linear) or `:log` (logarithmic; 0 Hz is left out)
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `title::String="default"`: plot title
- `leg::Bool=true`: if `true`, draw a legend (outside the plot, on the right)
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_study_psd(
    study::NeuroAnalyzer.STUDY;
    ch::Union{String, Vector{String}},
    factor::Symbol = :condition,
    levels::Union{Nothing, Vector{String}} = nothing,
    condition::Union{Nothing, String} = nothing,
    group::Union{Nothing, String} = nothing,
    db::Bool = true,
    band::Symbol = :se,
    subjects::Bool = false,
    diff::Bool = false,
    stats::Union{Nothing, NamedTuple} = nothing,
    flim::Union{Nothing, Tuple{Real, Real}} = nothing,
    frq::Symbol = :lin,
    xlabel::String = "default",
    ylabel::String = "default",
    title::String = "default",
    leg::Bool = true,
    mono::Bool = false,
)::GLMakie.Figure
    # validate
    _check_var(band, [:ci95, :se, :sd, :none], "band")
    _check_var(frq, [:lin, :log], "frq")
    m = _study_measure(study, :psd)
    d, idx = _study_measure_rows(study, m)
    chs = ch isa String ? [ch] : ch
    isempty(chs) && throw(ArgumentError("No channels selected."))
    chi = indexin(chs, m.labels)
    any(isnothing, chi) && throw(
        ArgumentError("Channels not in measure :psd: $(join(chs[isnothing.(chi)], ", "))."),
    )
    chi = Int64.(chi)
    length(unique(m.unit[chi])) == 1 ||
        throw(ArgumentError("All channels must have the same unit."))

    # frequencies
    f = m.f
    fsel = isnothing(flim) ? trues(length(f)) : (f .>= flim[1]) .& (f .<= flim[2])
    isnothing(flim) || flim[1] < flim[2] || throw(ArgumentError("flim[1] must be < flim[2]."))
    frq === :log && (fsel .&= f .> 0)
    count(fsel) >= 2 || throw(ArgumentError("flim must contain at least 2 frequencies."))
    fi = findall(fsel)
    f = f[fi]

    # levels of the factor: subject spectra (subject × frequency)
    lv, S, subj = _study_level_erps(m, d, idx, chi; factor, levels, condition, group, db)
    S = [s[:, fi] for s in S]
    diff && length(lv) != 2 &&
        throw(ArgumentError("diff requires 2 levels (levels: $(join(lv, ", ")))."))
    within = factor in (:condition, :session)
    ga = [_ga_band(s, band) for s in S]
    if diff
        if within
            common = intersect(subj[1], subj[2])
            length(common) >= 2 ||
                throw(ArgumentError("diff requires ≥ 2 subjects with data of both levels."))
            D = S[1][indexin(common, subj[1]), :] .- S[2][indexin(common, subj[2]), :]
            gd = _ga_band(D, band)
        else
            gd = _ga_band_welch(S[1], S[2], band)
        end
    end

    # significant frequencies at the plotted channel(s)
    sig = nothing
    if !isnothing(stats)
        (haskey(stats, :kind) && stats.kind === :psd && haskey(stats, :mask)) ||
            throw(ArgumentError("stats must be a result of study_stats() for :psd."))
        length(stats.f) == length(m.f) ||
            throw(ArgumentError("stats and the :psd measure have different frequencies."))
        si = indexin(chs, stats.labels)
        any(isnothing, si) && throw(
            ArgumentError("Channels not in stats: $(join(chs[isnothing.(si)], ", "))."),
        )
        sig = vec(any(stats.mask[Int64.(si), :]; dims = 1))[fi]
    end

    # prepare plot
    GLMakie.activate!(; title = "plot_study_psd()")
    fig = GLMakie.Figure(; size = (900, 450))
    unit = m.unit[chi[1]]
    fl = factor === :condition ? "condition" : string(factor)
    ch_str = length(chs) == 1 ? chs[1] : "mean of $(length(chs)) channels"
    xl, yl, tt = _set_defaults(
        xlabel, ylabel, title,
        "Frequency [Hz]", db ? "Power [dB $unit^2/Hz]" : "Power [$unit^2/Hz]",
        "Group PSD, $ch_str, by $fl" *
        (band === :none ? "" : band === :ci95 ? " (95% CI)" : band === :sd ? " (± SD)" : " (± SE)"),
    )
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel = xl,
        ylabel = yl,
        title = tt,
        xscale = frq === :lin ? identity : log10,
        xautolimitmargin = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )
    _style_axis!(ax)
    if frq === :log
        # ticks in Hz: 1, 2, 5, 10, 20, 50, ...
        tk = [m * 10.0^e for e = -2:4 for m in (1, 2, 5)]
        tk = tk[(tk .>= f[1]) .& (tk .<= f[end])]
        ax.xticks = (tk, [isinteger(v) ? string(Int64(v)) : string(v) for v in tk])
    end

    # significant frequencies: shaded (behind the spectra)
    if !isnothing(sig)
        runs = _true_runs(sig)
        for (r, (a, b)) in enumerate(runs)
            # half a frequency step around the first and last significant frequency
            x1 = a > 1 ? (f[a - 1] + f[a]) / 2 : f[a]
            x2 = b < length(f) ? (f[b] + f[b + 1]) / 2 : f[b]
            x1 == x2 && (x2 = x1 + eps(x1))
            GLMakie.vspan!(
                ax, x1, x2;
                color = (mono ? :black : :red, 0.12),
                label = r == 1 ? "significant" : nothing,
            )
        end
        isempty(runs) && _info("No significant frequencies at the plotted channel(s).")
    end

    colors = mono ? _trace_colors(length(lv); mono = true) :
             GLMakie.Makie.wong_colors()[mod1.(1:length(lv), 7)]
    if subjects
        for (l, s) in enumerate(S), k in axes(s, 1)
            all(isnan, @view(s[k, :])) && continue
            GLMakie.lines!(ax, f, s[k, :]; color = (colors[l], 0.25), linewidth = 0.75)
        end
    end
    for (l, g) in enumerate(ga)
        band === :none || GLMakie.band!(ax, f, g.ll, g.ul; color = (colors[l], 0.2))
    end
    diff && band !== :none && GLMakie.band!(ax, f, gd.ll, gd.ul; color = (:black, 0.12))
    for (l, g) in enumerate(ga)
        GLMakie.lines!(
            ax, f, g.m;
            color = colors[l], linewidth = 2, label = "$(lv[l]) (n = $(g.n))",
        )
    end
    if diff
        GLMakie.hlines!(ax, 0; color = :gray, linewidth = 0.5)
        GLMakie.lines!(
            ax, f, gd.m;
            color = :black, linewidth = 2, linestyle = :dash, label = "$(lv[1]) − $(lv[2])",
        )
    end
    GLMakie.xlims!(ax, f[1], f[end])

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
    db::Bool = false,
)
    lv, A, ids = _study_level_data(m, d, idx; factor, levels, condition, group, db)
    S = Matrix{Float64}[]
    subj = Vector{String}[]
    for l in eachindex(lv)
        # mean of the channels (bad channels: NaN, left out)
        R = fill(NaN, size(A[l], 1), size(A[l], 3))
        for k in axes(A[l], 1), j in axes(A[l], 3)
            v = Base.filter(!isnan, @view A[l][k, chi, j])
            isempty(v) || (R[k, j] = mean(v))
        end
        # subjects without data of this level at these channels
        has = [!all(isnan, @view R[k, :]) for k in axes(R, 1)]
        push!(S, R[has, :])
        push!(subj, ids[l][has])
    end
    return lv, S, subj
end

# subject averages of the levels of `factor`: level names, (subject × channel × points...)
# arrays (subjects without data of the level are left out) and the subject IDs of their rows;
# `db`: powers converted to dB before averaging
function _study_level_data(
    m, d::DataFrame, idx::Vector{Int64};
    factor::Symbol,
    levels::Union{Nothing, Vector{String}},
    condition::Union{Nothing, String},
    group::Union{Nothing, String},
    db::Bool = false,
)
    X = selectdim(m.data, 1, idx)
    db && (X = map(x -> x > 0 ? 10 * log10(x) : NaN, X))
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
    elseif !isnothing(condition)
        throw(ArgumentError("condition is not used with factor = :condition (see levels)."))
    end

    all_lv = factor === :condition ? m.conditions : unique(string.(d[gsel, factor]))
    lv = isnothing(levels) ? all_lv : levels
    miss = setdiff(lv, all_lv)
    isempty(miss) || throw(
        ArgumentError(
            "Levels not found: $(join(repr.(miss), ", ")) (levels: $(join(repr.(all_lv), ", "))).",
        ),
    )

    A = Array{Float64}[]
    subj = Vector{String}[]
    for l in lv
        if factor === :condition
            sel = gsel
            Xc = selectdim(X, 2, only(_study_conditions(m, l)))
        else
            sel = gsel .& (string.(d[!, factor]) .== l)
            Xc = selectdim(X, 2, c)
        end
        # (subject, channel, points...); with factor = :session each recording is one row
        S = _study_subject_avg(Xc, d, sel, factor === :session)
        rows = findall(sel)
        ids = factor === :session ? string.(d.subject[rows]) : unique(string.(d.subject[rows]))
        has = [!all(isnan, selectdim(S, 1, k)) for k in axes(S, 1)]
        push!(A, S[has, ntuple(_ -> Colon(), ndims(S) - 1)...])
        push!(subj, ids[has])
    end
    return lv, A, subj
end

# grand average over the rows of `s` (subjects) with a band: 95% CI (t distribution), SE or SD
function _ga_band(s::AbstractMatrix, band::Symbol)
    mu, sd, n = _nanstats(s)
    se = sd ./ sqrt.(max.(n, 1))
    h = band === :ci95 ? se .* [nn > 1 ? quantile(TDist(nn - 1), 0.975) : NaN for nn in n] :
        band === :se ? se : band === :sd ? sd : zero(sd)
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
        # (the difference of the means of independent groups has no subject-level SD)
        band === :sd && _info("band = :sd: the difference of independent groups is drawn without a band.")
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
            ax, _mask_padded(x, y, msk)...;
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

"""
    plot_study_topo(study, kind; <keyword arguments>)

Plot a topographic series of a study: topographical maps in consecutive time windows (`:erp`) or frequency bands (`:psd`), one row of maps for the grand average of each level of a factor (e.g. each group or condition), optionally a row for the difference of two levels and a row for the statistic of `study_stats()`. Channels significant at any point of a window (`mask` of `study_stats()`, e.g. the channels of the significant clusters) are marked in all maps of the window (as FieldTrip `ft_clusterplot`). Requires the measure (see `study_measure!()`) and channel locations.

The value of a channel in a window is the mean over the points of the window; the recordings of one subject are averaged first (unless `factor = :session`), so each subject counts once. Bad channels (`NaN`) are left out of the averages; channels without data or locations are not drawn.

# Arguments

- `study::NeuroAnalyzer.STUDY`
- `kind::Symbol=:erp`: measure: `:erp` (time windows) or `:psd` (frequency bands)
- `factor::Symbol=:condition`: factor whose levels are plotted: `:condition` or `:session` (within subjects), or a design column, e.g. `:group` (between subjects)
- `levels::Union{Nothing, Vector{String}}=nothing`: levels to plot, in this order; `nothing`: all levels
- `condition::Union{Nothing, String}=nothing`: condition (between-subject factors and `:session`); may be omitted if the measure has one condition
- `group::Union{Nothing, String}=nothing`: subjects of one group (within-subject factors); `nothing`: all subjects
- `windows::Union{Nothing, Vector{<:Tuple{Real, Real}}}=nothing`: time windows in seconds (`:erp`) or frequency bands in Hz (`:psd`); `nothing`: consecutive windows of `width` in `lim`
- `width::Union{Nothing, Real}=nothing`: width of the consecutive windows; `nothing`: 1/6 of `lim`
- `lim::Union{Nothing, Tuple{Real, Real}}=nothing`: range of the consecutive windows; `nothing`: from 0 s to the end of the epoch (`:erp`) or all frequencies (`:psd`)
- `diff::Bool=false`: if `true`, add a row for the difference of the two plotted levels (first − second): for within-subject factors the mean of the differences of the subjects, for between-subject factors the difference of the grand averages
- `stats::Union{Nothing, NamedTuple}=nothing`: result of `study_stats()` for the same measure: adds a row for its statistic (t or F, mean over the window) and marks the significant channels
- `db::Bool=true`: for `:psd`, average the power of each subject in dB
- `locs::Union{Nothing, DataFrame}=nothing`: channel locations; `nothing`: the locations of the first recording
- `title::String="default"`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_study_topo(
    study::NeuroAnalyzer.STUDY,
    kind::Symbol = :erp;
    factor::Symbol = :condition,
    levels::Union{Nothing, Vector{String}} = nothing,
    condition::Union{Nothing, String} = nothing,
    group::Union{Nothing, String} = nothing,
    windows::Union{Nothing, Vector{<:Tuple{Real, Real}}} = nothing,
    width::Union{Nothing, Real} = nothing,
    lim::Union{Nothing, Tuple{Real, Real}} = nothing,
    diff::Bool = false,
    stats::Union{Nothing, NamedTuple} = nothing,
    db::Bool = true,
    locs::Union{Nothing, DataFrame} = nothing,
    title::String = "default",
    mono::Bool = false,
)::GLMakie.Figure
    # validate
    _check_var(kind, [:erp, :psd], "kind")
    m = _study_measure(study, kind)
    d, idx = _study_measure_rows(study, m)
    db = db && kind === :psd
    x = kind === :erp ? m.t : m.f
    if !isnothing(stats)
        (haskey(stats, :kind) && haskey(stats, :mask) && haskey(stats, :stat)) ||
            throw(ArgumentError("stats must be a result of study_stats()."))
        stats.kind === kind ||
            throw(ArgumentError("stats must be a result of study_stats() for :$kind."))
        (stats.labels == m.labels && length(kind === :erp ? stats.t : stats.f) == length(x)) ||
            throw(ArgumentError("stats and the :$kind measure have different channels or points."))
    end

    # windows
    if isnothing(windows)
        isnothing(lim) && (lim = kind === :erp ? (max(0.0, x[1]), x[end]) : (x[1], x[end]))
        lim[1] < lim[2] || throw(ArgumentError("lim[1] must be < lim[2]."))
        w = isnothing(width) ? (lim[2] - lim[1]) / 6 : width
        w > 0 || throw(ArgumentError("width must be > 0."))
        nw = ceil(Int64, (lim[2] - lim[1]) / w - 1e-9)
        windows = [(lim[1] + (i - 1) * w, min(lim[1] + i * w, lim[2])) for i = 1:nw]
    else
        isnothing(width) && isnothing(lim) ||
            throw(ArgumentError("width and lim are used only without windows."))
        isempty(windows) && throw(ArgumentError("windows must not be empty."))
    end
    # points of the windows: a ≤ x < b (the last point of the range belongs to its window)
    pts = map(windows) do (a, b)
        a < b || throw(ArgumentError("Window ($a, $b): start must be < end."))
        p = findall(v -> a <= v < b || (v == b && v == x[end]), x)
        p = isempty(p) && b == x[end] ? [lastindex(x)] : p
        isempty(p) && throw(ArgumentError("Window ($a, $b) contains no points."))
        p
    end
    nw = length(windows)

    # channel locations: channels of the measure with locations
    if isnothing(locs)
        obj = _silent(() -> load_neuro(_study_path(study, d.file[1])))
        locs = obj.locs
    end
    DataFrames.nrow(locs) > 0 ||
        throw(ArgumentError("Channel locations are empty: use load_locs() or the locs argument."))
    locs = _locs_unit(locs)
    lrow = indexin(m.labels, locs[!, :label])
    chi = findall(!isnothing, lrow)
    length(chi) >= 3 || throw(ArgumentError("At least 3 channels with locations are needed."))
    locs = locs[Int64.(lrow[chi]), :]
    length(unique(zip(locs.loc_radius, locs.loc_theta))) >= 3 || throw(
        ArgumentError("Channel locations have no coordinates: use load_locs() or the locs argument."),
    )

    # values of the subjects in the windows: (subject, channel, window)
    lv, A, subj = _study_level_data(m, d, idx; factor, levels, condition, group, db)
    diff && length(lv) != 2 &&
        throw(ArgumentError("diff requires 2 levels (levels: $(join(lv, ", ")))."))
    function _win_means(a::AbstractArray)
        W = fill(NaN, size(a, 1), length(chi), nw)
        for k in axes(a, 1), (j, c) in enumerate(chi), w = 1:nw
            v = Base.filter(!isnan, @view a[k, c, pts[w]])
            isempty(v) || (W[k, j, w] = mean(v))
        end
        return W
    end
    W = [_win_means(a) for a in A]

    # maps: (name, values (channel, window), color range)
    unit = kind === :erp ? "[$(m.unit[chi[1]])]" : db ? "[dB]" : "[$(m.unit[chi[1]])²/Hz]"
    sym(v) = (vm = maximum(abs, Base.filter(isfinite, v); init = 0.0); vm > 0 ? (-vm, vm) : (-1.0, 1.0))
    rng(v) = (f = Base.filter(isfinite, v); isempty(f) ? (-1.0, 1.0) :
              (e = extrema(f); e[1] < e[2] ? e : (e[1] - 1, e[2] + 1)))
    ga = [_nanstats(w)[1] for w in W]
    # one color range for the levels
    cr = kind === :erp ? sym(reduce(vcat, vec.(ga))) : rng(reduce(vcat, vec.(ga)))
    rows = [("$(lv[l])\n(n = $(length(subj[l])))", ga[l], cr, unit) for l in eachindex(lv)]
    if diff
        if factor in (:condition, :session)
            common = intersect(subj[1], subj[2])
            length(common) >= 1 ||
                throw(ArgumentError("diff requires subjects with data of both levels."))
            D = W[1][Int64.(indexin(common, subj[1])), :, :] .-
                W[2][Int64.(indexin(common, subj[2])), :, :]
            dv = _nanstats(D)[1]
        else
            dv = ga[1] .- ga[2]
        end
        push!(rows, ("$(lv[1]) − $(lv[2])", dv, sym(dv), unit))
    end
    sig = [Int64[] for _ = 1:nw]
    if !isnothing(stats)
        is_t = stats.test in (:ttest, :ttest_paired)
        sv = fill(NaN, length(chi), nw)
        for (j, c) in enumerate(chi), w = 1:nw
            v = Base.filter(!isnan, @view stats.stat[c, pts[w]])
            isempty(v) || (sv[j, w] = mean(v))
            any(@view stats.mask[c, pts[w]]) && push!(sig[w], j)
        end
        sl = length(stats.levels) == 2 ? "$(stats.levels[1]) vs $(stats.levels[2])" :
             join(stats.levels, ", ")
        push!(rows, ("$(is_t ? "t" : "F")\n$sl", sv, is_t ? sym(sv) : (0.0, rng(sv)[2]), ""))
        all(isempty, sig) && _info("No significant channels in the windows.")
    end

    # window labels
    function _wlabel(a, b)
        kind === :erp && ((a, b) = (1000a, 1000b))
        fmt(v) = isinteger(round(v; digits = 1)) ? string(round(Int64, v)) : string(round(v; digits = 1))
        return "$(fmt(a))–$(fmt(b)) $(kind === :erp ? "ms" : "Hz")"
    end

    # prepare plot
    GLMakie.activate!(; title = "plot_study_topo()")
    ms = 130
    nr = length(rows)
    fig = GLMakie.Figure(; size = (140 + nw * (ms + 10) + 90, 60 + nr * (ms + 30)))
    fl = factor === :condition ? "condition" : string(factor)
    tt = title == "default" ?
         (kind === :erp ? "Grand average ERP" : "Grand average PSD") * " by $fl" *
         (isnothing(stats) ? "" : ", significant channels marked ($(stats.test), correction: $(stats.correction))") :
         title
    isempty(tt) || GLMakie.Label(fig[0, 1:(nw + 2)], tt; fontsize = 16, font = :bold, tellwidth = false)
    for (r, (name, v, crange, u)) in enumerate(rows)
        GLMakie.Label(fig[r, 1], name; fontsize = 12, rotation = 0, tellheight = false,
                      halign = :right, justification = :right)
        for w = 1:nw
            ok = findall(isfinite, @view v[:, w])
            ax = GLMakie.Axis(
                fig[r, w + 1];
                width = ms,
                height = ms,
                aspect = DataAspect(),
                title = r == 1 ? _wlabel(windows[w]...) : "",
                titlesize = 11,
                titlegap = 2,
                titlefont = :regular,
                backgroundcolor = :transparent,
                _AXIS_LOCK_KWARGS...,
            )
            hidedecorations!(ax)
            hidespines!(ax)
            length(ok) >= 3 || continue
            pl = plot_topo(
                v[ok, w];
                locs = locs[ok, :],
                cb = false,
                ps = :m,
                nmethod = :none,
                crange = crange,
                mono = mono,
                ch_labels = false,
                electrodes = false,
            )
            img = _figure_image(pl)
            GLMakie.image!(ax, rotr90(img))
            # electrodes drawn here (markers of plot_topo would be too small in the series):
            # significant channels as large white circles
            ex, ey = _topo_image_xy(pl, img, locs)
            sw = intersect(sig[w], ok)
            ns = setdiff(ok, sw)
            GLMakie.scatter!(ax, ex[ns], ey[ns]; markersize = 3, color = :black)
            isempty(sw) || GLMakie.scatter!(
                ax, ex[sw], ey[sw];
                markersize = 8, color = :white, strokecolor = :black, strokewidth = 1.25,
            )
        end
        GLMakie.Colorbar(
            fig[r, nw + 2];
            colormap = _palette(:bluesreds; mono = mono),
            limits = crange,
            label = u,
            labelsize = 10,
            ticklabelsize = 9,
            ticks = _cb_tickfun(; n = 3),
            height = round(Int64, 0.8 * ms),
            width = 8,
            tellheight = false,
        )
    end
    GLMakie.colgap!(fig.layout, 8)
    GLMakie.rowgap!(fig.layout, 8)
    resize_to_layout!(fig)

    return fig
end

# positions of the channels `locs` (unit sphere, polar projection as in plot_topo()) in the
# image `img` of the topographical map `pl` (pixels from the bottom-left corner)
function _topo_image_xy(pl::GLMakie.Figure, img::AbstractMatrix, locs::DataFrame)
    ax = only(Base.filter(c -> c isa GLMakie.Axis, pl.content))
    vp = ax.scene.viewport[]
    lim = ax.finallimits[]
    sc = size(img, 2) / pl.scene.viewport[].widths[1]
    xy = [pol2cart(r, t) for (r, t) in zip(locs.loc_radius, locs.loc_theta)]
    px = [sc * (vp.origin[1] + (p[1] - lim.origin[1]) / lim.widths[1] * vp.widths[1]) for p in xy]
    py = [sc * (vp.origin[2] + (p[2] - lim.origin[2]) / lim.widths[2] * vp.widths[2]) for p in xy]
    return px, py
end

"""
    plot_study_tf(study; <keyword arguments>)

Plot group time-frequency maps of a study: the mean time-frequency power (dB) of each level of a factor (e.g. each group or condition) and optionally the difference of two levels, for one channel or averaged over channels (e.g. a channel cluster). The significant points of `study_stats()` are outlined (cluster contours). Requires the `:tf` measure (see `study_measure!()`).

The recordings of one subject are averaged first (unless `factor = :session`), so each subject counts once. The power of each subject is converted to dB before averaging (`db = true`) and, with a baseline `bl`, expressed as the change from the mean of the baseline at each frequency. With several channels, the map of each subject is the mean of the channels; bad channels of a subject (`NaN`) are left out of its mean.

# Arguments

- `study::NeuroAnalyzer.STUDY`
- `ch::Union{String, Vector{String}}`: channel, or channels averaged (e.g. the channels of a cluster)
- `factor::Symbol=:condition`: factor whose levels are plotted: `:condition` or `:session` (within subjects), or a design column, e.g. `:group` (between subjects)
- `levels::Union{Nothing, Vector{String}}=nothing`: levels to plot, in this order; `nothing`: all levels
- `condition::Union{Nothing, String}=nothing`: condition (between-subject factors and `:session`); may be omitted if the measure has one condition
- `group::Union{Nothing, String}=nothing`: subjects of one group (within-subject factors); `nothing`: all subjects
- `db::Bool=true`: average the power of the subjects in dB
- `bl::Union{Nothing, Tuple{Real, Real}}=nothing`: baseline window in seconds; the mean of the baseline at each frequency is subtracted from the map of each subject (with `db = true`: change in dB)
- `diff::Bool=false`: if `true`, add the map of the difference of the two plotted levels (first − second): for within-subject factors the mean of the differences of the subjects with data of both levels, for between-subject factors the difference of the means
- `stats::Union{Nothing, NamedTuple}=nothing`: result of `study_stats()` for the `:tf` measure; its significant points (`mask`) at the plotted channel(s) (with several channels: at any of them) are outlined in all maps
- `flim::Union{Nothing, Tuple{Real, Real}}=nothing`: frequency limits; `nothing`: all frequencies
- `tlim::Union{Nothing, Tuple{Real, Real}}=nothing`: time limits; `nothing`: the whole epoch
- `title::String="default"`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_study_tf(
    study::NeuroAnalyzer.STUDY;
    ch::Union{String, Vector{String}},
    factor::Symbol = :condition,
    levels::Union{Nothing, Vector{String}} = nothing,
    condition::Union{Nothing, String} = nothing,
    group::Union{Nothing, String} = nothing,
    db::Bool = true,
    bl::Union{Nothing, Tuple{Real, Real}} = nothing,
    diff::Bool = false,
    stats::Union{Nothing, NamedTuple} = nothing,
    flim::Union{Nothing, Tuple{Real, Real}} = nothing,
    tlim::Union{Nothing, Tuple{Real, Real}} = nothing,
    title::String = "default",
    mono::Bool = false,
)::GLMakie.Figure
    # validate
    m = _study_measure(study, :tf)
    d, idx = _study_measure_rows(study, m)
    chs = ch isa String ? [ch] : ch
    isempty(chs) && throw(ArgumentError("No channels selected."))
    chi = indexin(chs, m.labels)
    any(isnothing, chi) && throw(
        ArgumentError("Channels not in measure :tf: $(join(chs[isnothing.(chi)], ", "))."),
    )
    chi = Int64.(chi)
    length(unique(m.unit[chi])) == 1 ||
        throw(ArgumentError("All channels must have the same unit."))
    function _sel(v, lim, name)
        isnothing(lim) && return collect(eachindex(v))
        lim[1] < lim[2] || throw(ArgumentError("$name[1] must be < $name[2]."))
        i = findall(x -> lim[1] <= x <= lim[2], v)
        length(i) >= 2 || throw(ArgumentError("$name must contain at least 2 points."))
        return i
    end
    fi = _sel(m.f, flim, "flim")
    ti = _sel(m.t, tlim, "tlim")
    bi = nothing
    if !isnothing(bl)
        bl[1] < bl[2] || throw(ArgumentError("bl[1] must be < bl[2]."))
        bi = findall(x -> bl[1] <= x <= bl[2], m.t)
        isempty(bi) && throw(ArgumentError("Baseline $bl contains no time points."))
    end

    # subject maps (subject × frequency × time): mean of the channels, baseline subtracted
    lv, A, subj = _study_level_data(m, d, idx; factor, levels, condition, group, db)
    diff && length(lv) != 2 &&
        throw(ArgumentError("diff requires 2 levels (levels: $(join(lv, ", ")))."))
    S = map(A) do a
        R = fill(NaN, size(a, 1), size(a, 3), size(a, 4))
        for k in axes(a, 1), i in axes(a, 3), j in axes(a, 4)
            v = Base.filter(!isnan, @view a[k, chi, i, j])
            isempty(v) || (R[k, i, j] = mean(v))
        end
        if !isnothing(bi)
            for k in axes(R, 1), i in axes(R, 2)
                v = Base.filter(!isnan, @view R[k, i, bi])
                R[k, i, :] .-= isempty(v) ? NaN : mean(v)
            end
        end
        R[:, fi, ti]
    end
    # subjects without data at these channels
    keep = [[!all(isnan, @view s[k, :, :]) for k in axes(s, 1)] for s in S]
    subj = [subj[l][keep[l]] for l in eachindex(S)]
    S = [S[l][keep[l], :, :] for l in eachindex(S)]
    any(isempty, S) && throw(ArgumentError("A level has no subjects with data at the channel(s)."))
    maps = [_nanstats(s)[1] for s in S]
    if diff
        if factor in (:condition, :session)
            common = intersect(subj[1], subj[2])
            isempty(common) &&
                throw(ArgumentError("diff requires subjects with data of both levels."))
            dm = _nanstats(S[1][Int64.(indexin(common, subj[1])), :, :] .-
                           S[2][Int64.(indexin(common, subj[2])), :, :])[1]
        else
            dm = maps[1] .- maps[2]
        end
    end

    # significant points at the plotted channel(s)
    sig = nothing
    if !isnothing(stats)
        (haskey(stats, :kind) && stats.kind === :tf && haskey(stats, :mask)) ||
            throw(ArgumentError("stats must be a result of study_stats() for :tf."))
        size(stats.mask)[2:3] == size(m.data)[4:5] ||
            throw(ArgumentError("stats and the :tf measure have different frequencies or time points."))
        si = indexin(chs, stats.labels)
        any(isnothing, si) && throw(
            ArgumentError("Channels not in stats: $(join(chs[isnothing.(si)], ", "))."),
        )
        sig = dropdims(any(stats.mask[Int64.(si), :, :]; dims = 1); dims = 1)[fi, ti]
        any(sig) || _info("No significant points at the plotted channel(s).")
    end

    # color ranges: one for the levels, one for the difference
    fin(v) = Base.filter(isfinite, v)
    function sym(v)
        vm = maximum(abs, fin(v); init = 0.0)
        return vm > 0 ? (-vm, vm) : (-1.0, 1.0)
    end
    function rng(v)
        f = fin(v)
        isempty(f) && return (0.0, 1.0)
        e = extrema(f)
        return e[1] < e[2] ? e : (e[1] - 1, e[2] + 1)
    end
    all_m = reduce(vcat, vec.(maps))
    # diverging scale for changes from the baseline
    cr = isnothing(bl) ? rng(all_m) : sym(all_m)
    pal = isnothing(bl) ? _palette(:viridis; mono = mono) : _palette(:bluesreds; mono = mono)
    unit = m.unit[chi[1]]
    cb_lbl = db ? (isnothing(bl) ? "Power [dB $unit^2/Hz]" : "Power change [dB]") :
             (isnothing(bl) ? "Power [$unit^2/Hz]" : "Power change [$unit^2/Hz]")

    # prepare plot
    GLMakie.activate!(; title = "plot_study_tf()")
    np = length(maps) + diff
    fig = GLMakie.Figure(; size = (80 + np * 300 + (diff ? 160 : 80), 420))
    fl = factor === :condition ? "condition" : string(factor)
    ch_str = length(chs) == 1 ? chs[1] : "mean of $(length(chs)) channels"
    tt = title == "default" ?
         "Group time-frequency power, $ch_str, by $fl" *
         (isnothing(stats) ? "" : "\nsignificant points outlined ($(stats.test), correction: $(stats.correction))") :
         title
    t = m.t[ti]
    f = m.f[fi]
    function _panel!(pos, M, name, crange, cmap, first)
        ax = GLMakie.Axis(
            pos;
            xlabel = "Time [s]",
            ylabel = first ? "Frequency [Hz]" : "",
            title = name,
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            _AXIS_LOCK_KWARGS...,
        )
        _style_axis!(ax)
        first || GLMakie.hideydecorations!(ax; grid = false, ticks = false)
        hm = GLMakie.heatmap!(ax, t, f, permutedims(M); colormap = cmap, colorrange = crange,
                              nan_color = :lightgray)
        if !isnothing(sig) && any(sig)
            GLMakie.contour!(ax, _mask_padded(t, f, permutedims(sig))...; levels = [0.5],
                             color = :black, linewidth = 1.5)
        end
        t[1] < 0 < t[end] &&
            GLMakie.vlines!(ax, 0; color = :black, linestyle = :dash, linewidth = 1)
        GLMakie.limits!(ax, _cell_lims(t)..., _cell_lims(f)...)
        return hm
    end
    hm = nothing
    for (l, M) in enumerate(maps)
        hm = _panel!(fig[1, l], M, "$(lv[l]) (n = $(length(subj[l])))", cr, pal, l == 1)
    end
    GLMakie.Colorbar(fig[1, length(maps) + 1], hm; label = cb_lbl, ticks = _cb_tickfun(; n = 5))
    if diff
        hd = _panel!(fig[1, length(maps) + 2], dm, "$(lv[1]) − $(lv[2])", sym(dm),
                     _palette(:bluesreds; mono = mono), false)
        GLMakie.Colorbar(fig[1, length(maps) + 3], hd;
                         label = db ? "Difference [dB]" : "Difference [$unit^2/Hz]",
                         ticks = _cb_tickfun(; n = 5))
    end
    # (after the panels: spans all columns)
    isempty(tt) || GLMakie.Label(fig[0, :], tt; fontsize = 16, font = :bold, tellwidth = false)

    return fig
end

# edges of the cells of the points `x` (half a step around the first and last point)
function _cell_lims(x::AbstractVector)::Tuple{Float64, Float64}
    d = length(x) > 1 ? (x[2] - x[1]) / 2 : 0.5
    return (x[1] - d, x[end] + d)
end

# mask `M` (x × y) padded with `false` outside its points, so that the contours of regions
# at the edges are closed (between the edge point and the padding, at the cell edge)
function _mask_padded(x::AbstractVector, y::AbstractVector, M::AbstractMatrix{Bool})
    pad(v) = (d = length(v) > 1 ? v[2] - v[1] : 1.0; [v[1] - d; collect(v); v[end] + d])
    P = falses(size(M) .+ 2)
    P[2:(end - 1), 2:(end - 1)] .= M
    return pad(x), pad(y), Float64.(P)
end

"""
    plot_study_clusters(st; <keyword arguments>)

Plot the clusters of a cluster-based permutation test of a study measure (`study_stats()` with `correction = :cluster`, see `cluster_summary()`):

- left: extent of the clusters in time (`:erp`), frequency (`:psd`) or both (`:tf`: bounding boxes), labelled with their number of channels and p-value; positive clusters red, negative blue (F-tests: red); non-significant clusters light
- right: the permutation null distribution of the maximal cluster mass (histogram), the observed cluster masses (absolute values; vertical lines labelled with the cluster number) and the 95th percentile of the null distribution (dashed)

# Arguments

- `st::NamedTuple`: result of `study_stats()` with `correction = :cluster`
- `n::Int64=10`: maximal number of clusters drawn (those with the smallest p-values)
- `title::String="default"`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_study_clusters(
    st::NamedTuple;
    n::Int64 = 10,
    title::String = "default",
    mono::Bool = false,
)::GLMakie.Figure
    cs = cluster_summary(st)
    haskey(st, :null) && !isempty(st.null) ||
        throw(ArgumentError("st has no null distribution (run study_stats() again)."))
    n >= 1 || throw(ArgumentError("n must be ≥ 1."))
    cs = first(cs, n)
    nc = DataFrames.nrow(cs)
    is_t = st.test in (:ttest, :ttest_paired)

    pos_c = mono ? :black : :red
    neg_c = mono ? :gray40 : :blue
    col(r) = (r.sign < 0 && is_t ? neg_c : pos_c, r.significant ? 0.8 : 0.25)
    pstr(p) = p < 0.001 ? "p < 0.001" : "p = $(round(p; digits = 3))"

    GLMakie.activate!(; title = "plot_study_clusters()")
    fig = GLMakie.Figure(; size = (1100, 450))
    lv = length(st.levels) == 2 ? "$(st.levels[1]) vs $(st.levels[2])" : join(st.levels, ", ")
    tt = title == "default" ?
         "Clusters, $lv ($(st.test), $(length(st.null)) permutations)" : title

    # extent of the clusters
    x_name = st.kind === :psd ? "Frequency [Hz]" : "Time [s]"
    ax1 = GLMakie.Axis(
        fig[1, 1];
        xlabel = x_name,
        ylabel = st.kind === :tf ? "Frequency [Hz]" : "Cluster",
        title = "Extent",
        _AXIS_LOCK_KWARGS...,
    )
    _style_axis!(ax1)
    if st.kind === :tf
        for r in eachrow(cs)
            dt = length(st.t) > 1 ? (st.t[2] - st.t[1]) / 2 : 0.5
            df = length(st.f) > 1 ? (st.f[2] - st.f[1]) / 2 : 0.5
            rect = GLMakie.Rect2f(r.t_start - dt, r.f_start - df,
                                  r.t_end - r.t_start + 2dt, r.f_end - r.f_start + 2df)
            GLMakie.poly!(ax1, rect; color = col(r), strokecolor = :black, strokewidth = 1)
            GLMakie.text!(ax1, r.t_start - dt, r.f_end + df;
                          text = "$(r.id): $(r.n_channels) ch, $(pstr(r.p))",
                          align = (:left, :bottom), fontsize = 10)
        end
        GLMakie.limits!(ax1, _cell_lims(st.t)..., _cell_lims(st.f)...)
    else
        x = st.kind === :psd ? st.f : st.t
        a, b = st.kind === :psd ? (:f_start, :f_end) : (:t_start, :t_end)
        dx = length(x) > 1 ? (x[2] - x[1]) / 2 : 0.5
        for (i, r) in enumerate(eachrow(cs))
            y = nc - i + 1
            GLMakie.lines!(ax1, [r[a] - dx, r[b] + dx], [y, y];
                           color = col(r), linewidth = 14, linecap = :butt)
            # (left of the bar if it ends in the right third of the axis)
            right = r[b] > x[1] + 2 * (x[end] - x[1]) / 3
            GLMakie.text!(ax1, right ? r[a] - dx : r[b] + dx, y;
                          text = right ? "$(r.n_channels) ch, $(pstr(r.p))  " :
                                 "  $(r.n_channels) ch, $(pstr(r.p))",
                          align = (right ? :right : :left, :center), fontsize = 10)
        end
        ax1.yticks = (collect(nc:-1:1), string.(cs.id))
        GLMakie.limits!(ax1, _cell_lims(x)..., 0.4, max(nc, 1) + 0.6)
    end
    nc == 0 && GLMakie.text!(ax1, 0.5, 0.5; text = "no clusters", space = :relative,
                             align = (:center, :center))

    # null distribution and observed masses
    ax2 = GLMakie.Axis(
        fig[1, 2];
        xlabel = "Maximal cluster mass (|$(is_t ? "t" : "F")| sum)",
        ylabel = "Permutations",
        title = "Null distribution",
        _AXIS_LOCK_KWARGS...,
    )
    _style_axis!(ax2)
    # histogram over the range of the null distribution and the observed masses
    lo, hi = extrema(vcat(st.null, abs.(cs.mass)))
    hi > lo || (hi = lo + 1)
    edges = range(lo, hi; length = 41)
    cnt = zeros(Int64, 40)
    for v in st.null
        cnt[clamp(searchsortedlast(edges, v), 1, 40)] += 1
    end
    GLMakie.barplot!(ax2, collect(edges[1:(end - 1)]) .+ step(edges) / 2, cnt;
                     width = step(edges), gap = 0, color = (:gray, 0.6),
                     strokecolor = :gray30, strokewidth = 0.5, label = "null")
    ytop = maximum(cnt)
    q = quantile(st.null, 0.95)
    GLMakie.vlines!(ax2, q; color = :black, linestyle = :dash, linewidth = 1.5,
                    label = "95th percentile")
    for (i, r) in enumerate(eachrow(cs))
        GLMakie.vlines!(ax2, abs(r.mass); color = (col(r)[1], r.significant ? 1.0 : 0.4),
                        linewidth = 2)
        # (labels in 3 staggered rows)
        GLMakie.text!(ax2, abs(r.mass), (1.03 + 0.07 * mod(i - 1, 3)) * ytop;
                      text = string(r.id), align = (:center, :bottom), fontsize = 10)
    end
    GLMakie.ylims!(ax2, 0, 1.3 * ytop)
    GLMakie.axislegend(ax2; position = :rc, labelsize = 10, backgroundcolor = (:white, 0.85))

    GLMakie.Label(fig[0, :], tt; fontsize = 16, font = :bold, tellwidth = false)

    return fig
end

# shared arguments of plot_study_values() and plot_study_covariate()
const _STUDY_VALUES_DOC = """
- `ch::Union{Nothing, String, Vector{String}}=nothing`: channel(s) averaged; `nothing`: all channels
- `tlim::Union{Nothing, Tuple{Real, Real}}=nothing`: time window in seconds (`:erp`, `:tf`); `nothing`: the whole epoch
- `flim::Union{Nothing, Tuple{Real, Real}}=nothing`: frequency band in Hz (`:psd`, `:tf`); `nothing`: all frequencies
- `stats::Union{Nothing, NamedTuple}=nothing`: result of `study_stats()` for the same measure: the mean over the points of a cluster (`correction = :cluster`, see `cluster`) or over the significant points, instead of `ch`, `tlim` and `flim`
- `cluster::Union{Nothing, Int64}=nothing`: cluster ID (see `cluster_summary()`); `nothing`: the cluster with the smallest p-value
- `db::Bool=true`: for `:psd` and `:tf`, average the power of each subject in dB"""

# y-axis label of subject-level values
function _study_values_label(kind::Symbol, unit::String, db::Bool)::String
    kind === :erp && return "Mean amplitude [$unit]"
    return db ? "Mean power [dB $unit^2/Hz]" : "Mean power [$unit^2/Hz]"
end

"""
    plot_study_values(study, kind; <keyword arguments>)

Plot subject-level values of a study measure (see `study_values()`): the mean of the measure over channels and a time window or frequency band, or over the points of a cluster of `study_stats()` (e.g. the mean amplitude of a cluster), one point per subject, by the levels of a factor (e.g. groups or conditions) as strip, box or violin plots. For within-subject factors the points of one subject are connected by lines.

The recordings of one subject are averaged first (unless `factor = :session`), so each subject counts once.

# Arguments

- `study::NeuroAnalyzer.STUDY`
- `kind::Symbol`: measure (`:erp`, `:psd` or `:tf`)
$(_STUDY_VALUES_DOC)
- `factor::Symbol=:condition`: factor: `:condition` or `:session` (within subjects), or a column of the design table, e.g. `:group` or a categorical covariate (between subjects)
- `levels::Union{Nothing, Vector{String}}=nothing`: levels to plot, in this order; `nothing`: all levels
- `condition::Union{Nothing, String}=nothing`: condition (between-subject factors and `:session`); may be omitted if the measure has one condition
- `group::Union{Nothing, String}=nothing`: subjects of one group (within-subject factors); `nothing`: all subjects
- `type::Symbol=:box`: `:strip` (points and the mean), `:box` (box plots and points) or `:violin` (violin plots and points)
- `lines::Bool=true`: for within-subject factors, connect the points of one subject
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `title::String="default"`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_study_values(
    study::NeuroAnalyzer.STUDY,
    kind::Symbol;
    ch::Union{Nothing, String, Vector{String}} = nothing,
    tlim::Union{Nothing, Tuple{Real, Real}} = nothing,
    flim::Union{Nothing, Tuple{Real, Real}} = nothing,
    stats::Union{Nothing, NamedTuple} = nothing,
    cluster::Union{Nothing, Int64} = nothing,
    db::Bool = true,
    factor::Symbol = :condition,
    levels::Union{Nothing, Vector{String}} = nothing,
    condition::Union{Nothing, String} = nothing,
    group::Union{Nothing, String} = nothing,
    type::Symbol = :box,
    lines::Bool = true,
    xlabel::String = "default",
    ylabel::String = "default",
    title::String = "default",
    mono::Bool = false,
)::GLMakie.Figure
    # validate
    _check_var(type, [:strip, :box, :violin], "type")
    m = _study_measure(study, kind)
    within = factor in (:condition, :session)
    factor === :subject && throw(ArgumentError("factor must not be :subject."))
    !within && !isnothing(group) && throw(
        ArgumentError("group selects subjects of within-subject factors; use levels for $factor."),
    )
    if factor === :condition
        isnothing(condition) ||
            throw(ArgumentError("condition is not used with factor = :condition (see levels)."))
    elseif isnothing(condition)
        length(m.conditions) == 1 || throw(
            ArgumentError(
                "condition must be given (conditions: $(join(repr.(m.conditions), ", "))).",
            ),
        )
        condition = m.conditions[1]
    end
    v = _study_values(study, kind; ch, tlim, flim, stats, cluster, condition,
                      db, sessions = factor === :session)
    df = v.df
    factor in propertynames(df) ||
        throw(ArgumentError("Design table has no column $factor."))
    if !isnothing(group)
        df = df[df.group .== group, :]
        DataFrames.nrow(df) > 0 || throw(ArgumentError("No subjects of group $(repr(group))."))
    end
    df = df[.!ismissing.(df[!, factor]), :]
    col = string.(df[!, factor])
    all_lv = factor === :condition ? Base.filter(in(col), m.conditions) : unique(col)
    lv = isnothing(levels) ? all_lv : levels
    miss = setdiff(lv, all_lv)
    isempty(miss) || throw(
        ArgumentError(
            "Levels not found: $(join(repr.(miss), ", ")) (levels: $(join(repr.(all_lv), ", "))).",
        ),
    )
    isempty(lv) && throw(ArgumentError("No levels to plot."))
    if !within
        # a subject must belong to one level
        for s in unique(df.subject)
            length(unique(col[df.subject .== s])) == 1 || throw(
                ArgumentError("Subject $s has several levels of $factor (a within-subject factor?)."),
            )
        end
    end

    # points: x offsets of the subjects (the same in all levels, so that lines are parallel)
    subj = sort(unique(df.subject))
    ns = length(subj)
    off = Dict(s => ns > 1 ? -0.15 + 0.3 * (k - 1) / (ns - 1) : 0.0 for (k, s) in enumerate(subj))
    vals = [df.value[col .== l] for l in lv]
    ids = [df.subject[col .== l] for l in lv]

    # prepare plot
    GLMakie.activate!(; title = "plot_study_values()")
    fig = GLMakie.Figure(; size = (max(700, 300 + 120 * length(lv)), 500))
    fl = factor === :condition ? "condition" : string(factor)
    kn = kind === :erp ? "ERP" : kind === :psd ? "PSD" : "TF power"
    xl, yl, tt = _set_defaults(
        xlabel, ylabel, title,
        factor === :condition ? "Condition" : uppercasefirst(string(factor)),
        _study_values_label(kind, v.unit, v.db),
        "$kn, $(v.desc), by $fl" * (factor === :condition ? "" : ", condition $condition"),
    )
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel = xl,
        ylabel = yl,
        title = tt,
        xticks = (collect(eachindex(lv)), ["$(l)\n(n = $(length(vals[i])))" for (i, l) in enumerate(lv)]),
        _AXIS_LOCK_KWARGS...,
    )
    _style_axis!(ax)
    kind === :erp && GLMakie.hlines!(ax, 0; color = :gray, linewidth = 0.5)
    colors = mono ? _trace_colors(length(lv); mono = true) :
             GLMakie.Makie.wong_colors()[mod1.(1:length(lv), 7)]

    # distributions
    for (l, y) in enumerate(vals)
        isempty(y) && continue
        if type === :box
            GLMakie.boxplot!(ax, fill(l, length(y)), y; width = 0.5, show_outliers = false,
                             color = (colors[l], 0.3), strokecolor = colors[l],
                             strokewidth = 1, mediancolor = colors[l], whiskercolor = colors[l])
        elseif type === :violin && length(y) >= 2
            GLMakie.violin!(ax, fill(l, length(y)), y; width = 0.7, color = (colors[l], 0.3),
                            strokecolor = colors[l], strokewidth = 1)
        end
    end

    # lines connecting the levels of a subject
    if within && lines && length(lv) >= 2
        for s in subj
            xs = Float64[]
            ys = Float64[]
            for l in eachindex(lv)
                k = findfirst(==(s), ids[l])
                isnothing(k) && continue
                push!(xs, l + off[s])
                push!(ys, vals[l][k])
            end
            length(xs) >= 2 && GLMakie.lines!(ax, xs, ys; color = (:gray, 0.6), linewidth = 1)
        end
    end

    # points and (strip) the means
    for (l, y) in enumerate(vals)
        isempty(y) && continue
        GLMakie.scatter!(ax, [l + off[s] for s in ids[l]], y; color = colors[l],
                         markersize = 9, strokecolor = :black, strokewidth = 0.5)
        type === :strip && GLMakie.lines!(ax, [l - 0.25, l + 0.25], fill(mean(y), 2);
                                          color = :black, linewidth = 2.5)
    end
    GLMakie.xlims!(ax, 0.4, length(lv) + 0.6)

    return fig
end

"""
    plot_study_covariate(study, kind; <keyword arguments>)

Plot subject-level values of a study measure (see `study_values()`) against a numeric covariate of the design table (e.g. age): one point per subject, colored by the levels of a design column (e.g. group), with the least-squares regression line of all plotted subjects and its Pearson correlation.

The recordings of one subject are averaged first; a covariate that differs between the sessions of a subject is averaged.

# Arguments

- `study::NeuroAnalyzer.STUDY`
- `kind::Symbol`: measure (`:erp`, `:psd` or `:tf`)
- `covariate::Symbol`: numeric column of the design table
$(_STUDY_VALUES_DOC)
- `condition::Union{Nothing, String}=nothing`: condition; may be omitted if the measure has one condition
- `group::Union{Nothing, String}=nothing`: subjects of one group; `nothing`: all subjects
- `by::Union{Nothing, Symbol}=:group`: design column whose levels color the points; `nothing`: one color
- `fit::Bool=true`: if `true`, draw the regression line (with its 95% confidence band) and the Pearson correlation (r, p) in the legend
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `title::String="default"`: plot title
- `leg::Bool=true`: if `true`, draw a legend (outside the plot, on the right)
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_study_covariate(
    study::NeuroAnalyzer.STUDY,
    kind::Symbol;
    covariate::Symbol,
    ch::Union{Nothing, String, Vector{String}} = nothing,
    tlim::Union{Nothing, Tuple{Real, Real}} = nothing,
    flim::Union{Nothing, Tuple{Real, Real}} = nothing,
    stats::Union{Nothing, NamedTuple} = nothing,
    cluster::Union{Nothing, Int64} = nothing,
    db::Bool = true,
    condition::Union{Nothing, String} = nothing,
    group::Union{Nothing, String} = nothing,
    by::Union{Nothing, Symbol} = :group,
    fit::Bool = true,
    xlabel::String = "default",
    ylabel::String = "default",
    title::String = "default",
    leg::Bool = true,
    mono::Bool = false,
)::GLMakie.Figure
    # validate
    m = _study_measure(study, kind)
    if isnothing(condition)
        length(m.conditions) == 1 || throw(
            ArgumentError(
                "condition must be given (conditions: $(join(repr.(m.conditions), ", "))).",
            ),
        )
        condition = m.conditions[1]
    end
    covariate in (:subject, :session, :group, :file, :include) &&
        throw(ArgumentError("covariate must be a covariate of the design table."))
    covariate in propertynames(study.design) ||
        throw(ArgumentError("Design table has no column $covariate."))
    eltype(study.design[!, covariate]) <: Union{Missing, Real} ||
        throw(ArgumentError("Covariate $covariate must be numeric."))
    isnothing(by) || by in propertynames(study.design) ||
        throw(ArgumentError("Design table has no column $by."))
    v = _study_values(study, kind; ch, tlim, flim, stats, cluster, condition, db,
                      sessions = false)
    df = v.df
    isnothing(group) || (df = df[df.group .== group, :])
    df = df[.!ismissing.(df[!, covariate]), :]
    isnothing(by) || (df = df[.!ismissing.(df[!, by]), :])
    df = df[isfinite.(Float64.(df[!, covariate])), :]
    DataFrames.nrow(df) > 0 || throw(ArgumentError("No subjects with values of $covariate."))
    x = Float64.(df[!, covariate])
    y = df.value
    lv = isnothing(by) ? ["all"] : unique(string.(df[!, by]))
    col = isnothing(by) ? fill("all", length(y)) : string.(df[!, by])

    # prepare plot
    GLMakie.activate!(; title = "plot_study_covariate()")
    fig = GLMakie.Figure(; size = (800, 500))
    kn = kind === :erp ? "ERP" : kind === :psd ? "PSD" : "TF power"
    xl, yl, tt = _set_defaults(
        xlabel, ylabel, title,
        string(covariate), _study_values_label(kind, v.unit, v.db),
        "$kn, $(v.desc), condition $condition, vs $covariate",
    )
    ax = GLMakie.Axis(fig[1, 1]; xlabel = xl, ylabel = yl, title = tt, _AXIS_LOCK_KWARGS...)
    _style_axis!(ax)
    colors = isnothing(by) ? [mono ? :black : GLMakie.Makie.wong_colors()[1]] :
             mono ? _trace_colors(length(lv); mono = true) :
             GLMakie.Makie.wong_colors()[mod1.(1:length(lv), 7)]

    # regression line and its 95% confidence band
    n = length(x)
    labelled = !isnothing(by)
    if fit
        if n >= 3 && std(x) > 0
            b = cov(x, y) / var(x)
            a = mean(y) - b * mean(x)
            r = cor(x, y)
            res = y .- (a .+ b .* x)
            s = sqrt(sum(abs2, res) / (n - 2))
            xs = collect(range(extrema(x)...; length = 100))
            h = quantile(TDist(n - 2), 0.975) .* s .*
                sqrt.(1 / n .+ (xs .- mean(x)) .^ 2 ./ sum(abs2, x .- mean(x)))
            yf = a .+ b .* xs
            tv = abs(r) < 1 ? r * sqrt((n - 2) / (1 - r^2)) : sign(r) * Inf
            p = 2 * ccdf(TDist(n - 2), abs(tv))
            pstr = p < 0.001 ? "p < 0.001" : "p = $(round(p; digits = 3))"
            GLMakie.band!(ax, xs, yf .- h, yf .+ h; color = (:gray, 0.2))
            GLMakie.lines!(ax, xs, yf; color = :black, linewidth = 2,
                           label = "r = $(round(r; digits = 2)), $pstr (n = $n)")
            labelled = true
        else
            _info("Fewer than 3 subjects or a constant covariate: regression line is not drawn.")
        end
    end
    for (l, name) in enumerate(lv)
        sel = col .== name
        GLMakie.scatter!(ax, x[sel], y[sel]; color = colors[l], markersize = 10,
                         strokecolor = :black, strokewidth = 0.5,
                         label = isnothing(by) ? nothing : "$name (n = $(count(sel)))")
    end

    leg && labelled && _side_legend!(fig, ax)

    return fig
end

"""
    plot_study_overview(study; <keyword arguments>)
    plot_study_overview(ov; <keyword arguments>)

Plot an overview of the study recordings (see `study_overview()`):

- left: data availability and the number of epochs (without the rejected ones) of each subject (rows) in each session and condition (columns); `–`: no recording, red numbers: no epochs of the condition, gray cells: recordings excluded from the analyses
- right: bad channels of each recording (rows) (`all_channels = false`: only the channels bad in at least one recording)

# Arguments

- `study::NeuroAnalyzer.STUDY`
- `ov::NamedTuple`: result of `study_overview()`
- `all::Bool=true`: if `true`, include the recordings excluded from the analyses (`study` only)
- `all_channels::Bool=false`: if `true`, plot all channels in the bad channels panel
- `title::String="default"`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_study_overview(
    study::NeuroAnalyzer.STUDY;
    all::Bool = true,
    all_channels::Bool = false,
    title::String = "default",
    mono::Bool = false,
)::GLMakie.Figure
    tt = title == "default" ? (isempty(study.name) ? "Study overview" : "Study $(study.name)") :
         title
    return plot_study_overview(study_overview(study; all); all_channels, title = tt, mono)
end

function plot_study_overview(
    ov::NamedTuple;
    all_channels::Bool = false,
    title::String = "default",
    mono::Bool = false,
)::GLMakie.Figure
    Base.all(k -> haskey(ov, k), (:recordings, :epochs, :conditions, :bad, :labels)) ||
        throw(ArgumentError("ov must be a result of study_overview()."))
    r = ov.recordings
    e = ov.epochs
    subj = unique(r.subject)
    ses = unique(r.session)
    cs = ov.conditions
    multi_ses = length(ses) > 1
    cname(c) = isempty(c) ? "(none)" : c

    # epochs: subject × (session, condition); NaN: no recording
    cols = [(s, c) for s in ses for c in cs]
    N = fill(NaN, length(subj), length(cols))
    X = falses(length(subj), length(cols))
    for x in eachrow(e)
        i = findfirst(==(x.subject), subj)
        j = findfirst(==((x.session, x.condition)), cols)
        N[i, j] = x.n_epochs
        X[i, j] = !x.include
    end
    # recordings that could not be loaded: no data
    for x in eachrow(r)
        x.ok && continue
        i = findfirst(==(x.subject), subj)
        for (j, (s, _)) in enumerate(cols)
            s == x.session && (N[i, j] = NaN)
        end
    end
    col_lbl = [multi_ses ? (length(cs) == 1 ? "ses $s" : "ses $s\n$(cname(c))") : cname(c)
               for (s, c) in cols]

    # bad channels: recording × channel
    rec_lbl = multi_ses ? r.subject .* "/" .* r.session : copy(r.subject)
    bi = all_channels ? collect(eachindex(ov.labels)) : findall(vec(any(ov.bad; dims = 1)))
    B = Float64.(ov.bad[:, bi])

    GLMakie.activate!(; title = "plot_study_overview()")
    nr = max(length(subj), length(rec_lbl))
    w = clamp(450 + 70 * length(cols) + 22 * max(length(bi), 8), 900, 1800)
    h = clamp(170 + 22 * nr, 400, 1400)
    fig = GLMakie.Figure(; size = (w, h))
    tt = title == "default" ? "Study overview" : title
    fs = nr > 40 ? 7 : 10

    # availability and epochs
    ax1 = GLMakie.Axis(
        fig[1, 1];
        title = "Epochs (subject × " * (multi_ses ? "session × condition)" : "condition)"),
        xticks = (collect(eachindex(cols)), col_lbl),
        yticks = (collect(eachindex(subj)), subj),
        yreversed = true,
        xticklabelsize = fs,
        yticklabelsize = fs,
        _AXIS_LOCK_KWARGS...,
    )
    _style_axis!(ax1)
    ax1.xgridvisible = false
    ax1.ygridvisible = false
    mx = Base.all(isnan, N) ? 1.0 : max(1.0, maximum(Base.filter(!isnan, N)))
    hm = GLMakie.heatmap!(
        ax1, 1:length(cols), 1:length(subj), permutedims(N);
        colormap = mono ? GLMakie.Reverse(:grays) : :Blues,
        colorrange = (0, mx),
        nan_color = :white,
    )
    any(X) && GLMakie.heatmap!(
        ax1, 1:length(cols), 1:length(subj), permutedims(ifelse.(X, 1.0, NaN));
        colormap = [(:gray, 0.6), (:gray, 0.6)], nan_color = :transparent,
    )
    if length(N) <= 2000
        for i in axes(N, 1), j in axes(N, 2)
            v = N[i, j]
            txt = isnan(v) ? "–" : string(Int64(v))
            c = isnan(v) ? :gray40 : v == 0 ? (mono ? :black : :red) :
                v > 0.6 * mx ? :white : :black
            GLMakie.text!(ax1, j, i; text = txt, align = (:center, :center), fontsize = fs,
                          color = c, font = v == 0 ? :bold : :regular)
        end
    end
    GLMakie.Colorbar(fig[1, 2], hm; label = "Epochs", width = 12)

    # bad channels
    ax2 = GLMakie.Axis(
        fig[1, 3];
        title = "Bad channels" * (all_channels ? "" : " ($(length(bi))/$(length(ov.labels)))"),
        xticks = (collect(eachindex(bi)), ov.labels[bi]),
        yticks = (collect(eachindex(rec_lbl)), rec_lbl .* " (" .* string.(r.n_bad) .* ")"),
        yreversed = true,
        xticklabelrotation = pi / 2,
        xticklabelsize = fs,
        yticklabelsize = fs,
        _AXIS_LOCK_KWARGS...,
    )
    _style_axis!(ax2)
    ax2.xgridvisible = false
    ax2.ygridvisible = false
    if isempty(bi)
        GLMakie.text!(ax2, 0.5, 0.5; text = "no bad channels", space = :relative,
                      align = (:center, :center))
        GLMakie.hidexdecorations!(ax2)
    else
        GLMakie.heatmap!(
            ax2, 1:length(bi), 1:length(rec_lbl), permutedims(B);
            colormap = [:white, mono ? :black : :red], colorrange = (0, 1),
        )
        ex = .!r.include
        any(ex) && GLMakie.heatmap!(
            ax2, 1:length(bi), 1:length(rec_lbl),
            permutedims(repeat(ifelse.(ex, 1.0, NaN), 1, length(bi)));
            colormap = [(:gray, 0.4), (:gray, 0.4)], nan_color = :transparent,
        )
    end
    GLMakie.colsize!(fig.layout, 1, GLMakie.Auto(max(length(cols), 2)))
    GLMakie.colsize!(fig.layout, 3, GLMakie.Auto(max(length(bi), 8) / 3))

    GLMakie.Label(fig[0, :], tt; fontsize = 16, font = :bold, tellwidth = false)

    return fig
end
