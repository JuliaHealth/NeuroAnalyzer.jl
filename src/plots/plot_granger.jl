export plot_granger

"""
    plot_granger(g; <keyword arguments>)

Plot Granger causality returned by `granger()` (matrix: rows = source, columns = target) or `granger_spectral()` (grid of spectra: panel in row i, column j = i → j).

# Arguments

- `g::NamedTuple`: output of `granger()` or `granger_spectral()`
- `labels::Union{Nothing, Vector{String}}=nothing`: channel labels; if `nothing`, `g.labels` (NEURO methods) or channel numbers
- `alpha::Float64=0.05`: significance level
- `sig::Bool=true`: if `true`, mark significant connections:
    - `granger()`: cells with FDR-adjusted permutation p-value (`p_perm_adj`) < `alpha`, or FDR-adjusted F-test p-value (`p_adj`) if the permutation test was not run (marked with `*`)
    - `granger_spectral()`: frequencies with p-value corrected over frequencies (`p_perm_max`) < `alpha` (thick line); requires the permutation test
- `values::Bool=true`: if `true`, show the values in the matrix cells (`granger()` only, up to 10 channels)
- `flim::Union{Nothing, Tuple{Real, Real}}=nothing`: frequency limits (`granger_spectral()` only); if `nothing`, 0 to Nyquist
- `title::String="Granger causality"`: plot title
- `cb_title::String="GC"`: colorbar title (`granger()` only)
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_granger(
    g::NamedTuple;
    labels::Union{Nothing, Vector{String}} = nothing,
    alpha::Float64 = 0.05,
    sig::Bool = true,
    values::Bool = true,
    flim::Union{Nothing, Tuple{Real, Real}} = nothing,
    title::String = "Granger causality",
    cb_title::String = "GC",
    mono::Bool = false,
)::GLMakie.Figure
    # validate
    hasproperty(g, :gc) ||
        throw(ArgumentError("g must be the output of granger() or granger_spectral()."))
    n = size(g.gc, 1)
    isnothing(labels) && (labels = hasproperty(g, :labels) ? g.labels : string.(1:n))
    length(labels) == n ||
        throw(ArgumentError("Length of labels ($(length(labels))) must be equal to the number of channels ($n)."))
    0 < alpha < 1 || throw(ArgumentError("alpha must be in (0, 1)."))

    if ndims(g.gc) == 2
        return _plot_granger_matrix(g, labels, alpha, sig, values, title, cb_title, mono)
    elseif ndims(g.gc) == 3 && hasproperty(g, :frq)
        isnothing(flim) && (flim = (0, g.frq[end]))
        _check_tuple(flim, extrema(g.frq), "flim")
        return _plot_granger_spectral(g, labels, alpha, sig, flim, title, mono)
    else
        throw(ArgumentError("g must be the output of granger() or granger_spectral()."))
    end
end

function _plot_granger_matrix(
    g::NamedTuple,
    labels::Vector{String},
    alpha::Float64,
    sig::Bool,
    values::Bool,
    title::String,
    cb_title::String,
    mono::Bool,
)::GLMakie.Figure
    n = size(g.gc, 1)
    m = copy(g.gc)
    # the diagonal is undefined
    for i = 1:n
        m[i, i] = NaN
    end
    # significance: permutation test if it was run, otherwise F-test
    pv = if hasproperty(g, :p_perm_adj) && !all(isnan, g.p_perm_adj)
        g.p_perm_adj
    elseif hasproperty(g, :p_adj)
        g.p_adj
    else
        nothing
    end

    pal = _palette(:viridis; mono = mono)

    # prepare plot
    GLMakie.activate!(; title = "plot_granger()")
    fig = GLMakie.Figure(; size = (800, 800))

    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel             = "To",
        ylabel             = "From",
        title              = title,
        xticks             = (1:n, labels),
        xticklabelrotation = n > 6 ? π / 2 : 0.0,
        xticksvisible      = false,
        yticks             = (1:n, labels),
        yticksvisible      = false,
        yreversed          = true,
        aspect             = GLMakie.DataAspect(),
        xautolimitmargin   = (0, 0),
        yautolimitmargin   = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )
    _style_axis!(ax)

    # m[i, j] (from i, to j) is drawn at x = j, y = i
    hm = GLMakie.heatmap!(ax, 1:n, 1:n, m'; colormap = pal, nan_color = :white)
    GLMakie.Colorbar(
        fig[1, 2],
        hm;
        label = cb_title,
        labelsize = 16,
        ticks = _cb_tickfun(),
        height = GLMakie.lift(vp -> vp.widths[2], ax.scene.viewport),
        tellheight = false,
    )

    # cell annotations: values and/or significance stars
    lo, hi = extrema(filter(!isnan, m))
    show_values = values && n <= 10
    for i = 1:n, j = 1:n
        i == j && continue
        star = sig && !isnothing(pv) && pv[i, j] < alpha ? "*" : ""
        s = show_values ? _cb_label(m[i, j], max(2, ceil(Int64, -log10(max(hi, eps()))) + 2)) * star : star
        isempty(s) && continue
        # text color contrasting with the cell
        light = hi > lo && (m[i, j] - lo) / (hi - lo) > 0.6
        GLMakie.text!(
            ax,
            j,
            i;
            text = s,
            align = (:center, :center),
            fontsize = n <= 6 ? 16 : 12,
            color = light ? :black : :white,
        )
    end
    if sig && !isnothing(pv)
        src = hasproperty(g, :p_perm_adj) && !all(isnan, g.p_perm_adj) ? "permutation test" : "F-test"
        GLMakie.Label(
            fig[2, 1],
            "* p < $alpha ($src, FDR-adjusted)";
            fontsize = 12,
            halign = :left,
            tellwidth = false,
        )
    end

    return fig
end

function _plot_granger_spectral(
    g::NamedTuple,
    labels::Vector{String},
    alpha::Float64,
    sig::Bool,
    flim::Tuple{Real, Real},
    title::String,
    mono::Bool,
)::GLMakie.Figure
    n = size(g.gc, 1)
    f = g.frq
    fidx = findall(x -> flim[1] <= x <= flim[2], f)
    has_perm = hasproperty(g, :p_perm_max) && !all(isnan, g.p_perm_max)

    # common y-axis limits
    ymax = maximum(g.gc[i, j, k] for i = 1:n, j = 1:n, k in fidx if i != j)
    ymax = ymax > 0 ? 1.05 * ymax : 1.0

    # prepare plot
    GLMakie.activate!(; title = "plot_granger()")
    ps = n <= 3 ? 300 : n <= 6 ? 200 : 120
    fig = GLMakie.Figure(; size = (max(600, n * ps + 100), max(400, n * ps + 100)))
    GLMakie.Label(fig[0, 1:n], title; fontsize = 18, font = :bold)

    for i = 1:n, j = 1:n
        if i == j
            # diagonal: channel label
            ax = GLMakie.Axis(fig[i, j]; _AXIS_LOCK_KWARGS...)
            GLMakie.hidedecorations!(ax)
            GLMakie.hidespines!(ax)
            GLMakie.text!(
                ax,
                0.5,
                0.5;
                text = labels[i],
                align = (:center, :center),
                fontsize = n <= 6 ? 20 : 14,
                font = :bold,
                space = :relative,
            )
            continue
        end
        ax = GLMakie.Axis(
            fig[i, j];
            title = n <= 6 ? "$(labels[i]) → $(labels[j])" : "",
            titlesize = 12,
            titlefont = :regular,
            xlabel = i == n || (i == n - 1 && j == n) ? "Frequency [Hz]" : "",
            ylabel = j == 1 || (j == 2 && i == 1) ? "GC" : "",
            xticks = GLMakie.LinearTicks(4),
            yticks = GLMakie.LinearTicks(3),
            xlabelsize = 12,
            ylabelsize = 12,
            xticklabelsize = 10,
            yticklabelsize = 10,
            xautolimitmargin = (0, 0),
            _AXIS_LOCK_KWARGS...,
        )
        GLMakie.xlims!(ax, flim)
        GLMakie.ylims!(ax, 0, ymax)
        GLMakie.lines!(ax, f[fidx], g.gc[i, j, fidx]; linewidth = 1.5, color = :black)
        if sig && has_perm
            y = [g.p_perm_max[i, j, k] < alpha ? g.gc[i, j, k] : NaN for k in fidx]
            GLMakie.lines!(ax, f[fidx], y; linewidth = 4, color = mono ? :gray40 : :red)
        end
        # tick labels only on the outer panels
        i == n || (i == n - 1 && j == n) || GLMakie.hidexdecorations!(ax; grid = false, ticks = false)
        j == 1 || (j == 2 && i == 1) || GLMakie.hideydecorations!(ax; grid = false, ticks = false)
    end
    if sig && has_perm
        GLMakie.Label(
            fig[n + 1, 1:n],
            "thick line: p < $alpha (permutation test, corrected over frequencies)";
            fontsize = 12,
            halign = :left,
            tellwidth = false,
        )
    end
    GLMakie.rowgap!(fig.layout, 8)
    GLMakie.colgap!(fig.layout, 8)

    return fig
end
