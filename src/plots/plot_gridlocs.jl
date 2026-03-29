export plot_gridlocs

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

const _GRIDLOCS_X = Float64[
    -0.5,  0.0,  0.5,
    -1.0, -0.5,  0.0,  0.5,  1.0,
    -1.0, -0.5,  0.0,  0.5,  1.0,
    -1.0, -0.5,  0.0,  0.5,  1.0,
    -0.5,  0.0,  0.5,
]
const _GRIDLOCS_Y = Float64[
     1.0,  1.0,  1.0,
     0.5,  0.5,  0.5,  0.5,  0.5,
     0.0,  0.0,  0.0,  0.0,  0.0,
    -0.5, -0.5, -0.5, -0.5, -0.5,
    -1.0, -1.0, -1.0,
]
const _GRIDLOCS_LABELS = [
    "Fp1", "Fpz", "Fp2",
    "F7",  "F3",  "Fz",  "F4",  "F8",
    "T3",  "C3",  "Cz",  "C4",  "T4",
    "T5",  "P3",  "Pz",  "P4",  "T6",
    "O1",  "Oz",  "O2",
]

# ---------------------------------------------------------------------------


"""
    plot_gridlocs()

Plot a simplified grid layout of standard 10-20 EEG channels for reference.

# Arguments

- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_gridlocs(; mono::Bool = false)::GLMakie.Figure
    # set color palette
    pal = mono ? :grays : :darktest

    ch_n = length(_GRIDLOCS_LABELS)
    cmap = GLMakie.resample_cmap(pal, ch_n)

    GLMakie.activate!(; title = "plot_gridlocs()")
    fig = GLMakie.Figure(; size = (800, 800), figure_padding = 0)
    ax  = GLMakie.Axis(
        fig[1, 1];
        aspect = 1,
        xlabel = "",
        ylabel = "",
        title  = "",
        xzoomlock  = true,
        yzoomlock  = true,
        xpanlock   = true,
        ypanlock   = true,
        xrectzoom  = false,
        yrectzoom  = false,
    )
    hidedecorations!(ax; grid = true)
    hidespines!(ax)
    GLMakie.xlims!(ax, (-1.2, 1.2))
    GLMakie.ylims!(ax, (-1.2, 1.2))

    # outer border — thin lines
    for (xs, ys) in [
        ([-1.0,  1.0], [-1.0, -1.0]),
        ([-1.0,  1.0], [ 1.0,  1.0]),
        ([-1.0, -1.0], [-1.0,  1.0]),
        ([ 1.0,  1.0], [-1.0,  1.0]),
    ]
        GLMakie.lines!(ax, xs, ys; color = :black, linewidth = 0.2)
    end

    # diagonal corner cuts
    for (xs, ys) in [
        ([-1.0, -0.5], [ 0.5,  1.0]),
        ([ 0.5,  1.0], [ 1.0,  0.5]),
        ([-1.0, -0.5], [-0.5, -1.0]),
        ([ 0.5,  1.0], [-1.0, -0.5]),
    ]
        GLMakie.lines!(ax, xs, ys; color = :black, linewidth = 0.5)
    end

    # horizontal grid lines
    for y in [-0.5, 0.0, 0.5]
        GLMakie.lines!(ax, [-1.0, 1.0], [y, y]; color = :black, linewidth = 0.5)
    end
    # partial horizontal lines at top and bottom (clipped by corners)
    GLMakie.lines!(ax, [-0.5,  0.5], [ 1.0,  1.0]; color = :black, linewidth = 0.5)
    GLMakie.lines!(ax, [-0.5,  0.5], [-1.0, -1.0]; color = :black, linewidth = 0.5)

    # vertical grid lines
    for x in [-0.5, 0.0, 0.5]
        GLMakie.lines!(ax, [x, x], [-1.0, 1.0]; color = :black, linewidth = 0.5)
    end
    # partial vertical lines at sides (clipped by corners)
    GLMakie.lines!(ax, [-1.0, -1.0], [-0.5,  0.5]; color = :black, linewidth = 0.5)
    GLMakie.lines!(ax, [ 1.0,  1.0], [-0.5,  0.5]; color = :black, linewidth = 0.5)

    # draw all channel markers in one call
    GLMakie.scatter!(
        ax, _GRIDLOCS_X, _GRIDLOCS_Y;
        color      = cmap,
        colormap   = pal,
        colorrange = 1:ch_n,
        markersize  = 16.0,
        strokewidth = 2,
        strokecolor = :black,
    )

    # draw all labels in one call
    GLMakie.text!(
        ax,
        _GRIDLOCS_X .+ 0.015,
        _GRIDLOCS_Y .+ 0.015;
        text     = _GRIDLOCS_LABELS,
        fontsize = 16,
    )

    return fig
end