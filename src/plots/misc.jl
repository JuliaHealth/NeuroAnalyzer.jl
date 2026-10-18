export plot_compose
export plot_empty
export add_inset

"""
    plot_compose(vfig; <keyword arguments>)

Compose a grid of `GLMakie.Figure` plots into a single figure using the specified row × column layout.

Empty plots are added automatically when `length(vfig) < layout[1]*layout[2]` so the grid is always fully populated.

# Arguments

- `vfig::Vector{GLMakie.Figure}`: plots to compose
- `layout::Tuple{Int64, Int64}`: grid dimensions as `(rows, columns)`

# Returns

- `GLMakie.Figure`: composite figure
"""
function plot_compose(
    vfig::Vector{GLMakie.Figure};
    layout::Tuple{Int64, Int64},
)::GLMakie.Figure
    # validate that the layout can accommodate all provided plots.
    layout[1] * layout[2] >= length(vfig) ||
        throw(
            ArgumentError(
                "Layout ($(layout[1]) × $(layout[2])) must be ≥ number of plots ($(length(vfig))).",
            ),
        )

    plot_size = (0, 0)
    for idx in eachindex(vfig)
        s = size(vfig[idx].scene)
        prod(s) > prod(plot_size) && (plot_size = s)
    end

    # warn when plots have mismatched sizes — compositing works but may look uneven.
    for idx in eachindex(vfig)
        size(vfig[idx].scene) != plot_size &&
            _warn("For best results all plots should be $(plot_size[1])×$(plot_size[2]).")
    end

    # total canvas size: columns determine width, rows determine height
    canvas_size = (plot_size[1] * layout[2], plot_size[2] * layout[1])

    # pad (a copy of) the vector with empty plots to fill every grid cell
    n_empty = layout[1] * layout[2] - length(vfig)
    vfig = vcat(vfig, [plot_empty() for _ = 1:n_empty])

    # build the composite figure using a GridLayout so axes are properly nested
    GLMakie.activate!(; title = "plot_compose()")
    pc = GLMakie.Figure(; size = canvas_size)
    gl = pc[1, 1] = GridLayout(layout[1], layout[2])

    p_idx = 1
    for idx1 = 1:layout[1], idx2 = 1:layout[2]
        # render each sub-figure to a temporary PNG, load it as a raster image,
        # then display it in a decoration-free axis.
        fname = tempname() * ".png"
        try
            GLMakie.save(fname, vfig[p_idx])
            pp = FileIO.load(fname)
            # place the axis inside the GridLayout, not the Figure directly.
            ax = GLMakie.Axis(
                gl[idx1, idx2];
                aspect = DataAspect(),
                xzoomlock = true,
                yzoomlock = true,
                xpanlock = true,
                ypanlock = true,
                xrectzoom = false,
                yrectzoom = false,
            )
            # rotr90 aligns the image orientation with GLMakie's y-up convention.
            GLMakie.image!(ax, rotr90(pp))
            hidedecorations!(ax)
            hidespines!(ax)
        finally
            isfile(fname) && rm(fname)
        end
        p_idx += 1
    end

    return pc
end

"""
    plot_empty()

Return an empty `GLMakie.Figure`, useful for padding a grid of plots.

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_empty()::GLMakie.Figure
    return GLMakie.Figure()
end

"""
    add_inset(fig, pl; <keyword arguments>)

Add plot `pl` (typically a locations plot, e.g. `plot_locs(obj; ps = :s)`) as an inset of `fig`. `pl` is rendered to an image whose white background is made transparent, so only its content (head outline, electrodes) is drawn.

# Arguments

- `fig::GLMakie.Figure`: primary plot (axis at `fig[1, 1]`)
- `pl::GLMakie.Figure`: plot to add
- `pos::Symbol=:rt`: inset position:
    - `:outside`: to the right of the plot (top-aligned)
    - `:rt`, `:lt`, `:rb`, `:lb`: inside the plot, in the right-top, left-top, right-bottom or left-bottom corner
- `width::Union{Nothing, Int64}=nothing`: inset width in pixels (aspect ratio is kept); `nothing`: the size of `pl`

# Returns

- `GLMakie.Figure`: `fig` with the inset (modified in place)
"""
function add_inset(
    fig::GLMakie.Figure,
    pl::GLMakie.Figure;
    pos::Symbol = :rt,
    width::Union{Nothing, Int64} = nothing,
)::GLMakie.Figure
    _figure_inset!(fig, pl; pos = pos, width = width)

    return fig
end
