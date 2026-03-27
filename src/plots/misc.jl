export plot_compose
export plot_empty
export add_pl

"""
    plot_compose(vfig; <keyword arguments>)

Compose a grid of `GLMakie.Figure` plots into a single figure using the specified row × column layout.

Empty plots are added automatically when `length(vfig) < layout[1]*layout[2]` so the grid is always fully populated.


# Arguments

- `vfig::Vector{GLMakie.Figure}`: plots to compose
- `layout::Tuple{Int64, Int64}`: grid dimensions as `(rows, columns)`

# Returns

- `GLMakie.Figure`: the plotted figure: composite figure
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

    # pad the vector with empty plots to fill every grid cell
    n_empty = layout[1] * layout[2] - length(vfig)
    for _ in 1:n_empty
        push!(vfig, plot_empty())
    end

    # build the composite figure using a GridLayout so axes are properly nested
    GLMakie.activate!(; title = "plot_compose()")
    pc = GLMakie.Figure(; size = canvas_size)
    gl = pc[1, 1] = GridLayout(layout[1], layout[2])

    p_idx = 1
    for idx1 in 1:layout[1], idx2 in 1:layout[2]
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
    add_pl(fig, pl; <keyword arguments>)

Overlay a locations plot `pl` onto the top-right corner of `fig`, making the white background of `pl` transparent so only the electrode markers are visible.

# Arguments

- `p1::GLMakie.Figure`: primary plot
- `p2::GLMakie.Figure`: locations plot to overlay

# Returns

- `GLMakie.Figure`: the plotted figure: `fig` with the locations overlay applied in-place
"""
function add_pl(fig::GLMakie.Figure, pl::GLMakie.Figure)::GLMakie.Figure

    # render the locations figure to an in-memory PNG stream
    io = IOBuffer()
    show(io, MIME"image/png"(), pl)

    # seek back to the beginning then load with an explicit PNG format hint
    seekstart(io)
    pp = FileIO.load(FileIO.Stream{FileIO.format"PNG"}(io))

    # make the white background transparent so only the electrode markers are composited onto the primary figure
    # three near-white values are handled to account for sub-pixel anti-aliasing on the background fill
    transparent_pp = map(c -> RGBA(color(c), 1.0), pp)
    for near_white in (
        RGBA(1.0, 1.0, 1.0, 1.0),
        RGBA(0.999, 0.999, 0.999, 1.0),
        RGBA(0.998, 0.998, 0.998, 1.0),
    )
        transparent_pp[transparent_pp .== near_white] .=
            RGBA(near_white.r, near_white.g, near_white.b, 0.0)
    end

    # determine the top-right corner of the primary axis in data coordinates
    ax = contents(fig[1, 1])[1]
    pos_x = (ax.targetlimits[].origin .+ ax.targetlimits[].widths)[1]
    pos_y = ax.targetlimits[].origin[2]

    half_size = GLMakie.Vec2f(size(transparent_pp) ./ -2)

    # scatter a single point with the locations image as its marker,
    # offset so the marker is centered on pos_x, pos_y.
    GLMakie.scatter!(
        fig[1, 1],
        pos_x,
        pos_y;
        marker_offset = half_size,
        marker = transparent_pp,
        markersize = size(transparent_pp),
        markerspace = :pixel,
    )

    return fig
end
