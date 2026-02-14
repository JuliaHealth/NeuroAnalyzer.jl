export plot_gridlocs

"""
    plot_gridlocs()

Plot a simplified plot of 10-20 EEG channels on a grid.

# Arguments

- `mono::Bool=false`: use color or gray palette

# Returns

- `p::GLMakie.Figure`
"""
function plot_gridlocs(; mono::Bool=false)::GLMakie.Figure
    pal = mono ? :grays : :darktest

    plot_size=(800, 800)
    p = GLMakie.Figure(; size=plot_size, figure_padding=0)
    ax = GLMakie.Axis(
        p[1, 1];
        aspect=1,
        xlabel="",
        ylabel="",
        title="",
        xautolimitmargin=(0, 0),
        yautolimitmargin=(0, 0),
        xzoomlock=true,
        yzoomlock=true,
        xpanlock=true,
        ypanlock=true,
        xrectzoom=false,
        yrectzoom=false,
    )
    hidedecorations!(ax; grid=true)
    hidespines!(ax)
    GLMakie.xlims!(ax, (-1.2, 1.2))
    GLMakie.ylims!(ax, (-1.2, 1.2))

    GLMakie.lines!([-1, 1], [-1, -1]; color=:black, linewidth=0.2)
    GLMakie.lines!([-1, 1], [1, 1]; color=:black, linewidth=0.2)

    GLMakie.lines!([-1, -1], [-1, 1]; color=:black, linewidth=0.2)
    GLMakie.lines!([1, 1], [-1, 1]; color=:black, linewidth=0.2)

    GLMakie.lines!([-1, -0.5], [0.5, 1]; color=:black, linewidth=0.5)
    GLMakie.lines!([0.5, 1], [1, 0.5]; color=:black, linewidth=0.5)
    GLMakie.lines!([-1, -0.5], [-0.5, -1]; color=:black, linewidth=0.5)
    GLMakie.lines!([0.5, 1], [-1, -0.5]; color=:black, linewidth=0.5)

    GLMakie.lines!([-0.5, 0.5], [-1, -1]; color=:black, linewidth=0.5)
    GLMakie.lines!([-1, 1], [-0.5, -0.5]; color=:black, linewidth=0.5)
    GLMakie.lines!([-1, 1], [0, 0]; color=:black, linewidth=0.5)
    GLMakie.lines!([-1, 1], [0.5, 0.5]; color=:black, linewidth=0.5)
    GLMakie.lines!([-0.5, 0.5], [1, 1]; color=:black, linewidth=0.5)

    GLMakie.lines!([-1, -1], [-0.5, 0.5]; color=:black, linewidth=0.5)
    GLMakie.lines!([-0.5, -0.5], [-1, 1]; color=:black, linewidth=0.5)
    GLMakie.lines!([0, 0], [-1, 1]; color=:black, linewidth=0.5)
    GLMakie.lines!([0.5, 0.5], [-1, 1]; color=:black, linewidth=0.5)
    GLMakie.lines!([1, 1], [-0.5, 0.5]; color=:black, linewidth=0.5)

    loc_x = [
        -0.5,
        0,
        0.5,
        -1,
        -0.5,
        0,
        0.5,
        1,
        -1,
        -0.5,
        0,
        0.5,
        1,
        -1,
        -0.5,
        0,
        0.5,
        1,
        -0.5,
        0,
        0.5,
    ]
    loc_y = [
        1,
        1,
        1,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0,
        0,
        0,
        0,
        0,
        -0.5,
        -0.5,
        -0.5,
        -0.5,
        -0.5,
        -1,
        -1,
        -1,
    ]
    loc_lab = [
        "Fp1",
        "Fpz",
        "Fp2",
        "F7",
        "F3",
        "Fz",
        "F4",
        "F8",
        "T3",
        "C3",
        "Cz",
        "C4",
        "T4",
        "T5",
        "P3",
        "Pz",
        "P4",
        "T6",
        "O1",
        "Oz",
        "O2",
    ]
    font_size = 16
    label_offset_x = 0.015
    label_offset_y = 0.015

    ch_n = length(loc_lab)
    cmap = GLMakie.resample_cmap(pal, ch_n)

    for idx in eachindex(loc_x)
        GLMakie.scatter!(
            loc_x[idx],
            loc_y[idx];
            colormap=pal,
            color=cmap[idx],
            colorrange=1:ch_n,
            markersize=16.0,
            strokewidth=2,
            strokecolor=:black,
        )
        GLMakie.text!(
            loc_x[idx] + label_offset_x,
            loc_y[idx] + label_offset_y;
            text=loc_lab[idx],
            fontsize=font_size,
        )
    end

    return p
end
