export plot_dipole2d

"""
    plot_dipole2d(d; <keyword arguments>)

Plot a 3D dipole in 2D views (top, side, front).

# Arguments

- `d::NeuroAnalyzer.DIPOLE`: dipole object with `pos` (position) and `mag` (magnitude) fields
    - `pos::Tuple{Float64, Float64, Float64}`: dipole position (x, y, z) in brain volume (range: -1.0 to +1.0)
    - `mag::Tuple{Float64, Float64, Float64}`: dipole magnitude vector (mx, my, mz)

# Returns

- `GLMakie.Figure`: the plotted figure

# Notes

- Brain volume is within -1.0 to +1.0 for all axes (x, y, z).
- The dipole position is marked with a red dot, and its magnitude is shown as a red line.
"""
function plot_dipole2d(d::NeuroAnalyzer.DIPOLE)::GLMakie.Figure

    # validate
    all(-1.0 .≤ d.pos .≤ 1.0) ||
        throw(ArgumentError("Position must be within [-1.0, 1.0]."))
    all(-1.0 .≤ d.mag .≤ 1.0) ||
        throw(ArgumentError("Magnitude must be within [-1.0, 1.0]."))

    # extract position and magnitude
    x, y, z = d.pos
    mx, my, mz = d.mag

    # prepare plot
    GLMakie.activate!(; title = "plot_dipole_2d()")
    plot_size = (1200, 400)
    fig = GLMakie.Figure(; size = plot_size)

    # create three subplots: top (xy), side (yz), front (xz)

    ax_xy = GLMakie.Axis(
        fig[1, 1];
        aspect = DataAspect(),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
        title = "Top view",
    )
    hidedecorations!(ax_xy)
    hidespines!(ax_xy)
    GLMakie.xlims!(ax_xy, -1.2, 1.2)
    GLMakie.ylims!(ax_xy, -1.2, 1.2)

    ax_yz = GLMakie.Axis(
        fig[1, 2];
        aspect = DataAspect(),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
        title = "Side view",
    )
    hidedecorations!(ax_yz)
    hidespines!(ax_yz)
    GLMakie.xlims!(ax_yz, -1.2, 1.2)
    GLMakie.ylims!(ax_yz, -1.2, 1.2)

    ax_xz = GLMakie.Axis(
        fig[1, 3];
        aspect = DataAspect(),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
        title = "Front view",
    )
    hidedecorations!(ax_xz)
    hidespines!(ax_xz)
    GLMakie.xlims!(ax_xz, -1.2, 1.2)
    GLMakie.ylims!(ax_xz, -1.2, 1.2)

    # draw head outline (nose, ears, head)
    lw = 2
    # nose
    GLMakie.lines!(ax_xy, [-0.2, 0], [0.98, 1.08]; linewidth = lw, color = :black)
    GLMakie.lines!(ax_xy, [0.2, 0], [0.98, 1.08]; linewidth = lw, color = :black)
    # left ear
    GLMakie.lines!(ax_xy, [-0.995, -1.03], [0.1, 0.15]; linewidth = lw, color = :black)
    GLMakie.lines!(ax_xy, [-1.03, -1.06], [0.15, 0.16]; linewidth = lw, color = :black)
    GLMakie.lines!(ax_xy, [-1.06, -1.1], [0.16, 0.14]; linewidth = lw, color = :black)
    GLMakie.lines!(ax_xy, [-1.1, -1.12], [0.14, 0.05]; linewidth = lw, color = :black)
    GLMakie.lines!(ax_xy, [-1.12, -1.1], [0.05, -0.1]; linewidth = lw, color = :black)
    GLMakie.lines!(ax_xy, [-1.1, -1.13], [-0.1, -0.3]; linewidth = lw, color = :black)
    GLMakie.lines!(ax_xy, [-1.13, -1.09], [-0.3, -0.37]; linewidth = lw, color = :black)
    GLMakie.lines!(ax_xy, [-1.09, -1.02], [-0.37, -0.39]; linewidth = lw, color = :black)
    GLMakie.lines!(ax_xy, [-1.02, -0.98], [-0.39, -0.33]; linewidth = lw, color = :black)
    GLMakie.lines!(ax_xy, [-0.98, -0.975], [-0.33, -0.22]; linewidth = lw, color = :black)
    # right ear
    GLMakie.lines!(ax_xy, [0.995, 1.03], [0.1, 0.15]; linewidth = lw, color = :black)
    GLMakie.lines!(ax_xy, [1.03, 1.06], [0.15, 0.16]; linewidth = lw, color = :black)
    GLMakie.lines!(ax_xy, [1.06, 1.1], [0.16, 0.14]; linewidth = lw, color = :black)
    GLMakie.lines!(ax_xy, [1.1, 1.12], [0.14, 0.05]; linewidth = lw, color = :black)
    GLMakie.lines!(ax_xy, [1.12, 1.1], [0.05, -0.1]; linewidth = lw, color = :black)
    GLMakie.lines!(ax_xy, [1.1, 1.13], [-0.1, -0.3]; linewidth = lw, color = :black)
    GLMakie.lines!(ax_xy, [1.13, 1.09], [-0.3, -0.37]; linewidth = lw, color = :black)
    GLMakie.lines!(ax_xy, [1.09, 1.02], [-0.37, -0.39]; linewidth = lw, color = :black)
    GLMakie.lines!(ax_xy, [1.02, 0.98], [-0.39, -0.33]; linewidth = lw, color = :black)
    GLMakie.lines!(ax_xy, [0.98, 0.975], [-0.33, -0.22]; linewidth = lw, color = :black)
    # head outline
    GLMakie.arc!(ax_xy, (0, 0), 1, 0, 2pi; linewidth = lw, color = :black)

    # head outline
    GLMakie.arc!(ax_yz, (0, 0), 1, 0, pi; linewidth = lw, color = :black)

    # head outline
    GLMakie.arc!(ax_xz, (0, 0), 1, 0, pi; linewidth = lw, color = :black)

    # draw dipole position and magnitude
    dipole_size = sqrt(sum(d.mag .^ 2)) * 20
    GLMakie.scatter!(ax_xy, x, y; markersize = dipole_size, color = :red)
    GLMakie.scatter!(ax_yz, y, z; markersize = dipole_size, color = :red)
    GLMakie.scatter!(ax_xz, x, z; markersize = dipole_size, color = :red)

    # draw dipole magnitude vectors
    GLMakie.lines!(ax_xy, [x, x + mx], [y, y + my]; color = :red)
    GLMakie.lines!(ax_yz, [y, y + my], [z, z + mz]; color = :red)
    GLMakie.lines!(ax_xz, [x, x + mx], [z, z + mz]; color = :red)

    return fig
end
