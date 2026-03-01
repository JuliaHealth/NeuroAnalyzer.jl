export plot_dipole2d

"""
    plot_dipole2d(d; <keyword arguments>)

Plot dipole in 2D.

# Arguments

  - `d::NeuroAnalyzer.DIPOLE`

# Returns

  - `p::GLMakie.Figure`

# Notes

Brain volume is within -1.0 to +1.0 (x-, y- and z-axis)
"""
function plot_dipole2d(d::NeuroAnalyzer.DIPOLE)::GLMakie.Figure

    # get dipole position
    x = d.pos[1]
    y = d.pos[2]
    z = d.pos[3]
    pos = (x, y, z)

    # get dipole magnitude
    mx = d.mag[1]
    my = d.mag[2]
    mz = d.mag[3]
    mag = (mx, my, mz)

    # prepare plot
    GLMakie.activate!(title = "plot_dipole_2d()")
    plot_size = (1200, 400)
    p = GLMakie.Figure(size = plot_size)
    ax_xy = GLMakie.Axis(
        p[1, 1];
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
        p[1, 2];
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
        p[1, 3];
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

    # draw head
    lw = 2
    # nose
    GLMakie.lines!(ax_xy, [-0.2, 0], [0.98, 1.08]; linewidth = lw, color = :black)
    GLMakie.lines!(ax_xy, [0.2, 0], [0.98, 1.08]; linewidth = lw, color = :black)
    # ears
    # left
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
    # right
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
    # head
    GLMakie.arc!(ax_xy, (0, 0), 1, 0, 2pi; linewidth = lw, color = :black)

    # head
    GLMakie.arc!(ax_yz, (0, 0), 1, 0, pi; linewidth = lw, color = :black)

    # head
    GLMakie.arc!(ax_xz, (0, 0), 1, 0, pi; linewidth = lw, color = :black)

    # draw dipole position
    GLMakie.scatter!(ax_xy, pos[1], pos[2]; markersize = sqrt(sum(mag .^ 2)) * 20, color = :red)
    GLMakie.scatter!(ax_yz, pos[2], pos[3]; markersize = sqrt(sum(mag .^ 2)) * 20, color = :red)
    GLMakie.scatter!(ax_xz, pos[1], pos[3]; markersize = sqrt(sum(mag .^ 2)) * 20, color = :red)

    # draw magnitude
    GLMakie.lines!(ax_xy, [pos[1], pos[1] + mag[1]], [pos[2], pos[2] + mag[2]]; color = :red)
    GLMakie.lines!(ax_yz, [pos[2], pos[2] + mag[2]], [pos[3], pos[3] + mag[3]]; color = :red)
    GLMakie.lines!(ax_xz, [pos[1], pos[1] + mag[1]], [pos[3], pos[3] + mag[3]]; color = :red)

    return p

end
