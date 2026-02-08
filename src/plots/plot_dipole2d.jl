export plot_dipole2d

"""
    plot_dipole2d(d; <keyword arguments>)

Plot dipole in 2D.

# Arguments

- `d::NeuroAnalyzer.DIPOLE`

# Returns

- `p::GLMakie.Figure`

# Notes

Brain volume is within -1.0 to +1.0 (X-, Y- and Z-axis)
"""
function plot_dipole2d(d::NeuroAnalyzer.DIPOLE)::GLMakie.Figure

    _wip()

    # load textures
    head_top_texture = rotr90(FileIO.load(joinpath(res_path, "head_t.png")))
    head_side_texture = rotr90(FileIO.load(joinpath(res_path, "head_s.png")))
    head_front_texture = rotr90(FileIO.load(joinpath(res_path, "head_f.png")))

    # scale according to textures
    xy_s1 = size(head_top_texture, 2)
    xy_s2 = size(head_top_texture, 1)
    yz_s1 = size(head_side_texture, 2)
    yz_s2 = size(head_side_texture, 1)
    xz_s1 = size(head_front_texture, 2)
    xz_s2 = size(head_front_texture, 1)

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
    plot_size = (600, 200)
    p = GLMakie.Figure(size=plot_size)
    ax_xy = GLMakie.Axis(p[1, 1],
                         aspect=DataAspect(),
                         xzoomlock=true,
                         yzoomlock=true,
                         xpanlock=true,
                         ypanlock=true,
                         xrectzoom=false,
                         yrectzoom=false)
    hidedecorations!(ax_xy)
    hidespines!(ax_xy)
    GLMakie.image!(ax_xy,
                   head_top_texture)

    ax_yz = GLMakie.Axis(p[1, 2],
                         aspect=DataAspect(),
                         xzoomlock=true,
                         yzoomlock=true,
                         xpanlock=true,
                         ypanlock=true,
                         xrectzoom=false,
                         yrectzoom=false)
    hidedecorations!(ax_yz)
    hidespines!(ax_yz)
    GLMakie.image!(ax_yz,
                   head_side_texture)

    ax_xz = GLMakie.Axis(p[1, 3],
                         aspect=DataAspect(),
                         xzoomlock=true,
                         yzoomlock=true,
                         xpanlock=true,
                         ypanlock=true,
                         xrectzoom=false,
                         yrectzoom=false)
    hidedecorations!(ax_xz)
    hidespines!(ax_xz)
    GLMakie.image!(ax_xz,
                   head_front_texture)

    # draw dipole position
    GLMakie.scatter!(ax_xy,
                     (xy_s2 / 2) + (pos[1] * (xy_s2 / 2)),
                     (xy_s1 / 2) + (pos[2] * (xy_s1 / 2)),
                     color=:red)
    GLMakie.scatter!(ax_yz,
                     (yz_s2 / 2) + (pos[2] * (yz_s2 / 2)),
                     (yz_s1 / 2) + (pos[3] * (yz_s1 / 2)),
                     color=:red)
    GLMakie.scatter!(ax_xz,
                     (xz_s2 / 2) + (pos[1] * (xz_s2 / 2)),
                     (xz_s1 / 2) + (pos[3] * (xz_s1 / 2)),
                     color=:red)
    # draw magnitude
    GLMakie.lines!(ax_xy,
                   [(xy_s2 / 2) + (pos[1] * (xy_s2 / 2)), (xy_s2 / 2) + (pos[1] * (xy_s2 / 2)) + (mag[1] * (xy_s2 / 2))], 
                   [(xy_s1 / 2) + (pos[2] * (xy_s1 / 2)), (xy_s1 / 2) + (pos[2] * (xy_s1 / 2)) + (mag[2] * (xy_s1 / 2))],
                   color=:red)
    GLMakie.lines!(ax_yz,
                   [(yz_s2 / 2) + (pos[2] * (yz_s2 / 2)), (yz_s2 / 2) + (pos[2] * (yz_s2 / 2)) + (mag[2] * (yz_s2 / 2))], 
                   [(yz_s1 / 2) + (pos[3] * (yz_s1 / 2)), (yz_s1 / 2) + (pos[3] * (yz_s1 / 2)) + (mag[3] * (yz_s1 / 2))],
                   color=:red)
    GLMakie.lines!(ax_xz,
                   [(xz_s2 / 2) + (pos[1] * (xz_s2 / 2)), (xz_s2 / 2) + (pos[1] * (xz_s2 / 2)) + (mag[1] * (xz_s2 / 2))], 
                   [(xz_s1 / 2) + (pos[3] * (xz_s1 / 2)), (xz_s1 / 2) + (pos[3] * (xz_s1 / 2)) + (mag[3] * (xz_s1 / 2))],
                   color=:red)

    return p

end
